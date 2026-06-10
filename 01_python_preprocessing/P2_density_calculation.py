#!/usr/bin/env python3
"""
P2_density_calculation.py
=========================
Compute 5 per-cell spatial density estimators for every sample,
grouped by dataset: 3 KNN-based densities with data-driven K
selection (2nd-difference / piecewise-regression / max-distance knee
detection), plus Voronoi and Delaunay densities. Inputs are the P1
cell-metadata coordinates; outputs are per-sample density tables,
diagnostics and QC summaries.

Paths: DATA_DIR (registry) and RESULTS_DIR (P1 input / P2 output);
see config/paths.py.

重构版 P2：三方法密度计算 + 多K值策略

功能：
  阶段0：按数据集分组（从 sample_registry.json）
  阶段1：对每个数据集，用三种拐点检测方法确定三个K值
         - 二阶差分法 (aggressive)
         - 分段回归法 (main/primary)
         - 最大距离法 (conservative)
  阶段2：逐样本计算5种密度
         - density_knn_aggr_2nd_diff
         - density_knn_main_piecewise
         - density_knn_cons_max_dist
         - density_voronoi
         - density_delaunay
  阶段3：汇总QC报告

输入：
  - sample_registry.json（全局配置）
  - P1_Results/{sample}/cell_metadata.csv（每样本的细胞坐标）

输出：
  - P2_Results/{dataset}/KNN_Optimization/（K值选择诊断）
  - P2_Results/{sample}/cell_density.csv（5列密度值）
  - P2_Results/{sample}/density_diagnostics.png（空间分布图）
  - P2_Results/{sample}/density_qc.csv（样本级QC）
  - P2_Results/ALL_SAMPLES_P2_QC.csv（汇总QC）
  - P2_Results/ALL_DATASETS_K_SELECTION.csv（K值选择汇总）

运行方式：
  cd 01_python_preprocessing
  nohup python3 P2_density_calculation.py > P2_run.log 2>&1 &
  tail -f P2_run.log
"""

import os
import json
import time
import warnings
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, Delaunay
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================
# 全局配置
# =============================================
# Paths from environment variables (see config/paths.py); no hardcoded home dirs.
REGISTRY_PATH = os.path.join(os.environ.get("DATA_DIR", "./data"), "sample_registry.json")
P1_RESULTS_DIR = os.path.join(os.environ.get("RESULTS_DIR", "./results"), "P1_Results")
P2_RESULTS_DIR = os.path.join(os.environ.get("RESULTS_DIR", "./results"), "P2_Results")

# 候选K值：前半段步长2（精细扫描拐点区），后半段步长5（粗扫平稳区）
K_CANDIDATES = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 30, 35, 40, 45, 50]

# 生物学参考尺度（仅用于后验解读，不用于约束K值选择）
BIO_SCALES = [
    (0, 10, "sub-cellular"),
    (10, 30, "direct contact (juxtacrine)"),
    (30, 200, "paracrine signaling"),
    (200, 1000, "tissue-level architecture"),
    (1000, float("inf"), "macro-scale"),
]

# =============================================
# 工具函数
# =============================================

def classify_bio_scale(distance_um):
    """将物理距离映射到生物学尺度"""
    for lo, hi, label in BIO_SCALES:
        if lo <= distance_um < hi:
            return label
    return "unknown"


def find_k_2nd_diff(cv_values, k_candidates):
    """方法1：二阶差分法（aggressive）
    找CV曲线弯曲最剧烈的位置 = 二阶差分绝对值最大处
    """
    if len(cv_values) < 3:
        return k_candidates[0]
    first_diff = np.diff(cv_values)
    second_diff = np.diff(first_diff)
    idx = np.argmax(np.abs(second_diff))
    return k_candidates[idx + 1]  # +1 因为二阶差分损失了两个元素


def find_k_piecewise(cv_values, k_candidates):
    """方法2：分段线性回归法（main）
    假设CV曲线由两段直线组成（急降+平稳），找总拟合误差最小的断点
    """
    k_arr = np.array(k_candidates, dtype=float)
    cv_arr = np.array(cv_values, dtype=float)
    n = len(k_arr)

    if n < 4:
        return k_candidates[n // 2]

    best_break = 2
    best_error = float("inf")

    # 遍历所有可能的断点位置（至少左右各2个点）
    for bp in range(2, n - 1):
        # 左段线性回归
        x_left = k_arr[:bp + 1]
        y_left = cv_arr[:bp + 1]
        if len(x_left) >= 2:
            coef_left = np.polyfit(x_left, y_left, 1)
            pred_left = np.polyval(coef_left, x_left)
            err_left = np.sum((y_left - pred_left) ** 2)
        else:
            err_left = 0

        # 右段线性回归
        x_right = k_arr[bp:]
        y_right = cv_arr[bp:]
        if len(x_right) >= 2:
            coef_right = np.polyfit(x_right, y_right, 1)
            pred_right = np.polyval(coef_right, x_right)
            err_right = np.sum((y_right - pred_right) ** 2)
        else:
            err_right = 0

        total_error = err_left + err_right
        if total_error < best_error:
            best_error = total_error
            best_break = bp

    return k_candidates[best_break]


def find_k_max_distance(cv_values, k_candidates):
    """方法3：最大距离法 / Kneedle（conservative）
    将曲线归一化到[0,1]x[0,1]，找到曲线到对角线距离最大的点
    """
    k_arr = np.array(k_candidates, dtype=float)
    cv_arr = np.array(cv_values, dtype=float)

    # 归一化到 [0, 1]
    k_norm = (k_arr - k_arr.min()) / (k_arr.max() - k_arr.min() + 1e-12)
    cv_norm = (cv_arr - cv_arr.min()) / (cv_arr.max() - cv_arr.min() + 1e-12)

    # 对角线从 (0, 1) 到 (1, 0)：方程 x + y = 1，即 y = 1 - x
    # 点到直线 x + y - 1 = 0 的距离 = |x_i + y_i - 1| / sqrt(2)
    distances = np.abs(k_norm + cv_norm - 1.0) / np.sqrt(2.0)

    # CV曲线是下降的，归一化后曲线在对角线上方，找最大距离
    idx = np.argmax(distances)
    return k_candidates[idx]


def compute_voronoi_density(coords):
    """计算Voronoi密度
    密度 = 1 / Voronoi多边形面积
    边界处开放多边形标记为NaN
    """
    n = len(coords)
    density = np.full(n, np.nan)

    if n < 4:
        return density

    try:
        vor = Voronoi(coords)
    except Exception:
        return density

    for i in range(n):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]

        # 开放多边形（含-1顶点）→ NaN
        if -1 in region or len(region) == 0:
            continue

        # 获取多边形顶点
        vertices = vor.vertices[region]

        # Shoelace公式计算面积
        x = vertices[:, 0]
        y = vertices[:, 1]
        area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

        if area > 0:
            density[i] = 1.0 / area

    # 异常大面积（极小密度）→ NaN（99百分位截断）
    valid = density[~np.isnan(density)]
    if len(valid) > 0:
        p1 = np.percentile(valid, 1)
        density[density < p1] = np.nan

    return density


def compute_delaunay_density(coords):
    """计算Delaunay密度
    密度 = 1 / 该细胞所有Delaunay邻居的平均边长
    """
    n = len(coords)
    density = np.full(n, np.nan)

    if n < 4:
        return density

    try:
        tri = Delaunay(coords)
    except Exception:
        return density

    # 构建邻接表：每个点的所有Delaunay邻居
    neighbor_distances = [[] for _ in range(n)]
    simplices = tri.simplices  # 每行三个顶点索引

    for simplex in simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                pi, pj = simplex[i], simplex[j]
                dist = np.linalg.norm(coords[pi] - coords[pj])
                neighbor_distances[pi].append(dist)
                neighbor_distances[pj].append(dist)

    for i in range(n):
        if len(neighbor_distances[i]) > 0:
            mean_dist = np.mean(neighbor_distances[i])
            if mean_dist > 0:
                density[i] = 1.0 / mean_dist

    # 异常长边（极小密度）→ NaN（1百分位截断）
    valid = density[~np.isnan(density)]
    if len(valid) > 0:
        p1 = np.percentile(valid, 1)
        density[density < p1] = np.nan

    return density


def plot_k_optimization(dataset_name, k_candidates, sample_cv_data,
                        mean_cv, k_aggr, k_main, k_cons, save_path):
    """绘制K值优化诊断图（四面板）"""
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # --- 面板1：CV vs K 曲线 + 三种方法标记 ---
    ax1 = fig.add_subplot(gs[0, 0])
    for sample_name, cv_vals in sample_cv_data.items():
        ax1.plot(k_candidates, cv_vals, "o-", alpha=0.3, markersize=3,
                 label=sample_name if len(sample_cv_data) <= 7 else None)
    ax1.plot(k_candidates, mean_cv, "k-", linewidth=2.5, label="Mean CV", zorder=10)
    ax1.axvline(k_aggr, color="red", linestyle="--", linewidth=1.5,
                label=f"2nd_diff K={k_aggr}")
    ax1.axvline(k_main, color="blue", linestyle="-", linewidth=2,
                label=f"Piecewise K={k_main}")
    ax1.axvline(k_cons, color="green", linestyle="--", linewidth=1.5,
                label=f"Max_dist K={k_cons}")
    ax1.set_xlabel("K value", fontsize=12)
    ax1.set_ylabel("CV (Coefficient of Variation)", fontsize=12)
    ax1.set_title(f"{dataset_name}: CV vs K", fontsize=14)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, alpha=0.3)

    # --- 面板2：中位距离 vs K ---
    ax2 = fig.add_subplot(gs[0, 1])
    # 中位距离从sample_cv_data的附加数据中获取（在调用时传入）
    ax2.set_xlabel("K value", fontsize=12)
    ax2.set_ylabel("Median k-th NN distance (µm)", fontsize=12)
    ax2.set_title(f"{dataset_name}: Physical Distance vs K", fontsize=14)
    ax2.text(0.5, 0.5, "See distance data in\nk_decision_table.csv",
             transform=ax2.transAxes, ha="center", va="center", fontsize=12, alpha=0.5)
    ax2.grid(True, alpha=0.3)

    # --- 面板3：一阶差分 ---
    ax3 = fig.add_subplot(gs[1, 0])
    first_diff = np.diff(mean_cv)
    k_mid = [(k_candidates[i] + k_candidates[i + 1]) / 2 for i in range(len(first_diff))]
    colors3 = ["red" if k_candidates[i + 1] == k_aggr else
               "blue" if k_candidates[i + 1] == k_main else
               "green" if k_candidates[i + 1] == k_cons else "gray"
               for i in range(len(first_diff))]
    ax3.bar(k_mid, first_diff, width=1.5, color=colors3, alpha=0.7)
    ax3.set_xlabel("K value (midpoint)", fontsize=12)
    ax3.set_ylabel("1st order difference (ΔCV)", fontsize=12)
    ax3.set_title(f"{dataset_name}: First Derivative of CV", fontsize=14)
    ax3.axhline(0, color="black", linewidth=0.5)
    ax3.grid(True, alpha=0.3)

    # --- 面板4：二阶差分 ---
    ax4 = fig.add_subplot(gs[1, 1])
    second_diff = np.diff(first_diff)
    k_mid2 = [(k_candidates[i + 1] + k_candidates[i + 2]) / 2
              for i in range(len(second_diff))]
    colors4 = ["red" if k_candidates[i + 1] == k_aggr else
               "blue" if k_candidates[i + 1] == k_main else
               "green" if k_candidates[i + 1] == k_cons else "gray"
               for i in range(len(second_diff))]
    ax4.bar(k_mid2, np.abs(second_diff), width=1.5, color=colors4, alpha=0.7)
    ax4.set_xlabel("K value (midpoint)", fontsize=12)
    ax4.set_ylabel("|2nd order difference|", fontsize=12)
    ax4.set_title(f"{dataset_name}: Second Derivative of CV (absolute)", fontsize=14)
    ax4.grid(True, alpha=0.3)

    plt.suptitle(f"K-value Optimization: {dataset_name}", fontsize=16, fontweight="bold", y=0.98)
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_density_diagnostics(sample_name, coords, densities, method_names,
                             corr_matrix, save_path):
    """绘制样本密度诊断图（上排：5个空间分布图，下排：5个直方图 + 相关矩阵）"""
    n_methods = len(method_names)
    fig = plt.figure(figsize=(4 * n_methods + 4, 10))
    gs = GridSpec(2, n_methods + 1, figure=fig, hspace=0.35, wspace=0.3)

    # 上排：空间分布图
    for i, (name, dens) in enumerate(zip(method_names, densities)):
        ax = fig.add_subplot(gs[0, i])
        valid = ~np.isnan(dens)
        if valid.sum() > 0:
            # 对密度取log以增强可视化对比度
            log_dens = np.log10(dens[valid] + 1e-20)
            vmin, vmax = np.percentile(log_dens, [2, 98])
            sc = ax.scatter(coords[valid, 0], coords[valid, 1],
                            c=log_dens, cmap="viridis", s=0.1, alpha=0.5,
                            vmin=vmin, vmax=vmax, rasterized=True)
            plt.colorbar(sc, ax=ax, shrink=0.6, label="log10(density)")
        ax.set_title(name, fontsize=9)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=6)

    # 下排：直方图
    for i, (name, dens) in enumerate(zip(method_names, densities)):
        ax = fig.add_subplot(gs[1, i])
        valid = dens[~np.isnan(dens)]
        if len(valid) > 0:
            log_valid = np.log10(valid + 1e-20)
            ax.hist(log_valid, bins=50, color="steelblue", alpha=0.7, edgecolor="none")
        ax.set_xlabel("log10(density)", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.set_title(f"{name}\nN_valid={np.sum(~np.isnan(dens))}", fontsize=8)
        ax.tick_params(labelsize=6)

    # 最后一格：相关矩阵热力图
    ax_corr = fig.add_subplot(gs[1, n_methods])
    im = ax_corr.imshow(corr_matrix, cmap="RdYlGn", vmin=0.5, vmax=1.0)
    ax_corr.set_xticks(range(n_methods))
    ax_corr.set_yticks(range(n_methods))
    short_names = [n.replace("density_", "").replace("knn_", "K_") for n in method_names]
    ax_corr.set_xticklabels(short_names, fontsize=6, rotation=45, ha="right")
    ax_corr.set_yticklabels(short_names, fontsize=6)
    for ii in range(n_methods):
        for jj in range(n_methods):
            ax_corr.text(jj, ii, f"{corr_matrix[ii, jj]:.3f}",
                         ha="center", va="center", fontsize=6,
                         color="black" if corr_matrix[ii, jj] > 0.7 else "white")
    plt.colorbar(im, ax=ax_corr, shrink=0.6, label="Spearman ρ")
    ax_corr.set_title("Method correlation", fontsize=9)

    plt.suptitle(f"Density Diagnostics: {sample_name}", fontsize=13, fontweight="bold", y=0.99)
    plt.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close()


# =============================================
# 主流程
# =============================================

def main():
    start_time = time.time()
    print("=" * 70)
    print("P2 密度计算（重构版 - 三K值策略）")
    print("=" * 70)

    # 读取 sample_registry
    with open(REGISTRY_PATH, "r") as f:
        registry = json.load(f)
    print(f"读取 sample_registry: {len(registry)} 个样本")

    os.makedirs(P2_RESULTS_DIR, exist_ok=True)

    # ==========================================
    # 阶段 0：按数据集分组
    # ==========================================
    print("\n" + "=" * 70)
    print("阶段 0：按数据集分组")
    print("=" * 70)

    dataset_groups = {}
    for sample_name, info in registry.items():
        dataset = sample_name.split("/")[0]
        if dataset not in dataset_groups:
            dataset_groups[dataset] = []
        dataset_groups[dataset].append(sample_name)

    for ds, samples in dataset_groups.items():
        print(f"  {ds}: {len(samples)} 个样本")

    # ==========================================
    # 阶段 1：为每个数据集确定三个K值
    # ==========================================
    print("\n" + "=" * 70)
    print("阶段 1：K值优化（三种拐点检测方法）")
    print("=" * 70)

    # 存储每个数据集的K值决定
    dataset_k_decisions = {}
    # 存储所有样本的CV数据（供诊断图使用）
    all_cv_data = {}

    for ds_name, sample_list in dataset_groups.items():
        print(f"\n--- 数据集: {ds_name} ({len(sample_list)} 个样本) ---")

        ds_opt_dir = os.path.join(P2_RESULTS_DIR, ds_name, "KNN_Optimization")
        os.makedirs(ds_opt_dir, exist_ok=True)

        # 对每个样本扫描所有候选K的CV和中位距离
        sample_cv_dict = {}   # {sample_name: [cv_for_each_k]}
        sample_dist_dict = {} # {sample_name: [median_dist_for_each_k]}

        for sample_name in sample_list:
            # 读取坐标
            meta_path = os.path.join(P1_RESULTS_DIR,
                                     sample_name.split("/")[0],
                                     sample_name.split("/")[1],
                                     "cell_metadata.csv")
            df = pd.read_csv(meta_path)
            coords = df[["x_centroid", "y_centroid"]].values
            n_cells = len(coords)

            # BallTree 一次性算到 K_max
            k_max = max(K_CANDIDATES)
            nn = NearestNeighbors(n_neighbors=k_max + 1, algorithm="ball_tree")
            nn.fit(coords)
            distances, _ = nn.kneighbors(coords)
            # distances[:, 0] 是自己到自己的距离(0)，第k近邻在 distances[:, k]

            cv_list = []
            dist_list = []
            for k in K_CANDIDATES:
                dist_k = distances[:, k]  # 第k近邻距离
                # 处理距离=0的极端情况
                dist_k_safe = dist_k.copy()
                dist_k_safe[dist_k_safe == 0] = np.finfo(float).eps
                density_k = 1.0 / dist_k_safe

                cv = np.std(density_k) / np.mean(density_k) if np.mean(density_k) > 0 else 0
                median_dist = np.median(dist_k)

                cv_list.append(cv)
                dist_list.append(median_dist)

            sample_cv_dict[sample_name] = cv_list
            sample_dist_dict[sample_name] = dist_list
            print(f"  扫描完成: {sample_name} ({n_cells} cells)")

        # 计算跨样本平均CV
        cv_matrix = np.array(list(sample_cv_dict.values()))  # (n_samples, n_k)
        mean_cv = np.mean(cv_matrix, axis=0)

        # 计算跨样本平均中位距离
        dist_matrix = np.array(list(sample_dist_dict.values()))
        mean_dist = np.mean(dist_matrix, axis=0)

        # 三种方法选K
        k_aggr = find_k_2nd_diff(mean_cv, K_CANDIDATES)
        k_main = find_k_piecewise(mean_cv, K_CANDIDATES)
        k_cons = find_k_max_distance(mean_cv, K_CANDIDATES)

        # 确保三个K值的排序：aggr <= main <= cons
        # 如果方法输出不满足这个预期，调整
        k_sorted = sorted([k_aggr, k_main, k_cons])
        k_aggr, k_main, k_cons = k_sorted[0], k_sorted[1], k_sorted[2]

        # 获取三个K对应的平均中位距离
        idx_aggr = K_CANDIDATES.index(k_aggr)
        idx_main = K_CANDIDATES.index(k_main)
        idx_cons = K_CANDIDATES.index(k_cons)
        dist_aggr = mean_dist[idx_aggr]
        dist_main = mean_dist[idx_main]
        dist_cons = mean_dist[idx_cons]

        # 生物学尺度分类
        bio_aggr = classify_bio_scale(dist_aggr)
        bio_main = classify_bio_scale(dist_main)
        bio_cons = classify_bio_scale(dist_cons)

        dataset_k_decisions[ds_name] = {
            "k_aggr_2nd_diff": k_aggr,
            "k_main_piecewise": k_main,
            "k_cons_max_dist": k_cons,
            "dist_aggr_um": round(dist_aggr, 1),
            "dist_main_um": round(dist_main, 1),
            "dist_cons_um": round(dist_cons, 1),
            "bio_scale_aggr": bio_aggr,
            "bio_scale_main": bio_main,
            "bio_scale_cons": bio_cons,
            "n_samples": len(sample_list),
        }

        print(f"\n  K值选择结果:")
        print(f"    Aggressive (2nd_diff):   K={k_aggr:3d}  →  {dist_aggr:.1f} µm  [{bio_aggr}]")
        print(f"    Main (piecewise):        K={k_main:3d}  →  {dist_main:.1f} µm  [{bio_main}]")
        print(f"    Conservative (max_dist): K={k_cons:3d}  →  {dist_cons:.1f} µm  [{bio_cons}]")

        # 保存K值决策表
        k_table_rows = []
        for i, k in enumerate(K_CANDIDATES):
            row = {
                "K": k,
                "mean_CV": round(mean_cv[i], 6),
                "std_CV": round(np.std(cv_matrix[:, i]), 6),
                "mean_median_distance_um": round(mean_dist[i], 2),
                "is_aggr_2nd_diff": k == k_aggr,
                "is_main_piecewise": k == k_main,
                "is_cons_max_dist": k == k_cons,
            }
            k_table_rows.append(row)
        pd.DataFrame(k_table_rows).to_csv(
            os.path.join(ds_opt_dir, "k_decision_table.csv"), index=False
        )

        # 保存所有样本CV数据（供后续分析和诊断）
        cv_save = {"K": K_CANDIDATES}
        for sn, cv_vals in sample_cv_dict.items():
            cv_save[sn] = cv_vals
        cv_save["mean"] = mean_cv.tolist()
        pd.DataFrame(cv_save).to_csv(
            os.path.join(ds_opt_dir, "all_samples_knn_cv.csv"), index=False
        )

        # 绘制诊断图
        plot_k_optimization(
            ds_name, K_CANDIDATES, sample_cv_dict,
            mean_cv, k_aggr, k_main, k_cons,
            os.path.join(ds_opt_dir, "k_optimization_diagnostic.png")
        )

        # 存储CV数据供全局使用
        all_cv_data[ds_name] = {
            "sample_cv": sample_cv_dict,
            "mean_cv": mean_cv,
        }

    # 保存全局K值选择汇总
    k_summary_rows = []
    for ds_name, info in dataset_k_decisions.items():
        row = {"dataset": ds_name}
        row.update(info)
        k_summary_rows.append(row)
    k_summary_df = pd.DataFrame(k_summary_rows)
    k_summary_df.to_csv(
        os.path.join(P2_RESULTS_DIR, "ALL_DATASETS_K_SELECTION.csv"), index=False
    )

    print("\n" + "=" * 70)
    print("K值选择汇总表:")
    print("=" * 70)
    print(f"{'Dataset':<25} {'K_aggr':>7} {'K_main':>7} {'K_cons':>7}  "
          f"{'Dist_aggr':>10} {'Dist_main':>10} {'Dist_cons':>10}")
    print("-" * 95)
    for _, row in k_summary_df.iterrows():
        print(f"{row['dataset']:<25} {row['k_aggr_2nd_diff']:>7} "
              f"{row['k_main_piecewise']:>7} {row['k_cons_max_dist']:>7}  "
              f"{row['dist_aggr_um']:>9.1f}µm {row['dist_main_um']:>9.1f}µm "
              f"{row['dist_cons_um']:>9.1f}µm")

    # ==========================================
    # 阶段 2：逐样本计算5种密度
    # ==========================================
    print("\n" + "=" * 70)
    print("阶段 2：逐样本密度计算（5种方法）")
    print("=" * 70)

    all_qc_rows = []
    method_names = [
        "density_knn_aggr_2nd_diff",
        "density_knn_main_piecewise",
        "density_knn_cons_max_dist",
        "density_voronoi",
        "density_delaunay",
    ]

    sample_count = 0
    total_samples = len(registry)

    for sample_name, sample_info in registry.items():
        sample_count += 1
        dataset = sample_name.split("/")[0]
        sub_name = sample_name.split("/")[1]
        k_info = dataset_k_decisions[dataset]

        k_aggr = k_info["k_aggr_2nd_diff"]
        k_main = k_info["k_main_piecewise"]
        k_cons = k_info["k_cons_max_dist"]
        k_max = max(k_aggr, k_main, k_cons)

        print(f"\n[{sample_count}/{total_samples}] {sample_name}")
        print(f"  K values: aggr={k_aggr}, main={k_main}, cons={k_cons}")

        # 创建输出目录
        sample_out_dir = os.path.join(P2_RESULTS_DIR, dataset, sub_name)
        os.makedirs(sample_out_dir, exist_ok=True)

        # 读取坐标
        meta_path = os.path.join(P1_RESULTS_DIR, dataset, sub_name, "cell_metadata.csv")
        df = pd.read_csv(meta_path)
        coords = df[["x_centroid", "y_centroid"]].values
        n_cells = len(coords)
        t0 = time.time()

        # ---- KNN密度（三个K值，一次BallTree调用） ----
        nn = NearestNeighbors(n_neighbors=k_max + 1, algorithm="ball_tree")
        nn.fit(coords)
        distances, _ = nn.kneighbors(coords)

        knn_densities = {}
        for k_val, col_name in [(k_aggr, "density_knn_aggr_2nd_diff"),
                                 (k_main, "density_knn_main_piecewise"),
                                 (k_cons, "density_knn_cons_max_dist")]:
            dist_k = distances[:, k_val].copy()
            dist_k[dist_k == 0] = np.finfo(float).eps
            knn_densities[col_name] = 1.0 / dist_k

        t_knn = time.time() - t0
        print(f"  KNN完成: {t_knn:.1f}s")

        # ---- Voronoi密度 ----
        t0 = time.time()
        density_voronoi = compute_voronoi_density(coords)
        t_vor = time.time() - t0
        n_valid_vor = np.sum(~np.isnan(density_voronoi))
        print(f"  Voronoi完成: {t_vor:.1f}s, 有效值: {n_valid_vor}/{n_cells} "
              f"({100 * n_valid_vor / n_cells:.1f}%)")

        # ---- Delaunay密度 ----
        t0 = time.time()
        density_delaunay = compute_delaunay_density(coords)
        t_del = time.time() - t0
        n_valid_del = np.sum(~np.isnan(density_delaunay))
        print(f"  Delaunay完成: {t_del:.1f}s, 有效值: {n_valid_del}/{n_cells} "
              f"({100 * n_valid_del / n_cells:.1f}%)")

        # ---- 构建输出DataFrame ----
        out_df = pd.DataFrame({
            "cell_id": df["cell_id"].values,
            "x_centroid": df["x_centroid"].values,
            "y_centroid": df["y_centroid"].values,
            "density_knn_aggr_2nd_diff": knn_densities["density_knn_aggr_2nd_diff"],
            "density_knn_main_piecewise": knn_densities["density_knn_main_piecewise"],
            "density_knn_cons_max_dist": knn_densities["density_knn_cons_max_dist"],
            "density_voronoi": density_voronoi,
            "density_delaunay": density_delaunay,
        })

        # 保存密度CSV
        out_df.to_csv(os.path.join(sample_out_dir, "cell_density.csv"), index=False)

        # ---- 方法间相关性（5x5 Spearman矩阵） ----
        all_densities = [
            knn_densities["density_knn_aggr_2nd_diff"],
            knn_densities["density_knn_main_piecewise"],
            knn_densities["density_knn_cons_max_dist"],
            density_voronoi,
            density_delaunay,
        ]

        n_methods = len(method_names)
        corr_matrix = np.eye(n_methods)
        for i in range(n_methods):
            for j in range(i + 1, n_methods):
                # 只用两种方法都有有效值的细胞
                valid = ~np.isnan(all_densities[i]) & ~np.isnan(all_densities[j])
                if valid.sum() > 10:
                    rho, _ = spearmanr(all_densities[i][valid], all_densities[j][valid])
                    corr_matrix[i, j] = rho
                    corr_matrix[j, i] = rho
                else:
                    corr_matrix[i, j] = np.nan
                    corr_matrix[j, i] = np.nan

        # 打印关键相关系数
        # KNN_main vs Voronoi, KNN_main vs Delaunay
        rho_mv = corr_matrix[1, 3]  # main vs voronoi
        rho_md = corr_matrix[1, 4]  # main vs delaunay
        rho_aa = corr_matrix[0, 1]  # aggr vs main (KNN内部一致性)
        rho_mc = corr_matrix[1, 2]  # main vs cons (KNN内部一致性)
        print(f"  相关性: KNN_main↔Vor={rho_mv:.3f}, KNN_main↔Del={rho_md:.3f}, "
              f"KNN_aggr↔main={rho_aa:.3f}, KNN_main↔cons={rho_mc:.3f}")

        # ---- 诊断图 ----
        plot_density_diagnostics(
            sample_name, coords, all_densities, method_names,
            corr_matrix,
            os.path.join(sample_out_dir, "density_diagnostics.png")
        )

        # ---- 样本级QC ----
        qc_row = {
            "sample": sample_name,
            "dataset": dataset,
            "species": sample_info.get("species", ""),
            "condition": sample_info.get("condition", ""),
            "n_cells": n_cells,
            "k_aggr": k_aggr,
            "k_main": k_main,
            "k_cons": k_cons,
            "median_density_knn_main": float(np.median(knn_densities["density_knn_main_piecewise"])),
            "n_valid_voronoi": int(n_valid_vor),
            "n_valid_delaunay": int(n_valid_del),
            "pct_valid_voronoi": round(100 * n_valid_vor / n_cells, 1),
            "pct_valid_delaunay": round(100 * n_valid_del / n_cells, 1),
            "corr_knn_main_voronoi": round(rho_mv, 4),
            "corr_knn_main_delaunay": round(rho_md, 4),
            "corr_knn_aggr_main": round(rho_aa, 4),
            "corr_knn_main_cons": round(rho_mc, 4),
            "time_knn_s": round(t_knn, 1),
            "time_voronoi_s": round(t_vor, 1),
            "time_delaunay_s": round(t_del, 1),
        }
        qc_row_df = pd.DataFrame([qc_row])
        qc_row_df.to_csv(os.path.join(sample_out_dir, "density_qc.csv"), index=False)
        all_qc_rows.append(qc_row)

    # ==========================================
    # 阶段 3：汇总QC
    # ==========================================
    print("\n" + "=" * 70)
    print("阶段 3：汇总QC报告")
    print("=" * 70)

    qc_df = pd.DataFrame(all_qc_rows)
    qc_df.to_csv(os.path.join(P2_RESULTS_DIR, "ALL_SAMPLES_P2_QC.csv"), index=False)

    # 打印汇总表
    print(f"\n{'Sample':<40} {'Cells':>8} {'K_a':>4} {'K_m':>4} {'K_c':>4}  "
          f"{'Vor%':>5} {'Del%':>5}  "
          f"{'ρ_m↔V':>6} {'ρ_m↔D':>6} {'ρ_a↔m':>6} {'ρ_m↔c':>6}")
    print("-" * 120)
    for _, r in qc_df.iterrows():
        print(f"{r['sample']:<40} {r['n_cells']:>8} "
              f"{r['k_aggr']:>4} {r['k_main']:>4} {r['k_cons']:>4}  "
              f"{r['pct_valid_voronoi']:>5.1f} {r['pct_valid_delaunay']:>5.1f}  "
              f"{r['corr_knn_main_voronoi']:>6.3f} {r['corr_knn_main_delaunay']:>6.3f} "
              f"{r['corr_knn_aggr_main']:>6.3f} {r['corr_knn_main_cons']:>6.3f}")

    # 打印总体统计
    print(f"\n总计: {len(qc_df)} 个样本, {qc_df['n_cells'].sum():,} 个细胞")
    print(f"KNN_main ↔ Voronoi  ρ: "
          f"mean={qc_df['corr_knn_main_voronoi'].mean():.3f}, "
          f"range=[{qc_df['corr_knn_main_voronoi'].min():.3f}, "
          f"{qc_df['corr_knn_main_voronoi'].max():.3f}]")
    print(f"KNN_main ↔ Delaunay ρ: "
          f"mean={qc_df['corr_knn_main_delaunay'].mean():.3f}, "
          f"range=[{qc_df['corr_knn_main_delaunay'].min():.3f}, "
          f"{qc_df['corr_knn_main_delaunay'].max():.3f}]")
    print(f"KNN_aggr ↔ KNN_main ρ: "
          f"mean={qc_df['corr_knn_aggr_main'].mean():.3f}, "
          f"range=[{qc_df['corr_knn_aggr_main'].min():.3f}, "
          f"{qc_df['corr_knn_aggr_main'].max():.3f}]")
    print(f"KNN_main ↔ KNN_cons ρ: "
          f"mean={qc_df['corr_knn_main_cons'].mean():.3f}, "
          f"range=[{qc_df['corr_knn_main_cons'].min():.3f}, "
          f"{qc_df['corr_knn_main_cons'].max():.3f}]")

    elapsed = time.time() - start_time
    print(f"\n{'=' * 70}")
    print(f"P2 全部完成! 总耗时: {elapsed / 60:.1f} 分钟")
    print(f"输出目录: {P2_RESULTS_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
