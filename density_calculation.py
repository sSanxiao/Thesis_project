#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
DENSITY CALCULATION PIPELINE  v2
三方法密度计算: KNN (自动最优K) + Voronoi + Delaunay
============================================================

已确认 (基于 inspect_result):
  - 坐标列名: x_centroid, y_centroid (全部 21 样本一致)
  - 坐标单位: 微米 (μm), 无缺失值
  - ATRT 7 样本: cells.csv, cell_id 为 string
  - 其余 14 样本: cells.parquet, 部分 cell_id 为 bytes 需 decode
  - 细胞数范围: 13,178 ~ 454,004

运行方式:
    cd /home/disk/wangqilu/Density_Caculation/Scripts/
    nohup python3 density_calculation.py > density_run.log 2>&1 &
    tail -f density_run.log
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, Delaunay
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =============================================================
# 1. 全局配置
# =============================================================

RESULT_DIR = "/home/disk/wangqilu/Density_Caculation/Results/"

# KNN 候选 K 值
K_CANDIDATES = list(range(5, 55, 5))  # [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# 脑组织局部微环境合理物理半径 (μm)
BIO_RANGE_MIN = 30
BIO_RANGE_MAX = 250

# Voronoi / Delaunay 边界裁剪百分位
BOUNDARY_PERCENTILE = 99

# =============================================================
# 2. 样本注册表 (file_type: "parquet" 或 "csv")
# =============================================================

SAMPLES = {
    # ---- Human (14 samples) ----
    # ATRT — cells.csv
    "ATRT_Human/28": {"path": "/home/disk/wangqilu/ATRT_Human/28/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/29": {"path": "/home/disk/wangqilu/ATRT_Human/29/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/30": {"path": "/home/disk/wangqilu/ATRT_Human/30/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/31": {"path": "/home/disk/wangqilu/ATRT_Human/31/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/32": {"path": "/home/disk/wangqilu/ATRT_Human/32/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/33": {"path": "/home/disk/wangqilu/ATRT_Human/33/data/", "species": "human", "file_type": "csv"},
    "ATRT_Human/34": {"path": "/home/disk/wangqilu/ATRT_Human/34/data/", "species": "human", "file_type": "csv"},
    # Medulloblastoma — cells.parquet
    "Medulloblastoma_Human/MB266": {"path": "/home/disk/wangqilu/Medulloblastoma_Human/GSM8840047_MB266/data/", "species": "human", "file_type": "parquet"},
    "Medulloblastoma_Human/MB299": {"path": "/home/disk/wangqilu/Medulloblastoma_Human/GSM8840049_MB299/data/", "species": "human", "file_type": "parquet"},
    "Medulloblastoma_Human/MB295": {"path": "/home/disk/wangqilu/Medulloblastoma_Human/GSM8840048_MB295/data/", "species": "human", "file_type": "parquet"},
    "Medulloblastoma_Human/MB263": {"path": "/home/disk/wangqilu/Medulloblastoma_Human/GSM8840046_MB263/data/", "species": "human", "file_type": "parquet"},
    # Brain_Human — cells.parquet
    "Brain_Human/Healthy": {"path": "/home/disk/wangqilu/Brain_Human/Healthy/data/", "species": "human", "file_type": "parquet"},
    "Brain_Human/Alz":     {"path": "/home/disk/wangqilu/Brain_Human/Alz/data/",     "species": "human", "file_type": "parquet"},
    "Brain_Human/Glio":    {"path": "/home/disk/wangqilu/Brain_Human/Gilo/data/",     "species": "human", "file_type": "parquet"},

    # ---- Mouse (7 samples) ----
    "Brain_Mouse/Normal": {"path": "/home/disk/wangqilu/Brain_Mouse/data/", "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/TgCRND8_17_9": {"path": "/home/disk/wangqilu/Alzheimer_Mouse/TgCRND8_17_9/data/", "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/TgCRND8_2_5":  {"path": "/home/disk/wangqilu/Alzheimer_Mouse/TgCRND8_2_5/data/",  "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/TgCRND8_5_7":  {"path": "/home/disk/wangqilu/Alzheimer_Mouse/TgCRND8_5_7/data/",  "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/Wild_2_5":     {"path": "/home/disk/wangqilu/Alzheimer_Mouse/Wild_2_5/data/",      "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/Wild_5_7":     {"path": "/home/disk/wangqilu/Alzheimer_Mouse/Wild_5_7/data/",      "species": "mouse", "file_type": "parquet"},
    "Alzheimer_Mouse/Wild_13_4":    {"path": "/home/disk/wangqilu/Alzheimer_Mouse/Wild_13_4/data/",     "species": "mouse", "file_type": "parquet"},
}


# =============================================================
# 3. 数据读取
# =============================================================

def load_coordinates(data_dir, file_type):
    """
    读取 cell_id, x_centroid, y_centroid
    自动适配 csv / parquet, 处理 bytes 类型 cell_id
    """
    if file_type == "csv":
        filepath = os.path.join(data_dir, "cells.csv")
        df = pd.read_csv(filepath, usecols=["cell_id", "x_centroid", "y_centroid"])
    else:
        filepath = os.path.join(data_dir, "cells.parquet")
        df = pd.read_parquet(filepath, columns=["cell_id", "x_centroid", "y_centroid"])

    # 修复 bytes 类型 cell_id (Brain_Human, Alzheimer_Mouse 部分样本)
    if len(df) > 0 and isinstance(df["cell_id"].iloc[0], bytes):
        df["cell_id"] = df["cell_id"].apply(
            lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x)
        )

    df = df.dropna(subset=["x_centroid", "y_centroid"]).reset_index(drop=True)
    return df


# =============================================================
# 4. 三种 Density 计算方法
# =============================================================

def compute_knn_density(coords, k):
    """KNN 密度: 第 k 近邻距离的倒数"""
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="ball_tree")
    nn.fit(coords)
    distances, _ = nn.kneighbors(coords)
    dist_k = distances[:, k].copy()
    dist_k[dist_k == 0] = np.finfo(float).eps
    return 1.0 / dist_k


def knn_scan(coords, k_candidates):
    """一次性算到 max(k), 返回距离统计 + CV + 原始距离矩阵"""
    k_max = max(k_candidates)
    nn = NearestNeighbors(n_neighbors=k_max + 1, algorithm="ball_tree")
    nn.fit(coords)
    distances, _ = nn.kneighbors(coords)

    dist_records = []
    cv_records = []
    for k in k_candidates:
        dk = distances[:, k]
        dk_safe = dk.copy()
        dk_safe[dk_safe == 0] = np.finfo(float).eps
        density = 1.0 / dk_safe

        dist_records.append({
            "k": k,
            "median_dist_um": np.median(dk),
            "mean_dist_um":   np.mean(dk),
            "q25_dist_um":    np.percentile(dk, 25),
            "q75_dist_um":    np.percentile(dk, 75),
        })
        cv_records.append({
            "k": k,
            "cv": np.std(density) / np.mean(density),
        })

    return pd.DataFrame(dist_records), pd.DataFrame(cv_records), distances


def compute_voronoi_density(coords, percentile_cap=99):
    """Voronoi 密度: 多边形面积的倒数, 边界裁剪为 NaN"""
    n = len(coords)
    vor = Voronoi(coords)

    areas = np.full(n, np.nan)
    for i in range(n):
        region_idx = vor.point_region[i]
        verts = vor.regions[region_idx]
        if -1 in verts or len(verts) == 0:
            continue
        polygon = vor.vertices[verts]
        xp, yp = polygon[:, 0], polygon[:, 1]
        areas[i] = 0.5 * np.abs(np.dot(xp, np.roll(yp, 1)) - np.dot(yp, np.roll(xp, 1)))

    valid = ~np.isnan(areas)
    if valid.sum() > 0:
        cap = np.percentile(areas[valid], percentile_cap)
        areas[valid & (areas > cap)] = np.nan

    density = np.full(n, np.nan)
    ok = ~np.isnan(areas) & (areas > 0)
    density[ok] = 1.0 / areas[ok]
    return density


def compute_delaunay_density(coords, percentile_cap=99):
    """Delaunay 密度: 平均自然邻居边长的倒数, 边界裁剪为 NaN"""
    n = len(coords)
    tri = Delaunay(coords)

    edge_lengths = defaultdict(list)
    for simplex in tri.simplices:
        for i_loc in range(3):
            for j_loc in range(i_loc + 1, 3):
                pi, pj = simplex[i_loc], simplex[j_loc]
                d = np.sqrt((coords[pi, 0] - coords[pj, 0])**2 +
                            (coords[pi, 1] - coords[pj, 1])**2)
                edge_lengths[pi].append(d)
                edge_lengths[pj].append(d)

    mean_edge = np.full(n, np.nan)
    for i in range(n):
        if edge_lengths[i]:
            mean_edge[i] = np.mean(edge_lengths[i])

    valid = ~np.isnan(mean_edge)
    if valid.sum() > 0:
        cap = np.percentile(mean_edge[valid], percentile_cap)
        mean_edge[valid & (mean_edge > cap)] = np.nan

    density = np.full(n, np.nan)
    ok = ~np.isnan(mean_edge) & (mean_edge > 0)
    density[ok] = 1.0 / mean_edge[ok]
    return density


# =============================================================
# 5. KNN 最优 K 值确定 (按物种)
# =============================================================

def determine_optimal_k(species_samples, species_name, result_dir):
    """
    阶段一 (物理尺度锚定) + 阶段二 (CV 跨样本稳定性)
    → 在生物学合理范围内选 CV 跨样本变异最小的 k
    """
    species_dir = os.path.join(result_dir, f"KNN_Optimization_{species_name}")
    os.makedirs(species_dir, exist_ok=True)

    all_dist = []
    all_cv = []

    for name, info in species_samples.items():
        print(f"    [K-scan] {name} ...", flush=True)
        df = load_coordinates(info["path"], info["file_type"])
        coords = df[["x_centroid", "y_centroid"]].values

        dist_df, cv_df, _ = knn_scan(coords, K_CANDIDATES)
        dist_df["sample"] = name
        cv_df["sample"] = name
        all_dist.append(dist_df)
        all_cv.append(cv_df)

    dist_all = pd.concat(all_dist, ignore_index=True)
    cv_all = pd.concat(all_cv, ignore_index=True)

    # 阶段一: 每个 k 跨样本的平均中位距离
    k_phys = dist_all.groupby("k")["median_dist_um"].mean().reset_index()
    k_phys.columns = ["k", "avg_median_dist"]
    k_phys["in_bio_range"] = (
        (k_phys["avg_median_dist"] >= BIO_RANGE_MIN) &
        (k_phys["avg_median_dist"] <= BIO_RANGE_MAX)
    )

    # 阶段二: CV 跨样本标准差
    cv_stab = cv_all.groupby("k")["cv"].agg(["mean", "std"]).reset_index()
    cv_stab.columns = ["k", "cv_mean", "cv_std"]

    decision = pd.merge(k_phys, cv_stab, on="k")
    decision.to_csv(os.path.join(species_dir, "k_decision_table.csv"), index=False)

    candidates = decision[decision["in_bio_range"]]
    if len(candidates) > 0:
        best_k = int(candidates.loc[candidates["cv_std"].idxmin(), "k"])
    else:
        mid = (BIO_RANGE_MIN + BIO_RANGE_MAX) / 2
        decision["gap"] = abs(decision["avg_median_dist"] - mid)
        best_k = int(decision.loc[decision["gap"].idxmin(), "k"])
        print(f"    [WARNING] No k in bio range for {species_name}, fallback k={best_k}")

    # ---- 诊断图 ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax1 = axes[0]
    for s in dist_all["sample"].unique():
        sub = dist_all[dist_all["sample"] == s]
        ax1.plot(sub["k"], sub["median_dist_um"], "o-", markersize=4, alpha=0.7,
                 label=s.split("/")[-1])
    ax1.axhspan(BIO_RANGE_MIN, BIO_RANGE_MAX, alpha=0.12, color="red",
                label=f"Bio range ({BIO_RANGE_MIN}–{BIO_RANGE_MAX} μm)")
    ax1.axvline(best_k, color="black", ls="--", lw=2, label=f"Selected k={best_k}")
    ax1.set_xlabel("k", fontsize=12)
    ax1.set_ylabel("Median k-th Neighbor Distance (μm)", fontsize=12)
    ax1.set_title(f"{species_name}: Physical Scale Anchoring", fontsize=13)
    ax1.legend(fontsize=7, loc="upper left")
    ax1.grid(True, alpha=0.3)

    ax2 = axes[1]
    for s in cv_all["sample"].unique():
        sub = cv_all[cv_all["sample"] == s]
        ax2.plot(sub["k"], sub["cv"], "s-", markersize=4, alpha=0.7,
                 label=s.split("/")[-1])
    ax2.axvline(best_k, color="black", ls="--", lw=2, label=f"Selected k={best_k}")
    ax2.set_xlabel("k", fontsize=12)
    ax2.set_ylabel("Coefficient of Variation (CV)", fontsize=12)
    ax2.set_title(f"{species_name}: Density CV Across Samples", fontsize=13)
    ax2.legend(fontsize=7, loc="upper right")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(species_dir, "k_optimization_diagnostic.png"), dpi=200)
    plt.close()

    dist_all.to_csv(os.path.join(species_dir, "all_samples_knn_distances.csv"), index=False)
    cv_all.to_csv(os.path.join(species_dir, "all_samples_knn_cv.csv"), index=False)

    avg_dist = k_phys.loc[k_phys["k"] == best_k, "avg_median_dist"].values[0]
    print(f"    >>> {species_name} optimal k = {best_k}  "
          f"(avg median dist ≈ {avg_dist:.1f} μm)", flush=True)

    return best_k


# =============================================================
# 6. 单样本: 三方法 Density + 诊断图
# =============================================================

def process_sample(sample_name, sample_info, k_opt, result_dir):
    """计算三种 density, 输出 CSV + 诊断图"""

    sample_out = os.path.join(result_dir, sample_name)
    os.makedirs(sample_out, exist_ok=True)

    df = load_coordinates(sample_info["path"], sample_info["file_type"])
    coords = df[["x_centroid", "y_centroid"]].values
    n = len(df)
    print(f"    Cells: {n:,}", flush=True)

    # 三种 density
    t0 = time.time()
    print(f"    KNN (k={k_opt}) ...", end="", flush=True)
    d_knn = compute_knn_density(coords, k_opt)
    print(f" {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    print(f"    Voronoi ...", end="", flush=True)
    d_vor = compute_voronoi_density(coords, BOUNDARY_PERCENTILE)
    print(f" {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    print(f"    Delaunay ...", end="", flush=True)
    d_del = compute_delaunay_density(coords, BOUNDARY_PERCENTILE)
    print(f" {time.time()-t0:.1f}s", flush=True)

    # 核心输出 CSV
    result = pd.DataFrame({
        "cell_id":          df["cell_id"].values,
        "x_centroid":       df["x_centroid"].values,
        "y_centroid":       df["y_centroid"].values,
        "density_knn":      d_knn,
        "density_voronoi":  d_vor,
        "density_delaunay": d_del,
    })
    result.to_csv(os.path.join(sample_out, "cell_density_three_methods.csv"), index=False)

    # 方法间 Spearman
    mask = ~np.isnan(d_knn) & ~np.isnan(d_vor) & ~np.isnan(d_del)
    rho_kv, _ = spearmanr(d_knn[mask], d_vor[mask])
    rho_kd, _ = spearmanr(d_knn[mask], d_del[mask])
    rho_vd, _ = spearmanr(d_vor[mask], d_del[mask])

    cor_df = pd.DataFrame({
        "pair":         ["KNN_vs_Voronoi", "KNN_vs_Delaunay", "Voronoi_vs_Delaunay"],
        "spearman_rho": [rho_kv, rho_kd, rho_vd],
    })
    cor_df.to_csv(os.path.join(sample_out, "density_method_correlation.csv"), index=False)

    # QC
    qc = {
        "sample":  sample_name,
        "species": sample_info["species"],
        "n_cells": n,
        "k_used":  k_opt,
        "knn_valid_pct":      round((~np.isnan(d_knn)).mean() * 100, 1),
        "voronoi_valid_pct":  round((~np.isnan(d_vor)).mean() * 100, 1),
        "delaunay_valid_pct": round((~np.isnan(d_del)).mean() * 100, 1),
        "knn_median":      np.nanmedian(d_knn),
        "voronoi_median":  np.nanmedian(d_vor),
        "delaunay_median": np.nanmedian(d_del),
        "rho_knn_voronoi":  rho_kv,
        "rho_knn_delaunay": rho_kd,
        "rho_vor_delaunay": rho_vd,
    }
    pd.DataFrame([qc]).to_csv(os.path.join(sample_out, "density_qc.csv"), index=False)

    # 诊断图 1: 空间分布
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    triplets = [
        (f"KNN (k={k_opt})",   d_knn, None),
        (f"Voronoi",           d_vor, rho_kv),
        (f"Delaunay",          d_del, rho_kd),
    ]
    for ax, (title, dens, rho) in zip(axes, triplets):
        v = ~np.isnan(dens)
        if v.sum() == 0:
            ax.set_title(f"{title}: no valid values")
            continue
        vmin, vmax = np.percentile(dens[v], [1, 99])
        sc = ax.scatter(coords[v, 0], coords[v, 1], c=dens[v], s=0.1,
                        alpha=0.5, cmap="inferno", vmin=vmin, vmax=vmax, rasterized=True)
        ax.set_aspect("equal")
        sub = f"\nρ(vs KNN)={rho:.3f}" if rho is not None else ""
        ax.set_title(f"{title}{sub}", fontsize=11)
        ax.set_xlabel("X (μm)")
        ax.set_ylabel("Y (μm)")
        plt.colorbar(sc, ax=ax, shrink=0.8, label="Density")

    plt.suptitle(f"{sample_name}  |  n={n:,} cells", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(sample_out, "density_spatial_comparison.png"), dpi=150)
    plt.close()

    # 诊断图 2: 直方图
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = ["#2c7bb6", "#d7191c", "#1a9850"]
    for ax, (title, dens, _), clr in zip(axes, triplets, colors):
        vd = dens[~np.isnan(dens)]
        if len(vd) == 0:
            continue
        clip = np.percentile(vd, 99)
        ax.hist(vd[vd <= clip], bins=80, color=clr, alpha=0.85, edgecolor="white", linewidth=0.3)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Density")
        ax.set_ylabel("Count")

    plt.suptitle(f"{sample_name}  |  Density Distributions", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(sample_out, "density_histograms.png"), dpi=150)
    plt.close()

    return qc


# =============================================================
# 7. 主流程
# =============================================================

def main():
    t_start = time.time()

    print("=" * 60)
    print("  DENSITY CALCULATION PIPELINE  v2")
    print("  Methods: KNN + Voronoi + Delaunay")
    print("  Samples: 21  |  Species: Human (14) + Mouse (7)")
    print("=" * 60)
    print(flush=True)

    os.makedirs(RESULT_DIR, exist_ok=True)

    human = {k: v for k, v in SAMPLES.items() if v["species"] == "human"}
    mouse = {k: v for k, v in SAMPLES.items() if v["species"] == "mouse"}

    # ---- Step 1: 确定最优 K ----
    print("=" * 60)
    print("  STEP 1: KNN OPTIMAL K DETERMINATION")
    print("=" * 60, flush=True)

    print(f"\n  [Human] {len(human)} samples", flush=True)
    k_human = determine_optimal_k(human, "Human", RESULT_DIR)

    print(f"\n  [Mouse] {len(mouse)} samples", flush=True)
    k_mouse = determine_optimal_k(mouse, "Mouse", RESULT_DIR)

    k_report = pd.DataFrame([
        {"species": "human", "optimal_k": k_human},
        {"species": "mouse", "optimal_k": k_mouse},
    ])
    k_report.to_csv(os.path.join(RESULT_DIR, "OPTIMAL_K_VALUES.csv"), index=False)

    print(f"\n  >>> Human k = {k_human}")
    print(f"  >>> Mouse k = {k_mouse}")
    print(flush=True)

    # ---- Step 2: 逐样本三方法 density ----
    print("=" * 60)
    print("  STEP 2: THREE-METHOD DENSITY CALCULATION")
    print("=" * 60, flush=True)

    all_qc = []
    for i, (name, info) in enumerate(SAMPLES.items(), 1):
        k_opt = k_human if info["species"] == "human" else k_mouse
        print(f"\n  [{i}/{len(SAMPLES)}] {name}  ({info['species']}, k={k_opt})", flush=True)

        try:
            qc = process_sample(name, info, k_opt, RESULT_DIR)
            all_qc.append(qc)
            print(f"    OK.", flush=True)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)
            all_qc.append({"sample": name, "species": info["species"], "error": str(e)})

    # ---- 汇总 ----
    qc_df = pd.DataFrame(all_qc)
    qc_df.to_csv(os.path.join(RESULT_DIR, "SUMMARY_all_samples_qc.csv"), index=False)

    # 汇总图: 方法间相关性
    qc_valid = qc_df.dropna(subset=["rho_knn_voronoi"])
    if len(qc_valid) > 0:
        fig, ax = plt.subplots(figsize=(12, 6))
        x_pos = np.arange(len(qc_valid))
        w = 0.25
        ax.bar(x_pos - w, qc_valid["rho_knn_voronoi"],  w, label="KNN vs Voronoi",  color="#2c7bb6")
        ax.bar(x_pos,     qc_valid["rho_knn_delaunay"],  w, label="KNN vs Delaunay", color="#d7191c")
        ax.bar(x_pos + w, qc_valid["rho_vor_delaunay"],  w, label="Vor vs Delaunay", color="#1a9850")
        ax.set_xticks(x_pos)
        labels = [s.split("/")[-1] for s in qc_valid["sample"]]
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Spearman ρ")
        ax.set_title("Cross-Method Density Correlation per Sample")
        ax.legend()
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULT_DIR, "SUMMARY_method_correlation.png"), dpi=200)
        plt.close()

    elapsed = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"  ALL DONE  ({elapsed/60:.1f} min)")
    print(f"  Results: {RESULT_DIR}")
    print(f"  Human k = {k_human},  Mouse k = {k_mouse}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
