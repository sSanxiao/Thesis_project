#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
P1_data_loading.py

Load each sample's official cell_feature_matrix.h5, drop control
features and empty cells, and emit a standard 10x-format h5 plus a
cell-metadata csv for the R pipeline.

Input : sample_registry.json + each sample's h5 and cells.parquet
Output: <RESULTS_DIR>/P1_Results/<sample>/
          - filtered_matrix.h5    (standard 10x format, read by Seurat::Read10X_h5)
          - cell_metadata.csv     (cell_id + coords + QC columns, read by R fread)
        + ALL_SAMPLES_P1_QC.csv   (22-row summary)
Run   : python3 P1_data_loading.py

功能: 从官方 cell_feature_matrix.h5 加载数据，过滤对照feature
      和空细胞，输出标准10x格式h5 + 元数据csv
============================================================
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix

# =============================================================
# 配置
# =============================================================

# Paths are resolved from this script's location and from RESULTS_DIR
# (see config/paths.py); no hardcoded home directories.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REGISTRY_PATH = os.path.join(SCRIPT_DIR, "sample_registry.json")
RESULT_DIR = os.path.join(os.environ.get("RESULTS_DIR", "./results"), "P1_Results")
MIN_TRANSCRIPTS = 10

# 对照 feature 的名称前缀（不区分大小写匹配）
CONTROL_PREFIXES = [
    "NegControlProbe",
    "NegControlCodeword",
    "BLANK",
    "DeprecatedCodeword",
    "UnassignedCodeword",
    "antisense",
]

# =============================================================
# 核心函数
# =============================================================

def load_h5(h5_path):
    """
    读取 cell_feature_matrix.h5，自动适配 matrix/ 或 unknown/ 路径。
    返回: gene_names(list[str]), barcodes(list[str]),
          sparse_matrix(scipy csc_matrix, shape=(n_genes, n_cells)),
          path_type(str)
    """
    with h5py.File(h5_path, "r") as f:
        top_keys = list(f.keys())

        if "matrix" in top_keys:
            path_type = "matrix"
            genes_path = "matrix/features/name"
            barcodes_path = "matrix/barcodes"
            data_path = "matrix/data"
            indices_path = "matrix/indices"
            indptr_path = "matrix/indptr"
            shape_path = "matrix/shape"
        elif "unknown" in top_keys:
            path_type = "unknown"
            genes_path = "unknown/gene_names"
            barcodes_path = "unknown/barcodes"
            data_path = "unknown/data"
            indices_path = "unknown/indices"
            indptr_path = "unknown/indptr"
            shape_path = "unknown/shape"
        else:
            raise ValueError(f"无法识别的 h5 结构: 顶层 keys = {top_keys}")

        # 读取数据
        raw_genes = f[genes_path][:]
        raw_barcodes = f[barcodes_path][:]
        data = f[data_path][:]
        indices = f[indices_path][:]
        indptr = f[indptr_path][:]
        shape = tuple(f[shape_path][:])

    # decode bytes → str
    gene_names = [g.decode("utf-8") if isinstance(g, bytes) else str(g)
                  for g in raw_genes]
    barcodes = [b.decode("utf-8") if isinstance(b, bytes) else str(b)
                for b in raw_barcodes]

    # 重建稀疏矩阵
    mat = csc_matrix((data, indices, indptr), shape=shape)

    return gene_names, barcodes, mat, path_type


def identify_controls(gene_names):
    """
    识别对照 feature。返回两个列表:
    - gene_mask: bool数组，True=真基因，False=对照
    - control_names: 被识别为对照的 feature 名称列表
    """
    gene_mask = np.ones(len(gene_names), dtype=bool)
    control_names = []

    for i, name in enumerate(gene_names):
        name_lower = name.lower()
        for prefix in CONTROL_PREFIXES:
            if name_lower.startswith(prefix.lower()):
                gene_mask[i] = False
                control_names.append(name)
                break

    return gene_mask, control_names


def filter_matrix_rows(mat, mask):
    """从 CSC 稀疏矩阵中删除指定行（基因维度）"""
    # CSC 按列存储，删行需要转为 CSR 再切片再转回
    return mat.tocsr()[mask, :].tocsc()


def filter_matrix_cols(mat, mask):
    """从 CSC 稀疏矩阵中删除指定列（细胞维度）"""
    return mat[:, mask]


def load_cells(parquet_path):
    """
    读取 cells.parquet，提取需要的列。
    缺失的列会被跳过并警告。
    """
    df = pd.read_parquet(parquet_path)

    # 必需列
    required = ["cell_id", "x_centroid", "y_centroid"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"cells.parquet 缺少必需列: {col}")

    # 可选列
    optional = ["transcript_counts", "cell_area", "nucleus_area",
                "control_probe_counts", "control_codeword_counts"]
    cols_to_keep = required.copy()
    missing_optional = []
    for col in optional:
        if col in df.columns:
            cols_to_keep.append(col)
        else:
            missing_optional.append(col)

    df = df[cols_to_keep].copy()

    # cell_id 应该已经是 str（我们在上传前统一过），但防御性检查
    if len(df) > 0 and isinstance(df["cell_id"].iloc[0], bytes):
        df["cell_id"] = df["cell_id"].apply(
            lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x)
        )

    return df, missing_optional


def align_and_filter(gene_names, barcodes, mat, cells_df, min_transcripts):
    """
    对齐 h5 barcodes 和 cells cell_id，然后过滤空细胞。
    返回过滤后的 barcodes, mat, cells_df, 以及过滤统计。
    """
    h5_set = set(barcodes)
    cells_set = set(cells_df["cell_id"])

    common = h5_set & cells_set
    only_h5 = len(h5_set - cells_set)
    only_cells = len(cells_set - h5_set)

    if only_h5 > 0 or only_cells > 0:
        print(f"      ⚠ cell_id 不完全一致: "
              f"共有={len(common)}, h5独有={only_h5}, cells独有={only_cells}")

    # 按 h5 barcodes 的顺序对齐
    barcode_to_idx = {b: i for i, b in enumerate(barcodes)}
    common_ordered = [b for b in barcodes if b in common]
    col_indices = [barcode_to_idx[b] for b in common_ordered]

    # 切片矩阵列
    mat_aligned = mat[:, col_indices]

    # 对齐 cells_df
    cells_df = cells_df.set_index("cell_id").loc[common_ordered].reset_index()

    n_before = len(cells_df)

    # 过滤 transcript_counts < min_transcripts
    if "transcript_counts" in cells_df.columns:
        keep_mask = cells_df["transcript_counts"].values >= min_transcripts
        cells_df = cells_df[keep_mask].reset_index(drop=True)
        mat_aligned = filter_matrix_cols(mat_aligned, keep_mask)
        barcodes_filtered = [common_ordered[i] for i, k in enumerate(keep_mask) if k]
    else:
        print(f"      ⚠ transcript_counts 列不存在，跳过空细胞过滤")
        barcodes_filtered = common_ordered

    n_after = len(cells_df)
    n_removed = n_before - n_after

    return barcodes_filtered, mat_aligned, cells_df, n_before, n_after, n_removed


def write_10x_h5(h5_path, gene_names, barcodes, mat):
    """
    输出标准 10x 格式 h5 文件。
    Seurat::Read10X_h5() 可直接读取。
    """
    # 确保是 CSC 格式
    if not isinstance(mat, csc_matrix):
        mat = mat.tocsc()

    # 编码字符串为 bytes（h5 标准）
    gene_bytes = np.array(gene_names, dtype="S")
    barcode_bytes = np.array(barcodes, dtype="S")
    feature_type = np.array(["Gene Expression"] * len(gene_names), dtype="S")

    with h5py.File(h5_path, "w") as f:
        g = f.create_group("matrix")

        # 稀疏矩阵数据
        g.create_dataset("data", data=mat.data.astype(np.int32))
        g.create_dataset("indices", data=mat.indices.astype(np.int64))
        g.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(mat.shape, dtype=np.int32))

        # barcodes
        g.create_dataset("barcodes", data=barcode_bytes)

        # features
        feat = g.create_group("features")
        feat.create_dataset("name", data=gene_bytes)
        feat.create_dataset("id", data=gene_bytes)
        feat.create_dataset("feature_type", data=feature_type)
        feat.create_dataset("genome", data=np.array(["unknown"] * len(gene_names), dtype="S"))

    return os.path.getsize(h5_path)


# =============================================================
# 主流程
# =============================================================

def main():
    t_global = time.time()

    print("=" * 65)
    print("  P1_data_loading.py — 数据加载、过滤、标准化输出")
    print("=" * 65)
    print(f"  Registry:  {REGISTRY_PATH}")
    print(f"  Output:    {RESULT_DIR}")
    print(f"  Min transcripts: {MIN_TRANSCRIPTS}")
    print()

    # 读取 registry
    with open(REGISTRY_PATH, "r") as f:
        registry = json.load(f)
    print(f"  共 {len(registry)} 个样本\n")

    os.makedirs(RESULT_DIR, exist_ok=True)
    all_qc = []

    # ===========================================================
    # 阶段 0: 探查所有 h5 中的对照 feature（首次运行时查看）
    # ===========================================================
    print("=" * 65)
    print("  阶段 0: 探查对照 feature")
    print("=" * 65)

    for sample_name, info in registry.items():
        h5_path = os.path.join(info["path"], "cell_feature_matrix.h5")
        gene_names, _, _, path_type = load_h5(h5_path)
        gene_mask, control_names = identify_controls(gene_names)
        n_ctrl = len(control_names)
        n_gene = gene_mask.sum()

        if n_ctrl > 0:
            # 统计各前缀的数量
            prefix_counts = {}
            for cn in control_names:
                for pfx in CONTROL_PREFIXES:
                    if cn.lower().startswith(pfx.lower()):
                        prefix_counts[pfx] = prefix_counts.get(pfx, 0) + 1
                        break
            prefix_str = ", ".join(f"{k}={v}" for k, v in sorted(prefix_counts.items()))
            print(f"  {sample_name}: {len(gene_names)} → {n_gene} genes "
                  f"+ {n_ctrl} controls ({prefix_str})")
        else:
            print(f"  {sample_name}: {len(gene_names)} → {n_gene} genes "
                  f"+ 0 controls (已被数据源预过滤)")

        # 检查是否有未被识别的可疑 feature（不像基因名的条目）
        for gn in gene_names:
            if gene_mask[gene_names.index(gn)]:
                # 已标记为基因，检查是否有可疑前缀
                gn_lower = gn.lower()
                if any(gn_lower.startswith(s) for s in
                       ["control", "neg", "blank", "deprecated", "unassigned"]):
                    print(f"    ⚠ 可疑基因名（未被过滤）: {gn}")

    print()

    # ===========================================================
    # 阶段 1: 逐样本处理
    # ===========================================================
    print("=" * 65)
    print("  阶段 1: 逐样本加载、过滤、输出")
    print("=" * 65)

    for i, (sample_name, info) in enumerate(registry.items(), 1):
        t0 = time.time()
        print(f"\n  [{i}/{len(registry)}] {sample_name}")
        print(f"  {'─' * 55}")

        sample_out = os.path.join(RESULT_DIR, sample_name)
        os.makedirs(sample_out, exist_ok=True)

        # ---- 步骤 1: 读 h5 ----
        h5_path = os.path.join(info["path"], "cell_feature_matrix.h5")
        print(f"    读取 h5: {h5_path}")
        gene_names, barcodes, mat, path_type = load_h5(h5_path)
        n_features_raw = len(gene_names)
        n_cells_raw = len(barcodes)
        print(f"    h5 路径类型: {path_type}")
        print(f"    原始维度: {n_features_raw} features × {n_cells_raw} cells")

        # ---- 步骤 2: 过滤对照 ----
        gene_mask, control_names = identify_controls(gene_names)
        n_controls = len(control_names)
        gene_names_filtered = [g for g, m in zip(gene_names, gene_mask) if m]
        mat = filter_matrix_rows(mat, gene_mask)
        n_genes_final = len(gene_names_filtered)
        print(f"    过滤对照: {n_features_raw} → {n_genes_final} genes "
              f"(删除 {n_controls} controls)")

        # ---- 步骤 3: 读 cells.parquet ----
        cells_path = os.path.join(info["path"], "cells.parquet")
        print(f"    读取 cells: {cells_path}")
        cells_df, missing_cols = load_cells(cells_path)
        if missing_cols:
            print(f"      可选列缺失（已跳过）: {missing_cols}")

        # ---- 步骤 4: 对齐 + 过滤空细胞 ----
        barcodes_final, mat_final, cells_final, n_before, n_after, n_removed = \
            align_and_filter(gene_names_filtered, barcodes, mat, cells_df,
                             MIN_TRANSCRIPTS)
        print(f"    空细胞过滤 (transcript_counts < {MIN_TRANSCRIPTS}): "
              f"{n_before} → {n_after} cells (删除 {n_removed})")

        # ---- 步骤 5: 输出 filtered_matrix.h5 ----
        h5_out = os.path.join(sample_out, "filtered_matrix.h5")
        h5_size = write_10x_h5(h5_out, gene_names_filtered, barcodes_final,
                               mat_final)
        print(f"    输出 h5: {h5_out} ({h5_size / 1024 / 1024:.2f} MB)")

        # ---- 步骤 6: 输出 cell_metadata.csv ----
        csv_out = os.path.join(sample_out, "cell_metadata.csv")
        cells_final.to_csv(csv_out, index=False)
        print(f"    输出 csv: {csv_out} ({os.path.getsize(csv_out) / 1024 / 1024:.2f} MB)")

        # ---- 步骤 7: QC 统计 ----
        nnz = mat_final.nnz
        total_elements = mat_final.shape[0] * mat_final.shape[1]
        nonzero_frac = nnz / total_elements if total_elements > 0 else 0

        qc = {
            "sample_name": sample_name,
            "species": info["species"],
            "condition": info["condition"],
            "preservation": info["preservation"],
            "data_quality_tier": info["data_quality_tier"],
            "h5_path_type": path_type,
            "n_features_raw": n_features_raw,
            "n_controls_removed": n_controls,
            "n_genes_final": n_genes_final,
            "n_cells_raw": n_cells_raw,
            "n_cells_removed": n_removed,
            "n_cells_final": n_after,
            "nonzero_elements": nnz,
            "nonzero_fraction": round(nonzero_frac, 4),
            "sparsity_pct": round((1 - nonzero_frac) * 100, 2),
        }

        # 如果有 transcript_counts 列，加入中位数统计
        if "transcript_counts" in cells_final.columns:
            qc["median_transcripts"] = int(cells_final["transcript_counts"].median())
            qc["mean_transcripts"] = round(cells_final["transcript_counts"].mean(), 1)
        if "cell_area" in cells_final.columns:
            qc["median_cell_area"] = round(cells_final["cell_area"].median(), 2)
        if "nucleus_area" in cells_final.columns:
            qc["median_nucleus_area"] = round(cells_final["nucleus_area"].median(), 2)

        all_qc.append(qc)

        elapsed = time.time() - t0
        print(f"    最终维度: {n_genes_final} genes × {n_after} cells, "
              f"稀疏度 {qc['sparsity_pct']}%, 耗时 {elapsed:.1f}s")

    # ===========================================================
    # 阶段 2: 汇总 QC
    # ===========================================================
    print(f"\n{'=' * 65}")
    print("  阶段 2: 汇总 QC")
    print(f"{'=' * 65}")

    qc_df = pd.DataFrame(all_qc)
    qc_path = os.path.join(RESULT_DIR, "ALL_SAMPLES_P1_QC.csv")
    qc_df.to_csv(qc_path, index=False)
    print(f"\n  QC 汇总已保存: {qc_path}")

    # 打印汇总表
    print(f"\n  {'样本':<40} {'基因':>6} {'细胞(原)':>10} {'细胞(后)':>10} "
          f"{'删除':>6} {'稀疏度':>7}")
    print("  " + "─" * 85)
    for _, row in qc_df.iterrows():
        print(f"  {row['sample_name']:<40} {row['n_genes_final']:>6} "
              f"{row['n_cells_raw']:>10,} {row['n_cells_final']:>10,} "
              f"{row['n_cells_removed']:>6} {row['sparsity_pct']:>6.1f}%")

    total_cells = qc_df["n_cells_final"].sum()
    total_elapsed = time.time() - t_global

    print(f"\n  {'─' * 85}")
    print(f"  总计: {len(registry)} 个样本, {total_cells:,} 个细胞")
    print(f"  总耗时: {total_elapsed:.1f}s ({total_elapsed/60:.1f}min)")
    print(f"\n{'=' * 65}")
    print("  P1 完成")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
