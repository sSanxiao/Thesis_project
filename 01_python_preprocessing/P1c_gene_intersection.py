#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
P1c_gene_intersection.py

Cross-dataset panel gene-intersection analysis: for each species,
compute pairwise panel overlaps, common gene lists, and per-dataset
unique genes (handling human/mouse case mismatch).

Input : the 22 filtered_matrix.h5 files produced by P1
Output: <RESULTS_DIR>/P1_Results/Gene_Intersection/
Run   : python3 P1c_gene_intersection.py

功能: 跨数据集 panel 基因交集分析
旧输出目录: Results_New/P1_Results/Gene_Intersection/
      - pairwise_intersection_matrix_mouse.csv   (鼠 两两交集矩阵)
      - pairwise_intersection_matrix_human.csv   (人 两两交集矩阵)
      - pairwise_intersection_percent_human.csv  (人 两两交集占比矩阵)
      - common_genes_mouse.txt                   (鼠 公共基因列表)
      - common_genes_human.txt                   (人 公共基因列表)
      - common_genes_cross_species.txt           (跨物种 大写统一后公共基因)
      - unique_genes_per_dataset.csv             (各数据集独有基因)
      - preview_internal_comparison.csv           (Preview 三样本内部比较)
      - full_intersection_report.txt             (完整文字报告)
============================================================
"""

import os
import json
import h5py
import numpy as np
import pandas as pd
from itertools import combinations

# =============================================================
# 配置
# =============================================================

# Paths are resolved from this script's location and from RESULTS_DIR
# (see config/paths.py); no hardcoded home directories.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REGISTRY_PATH = os.path.join(SCRIPT_DIR, "sample_registry.json")
P1_RESULT_DIR = os.path.join(os.environ.get("RESULTS_DIR", "./results"), "P1_Results")
OUT_DIR = os.path.join(P1_RESULT_DIR, "Gene_Intersection")

# =============================================================
# 辅助函数
# =============================================================

def read_genes_from_h5(h5_path):
    """从 P1 输出的 filtered_matrix.h5 中读取基因名列表"""
    with h5py.File(h5_path, "r") as f:
        # P1 输出的 h5 统一使用 matrix/ 路径
        raw = f["matrix/features/name"][:]
    return [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in raw]


def intersection_matrix(gene_dict):
    """
    计算两两交集矩阵。
    gene_dict: {name: set(genes), ...}
    返回: (count_df, percent_df)
    """
    names = list(gene_dict.keys())
    n = len(names)
    counts = np.zeros((n, n), dtype=int)
    percents = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            inter = len(gene_dict[names[i]] & gene_dict[names[j]])
            counts[i, j] = inter
            # 占比 = 交集 / min(两者基因数)，衡量小panel被大panel覆盖的程度
            min_size = min(len(gene_dict[names[i]]), len(gene_dict[names[j]]))
            percents[i, j] = round(inter / min_size * 100, 1) if min_size > 0 else 0

    count_df = pd.DataFrame(counts, index=names, columns=names)
    percent_df = pd.DataFrame(percents, index=names, columns=names)
    return count_df, percent_df


# =============================================================
# 主流程
# =============================================================

def main():
    print("=" * 65)
    print("  P1c_gene_intersection.py — 跨数据集基因交集分析")
    print("=" * 65)

    os.makedirs(OUT_DIR, exist_ok=True)

    # 读取 registry
    with open(REGISTRY_PATH, "r") as f:
        registry = json.load(f)

    # =========================================================
    # 第一步：从所有 P1 输出的 h5 中提取基因列表
    # =========================================================
    print("\n  提取基因列表...")

    sample_genes = {}  # sample_name → set(genes)
    dataset_genes = {}  # dataset_name → set(genes)  (同一数据集的样本合并)
    sample_to_dataset = {}
    sample_species = {}

    for sample_name, info in registry.items():
        h5_path = os.path.join(P1_RESULT_DIR, sample_name, "filtered_matrix.h5")
        if not os.path.exists(h5_path):
            print(f"    ⚠ {sample_name}: filtered_matrix.h5 不存在，跳过")
            continue

        genes = read_genes_from_h5(h5_path)
        gene_set = set(genes)
        sample_genes[sample_name] = gene_set

        # 提取数据集名 (sample_name 的第一段，如 "Alzheimer_Mouse")
        dataset = sample_name.split("/")[0]
        sample_to_dataset[sample_name] = dataset
        sample_species[sample_name] = info["species"]

        if dataset not in dataset_genes:
            dataset_genes[dataset] = gene_set.copy()
        else:
            # 同一数据集内不同样本可能有不同的基因列表（如 Preview）
            # 用并集记录该数据集出现过的所有基因
            dataset_genes[dataset] = dataset_genes[dataset] | gene_set

        print(f"    {sample_name}: {len(genes)} genes")

    # =========================================================
    # 第二步：检查同一数据集内样本间基因列表是否一致
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  同一数据集内样本间基因一致性检查")
    print(f"{'=' * 65}")

    # 按数据集分组
    dataset_sample_map = {}
    for sn in sample_genes:
        ds = sample_to_dataset[sn]
        dataset_sample_map.setdefault(ds, []).append(sn)

    for ds, samples in sorted(dataset_sample_map.items()):
        if len(samples) == 1:
            print(f"\n  {ds}: 单样本，无需检查")
            continue

        gene_sets = [sample_genes[s] for s in samples]
        all_same = all(gs == gene_sets[0] for gs in gene_sets)

        if all_same:
            print(f"\n  {ds}: {len(samples)} 个样本基因列表完全一致 ({len(gene_sets[0])} genes)")
        else:
            print(f"\n  {ds}: ⚠ 样本间基因列表不一致!")
            union = set().union(*gene_sets)
            inter = set.intersection(*gene_sets)
            print(f"    并集: {len(union)} genes, 交集: {len(inter)} genes")
            for s, gs in zip(samples, gene_sets):
                unique_to_this = gs - inter
                short_name = s.split("/")[-1]
                if unique_to_this:
                    print(f"    {short_name}: {len(gs)} genes, "
                          f"独有 {len(unique_to_this)} genes")
                    if len(unique_to_this) <= 50:
                        print(f"      独有基因: {sorted(unique_to_this)}")
                else:
                    print(f"    {short_name}: {len(gs)} genes, 无独有基因")

    # =========================================================
    # 第三步：Brain_Human_Preview 内部详细比较
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  Brain_Human_Preview 内部比较")
    print(f"{'=' * 65}")

    preview_samples = [s for s in sample_genes if "Brain_Human_Preview" in s]
    if len(preview_samples) >= 2:
        preview_dict = {s.split("/")[-1]: sample_genes[s] for s in preview_samples}
        p_count, p_pct = intersection_matrix(preview_dict)
        print(f"\n  交集数量矩阵:")
        print(p_count.to_string())
        print(f"\n  交集占比矩阵 (% of min):")
        print(p_pct.to_string())

        p_count.to_csv(os.path.join(OUT_DIR, "preview_internal_comparison.csv"))

        # Preview 三样本公共基因
        preview_common = set.intersection(*[sample_genes[s] for s in preview_samples])
        print(f"\n  Preview 三样本公共基因: {len(preview_common)}")

    # =========================================================
    # 第四步：按物种分组，数据集级别两两交集
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  第一层 + 第二层：同物种内数据集两两交集 + 公共基因")
    print(f"{'=' * 65}")

    # 确定每个数据集的物种（取第一个样本的species）
    dataset_species = {}
    for sn, sp in sample_species.items():
        ds = sample_to_dataset[sn]
        dataset_species[ds] = sp

    # 为了交集分析，每个数据集用"所有样本的交集"作为该数据集的基因列表
    # （而非并集，因为跨样本比较需要所有样本都有的基因）
    dataset_common_genes = {}
    for ds, samples in dataset_sample_map.items():
        gene_sets = [sample_genes[s] for s in samples]
        dataset_common_genes[ds] = set.intersection(*gene_sets)

    # Preview 特殊处理：三个样本基因列表不同
    # 拆分为两组显示
    # 但在人类交集分析中，用 Preview 三样本的公共基因代表该数据集
    if "Brain_Human_Preview" in dataset_common_genes:
        preview_common_count = len(dataset_common_genes["Brain_Human_Preview"])
        print(f"\n  注意: Brain_Human_Preview 内部不一致，"
              f"用三样本公共基因 ({preview_common_count}) 代表该数据集")

    # --- 鼠 ---
    mouse_datasets = {ds: dataset_common_genes[ds]
                      for ds in dataset_common_genes if dataset_species.get(ds) == "mouse"}
    print(f"\n  === 鼠 ({len(mouse_datasets)} 个数据集) ===")
    for ds, gs in sorted(mouse_datasets.items()):
        print(f"    {ds}: {len(gs)} genes")

    if len(mouse_datasets) >= 2:
        m_count, m_pct = intersection_matrix(mouse_datasets)
        print(f"\n  鼠 两两交集数量:")
        print("  " + m_count.to_string().replace("\n", "\n  "))
        print(f"\n  鼠 两两交集占比 (% of min):")
        print("  " + m_pct.to_string().replace("\n", "\n  "))
        m_count.to_csv(os.path.join(OUT_DIR, "pairwise_intersection_matrix_mouse.csv"))

        mouse_common = set.intersection(*mouse_datasets.values())
        print(f"\n  鼠 全数据集公共基因: {len(mouse_common)}")
        with open(os.path.join(OUT_DIR, "common_genes_mouse.txt"), "w") as f:
            for g in sorted(mouse_common):
                f.write(g + "\n")
    elif len(mouse_datasets) == 1:
        ds_name = list(mouse_datasets.keys())[0]
        mouse_common = mouse_datasets[ds_name]
        print(f"  只有一个鼠数据集，公共基因 = 该数据集的 {len(mouse_common)} genes")

    # --- 人 ---
    human_datasets = {ds: dataset_common_genes[ds]
                      for ds in dataset_common_genes if dataset_species.get(ds) == "human"}
    print(f"\n  === 人 ({len(human_datasets)} 个数据集) ===")
    for ds, gs in sorted(human_datasets.items()):
        print(f"    {ds}: {len(gs)} genes")

    if len(human_datasets) >= 2:
        h_count, h_pct = intersection_matrix(human_datasets)
        print(f"\n  人 两两交集数量:")
        print("  " + h_count.to_string().replace("\n", "\n  "))
        print(f"\n  人 两两交集占比 (% of min):")
        print("  " + h_pct.to_string().replace("\n", "\n  "))
        h_count.to_csv(os.path.join(OUT_DIR, "pairwise_intersection_matrix_human.csv"))
        h_pct.to_csv(os.path.join(OUT_DIR, "pairwise_intersection_percent_human.csv"))

        human_common = set.intersection(*human_datasets.values())
        print(f"\n  人 全数据集公共基因: {len(human_common)}")
        with open(os.path.join(OUT_DIR, "common_genes_human.txt"), "w") as f:
            for g in sorted(human_common):
                f.write(g + "\n")

    # =========================================================
    # 第五步：各数据集独有基因
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  各数据集独有基因（只在该数据集中出现的基因）")
    print(f"{'=' * 65}")

    unique_records = []

    # 鼠独有
    all_mouse_genes = set().union(*mouse_datasets.values()) if mouse_datasets else set()
    for ds, gs in sorted(mouse_datasets.items()):
        others = all_mouse_genes - gs
        unique = gs - set().union(*[v for k, v in mouse_datasets.items() if k != ds])
        print(f"  {ds}: {len(unique)} 独有基因 (在其他鼠数据集中没有)")
        for g in sorted(unique):
            unique_records.append({"dataset": ds, "species": "mouse", "gene": g})

    # 人独有
    all_human_genes = set().union(*human_datasets.values()) if human_datasets else set()
    for ds, gs in sorted(human_datasets.items()):
        unique = gs - set().union(*[v for k, v in human_datasets.items() if k != ds])
        print(f"  {ds}: {len(unique)} 独有基因 (在其他人数据集中没有)")
        for g in sorted(unique):
            unique_records.append({"dataset": ds, "species": "human", "gene": g})

    if unique_records:
        pd.DataFrame(unique_records).to_csv(
            os.path.join(OUT_DIR, "unique_genes_per_dataset.csv"), index=False)

    # =========================================================
    # 第六步：跨物种交集（大写统一法）
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  第四层：跨物种基因交集（大写统一法）")
    print(f"{'=' * 65}")

    # 鼠公共基因转大写
    mouse_common_upper = {g.upper() for g in mouse_common} if mouse_datasets else set()
    # 人公共基因转大写
    human_common_upper = {g.upper() for g in human_common} if len(human_datasets) >= 2 else set()

    cross_species_common = mouse_common_upper & human_common_upper

    print(f"\n  鼠公共基因 (原始): {len(mouse_common)} → 大写: {len(mouse_common_upper)}")
    print(f"  人公共基因 (原始): {len(human_common)} → 大写: {len(human_common_upper)}")
    print(f"  跨物种公共基因 (大写统一): {len(cross_species_common)}")

    if cross_species_common:
        with open(os.path.join(OUT_DIR, "common_genes_cross_species.txt"), "w") as f:
            for g in sorted(cross_species_common):
                f.write(g + "\n")

        # 展示前20个
        preview_list = sorted(cross_species_common)[:20]
        print(f"  前20个: {preview_list}")
        if len(cross_species_common) > 20:
            print(f"  ... 共 {len(cross_species_common)} 个")

    # 也做每个鼠数据集 vs 每个人数据集的两两跨物种交集
    print(f"\n  各鼠数据集 vs 各人数据集（大写统一后两两交集）:")
    cross_records = []
    for m_ds, m_genes in sorted(mouse_datasets.items()):
        m_upper = {g.upper() for g in m_genes}
        for h_ds, h_genes in sorted(human_datasets.items()):
            h_upper = {g.upper() for g in h_genes}
            inter = len(m_upper & h_upper)
            pct_of_min = round(inter / min(len(m_upper), len(h_upper)) * 100, 1)
            print(f"    {m_ds} ({len(m_upper)}) ∩ {h_ds} ({len(h_upper)}) "
                  f"= {inter} ({pct_of_min}%)")
            cross_records.append({
                "mouse_dataset": m_ds, "mouse_genes": len(m_upper),
                "human_dataset": h_ds, "human_genes": len(h_upper),
                "intersection": inter, "pct_of_min": pct_of_min
            })

    if cross_records:
        pd.DataFrame(cross_records).to_csv(
            os.path.join(OUT_DIR, "cross_species_pairwise.csv"), index=False)

    # =========================================================
    # 第七步：生成完整文字报告
    # =========================================================
    report_path = os.path.join(OUT_DIR, "full_intersection_report.txt")
    with open(report_path, "w", encoding="utf-8") as rpt:
        rpt.write("跨数据集基因交集分析报告\n")
        rpt.write(f"生成时间: P1c_gene_intersection.py\n")
        rpt.write(f"样本数: {len(sample_genes)}\n")
        rpt.write("=" * 60 + "\n\n")

        rpt.write("一、各数据集基因数\n")
        for ds in sorted(dataset_common_genes.keys()):
            sp = dataset_species.get(ds, "?")
            n_samples = len(dataset_sample_map.get(ds, []))
            n_genes = len(dataset_common_genes[ds])
            rpt.write(f"  {ds}: {n_genes} genes, {n_samples} samples, {sp}\n")

        rpt.write(f"\n二、鼠公共基因: {len(mouse_common)}\n")
        if len(human_datasets) >= 2:
            rpt.write(f"三、人公共基因: {len(human_common)}\n")
        rpt.write(f"四、跨物种公共基因 (大写统一): {len(cross_species_common)}\n")

        rpt.write(f"\n五、人类两两交集矩阵\n")
        if len(human_datasets) >= 2:
            rpt.write(h_count.to_string() + "\n")

        rpt.write(f"\n六、独有基因统计\n")
        for ds in sorted(dataset_common_genes.keys()):
            sp = dataset_species.get(ds, "?")
            pool = mouse_datasets if sp == "mouse" else human_datasets
            unique = dataset_common_genes[ds] - set().union(
                *[v for k, v in pool.items() if k != ds])
            rpt.write(f"  {ds}: {len(unique)} 独有基因\n")

    print(f"\n  完整报告已保存: {report_path}")

    # =========================================================
    # 汇总
    # =========================================================
    print(f"\n{'=' * 65}")
    print("  分析完成")
    print(f"{'=' * 65}")
    print(f"  输出目录: {OUT_DIR}")
    print(f"  鼠公共基因: {len(mouse_common)}")
    if len(human_datasets) >= 2:
        print(f"  人公共基因: {len(human_common)}")
    print(f"  跨物种公共基因: {len(cross_species_common)}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
