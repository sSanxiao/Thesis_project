#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
P1a_generate_registry.py

Generate sample_registry.json: a unified configuration for all 22
Xenium samples (6 datasets), recording per-sample path, species,
condition, panel, segmentation and data-quality metadata. This
registry is the single source of truth consumed by the downstream
P1/P2 preprocessing and the R pipeline.

Input : raw Xenium sample folders under DATA_DIR
Output: sample_registry.json (written next to this script)
Run   : python3 P1a_generate_registry.py

功能: 生成 sample_registry.json (22个样本的统一配置文件)
============================================================
"""

import json
import os

# Root directory holding the raw Xenium sample folders.
# Configure via the DATA_DIR environment variable (see config/paths.py).
BASE = os.environ.get("DATA_DIR", "./data")

registry = {
    # ================================================================
    # Alzheimer_Mouse — 6 samples
    # ================================================================
    "Alzheimer_Mouse/TgCRND8_17_9": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "TgCRND8_17_9"),
        "species": "mouse",
        "condition": "TgCRND8_AD",
        "age_months": 17.9,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },
    "Alzheimer_Mouse/TgCRND8_2_5": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "TgCRND8_2_5"),
        "species": "mouse",
        "condition": "TgCRND8_AD",
        "age_months": 2.5,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },
    "Alzheimer_Mouse/TgCRND8_5_7": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "TgCRND8_5_7"),
        "species": "mouse",
        "condition": "TgCRND8_AD",
        "age_months": 5.7,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },
    "Alzheimer_Mouse/Wildtype_13_4": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "Wildtype_13_4"),
        "species": "mouse",
        "condition": "wild_type",
        "age_months": 13.2,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },
    "Alzheimer_Mouse/Wildtype_2_5": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "Wildtype_2_5"),
        "species": "mouse",
        "condition": "wild_type",
        "age_months": 2.5,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },
    "Alzheimer_Mouse/Wildtype_5_7": {
        "path": os.path.join(BASE, "Alzheimer_Mouse", "Wildtype_5_7"),
        "species": "mouse",
        "condition": "wild_type",
        "age_months": 5.7,
        "preservation": "FFPE",
        "panel_name": "Mouse_Brain_347",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.4.0",
        "data_quality_tier": "high",
        "data_source": "10x_official"
    },

    # ================================================================
    # ATRT_Human — 7 samples
    # ================================================================
    "ATRT_Human/GSM8672828": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672828"),
        "species": "human",
        "condition": "ATRT_SHH",
        "sample_name": "ATRT-05",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },
    "ATRT_Human/GSM8672829": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672829"),
        "species": "human",
        "condition": "ATRT_TYR",
        "sample_name": "ATRT-15-RV4",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },
    "ATRT_Human/GSM8672830": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672830"),
        "species": "human",
        "condition": "ATRT_SHH",
        "sample_name": "ATRT-173",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },
    "ATRT_Human/GSM8672831": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672831"),
        "species": "human",
        "condition": "ATRT_MYC",
        "sample_name": "ATRT-207",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832",
        "note": "empty_cell_rate_1.50%"
    },
    "ATRT_Human/GSM8672832": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672832"),
        "species": "human",
        "condition": "ATRT_MYC",
        "sample_name": "ATRT-243",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },
    "ATRT_Human/GSM8672833": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672833"),
        "species": "human",
        "condition": "ATRT_SHH",
        "sample_name": "ATRT-256",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },
    "ATRT_Human/GSM8672834": {
        "path": os.path.join(BASE, "ATRT_Human", "GSM8672834"),
        "species": "human",
        "condition": "ATRT_TYR",
        "sample_name": "ATRT-340",
        "preservation": "FFPE",
        "panel_name": "Human_Probe_475",
        "segmentation": "nuclei_expansion",
        "xoa_version": "unknown",
        "data_quality_tier": "high",
        "data_source": "GSE283832"
    },

    # ================================================================
    # Brain_Human_Preview — 3 samples
    # ================================================================
    "Brain_Human_Preview/Alzheimers": {
        "path": os.path.join(BASE, "Brain_Human_Preview", "Alzheimers"),
        "species": "human",
        "condition": "Alz_preview",
        "preservation": "FFPE",
        "panel_name": "Human_Brain_Preview_374",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.3.0",
        "data_quality_tier": "low",
        "data_source": "10x_preview",
        "note": "preview_data_dev_panel_100gene_addon"
    },
    "Brain_Human_Preview/GBM": {
        "path": os.path.join(BASE, "Brain_Human_Preview", "GBM"),
        "species": "human",
        "condition": "Glio_preview",
        "preservation": "FFPE",
        "panel_name": "Human_Brain_Preview_339",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.3.0",
        "data_quality_tier": "low",
        "data_source": "10x_preview",
        "note": "preview_data_dev_panel_65gene_addon"
    },
    "Brain_Human_Preview/Healthy": {
        "path": os.path.join(BASE, "Brain_Human_Preview", "Healthy"),
        "species": "human",
        "condition": "Healthy_preview",
        "preservation": "FFPE",
        "panel_name": "Human_Brain_Preview_339",
        "segmentation": "nuclei_expansion",
        "xoa_version": "v1.3.0",
        "data_quality_tier": "low",
        "data_source": "10x_preview",
        "note": "preview_data_dev_panel_65gene_addon_only_normal_brain"
    },

    # ================================================================
    # Brain_Mouse — 1 sample
    # ================================================================
    "Brain_Mouse/Single_Sample": {
        "path": os.path.join(BASE, "Brain_Mouse", "Single_Sample"),
        "species": "mouse",
        "condition": "WT_normal",
        "preservation": "FF",
        "panel_name": "Xenium_Prime_5K",
        "segmentation": "multimodal_membrane",
        "xoa_version": "v3.0.0",
        "data_quality_tier": "high",
        "data_source": "10x_official",
        "note": "highest_quality_reference"
    },

    # ================================================================
    # Glioblastoma_Human — 1 sample
    # ================================================================
    "Glioblastoma_Human/Single_Sample": {
        "path": os.path.join(BASE, "Glioblastoma_Human", "Single_Sample"),
        "species": "human",
        "condition": "GBM",
        "preservation": "FFPE",
        "panel_name": "Human_IO_500",
        "segmentation": "multimodal_membrane",
        "xoa_version": "v2.0.0",
        "data_quality_tier": "high",
        "data_source": "10x_official",
        "note": "IO_panel_immune_focused"
    },

    # ================================================================
    # Medulloblastoma_Human — 4 samples
    # ================================================================
    "Medulloblastoma_Human/GSM8840046": {
        "path": os.path.join(BASE, "Medulloblastoma_Human", "GSM8840046"),
        "species": "human",
        "condition": "Medulloblastoma",
        "sample_name": "MB263",
        "preservation": "FF",
        "panel_name": "Human_MB_379",
        "segmentation": "unknown",
        "xoa_version": "v1.7.6.0",
        "data_quality_tier": "high",
        "data_source": "GSE291688"
    },
    "Medulloblastoma_Human/GSM8840047": {
        "path": os.path.join(BASE, "Medulloblastoma_Human", "GSM8840047"),
        "species": "human",
        "condition": "Medulloblastoma",
        "sample_name": "MB266",
        "preservation": "FF",
        "panel_name": "Human_MB_379",
        "segmentation": "unknown",
        "xoa_version": "v1.7.6.0",
        "data_quality_tier": "high",
        "data_source": "GSE291688"
    },
    "Medulloblastoma_Human/GSM8840048": {
        "path": os.path.join(BASE, "Medulloblastoma_Human", "GSM8840048"),
        "species": "human",
        "condition": "Medulloblastoma",
        "sample_name": "MB295",
        "preservation": "FF",
        "panel_name": "Human_MB_379",
        "segmentation": "unknown",
        "xoa_version": "v1.7.6.0",
        "data_quality_tier": "high",
        "data_source": "GSE291688"
    },
    "Medulloblastoma_Human/GSM8840049": {
        "path": os.path.join(BASE, "Medulloblastoma_Human", "GSM8840049"),
        "species": "human",
        "condition": "Medulloblastoma",
        "sample_name": "MB299",
        "preservation": "FF",
        "panel_name": "Human_MB_379",
        "segmentation": "unknown",
        "xoa_version": "v1.7.6.0",
        "data_quality_tier": "high",
        "data_source": "GSE291688"
    },
}

# 输出: write the registry next to this script (use a relative path, no hardcoded home dir).
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_registry.json")
os.makedirs(os.path.dirname(out_path), exist_ok=True)

with open(out_path, "w", encoding="utf-8") as f:
    json.dump(registry, f, indent=2, ensure_ascii=False)

print(f"sample_registry.json 已生成: {out_path}")
print(f"共 {len(registry)} 个样本")

# 验证
for name, info in registry.items():
    p = info["path"]
    h5 = os.path.join(p, "cell_feature_matrix.h5")
    cells = os.path.join(p, "cells.parquet")
    if not os.path.exists(h5):
        print(f"  ⚠ {name}: h5 不存在 ({h5})")
    if not os.path.exists(cells):
        print(f"  ⚠ {name}: cells.parquet 不存在 ({cells})")
