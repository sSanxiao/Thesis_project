# Thesis_project

## Study of Cell Spatial Density in Cell Development and Cancer Oncogenesis

Computational pipeline for the MSc thesis *"Cell-Density-Coupled Gene
Signatures in Medulloblastoma with Extensive Nodularity"*
(Health Informatics, KI / Stockholm University).

This repository contains **analysis code only**. No patient data, raw
matrices, or intermediate results are included. All input datasets are
public and referenced by their GEO accession numbers.

---

## Project design (four stages, per Prof. Feng Zhang)

**Stage 1. Data Collection (~1 month):** Collect publicly available
downloadable Xenium datasets (PubMed / Google Scholar), focusing on brain
and cerebellum development and cancer.

**Stage 2. Data Analysis (~1 month):** Basic analysis using the provided
base scripts.

**Stage 3. Further Analysis:** Identify spatial-density-related genes.
Perform integrated analysis (statistics / machine learning) of these genes
against single-cell and bulk data to explore their relationship with
development and oncogenesis.

**Stage 4. Result Summary:** Three tiers of conclusions depending on
correlation strength — (I) strong: build and validate predictive models
(e.g. on TCGA); (II) moderate: descriptive summaries of patterns / pathways;
(III) weak: characterise a few example genes via literature review.

---

## Pipeline overview

| Stage | Scripts | Purpose |
|-------|---------|---------|
| Preprocessing (Python) | `01_python_preprocessing/P1_*`, `P2_*` | Build sample registry, load/filter Xenium data, gene intersection, compute 5 cell-density estimators (3 KNN-based + Voronoi + Delaunay) |
| Core analysis (R) | `02_R_core_pipeline/R1`–`R9` | Seurat object construction, SCTransform, density–gene correlation, filtering, visualization, cell-state coupling, sample integration, cross-dataset comparison, tier decision |
| Signature & validation (R) | `02_R_core_pipeline/R10`–`R16` | Cavalli 2017 bulk validation, signature expansion (sig_94), MBEN subtype validation, signature cleanup, feasibility checks |
| Spatial & pseudotime (R) | `02_R_core_pipeline/R17`–`R20` | Seurat assembly, signature scoring, pseudotime, spatial analyses, QC, multiple-testing correction |
| External validation (R) | `02_R_core_pipeline/R21*` | Ghasemi 2024 snRNA-seq + Aldinger 2021 fetal cerebellum validation |
| Setup | `setup/install_deps.R` | R package installation |
| Archive | `archive/` | Debugging / version-iteration scripts retained for provenance; not part of the clean pipeline |

---

## Datasets (all public)

| Dataset | Accession | Role |
|---------|-----------|------|
| Xenium spatial cohort (22 samples, 6 datasets) | various GEO (e.g. GSE283832) | Primary spatial data |
| Cavalli 2017 bulk | n=763 | Bulk validation |
| GSE124814 bulk | n=476 | Bulk validation |
| Aldinger 2021 fetal cerebellum | n=69,174 cells | Developmental reference |
| Ghasemi 2024 snRNA-seq | GSM6604617/18/20/21 | External single-nucleus validation |

---

## Reproducing the pipeline

### 1. Configure data paths

Scripts read input/output locations from environment variables (no hardcoded
paths). Set them to your local directories before running:

```bash
export DATA_DIR=/path/to/Xenium_datasets
export EXTDATA_DIR=/path/to/external_data
export RESULTS_DIR=/path/to/results
```

(On Windows PowerShell: `$env:DATA_DIR = "D:\..."`.)

See `config/paths.R` and `config/paths.py` for the full list.

### 2. Install dependencies

```bash
Rscript setup/install_deps.R
```

Key versions (compatibility-sensitive): Seurat 5.2.1, SCTransform,
ggplot2 4.0.2, harmony. A frozen `sessionInfo()` is recommended for exact
reproduction.

### 3. Run in order

```bash
# Preprocessing (run in this order)
python 01_python_preprocessing/P1a_generate_registry.py
python 01_python_preprocessing/P1b_data_loading.py
python 01_python_preprocessing/P1c_gene_intersection.py
python 01_python_preprocessing/P2_density_calculation.py

# R pipeline (run R01 → R21d in order)
Rscript 02_R_core_pipeline/R01_build_seurat.R
# ... etc

# Note: numeric prefixes are zero-padded (R01..R09) and the preprocessing
# steps use a/b/c suffixes so that the on-disk/file-listing order matches
# the intended execution order.
```

---

## Citation

If you use this code, please cite the associated thesis. Method development
builds on prior work (Wang et al., *Cellular Oncology*, 2024).
