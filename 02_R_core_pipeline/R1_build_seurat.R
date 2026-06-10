#!/usr/bin/env Rscript
# ============================================================================
# R1_build_seurat.R
# 功能：从P1和P2的输出构建Seurat对象
# 输入：
#   - P1: filtered_matrix.h5 + cell_metadata.csv (每样本)
#   - P2: cell_density.csv (每样本，含5列密度值)
#   - sample_registry.json (全局配置)
# 输出：
#   - 每样本一个 .rds 文件 (Seurat对象)
#   - ALL_SAMPLES_R1_QC.csv (22行汇总表)
# Run (EN): Rscript R1_build_seurat.R
#   Purpose: build one Seurat object per sample from the P1/P2 outputs,
#            attaching the 5 cell-density estimators as cell metadata.
#   Paths configured via env vars (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# ============================================================================
# 全局配置
# ============================================================================

# Configurable roots (see config/paths.R); override via environment variables.
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")

REGISTRY_PATH  <- file.path(DATA_DIR, "sample_registry.json")
P1_DIR         <- file.path(RESULTS_DIR, "P1_Results")
P2_DIR         <- file.path(RESULTS_DIR, "P2_Results")
OUTPUT_DIR     <- file.path(RESULTS_DIR, "R1_Results")

cat("============================================\n")
cat("R1: 构建 Seurat 对象\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 读取 sample_registry
# ============================================================================

registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
n_samples <- length(sample_names)
cat("从 sample_registry 读取到", n_samples, "个样本\n\n")

# ============================================================================
# 逐样本处理
# ============================================================================

qc_list <- list()

for (i in seq_along(sample_names)) {

  sample_name <- sample_names[i]
  sample_info <- registry[[sample_name]]

  cat("========================================\n")
  cat(sprintf("[%d/%d] %s\n", i, n_samples, sample_name))
  cat("========================================\n")

  # --- 确定文件路径 ---

  # 从sample_name中提取 dataset/sample 结构
  parts <- strsplit(sample_name, "/")[[1]]
  dataset_name <- parts[1]
  sample_subname <- parts[2]

  h5_path       <- file.path(P1_DIR, dataset_name, sample_subname, "filtered_matrix.h5")
  meta_path     <- file.path(P1_DIR, dataset_name, sample_subname, "cell_metadata.csv")
  density_path  <- file.path(P2_DIR, dataset_name, sample_subname, "cell_density.csv")

  # --- 检查文件是否存在 ---

  missing_files <- c()
  if (!file.exists(h5_path))      missing_files <- c(missing_files, "filtered_matrix.h5")
  if (!file.exists(meta_path))    missing_files <- c(missing_files, "cell_metadata.csv")
  if (!file.exists(density_path)) missing_files <- c(missing_files, "cell_density.csv")

  if (length(missing_files) > 0) {
    cat("  ✗ 缺少文件:", paste(missing_files, collapse = ", "), "\n")
    cat("  跳过该样本\n\n")
    next
  }

  # ---------------------------------------------------------------
  # 第一步：读 h5 构建 Seurat 对象
  # ---------------------------------------------------------------

  cat("  [1/6] 读取 filtered_matrix.h5 ...\n")
  counts <- Read10X_h5(h5_path)

  n_genes <- nrow(counts)
  n_cells <- ncol(counts)
  cat(sprintf("        %d 基因 × %d 细胞\n", n_genes, n_cells))

  # 创建Seurat对象（不做任何过滤，min.cells=0, min.features=0）
  seurat_obj <- CreateSeuratObject(
    counts    = counts,
    project   = sample_name,
    min.cells = 0,
    min.features = 0
  )

  cat(sprintf("        Seurat对象创建完成: %d 细胞\n", ncol(seurat_obj)))

  # ---------------------------------------------------------------
  # 第二步：合并细胞元数据
  # ---------------------------------------------------------------

  cat("  [2/6] 合并 cell_metadata.csv ...\n")
  cell_meta <- fread(meta_path)

  # 确保cell_id是字符型
  cell_meta[, cell_id := as.character(cell_id)]

  # 以Seurat对象中的细胞名为基准对齐
  seurat_cells <- colnames(seurat_obj)

  # 检查一致性
  meta_cells <- cell_meta$cell_id
  common_cells <- intersect(seurat_cells, meta_cells)
  n_only_seurat <- length(setdiff(seurat_cells, meta_cells))
  n_only_meta   <- length(setdiff(meta_cells, seurat_cells))

  if (n_only_seurat > 0 || n_only_meta > 0) {
    cat(sprintf("        ⚠ cell_id不完全一致: h5独有=%d, meta独有=%d, 交集=%d\n",
                n_only_seurat, n_only_meta, length(common_cells)))
  } else {
    cat(sprintf("        ✓ cell_id完全一致: %d 细胞\n", length(common_cells)))
  }

  # 按Seurat对象中的细胞顺序对齐
  cell_meta <- cell_meta[match(seurat_cells, cell_id), ]

  # 添加到 meta.data
  meta_cols <- c("x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area")
  for (col in meta_cols) {
    if (col %in% names(cell_meta)) {
      seurat_obj[[col]] <- cell_meta[[col]]
    } else {
      cat(sprintf("        ⚠ 列 %s 不存在于 cell_metadata.csv 中, 跳过\n", col))
    }
  }

  cat("        ✓ 元数据合并完成\n")

  # ---------------------------------------------------------------
  # 第三步：合并密度值
  # ---------------------------------------------------------------

  cat("  [3/6] 合并 cell_density.csv ...\n")
  density_data <- fread(density_path)
  density_data[, cell_id := as.character(cell_id)]

  # 按Seurat对象中的细胞顺序对齐
  density_data <- density_data[match(seurat_cells, cell_id), ]

  # 密度列名
  density_cols <- c("density_knn_aggr_2nd_diff", "density_knn_main_piecewise", "density_knn_cons_max_dist",
                    "density_voronoi", "density_delaunay")

  for (col in density_cols) {
    if (col %in% names(density_data)) {
      seurat_obj[[col]] <- density_data[[col]]
    } else {
      cat(sprintf("        ⚠ 列 %s 不存在于 cell_density.csv 中, 跳过\n", col))
    }
  }

  # 统计密度NA数
  n_na_vor <- sum(is.na(seurat_obj$density_voronoi))
  n_na_del <- sum(is.na(seurat_obj$density_delaunay))
  cat(sprintf("        ✓ 密度值合并完成 (Voronoi NA=%d, Delaunay NA=%d)\n", n_na_vor, n_na_del))

  # ---------------------------------------------------------------
  # 第四步：合并样本级元数据
  # ---------------------------------------------------------------

  cat("  [4/6] 合并样本级元数据 ...\n")

  sample_level_fields <- c("species", "condition", "preservation",
                           "panel_name", "segmentation", "data_quality_tier")

  for (field in sample_level_fields) {
    if (!is.null(sample_info[[field]])) {
      seurat_obj[[field]] <- sample_info[[field]]
    }
  }

  # 额外添加 dataset 和 sample_name
  seurat_obj[["dataset"]]     <- dataset_name
  seurat_obj[["sample_name"]] <- sample_name

  cat("        ✓ 样本级元数据合并完成\n")

  # ---------------------------------------------------------------
  # 第五步：验证
  # ---------------------------------------------------------------

  cat("  [5/6] 验证 ...\n")

  # 验证维度
  n_meta_rows <- nrow(seurat_obj@meta.data)
  n_count_cols <- ncol(seurat_obj)
  if (n_meta_rows != n_count_cols) {
    cat(sprintf("        ✗ 维度不一致! meta.data行数=%d, counts列数=%d\n",
                n_meta_rows, n_count_cols))
  } else {
    cat(sprintf("        ✓ 维度一致: %d 细胞\n", n_count_cols))
  }

  # 验证坐标无NA
  n_na_x <- sum(is.na(seurat_obj$x_centroid))
  n_na_y <- sum(is.na(seurat_obj$y_centroid))
  if (n_na_x > 0 || n_na_y > 0) {
    cat(sprintf("        ⚠ 坐标存在NA: x=%d, y=%d\n", n_na_x, n_na_y))
  } else {
    cat("        ✓ 坐标无NA\n")
  }

  # 验证KNN密度无NA（KNN不应有NA，Voronoi/Delaunay可以有）
  n_na_knn <- sum(is.na(seurat_obj$density_knn_main_piecewise))
  if (n_na_knn > 0) {
    cat(sprintf("        ⚠ KNN_main密度有 %d 个NA\n", n_na_knn))
  } else {
    cat("        ✓ KNN密度无NA\n")
  }

  # 验证nCount_RNA和transcript_counts的一致性
  ncount_cor <- cor(seurat_obj$nCount_RNA, seurat_obj$transcript_counts,
                    method = "spearman", use = "complete.obs")
  cat(sprintf("        nCount_RNA ↔ transcript_counts Spearman ρ = %.4f\n", ncount_cor))
  if (ncount_cor < 0.99) {
    cat("        ⚠ 相关性偏低，可能存在对照feature计数差异\n")
  }

  cat("        ✓ 验证完成\n")

  # ---------------------------------------------------------------
  # 第六步：输出
  # ---------------------------------------------------------------

  cat("  [6/6] 保存 .rds ...\n")

  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(out_dir, paste0(sample_subname, "_seurat.rds"))

  saveRDS(seurat_obj, rds_path)

  rds_size_mb <- round(file.size(rds_path) / 1024 / 1024, 1)
  cat(sprintf("        ✓ 保存完成: %s (%.1f MB)\n", rds_path, rds_size_mb))

  # ---------------------------------------------------------------
  # 收集QC指标
  # ---------------------------------------------------------------

  qc_list[[sample_name]] <- data.frame(
    sample_name        = sample_name,
    dataset            = dataset_name,
    species            = ifelse(!is.null(sample_info$species), sample_info$species, NA),
    condition          = ifelse(!is.null(sample_info$condition), sample_info$condition, NA),
    data_quality_tier  = ifelse(!is.null(sample_info$data_quality_tier), sample_info$data_quality_tier, NA),
    n_genes            = n_genes,
    n_cells            = ncol(seurat_obj),
    median_nCount      = median(seurat_obj$nCount_RNA),
    median_nFeature    = median(seurat_obj$nFeature_RNA),
    median_density_knn = median(seurat_obj$density_knn_main_piecewise, na.rm = TRUE),
    median_cell_area   = median(seurat_obj$cell_area, na.rm = TRUE),
    vor_na_count       = n_na_vor,
    del_na_count       = n_na_del,
    ncount_tc_cor      = round(ncount_cor, 4),
    rds_size_mb        = rds_size_mb,
    stringsAsFactors   = FALSE
  )

  cat("\n")
}

# ============================================================================
# 汇总QC
# ============================================================================

cat("============================================\n")
cat("汇总 QC\n")
cat("============================================\n")

qc_df <- do.call(rbind, qc_list)
rownames(qc_df) <- NULL

qc_path <- file.path(OUTPUT_DIR, "ALL_SAMPLES_R1_QC.csv")
fwrite(qc_df, qc_path)
cat("QC汇总已保存:", qc_path, "\n\n")

# 打印汇总表
cat("--- 22 样本 QC 汇总 ---\n")
cat(sprintf("%-45s %6s %8s %8s %8s %10s %8s %6s\n",
            "sample", "genes", "cells", "med_nCt", "med_nFt", "med_dens", "med_area", "rds_MB"))
cat(paste(rep("-", 130), collapse = ""), "\n")

for (j in 1:nrow(qc_df)) {
  row <- qc_df[j, ]
  cat(sprintf("%-45s %6d %8d %8.0f %8.0f %10.4f %8.1f %6.1f\n",
              row$sample_name,
              row$n_genes,
              row$n_cells,
              row$median_nCount,
              row$median_nFeature,
              row$median_density_knn,
              row$median_cell_area,
              row$rds_size_mb))
}

# 汇总统计
cat("\n--- 总计 ---\n")
cat("总样本数:", nrow(qc_df), "\n")
cat("总细胞数:", format(sum(qc_df$n_cells), big.mark = ","), "\n")
cat("总 .rds 大小:", round(sum(qc_df$rds_size_mb) / 1024, 2), "GB\n")

# 按数据集分组统计
cat("\n--- 按数据集分组 ---\n")
for (ds in unique(qc_df$dataset)) {
  sub <- qc_df[qc_df$dataset == ds, ]
  cat(sprintf("  %-30s %2d 样本, %10s 细胞, med_nCount=%.0f, med_nFeature=%.0f\n",
              ds,
              nrow(sub),
              format(sum(sub$n_cells), big.mark = ","),
              median(sub$median_nCount),
              median(sub$median_nFeature)))
}

# nCount_RNA vs transcript_counts 一致性
cat("\n--- nCount_RNA ↔ transcript_counts 相关性 ---\n")
low_cor <- qc_df[qc_df$ncount_tc_cor < 0.99, ]
if (nrow(low_cor) > 0) {
  cat("以下样本相关性 < 0.99:\n")
  for (j in 1:nrow(low_cor)) {
    cat(sprintf("  %s: ρ = %.4f\n", low_cor$sample_name[j], low_cor$ncount_tc_cor[j]))
  }
  cat("\n说明: transcript_counts (来自cells.parquet原始值) 可能包含了对照feature的计数,\n")
  cat("而 nCount_RNA (来自过滤后h5) 只包含靶标基因的计数。\n")
  cat("差异是正常的, 后续分析使用 nCount_RNA。\n")
} else {
  cat("所有样本 ρ ≥ 0.99, 一致性良好。\n")
}

cat("\n============================================\n")
cat("R1 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
