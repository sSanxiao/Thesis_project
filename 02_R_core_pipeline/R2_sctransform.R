#!/usr/bin/env Rscript
# ============================================================================
# R2_sctransform.R
# 功能：对每个样本做SCTransform标准化 + PCA + UMAP + 聚类
# 输入：R1输出的 .rds 文件（22个样本）
# 输出：更新后的 .rds 文件 + 诊断图 + QC汇总
# Run (EN): Rscript R2_sctransform.R
#   Purpose: per-sample SCTransform normalization + PCA + UMAP + clustering.
#   Paths configured via env vars (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

# ============================================================================
# 全局配置
# ============================================================================

# Configurable roots (see config/paths.R); override via environment variables.
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")

REGISTRY_PATH  <- file.path(DATA_DIR, "sample_registry.json")
R1_DIR         <- file.path(RESULTS_DIR, "R1_Results")
OUTPUT_DIR     <- file.path(RESULTS_DIR, "R2_Results")
CLUSTER_RES    <- 0.8
MAX_PCS_COMPUTE <- 50   # 先算50个PC，再自动选最优数量

cat("============================================\n")
cat("R2: SCTransform + PCA + UMAP + 聚类\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("聚类分辨率:", CLUSTER_RES, "\n\n")

# ============================================================================
# 读取 sample_registry
# ============================================================================

registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
n_samples <- length(sample_names)
cat("从 sample_registry 读取到", n_samples, "个样本\n\n")

# ============================================================================
# 自动选PC数的函数（二阶差分法）
# ============================================================================

auto_select_pcs <- function(seurat_obj, max_pcs = 50) {
  # 获取每个PC的标准差
  sdev <- seurat_obj[["pca"]]@stdev
  n_available <- length(sdev)
  
  if (n_available < 5) {
    return(list(n_pcs = n_available, sdev = sdev))
  }
  
  # 方差解释率（标准差的平方，归一化）
  var_explained <- sdev^2
  var_pct <- var_explained / sum(var_explained) * 100
  
  # 累计方差解释率
  cum_var_pct <- cumsum(var_pct)
  
  # 方法1：二阶差分找拐点
  if (length(var_pct) >= 4) {
    first_diff <- diff(var_pct)
    second_diff <- diff(first_diff)
    # 二阶差分最大的位置（曲线弯曲最剧烈处）
    elbow_idx <- which.max(abs(second_diff)) + 1  # +1 因为差分损失元素
    # 至少选5个PC，最多选max_pcs个
    n_pcs_elbow <- max(5, min(elbow_idx + 2, max_pcs, n_available))  # +2 给一点余量
  } else {
    n_pcs_elbow <- n_available
  }
  
  # 方法2：方差解释率低于1%的第一个PC
  below_threshold <- which(var_pct < 1)
  if (length(below_threshold) > 0) {
    n_pcs_threshold <- max(5, below_threshold[1] - 1)
  } else {
    n_pcs_threshold <- n_available
  }
  
  # 取两种方法的较大值（保守一些，宁可多保留几个有意义的PC）
  n_pcs <- min(max(n_pcs_elbow, n_pcs_threshold), n_available, max_pcs)
  
  return(list(
    n_pcs = n_pcs,
    n_pcs_elbow = n_pcs_elbow,
    n_pcs_threshold = n_pcs_threshold,
    sdev = sdev,
    var_pct = var_pct,
    cum_var_pct = cum_var_pct
  ))
}

# ============================================================================
# ElbowPlot 诊断图函数
# ============================================================================

plot_elbow <- function(pc_info, sample_name, out_path) {
  n <- length(pc_info$var_pct)
  df <- data.frame(
    PC = 1:n,
    VarPct = pc_info$var_pct,
    CumVarPct = pc_info$cum_var_pct
  )
  
  p1 <- ggplot(df, aes(x = PC, y = VarPct)) +
    geom_point(size = 1.5) +
    geom_line() +
    geom_vline(xintercept = pc_info$n_pcs, color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = pc_info$n_pcs + 1, y = max(df$VarPct) * 0.8,
             label = paste0("Selected: ", pc_info$n_pcs, " PCs"),
             color = "red", hjust = 0, size = 3.5) +
    labs(title = paste0(sample_name, " - ElbowPlot"),
         x = "PC", y = "Variance Explained (%)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  p2 <- ggplot(df, aes(x = PC, y = CumVarPct)) +
    geom_point(size = 1.5) +
    geom_line() +
    geom_vline(xintercept = pc_info$n_pcs, color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = df$CumVarPct[pc_info$n_pcs], color = "blue",
               linetype = "dotted", linewidth = 0.5) +
    annotate("text", x = pc_info$n_pcs + 1, y = df$CumVarPct[pc_info$n_pcs] - 5,
             label = paste0(round(df$CumVarPct[pc_info$n_pcs], 1), "% cumulative"),
             color = "blue", hjust = 0, size = 3.5) +
    labs(title = "Cumulative Variance",
         x = "PC", y = "Cumulative Variance Explained (%)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  # 合并两张图
  png(out_path, width = 1200, height = 500, res = 150)
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

# ============================================================================
# UMAP诊断图函数
# ============================================================================

plot_umap_diagnostics <- function(seurat_obj, sample_name, out_dir) {
  
  # 图1：UMAP按聚类着色
  p_cluster <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
                       label = TRUE, label.size = 3, pt.size = 0.1) +
    labs(title = paste0(sample_name, " - Clusters")) +
    theme(plot.title = element_text(size = 10),
          legend.position = "right")
  
  ggsave(file.path(out_dir, "umap_clusters.png"), p_cluster,
         width = 8, height = 6, dpi = 150)
  
  # 图2：UMAP按KNN密度着色
  if ("density_knn_main_piecewise" %in% colnames(seurat_obj@meta.data)) {
    p_density <- FeaturePlot(seurat_obj, reduction = "umap",
                             features = "density_knn_main_piecewise",
                             pt.size = 0.1) +
      scale_color_viridis_c() +
      labs(title = paste0(sample_name, " - KNN Density (main)")) +
      theme(plot.title = element_text(size = 10))
    
    ggsave(file.path(out_dir, "umap_density.png"), p_density,
           width = 8, height = 6, dpi = 150)
  }
  
  # 图3：UMAP按nCount_RNA着色
  p_ncount <- FeaturePlot(seurat_obj, reduction = "umap",
                          features = "nCount_RNA",
                          pt.size = 0.1) +
    scale_color_viridis_c() +
    labs(title = paste0(sample_name, " - nCount_RNA")) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(file.path(out_dir, "umap_ncount.png"), p_ncount,
         width = 8, height = 6, dpi = 150)
  
  # 图4：空间坐标按聚类着色
  if ("x_centroid" %in% colnames(seurat_obj@meta.data)) {
    df_spatial <- data.frame(
      x = seurat_obj$x_centroid,
      y = seurat_obj$y_centroid,
      cluster = seurat_obj$seurat_clusters
    )
    
    p_spatial <- ggplot(df_spatial, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.05, alpha = 0.3) +
      labs(title = paste0(sample_name, " - Spatial Clusters"),
           x = "x (µm)", y = "y (µm)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10)) +
      coord_fixed()
    
    ggsave(file.path(out_dir, "spatial_clusters.png"), p_spatial,
           width = 8, height = 8, dpi = 150)
  }
}

# ============================================================================
# 检查gridExtra是否可用
# ============================================================================

has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)
if (!has_gridExtra) {
  cat("⚠ gridExtra 未安装, ElbowPlot将只输出单面板\n")
  cat("  安装方法: install.packages('gridExtra')\n\n")
}

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
  
  t_start <- Sys.time()
  
  # --- 确定文件路径 ---
  parts <- strsplit(sample_name, "/")[[1]]
  dataset_name <- parts[1]
  sample_subname <- parts[2]
  
  rds_path <- file.path(R1_DIR, dataset_name, sample_subname,
                        paste0(sample_subname, "_seurat.rds"))
  
  if (!file.exists(rds_path)) {
    cat("  ✗ .rds 文件不存在:", rds_path, "\n\n")
    next
  }
  
  # --- 创建输出目录 ---
  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------
  # 第一步：读取 R1 的 Seurat 对象
  # ---------------------------------------------------------------
  
  cat("  [1/6] 读取 .rds ...\n")
  seurat_obj <- readRDS(rds_path)
  n_genes <- nrow(seurat_obj)
  n_cells <- ncol(seurat_obj)
  cat(sprintf("        %d 基因 × %d 细胞\n", n_genes, n_cells))
  
  # ---------------------------------------------------------------
  # 第二步：SCTransform
  # ---------------------------------------------------------------
  
  cat("  [2/6] SCTransform ...\n")
  
  # 确定PCA可用的最大PC数（不能超过min(基因数, 细胞数) - 1）
  max_pcs_possible <- min(n_genes, n_cells, MAX_PCS_COMPUTE) - 1
  npcs_to_compute <- min(MAX_PCS_COMPUTE, max_pcs_possible)
  
  seurat_obj <- SCTransform(
    seurat_obj,
    variable.features.n = min(3000, n_genes),
    return.only.var.genes = FALSE,
    verbose = FALSE
  )
  
  n_var_features <- length(VariableFeatures(seurat_obj))
  cat(sprintf("        ✓ SCTransform 完成, %d 个高变基因\n", n_var_features))
  
  # ---------------------------------------------------------------
  # 第三步：PCA
  # ---------------------------------------------------------------
  
  cat("  [3/6] PCA ...\n")
  
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs_to_compute, verbose = FALSE)
  
  # 自动选择PC数
  pc_info <- auto_select_pcs(seurat_obj, max_pcs = npcs_to_compute)
  n_pcs <- pc_info$n_pcs
  cum_var <- round(pc_info$cum_var_pct[n_pcs], 1)
  
  cat(sprintf("        ✓ PCA 完成, 自动选择 %d 个PC (累计方差 %.1f%%)\n", n_pcs, cum_var))
  cat(sprintf("          (二阶差分拐点=%d, 1%%阈值=%d)\n",
              pc_info$n_pcs_elbow, pc_info$n_pcs_threshold))
  
  # ElbowPlot诊断图
  if (has_gridExtra) {
    elbow_path <- file.path(out_dir, "elbow_plot.png")
    tryCatch({
      plot_elbow(pc_info, sample_name, elbow_path)
      cat("        ✓ ElbowPlot 已保存\n")
    }, error = function(e) {
      cat(sprintf("        ⚠ ElbowPlot 生成失败: %s\n", e$message))
    })
  }
  
  # 存储选择的PC数到meta.data（方便后续追溯）
  seurat_obj[["n_pcs_selected"]] <- n_pcs
  
  # ---------------------------------------------------------------
  # 第四步：UMAP
  # ---------------------------------------------------------------
  
  cat("  [4/6] UMAP ...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
  cat("        ✓ UMAP 完成\n")
  
  # ---------------------------------------------------------------
  # 第五步：聚类
  # ---------------------------------------------------------------
  
  cat("  [5/6] 聚类 (resolution=", CLUSTER_RES, ") ...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = CLUSTER_RES, verbose = FALSE)
  
  n_clusters <- length(levels(Idents(seurat_obj)))
  cat(sprintf("        ✓ 聚类完成, %d 个 clusters\n", n_clusters))
  
  # ---------------------------------------------------------------
  # 第六步：保存 + 诊断图
  # ---------------------------------------------------------------
  
  cat("  [6/6] 保存 + 诊断图 ...\n")
  
  # 保存更新后的Seurat对象
  out_rds <- file.path(out_dir, paste0(sample_subname, "_seurat_R2.rds"))
  saveRDS(seurat_obj, out_rds)
  rds_size_mb <- round(file.size(out_rds) / 1024 / 1024, 1)
  cat(sprintf("        ✓ .rds 保存完成 (%.1f MB)\n", rds_size_mb))
  
  # 诊断图
  tryCatch({
    plot_umap_diagnostics(seurat_obj, sample_name, out_dir)
    cat("        ✓ 诊断图已保存\n")
  }, error = function(e) {
    cat(sprintf("        ⚠ 诊断图生成失败: %s\n", e$message))
  })
  
  # --- 计时 ---
  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  cat(sprintf("  完成! 耗时 %.1f 秒\n\n", t_elapsed))
  
  # ---------------------------------------------------------------
  # 收集QC指标
  # ---------------------------------------------------------------
  
  # SCT assay的残差统计
  sct_data <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
  residual_mean <- round(mean(sct_data@x), 4)  # 非零元素均值
  residual_sd <- round(sd(sct_data@x), 4)
  
  qc_list[[sample_name]] <- data.frame(
    sample_name       = sample_name,
    dataset           = dataset_name,
    n_genes           = n_genes,
    n_cells           = n_cells,
    n_var_features    = n_var_features,
    n_pcs_selected    = n_pcs,
    n_pcs_elbow       = pc_info$n_pcs_elbow,
    n_pcs_threshold   = pc_info$n_pcs_threshold,
    cum_var_pct       = cum_var,
    n_clusters        = n_clusters,
    residual_mean     = residual_mean,
    residual_sd       = residual_sd,
    time_seconds      = t_elapsed,
    rds_size_mb       = rds_size_mb,
    stringsAsFactors   = FALSE
  )
  
  # 释放内存
  rm(seurat_obj, sct_data)
  gc(verbose = FALSE)
}

# ============================================================================
# 汇总QC
# ============================================================================

cat("============================================\n")
cat("汇总 QC\n")
cat("============================================\n")

qc_df <- do.call(rbind, qc_list)
rownames(qc_df) <- NULL

qc_path <- file.path(OUTPUT_DIR, "ALL_SAMPLES_R2_QC.csv")
fwrite(qc_df, qc_path)
cat("QC汇总已保存:", qc_path, "\n\n")

# 打印汇总表
cat("--- 22 样本 R2 QC 汇总 ---\n")
cat(sprintf("%-45s %6s %8s %6s %4s %4s %4s %7s %4s %8s %7s\n",
            "sample", "genes", "cells", "varFt", "PCs", "elb", "thr", "cumVar", "clst", "time_s", "rds_MB"))
cat(paste(rep("-", 140), collapse = ""), "\n")

for (j in 1:nrow(qc_df)) {
  row <- qc_df[j, ]
  cat(sprintf("%-45s %6d %8d %6d %4d %4d %4d %6.1f%% %4d %8.1f %7.1f\n",
              row$sample_name,
              row$n_genes,
              row$n_cells,
              row$n_var_features,
              row$n_pcs_selected,
              row$n_pcs_elbow,
              row$n_pcs_threshold,
              row$cum_var_pct,
              row$n_clusters,
              row$time_seconds,
              row$rds_size_mb))
}

# 汇总统计
cat("\n--- 总计 ---\n")
cat("总样本数:", nrow(qc_df), "\n")
cat("总细胞数:", format(sum(qc_df$n_cells), big.mark = ","), "\n")
cat("总耗时:", round(sum(qc_df$time_seconds) / 60, 1), "分钟\n")
cat("总 .rds 大小:", round(sum(qc_df$rds_size_mb) / 1024, 2), "GB\n")

# 按数据集分组
cat("\n--- 按数据集分组 ---\n")
for (ds in unique(qc_df$dataset)) {
  sub <- qc_df[qc_df$dataset == ds, ]
  cat(sprintf("  %-30s %2d 样本, PCs=%d-%d, clusters=%d-%d, 耗时=%.0f秒\n",
              ds,
              nrow(sub),
              min(sub$n_pcs_selected), max(sub$n_pcs_selected),
              min(sub$n_clusters), max(sub$n_clusters),
              sum(sub$time_seconds)))
}

cat("\n============================================\n")
cat("R2 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
