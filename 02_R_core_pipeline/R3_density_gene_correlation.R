#!/usr/bin/env Rscript
# ============================================================================
# R3_density_gene_correlation.R
# 功能：对每个样本的每个基因，计算 SCT残差 vs 密度 的 Spearman 相关
# 输入：R2 输出的 .rds 文件（22个样本）
# 输出：每样本一个 density_gene_correlations.csv + 全局汇总
# Run (EN): Rscript R3_density_gene_correlation.R
#   Purpose: per gene per sample, Spearman correlation of SCT residual vs density.
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
R2_DIR         <- file.path(RESULTS_DIR, "R2_Results")
OUTPUT_DIR     <- file.path(RESULTS_DIR, "R3_Results")

# 密度列名（必须和P2输出完全一致）
DENSITY_COLS <- c(
  "density_knn_aggr_2nd_diff",
  "density_knn_main_piecewise",
  "density_knn_cons_max_dist",
  "density_voronoi",
  "density_delaunay"
)

# 简称（用于输出列名）
DENSITY_SHORT <- c("knn_aggr", "knn_main", "knn_cons", "voronoi", "delaunay")

# 三个KNN列的索引（用于收束判断）
KNN_INDICES <- 1:3

cat("============================================\n")
cat("R3: Density-Gene Spearman 相关分析\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 矩阵化 Spearman 相关函数
# ============================================================================

# 对稀疏矩阵的每一行和一个密度向量算 Spearman 相关
# 返回：data.frame(rho, pvalue) 每行一个基因
matrix_spearman <- function(residual_matrix, density_vec) {
  
  # 找到密度非NA的细胞
  valid <- !is.na(density_vec)
  n_valid <- sum(valid)
  
  if (n_valid < 10) {
    # 有效细胞数太少，返回全NA
    return(data.frame(
      rho = rep(NA_real_, nrow(residual_matrix)),
      pvalue = rep(NA_real_, nrow(residual_matrix)),
      n_cells = rep(n_valid, nrow(residual_matrix))
    ))
  }
  
  # 子集化：只保留密度非NA的细胞
  sub_matrix <- residual_matrix[, valid, drop = FALSE]
  sub_density <- density_vec[valid]
  
  # 对密度向量做排名（一次性）
  rank_density <- rank(sub_density, ties.method = "average")
  
  # 对残差矩阵逐行做排名
  # 稀疏矩阵转稠密后逐行rank
  n_genes <- nrow(sub_matrix)
  rho_vec <- numeric(n_genes)
  
  # 预计算密度排名的均值和标准差（所有基因共用）
  mean_d <- mean(rank_density)
  sd_d <- sd(rank_density)
  
  # 分批处理以控制内存（每批100个基因）
  batch_size <- 100
  n_batches <- ceiling(n_genes / batch_size)
  
  for (b in 1:n_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, n_genes)
    batch_indices <- start_idx:end_idx
    
    # 提取这批基因的残差（转为稠密矩阵）
    batch_dense <- as.matrix(sub_matrix[batch_indices, , drop = FALSE])
    
    # 逐行排名
    batch_ranks <- t(apply(batch_dense, 1, rank, ties.method = "average"))
    
    # Pearson(排名, 排名) = Spearman
    # rho = cor(rank_gene, rank_density)
    # 手动算避免循环调用cor()
    mean_g <- rowMeans(batch_ranks)
    
    # 中心化
    centered_g <- batch_ranks - mean_g  # 每行减去该行均值
    centered_d <- rank_density - mean_d  # 密度排名中心化
    
    # 协方差 = rowSums(centered_g * centered_d) / (n-1)
    # 标准差 = sqrt(rowSums(centered_g^2) / (n-1))
    cov_gd <- as.numeric(centered_g %*% centered_d) / (n_valid - 1)
    sd_g <- sqrt(rowSums(centered_g^2) / (n_valid - 1))
    
    rho_vec[batch_indices] <- cov_gd / (sd_g * sd_d)
  }
  
  # 处理sd=0的情况（某基因在所有细胞中残差完全相同）
  rho_vec[is.nan(rho_vec)] <- 0
  
  # 用t分布近似计算p值
  # t = rho * sqrt((n-2) / (1-rho^2))
  # p = 2 * pt(-abs(t), df=n-2)
  t_stat <- rho_vec * sqrt((n_valid - 2) / (1 - rho_vec^2))
  t_stat[is.nan(t_stat) | is.infinite(t_stat)] <- 0
  p_vec <- 2 * pt(-abs(t_stat), df = n_valid - 2)
  
  return(data.frame(
    rho = rho_vec,
    pvalue = p_vec,
    n_cells = rep(n_valid, n_genes)
  ))
}

# ============================================================================
# 收束标签函数
# ============================================================================

assign_convergence <- function(result_df) {
  # q值列名
  q_cols_knn <- paste0("q_", DENSITY_SHORT[KNN_INDICES])
  q_cols_all <- paste0("q_", DENSITY_SHORT)
  rho_col_main <- paste0("rho_", DENSITY_SHORT[2])  # knn_main
  
  labels <- character(nrow(result_df))
  
  for (i in 1:nrow(result_df)) {
    # 三个KNN的q值
    q_knn <- as.numeric(result_df[i, q_cols_knn])
    # 所有5种方法的q值
    q_all <- as.numeric(result_df[i, q_cols_all])
    
    # 处理NA（Voronoi/Delaunay可能有NA）
    knn_sig <- sum(q_knn < 0.05, na.rm = TRUE)
    knn_total <- sum(!is.na(q_knn))
    all_sig <- sum(q_all < 0.05, na.rm = TRUE)
    all_total <- sum(!is.na(q_all))
    
    if (knn_total == 0) {
      labels[i] <- "no_data"
    } else if (knn_sig == knn_total && all_sig == all_total) {
      # 所有方法都显著
      labels[i] <- "method_robust"
    } else if (knn_sig == knn_total) {
      # 三个KNN都显著，但Voronoi/Delaunay不全显著
      labels[i] <- "high_confidence"
    } else if (knn_sig == 0) {
      # 三个KNN都不显著
      labels[i] <- "not_significant"
    } else {
      # 部分KNN显著
      labels[i] <- "K_sensitive"
    }
  }
  
  return(labels)
}

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

summary_list <- list()

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
  
  rds_path <- file.path(R2_DIR, dataset_name, sample_subname,
                        paste0(sample_subname, "_seurat_R2.rds"))
  
  if (!file.exists(rds_path)) {
    cat("  ✗ .rds 文件不存在:", rds_path, "\n\n")
    next
  }
  
  # --- 创建输出目录 ---
  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------
  # 第一步：读取 Seurat 对象
  # ---------------------------------------------------------------
  
  cat("  [1/5] 读取 .rds ...\n")
  seurat_obj <- readRDS(rds_path)
  n_genes <- nrow(seurat_obj)
  n_cells <- ncol(seurat_obj)
  cat(sprintf("        %d 基因 × %d 细胞\n", n_genes, n_cells))
  
  # ---------------------------------------------------------------
  # 第二步：提取残差矩阵和密度值
  # ---------------------------------------------------------------
  
  cat("  [2/5] 提取 SCT 残差 + 密度值 ...\n")
  
  # SCT 残差矩阵（稀疏，基因×细胞）
  residuals <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
  gene_names <- rownames(residuals)
  
  # 密度值
  density_list <- list()
  for (d in seq_along(DENSITY_COLS)) {
    col_name <- DENSITY_COLS[d]
    if (col_name %in% colnames(seurat_obj@meta.data)) {
      density_list[[d]] <- seurat_obj@meta.data[[col_name]]
    } else {
      cat(sprintf("        ⚠ 密度列 %s 不存在, 跳过\n", col_name))
      density_list[[d]] <- rep(NA_real_, n_cells)
    }
  }
  
  cat("        ✓ 数据提取完成\n")
  
  # 释放Seurat对象节省内存
  rm(seurat_obj)
  gc(verbose = FALSE)
  
  # ---------------------------------------------------------------
  # 第三步：矩阵化 Spearman 相关
  # ---------------------------------------------------------------
  
  cat("  [3/5] 计算 Spearman 相关 ...\n")
  
  # 初始化结果 data.frame
  result_df <- data.frame(gene = gene_names, stringsAsFactors = FALSE)
  
  for (d in seq_along(DENSITY_COLS)) {
    short <- DENSITY_SHORT[d]
    cat(sprintf("        [%d/5] %s ...", d, short))
    
    spear_result <- matrix_spearman(residuals, density_list[[d]])
    
    result_df[[paste0("rho_", short)]] <- spear_result$rho
    result_df[[paste0("p_", short)]] <- spear_result$pvalue
    result_df[[paste0("n_cells_", short)]] <- spear_result$n_cells
    
    # 报告有效细胞数
    n_valid <- spear_result$n_cells[1]
    n_sig_raw <- sum(spear_result$pvalue < 0.05, na.rm = TRUE)
    cat(sprintf(" n=%d, raw_sig=%d\n", n_valid, n_sig_raw))
  }
  
  cat("        ✓ Spearman 相关完成\n")
  
  # ---------------------------------------------------------------
  # 第四步：FDR 校正 + 收束标签
  # ---------------------------------------------------------------
  
  cat("  [4/5] FDR 校正 + 收束标签 ...\n")
  
  # 对每种密度方法独立做 FDR 校正
  for (d in seq_along(DENSITY_SHORT)) {
    short <- DENSITY_SHORT[d]
    p_col <- paste0("p_", short)
    q_col <- paste0("q_", short)
    result_df[[q_col]] <- p.adjust(result_df[[p_col]], method = "BH")
  }
  
  # 双阈值标记
  for (d in seq_along(DENSITY_SHORT)) {
    short <- DENSITY_SHORT[d]
    q_col <- paste0("q_", short)
    sig_col <- paste0("sig_", short)
    q_vals <- result_df[[q_col]]
    result_df[[sig_col]] <- ifelse(is.na(q_vals), "NA",
                              ifelse(q_vals < 0.01, "highly_significant",
                                ifelse(q_vals < 0.05, "significant",
                                  "not_significant")))
  }
  
  # 收束标签
  result_df$convergence <- assign_convergence(result_df)
  
  # 按 knn_main 的 |rho| 降序排列
  result_df <- result_df[order(-abs(result_df$rho_knn_main)), ]
  
  cat("        ✓ FDR 校正完成\n")
  
  # ---------------------------------------------------------------
  # 第五步：输出 + 汇总
  # ---------------------------------------------------------------
  
  cat("  [5/5] 输出 ...\n")
  
  # 保存完整结果
  out_path <- file.path(out_dir, "density_gene_correlations.csv")
  fwrite(result_df, out_path)
  
  # 汇总统计
  q_main <- result_df$q_knn_main
  rho_main <- result_df$rho_knn_main
  
  n_sig_005 <- sum(q_main < 0.05, na.rm = TRUE)
  n_sig_001 <- sum(q_main < 0.01, na.rm = TRUE)
  n_pos <- sum(q_main < 0.05 & rho_main > 0, na.rm = TRUE)
  n_neg <- sum(q_main < 0.05 & rho_main < 0, na.rm = TRUE)
  
  # 收束统计
  conv_table <- table(result_df$convergence)
  n_method_robust <- ifelse("method_robust" %in% names(conv_table), conv_table["method_robust"], 0)
  n_high_conf <- ifelse("high_confidence" %in% names(conv_table), conv_table["high_confidence"], 0)
  n_k_sensitive <- ifelse("K_sensitive" %in% names(conv_table), conv_table["K_sensitive"], 0)
  n_not_sig <- ifelse("not_significant" %in% names(conv_table), conv_table["not_significant"], 0)
  
  # Top 5 基因（按|rho_knn_main|）
  top5 <- head(result_df, 5)
  top5_str <- paste(sprintf("%s(%.3f)", top5$gene, top5$rho_knn_main), collapse = ", ")
  
  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  
  cat(sprintf("        KNN_main: %d/%d 基因显著 (q<0.05), 其中正相关=%d, 负相关=%d\n",
              n_sig_005, n_genes, n_pos, n_neg))
  cat(sprintf("        收束: method_robust=%d, high_conf=%d, K_sensitive=%d, not_sig=%d\n",
              n_method_robust, n_high_conf, n_k_sensitive, n_not_sig))
  cat(sprintf("        Top5: %s\n", top5_str))
  cat(sprintf("  完成! 耗时 %.1f 秒\n\n", t_elapsed))
  
  # 收集到汇总列表
  summary_list[[sample_name]] <- data.frame(
    sample_name = sample_name,
    dataset = dataset_name,
    species = ifelse(!is.null(sample_info$species), sample_info$species, NA),
    condition = ifelse(!is.null(sample_info$condition), sample_info$condition, NA),
    n_genes = n_genes,
    n_cells = n_cells,
    n_sig_q005 = n_sig_005,
    n_sig_q001 = n_sig_001,
    n_pos_q005 = n_pos,
    n_neg_q005 = n_neg,
    pct_sig = round(n_sig_005 / n_genes * 100, 1),
    median_abs_rho = round(median(abs(rho_main), na.rm = TRUE), 4),
    max_abs_rho = round(max(abs(rho_main), na.rm = TRUE), 4),
    n_method_robust = as.integer(n_method_robust),
    n_high_confidence = as.integer(n_high_conf),
    n_K_sensitive = as.integer(n_k_sensitive),
    n_not_significant = as.integer(n_not_sig),
    top1_gene = top5$gene[1],
    top1_rho = round(top5$rho_knn_main[1], 4),
    time_seconds = t_elapsed,
    stringsAsFactors = FALSE
  )
  
  # 释放内存
  rm(residuals, density_list, result_df)
  gc(verbose = FALSE)
}

# ============================================================================
# 全局汇总
# ============================================================================

cat("============================================\n")
cat("全局汇总\n")
cat("============================================\n")

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

summary_path <- file.path(OUTPUT_DIR, "ALL_SAMPLES_R3_SUMMARY.csv")
fwrite(summary_df, summary_path)
cat("汇总已保存:", summary_path, "\n\n")

# 打印汇总表
cat("--- 22 样本 R3 汇总 ---\n")
cat(sprintf("%-45s %6s %8s %6s %6s %5s %5s %6s %6s %6s %7s\n",
            "sample", "genes", "cells", "sig05", "sig01", "pos", "neg", "pct%",
            "robust", "hiConf", "time_s"))
cat(paste(rep("-", 140), collapse = ""), "\n")

for (j in 1:nrow(summary_df)) {
  row <- summary_df[j, ]
  cat(sprintf("%-45s %6d %8d %6d %6d %5d %5d %5.1f%% %6d %6d %7.1f\n",
              row$sample_name,
              row$n_genes,
              row$n_cells,
              row$n_sig_q005,
              row$n_sig_q001,
              row$n_pos_q005,
              row$n_neg_q005,
              row$pct_sig,
              row$n_method_robust,
              row$n_high_confidence,
              row$time_seconds))
}

# 按数据集分组
cat("\n--- 按数据集分组 ---\n")
for (ds in unique(summary_df$dataset)) {
  sub <- summary_df[summary_df$dataset == ds, ]
  cat(sprintf("  %-30s %2d 样本, 显著基因(q<0.05): %d-%d (%.1f%%-%.1f%%), method_robust: %d-%d\n",
              ds,
              nrow(sub),
              min(sub$n_sig_q005), max(sub$n_sig_q005),
              min(sub$pct_sig), max(sub$pct_sig),
              min(sub$n_method_robust), max(sub$n_method_robust)))
}

# 收束统计全局
cat("\n--- 收束标签全局统计 ---\n")
cat(sprintf("  method_robust (所有方法显著):   %d-%d 基因/样本\n",
            min(summary_df$n_method_robust), max(summary_df$n_method_robust)))
cat(sprintf("  high_confidence (3个KNN显著):   %d-%d 基因/样本\n",
            min(summary_df$n_high_confidence), max(summary_df$n_high_confidence)))
cat(sprintf("  K_sensitive (部分KNN显著):      %d-%d 基因/样本\n",
            min(summary_df$n_K_sensitive), max(summary_df$n_K_sensitive)))
cat(sprintf("  not_significant:                %d-%d 基因/样本\n",
            min(summary_df$n_not_significant), max(summary_df$n_not_significant)))

cat("\n============================================\n")
cat("R3 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
