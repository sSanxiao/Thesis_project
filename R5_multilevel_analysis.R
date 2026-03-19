############################################################
# R5_multilevel_analysis.R
# 功能: 三种分析策略捕捉不同层级的 density-gene 关联
#   策略A: 全局 Spearman (R3已完成, 此处读取结果)
#   策略B: Density 分箱差异表达 (Wilcoxon + FC)
#   策略C: Cluster 级别 Spearman
#   + PCA dims / resolution / signature 权重验证 (代表样本)
#
# 输入: seurat_umap.rds (R2) + density_results_KNN.csv (R3)
#       + cell_density_three_methods.csv (P2)
# 输出: 每样本 binned_DE_results.csv, cluster_density_results.csv
#       全局  multilevel_gene_classification.csv
#       参数验证报告
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R5_multilevel_analysis.R > R5_run.log 2>&1 &
#   tail -f R5_run.log
############################################################

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"
DENSITY_ROOT <- "/home/disk/wangqilu/Density_Caculation/Results"

FDR_THRESHOLD   <- 0.05
COR_THRESHOLD   <- 0.05   # R4 选定
FC_THRESHOLD    <- 0.25   # log2 FC, 对应约 19% 表达差异
CLUSTER_MIN_CELLS <- 50   # cluster 内最少细胞数

# P2 路径映射 (与 R3 一致)
P2_DIR_MAP <- list(
  "Alzheimer_Mouse/Wild_13_4"    = "Alzheimer_Mouse/Wild_13_4",
  "Alzheimer_Mouse/Wild_5_7"     = "Alzheimer_Mouse/Wild_5_7",
  "Alzheimer_Mouse/Wild_2_5"     = "Alzheimer_Mouse/Wild_2_5",
  "Alzheimer_Mouse/TgCRND8_17_9" = "Alzheimer_Mouse/TgCRND8_17_9",
  "Alzheimer_Mouse/TgCRND8_5_7"  = "Alzheimer_Mouse/TgCRND8_5_7",
  "Alzheimer_Mouse/TgCRND8_2_5"  = "Alzheimer_Mouse/TgCRND8_2_5",
  "Brain_Human/Alz"              = "Brain_Human/Alz",
  "Brain_Human/Gilo"             = "Brain_Human/Glio",
  "Brain_Human/Healthy"          = "Brain_Human/Healthy",
  "Brain_Mouse/single"           = "Brain_Mouse/Normal",
  "ATRT_Human/28"                = "ATRT_Human/28",
  "ATRT_Human/29"                = "ATRT_Human/29",
  "ATRT_Human/30"                = "ATRT_Human/30",
  "ATRT_Human/31"                = "ATRT_Human/31",
  "ATRT_Human/32"                = "ATRT_Human/32",
  "ATRT_Human/33"                = "ATRT_Human/33",
  "ATRT_Human/34"                = "ATRT_Human/34",
  "Medulloblastoma_Human/GSM8840046_MB263" = "Medulloblastoma_Human/MB263",
  "Medulloblastoma_Human/GSM8840047_MB266" = "Medulloblastoma_Human/MB266",
  "Medulloblastoma_Human/GSM8840048_MB295" = "Medulloblastoma_Human/MB295",
  "Medulloblastoma_Human/GSM8840049_MB299" = "Medulloblastoma_Human/MB299"
)

# ===========================================================
# 样本注册表
# ===========================================================

sample_list <- list(
  list(project="Alzheimer_Mouse", sample="Wild_13_4",     species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_5_7",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_2_5",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_17_9",  species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_5_7",   species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_2_5",   species="mouse", condition="TgCRND8_AD"),
  list(project="Brain_Human",     sample="Alz",           species="human", condition="alzheimer"),
  list(project="Brain_Human",     sample="Gilo",          species="human", condition="glioblastoma"),
  list(project="Brain_Human",     sample="Healthy",       species="human", condition="healthy_brain"),
  list(project="Brain_Mouse",     sample="single",        species="mouse", condition="normal_brain"),
  list(project="ATRT_Human",      sample="28",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="29",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="30",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="31",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="32",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="33",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="34",            species="human", condition="ATRT"),
  list(project="Medulloblastoma_Human", sample="GSM8840046_MB263", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840047_MB266", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840048_MB295", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840049_MB299", species="human", condition="medulloblastoma")
)

# ===========================================================
# 向量化 Spearman
# ===========================================================

compute_spearman_vectorized <- function(expr_mat, density_vec) {
  expr_dense <- as.matrix(expr_mat)
  n          <- ncol(expr_dense)
  rank_d     <- rank(density_vec, ties.method = "average")
  rd_c       <- rank_d - mean(rank_d)
  ss_d       <- sum(rd_c^2)
  rank_e     <- t(apply(expr_dense, 1, rank, ties.method = "average"))
  re_c       <- rank_e - rowMeans(rank_e)
  ss_e       <- rowSums(re_c^2)
  cross      <- as.numeric(re_c %*% rd_c)
  denom      <- sqrt(ss_e * ss_d)
  denom[denom == 0] <- 1e-30
  rho        <- pmax(pmin(cross / denom, 1), -1)
  tstat      <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-10))
  pval       <- 2 * pt(-abs(tstat), df = n - 2)
  data.frame(
    gene         = rownames(expr_mat),
    spearman_cor = rho,
    p_value      = pval,
    stringsAsFactors = FALSE
  )
}

# ===========================================================
# 策略B: Density 分箱差异表达
# ===========================================================

compute_binned_DE <- function(expr_mat, density_vec) {
  # 四分位分箱: Q1(低密度) vs Q4(高密度)
  q25 <- quantile(density_vec, 0.25, na.rm = TRUE)
  q75 <- quantile(density_vec, 0.75, na.rm = TRUE)

  low_idx  <- which(density_vec <= q25)
  high_idx <- which(density_vec >= q75)

  if (length(low_idx) < 20 || length(high_idx) < 20) {
    return(NULL)
  }

  expr_dense <- as.matrix(expr_mat)
  results <- data.frame(
    gene     = rownames(expr_mat),
    mean_high = NA_real_,
    mean_low  = NA_real_,
    log2FC    = NA_real_,
    p_value   = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(expr_dense))) {
    high_vals <- expr_dense[i, high_idx]
    low_vals  <- expr_dense[i, low_idx]

    results$mean_high[i] <- mean(high_vals)
    results$mean_low[i]  <- mean(low_vals)

    # log2 FC with pseudocount
    results$log2FC[i] <- log2((mean(high_vals) + 0.01) / (mean(low_vals) + 0.01))

    # Wilcoxon 检验
    wt <- suppressWarnings(wilcox.test(high_vals, low_vals))
    results$p_value[i] <- wt$p.value
  }

  results$FDR <- p.adjust(results$p_value, method = "BH")
  results <- results[order(-abs(results$log2FC)), ]
  return(results)
}

# ===========================================================
# 主函数: 单样本多层级分析
# ===========================================================

process_one_sample <- function(s) {

  project <- s$project
  sname   <- s$sample
  r_label <- paste0(project, "/", sname)
  sample_dir <- file.path(RESULTS_ROOT, project, sname)

  rds_path <- file.path(sample_dir, "seurat_umap.rds")
  p2_label <- P2_DIR_MAP[[r_label]]
  density_csv <- file.path(DENSITY_ROOT, p2_label, "cell_density_three_methods.csv")

  # 检查输入
  if (!file.exists(rds_path)) { cat("  [SKIP] no seurat_umap.rds\n"); return(NULL) }
  if (!file.exists(density_csv)) { cat("  [SKIP] no density CSV\n"); return(NULL) }

  t0 <- Sys.time()

  # --- 加载数据 ---
  cat("  Loading data...\n")
  so <- readRDS(rds_path)
  dens_df <- fread(density_csv)
  dens_df$cell_id <- as.character(dens_df$cell_id)

  # bytes 修复
  if (nrow(dens_df) > 0 && grepl("^b'", dens_df$cell_id[1])) {
    dens_df$cell_id <- gsub("^b'|'$", "", dens_df$cell_id)
  }

  # cell_id 对齐
  common <- intersect(colnames(so), dens_df$cell_id)
  if (length(common) < 100) { cat("  [FAIL] <100 shared cells\n"); return(NULL) }

  so_sub   <- so[, common]
  dens_sub <- dens_df[match(common, dens_df$cell_id), ]
  d_knn    <- dens_sub$density_knn

  # SCT 表达矩阵
  expr_mat <- GetAssayData(so_sub, assay = "SCT", layer = "data")
  gvar     <- apply(expr_mat, 1, var)
  valid    <- gvar > 0
  expr_v   <- expr_mat[valid, , drop = FALSE]

  cat("  Cells:", length(common), " Valid genes:", sum(valid), "\n")

  # =========================================================
  # 策略A: 全局 Spearman (读取 R3 结果)
  # =========================================================
  global_file <- file.path(sample_dir, "density_results_KNN.csv")
  if (file.exists(global_file)) {
    global_res <- read.csv(global_file, stringsAsFactors = FALSE)
    cat("  Strategy A (global Spearman): loaded from R3\n")
  } else {
    cat("  Strategy A: computing...\n")
    global_res <- compute_spearman_vectorized(expr_v, d_knn)
    global_res$FDR <- p.adjust(global_res$p_value, method = "BH")
    global_res <- global_res[order(-global_res$spearman_cor), ]
  }

  # =========================================================
  # 策略B: 分箱差异表达 (Q1 vs Q4)
  # =========================================================
  cat("  Strategy B (binned DE, Q1 vs Q4)...\n")
  binned_res <- compute_binned_DE(expr_v, d_knn)

  if (!is.null(binned_res)) {
    write.csv(binned_res,
              file.path(sample_dir, "binned_DE_results.csv"),
              row.names = FALSE)

    n_de_sig <- sum(binned_res$FDR < FDR_THRESHOLD & abs(binned_res$log2FC) > FC_THRESHOLD,
                    na.rm = TRUE)
    cat("    Sig DE genes (FDR<0.05 & |log2FC|>", FC_THRESHOLD, "):", n_de_sig, "\n")

    # 强信号: |log2FC| > 1 (2倍差异)
    n_strong <- sum(binned_res$FDR < FDR_THRESHOLD & abs(binned_res$log2FC) > 1, na.rm = TRUE)
    cat("    Strong DE genes (|log2FC| > 1):", n_strong, "\n")
  } else {
    cat("    [WARN] Not enough cells for binning\n")
    n_de_sig <- 0
    n_strong <- 0
  }

  # =========================================================
  # 策略C: Cluster 级别 Spearman
  # =========================================================
  cat("  Strategy C (cluster-level Spearman)...\n")
  clusters <- Idents(so_sub)
  cl_results_all <- list()

  for (cl in unique(clusters)) {
    cells_cl <- names(clusters[clusters == cl])
    if (length(cells_cl) < CLUSTER_MIN_CELLS) next

    d_cl   <- d_knn[match(cells_cl, common)]
    e_cl   <- expr_v[, cells_cl, drop = FALSE]

    # 去除 cluster 内零方差基因
    cl_var <- apply(e_cl, 1, var)
    e_cl   <- e_cl[cl_var > 0, , drop = FALSE]
    if (nrow(e_cl) < 5) next

    cl_res <- compute_spearman_vectorized(e_cl, d_cl)
    cl_res$FDR     <- p.adjust(cl_res$p_value, method = "BH")
    cl_res$cluster <- cl
    cl_res$n_cells <- length(cells_cl)
    cl_results_all[[as.character(cl)]] <- cl_res
  }

  if (length(cl_results_all) > 0) {
    cl_all <- do.call(rbind, cl_results_all)
    write.csv(cl_all,
              file.path(sample_dir, "cluster_density_results.csv"),
              row.names = FALSE)

    # 找 cluster 级别的强信号
    cl_strong <- cl_all[cl_all$FDR < FDR_THRESHOLD & abs(cl_all$spearman_cor) > 0.15, ]
    n_cl_strong <- nrow(cl_strong)

    # 找在全局弱但 cluster 内强的基因 (Simpson 悖论恢复)
    if (nrow(cl_strong) > 0) {
      global_weak <- global_res$gene[is.na(global_res$FDR) |
                                       global_res$FDR >= FDR_THRESHOLD |
                                       abs(global_res$spearman_cor) <= COR_THRESHOLD]
      rescued <- cl_strong$gene[cl_strong$gene %in% global_weak]
      n_rescued <- length(unique(rescued))
    } else {
      n_rescued <- 0
    }

    cat("    Clusters analyzed:", length(cl_results_all), "\n")
    cat("    Strong cluster-level genes (|ρ|>0.15):", n_cl_strong, "\n")
    cat("    Rescued genes (weak global, strong in cluster):", n_rescued, "\n")
  } else {
    n_cl_strong <- 0
    n_rescued   <- 0
    cat("    [WARN] No clusters with enough cells\n")
  }

  # =========================================================
  # 基因分层分类: 强 / 中 / 弱
  # =========================================================
  all_genes <- unique(c(global_res$gene,
                        if (!is.null(binned_res)) binned_res$gene else character(0)))

  classify <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

  # 全局 Spearman
  idx_g <- match(classify$gene, global_res$gene)
  classify$global_rho <- global_res$spearman_cor[idx_g]
  classify$global_fdr <- global_res$FDR[idx_g]

  # 分箱 DE
  if (!is.null(binned_res)) {
    idx_b <- match(classify$gene, binned_res$gene)
    classify$binned_log2FC <- binned_res$log2FC[idx_b]
    classify$binned_fdr    <- binned_res$FDR[idx_b]
  } else {
    classify$binned_log2FC <- NA
    classify$binned_fdr    <- NA
  }

  # Cluster 最强信号
  if (length(cl_results_all) > 0) {
    cl_all_sorted <- cl_all[order(-abs(cl_all$spearman_cor)), ]
    cl_best <- cl_all_sorted[!duplicated(cl_all_sorted$gene), ]
    idx_c <- match(classify$gene, cl_best$gene)
    classify$cluster_best_rho     <- cl_best$spearman_cor[idx_c]
    classify$cluster_best_fdr     <- cl_best$FDR[idx_c]
    classify$cluster_best_cluster <- cl_best$cluster[idx_c]
  } else {
    classify$cluster_best_rho     <- NA
    classify$cluster_best_fdr     <- NA
    classify$cluster_best_cluster <- NA
  }

  # 分类逻辑
  classify$tier <- "none"

  # 强: 分箱 |log2FC| > 1 且 FDR<0.05,
  #     或 cluster 内 |ρ| > 0.3 且 FDR<0.05
  strong_binned  <- !is.na(classify$binned_fdr) &
                    classify$binned_fdr < FDR_THRESHOLD &
                    abs(classify$binned_log2FC) > 1
  strong_cluster <- !is.na(classify$cluster_best_fdr) &
                    classify$cluster_best_fdr < FDR_THRESHOLD &
                    abs(classify$cluster_best_rho) > 0.3
  classify$tier[strong_binned | strong_cluster] <- "strong"

  # 中: 分箱 |log2FC| > 0.25 且 FDR<0.05,
  #     或 cluster 内 |ρ| > 0.15,
  #     或 全局 |ρ| > 0.1
  #     (但不是 strong)
  moderate_binned  <- !is.na(classify$binned_fdr) &
                      classify$binned_fdr < FDR_THRESHOLD &
                      abs(classify$binned_log2FC) > FC_THRESHOLD
  moderate_cluster <- !is.na(classify$cluster_best_fdr) &
                      classify$cluster_best_fdr < FDR_THRESHOLD &
                      abs(classify$cluster_best_rho) > 0.15
  moderate_global  <- !is.na(classify$global_fdr) &
                      classify$global_fdr < FDR_THRESHOLD &
                      abs(classify$global_rho) > 0.1
  is_moderate <- (moderate_binned | moderate_cluster | moderate_global) &
                 classify$tier != "strong"
  classify$tier[is_moderate] <- "moderate"

  # 弱: 全局 FDR<0.05 且 |ρ| > 0.05 (但不是 strong 或 moderate)
  weak_global <- !is.na(classify$global_fdr) &
                 classify$global_fdr < FDR_THRESHOLD &
                 abs(classify$global_rho) > COR_THRESHOLD
  is_weak <- weak_global & classify$tier == "none"
  classify$tier[is_weak] <- "weak"

  classify <- classify[classify$tier != "none", ]
  classify <- classify[order(factor(classify$tier, levels = c("strong", "moderate", "weak")),
                             -abs(classify$global_rho)), ]

  write.csv(classify,
            file.path(sample_dir, "gene_tier_classification.csv"),
            row.names = FALSE)

  n_strong_tier   <- sum(classify$tier == "strong")
  n_moderate_tier <- sum(classify$tier == "moderate")
  n_weak_tier     <- sum(classify$tier == "weak")

  cat("  Gene tiers: strong=", n_strong_tier,
      " moderate=", n_moderate_tier,
      " weak=", n_weak_tier, "\n")

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  cat("  Done (", elapsed, " min)\n", sep = "")

  return(data.frame(
    sample         = r_label,
    species        = s$species,
    condition      = s$condition,
    n_cells        = length(common),
    n_genes_valid  = sum(valid),
    n_global_sig   = sum(global_res$FDR < FDR_THRESHOLD &
                           abs(global_res$spearman_cor) > COR_THRESHOLD, na.rm = TRUE),
    n_binned_sig   = n_de_sig,
    n_binned_strong = n_strong,
    n_cluster_strong = n_cl_strong,
    n_rescued      = n_rescued,
    tier_strong    = n_strong_tier,
    tier_moderate  = n_moderate_tier,
    tier_weak      = n_weak_tier,
    time_min       = elapsed,
    stringsAsFactors = FALSE
  ))
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R5: MULTI-LEVEL DENSITY-GENE ANALYSIS\n")
cat("  Strategy A: Global Spearman (from R3)\n")
cat("  Strategy B: Binned DE (Q1 vs Q4, Wilcoxon + log2FC)\n")
cat("  Strategy C: Cluster-level Spearman\n")
cat("  Samples:", length(sample_list), "\n")
cat("============================================================\n\n")

t_total <- Sys.time()
all_qc  <- list()

for (i in seq_along(sample_list)) {
  s     <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("\n[%d/%d] %s (%s, %s)\n", i, length(sample_list), label, s$species, s$condition))

  qc <- process_one_sample(s)
  if (!is.null(qc)) {
    all_qc[[length(all_qc) + 1]] <- qc
  }
}

# ===========================================================
# 全局汇总
# ===========================================================

qc_dir <- file.path(RESULTS_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

qc_all <- do.call(rbind, all_qc)
write.csv(qc_all,
          file.path(qc_dir, "R5_multilevel_qc.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("  TIER SUMMARY ACROSS ALL SAMPLES\n")
cat("============================================================\n\n")

cat("  Per-sample averages:\n")
cat("    Strong tier:  ", round(mean(qc_all$tier_strong), 1), " genes (mean)\n")
cat("    Moderate tier:", round(mean(qc_all$tier_moderate), 1), " genes (mean)\n")
cat("    Weak tier:    ", round(mean(qc_all$tier_weak), 1), " genes (mean)\n")
cat("    Rescued:      ", round(mean(qc_all$n_rescued), 1),
    " genes (weak globally, strong in cluster)\n\n")

# 按 condition 汇总
cond_summary <- qc_all %>%
  group_by(condition) %>%
  summarise(
    n_samples      = n(),
    mean_strong    = round(mean(tier_strong), 1),
    mean_moderate  = round(mean(tier_moderate), 1),
    mean_weak      = round(mean(tier_weak), 1),
    mean_rescued   = round(mean(n_rescued), 1),
    .groups = "drop"
  )

cat("  By condition:\n")
print(as.data.frame(cond_summary))

write.csv(as.data.frame(cond_summary),
          file.path(qc_dir, "R5_tier_by_condition.csv"),
          row.names = FALSE)

# ===========================================================
# 跨样本 tier 基因汇总: 哪些基因在多个样本中是 strong/moderate
# ===========================================================

cat("\n  Building cross-sample tier matrix...\n")

all_tier_files <- list()
for (s in sample_list) {
  f <- file.path(RESULTS_ROOT, s$project, s$sample, "gene_tier_classification.csv")
  if (file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$sample <- paste0(s$project, "/", s$sample)
    all_tier_files[[length(all_tier_files) + 1]] <- df
  }
}

if (length(all_tier_files) > 0) {
  tier_all <- do.call(rbind, all_tier_files)

  # 每个基因: 在多少样本中是 strong, moderate, weak
  gene_tier_counts <- tier_all %>%
    group_by(gene) %>%
    summarise(
      n_strong   = sum(tier == "strong"),
      n_moderate = sum(tier == "moderate"),
      n_weak     = sum(tier == "weak"),
      n_any      = n(),
      mean_global_rho = round(mean(global_rho, na.rm = TRUE), 4),
      mean_log2FC     = round(mean(binned_log2FC, na.rm = TRUE), 4),
      .groups = "drop"
    ) %>%
    arrange(desc(n_strong), desc(n_moderate), desc(n_any))

  write.csv(gene_tier_counts,
            file.path(qc_dir, "R5_cross_sample_gene_tiers.csv"),
            row.names = FALSE)

  n_cross_strong <- sum(gene_tier_counts$n_strong >= 3)
  n_cross_moderate <- sum(gene_tier_counts$n_moderate >= 3)

  cat("  Genes 'strong' in >=3 samples:", n_cross_strong, "\n")
  cat("  Genes 'moderate' in >=3 samples:", n_cross_moderate, "\n")
  cat("  Top 10 strongest genes:\n")
  print(head(gene_tier_counts, 10))
}

total_time <- round(difftime(Sys.time(), t_total, units = "mins"), 1)

cat("\n============================================================\n")
cat("  R5 COMPLETE\n")
cat("  Total time:", total_time, "min\n")
cat("============================================================\n")
