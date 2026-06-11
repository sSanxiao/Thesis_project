#!/usr/bin/env Rscript
# ============================================================================
# R09_tier_decision.R
# 功能: 基于 R3/R4/R7/R8 的结果 + 重新扫描 22 个 R2 .rds 的 signature AUC,
#       对 3 个 "局部升舱位 Tier I" 候选方向做可行性评估, 输出 Tier 判断报告
#
# 输入 (under RESULTS_DIR):
#   - R3_Results/*/density_gene_correlations.csv
#   - R4_Results/*/filtered_density_genes.csv  + ALL_SAMPLES_R4_SUMMARY.csv (如有)
#   - R7_Results/ALL_SAMPLES_R7_PROFILE.csv, ALL_DATASETS_R7_CONSISTENCY.csv
#   - R8_Results/ALL_COMPARISONS_R8_SUMMARY.csv (如有)
#   - R2_Results/*/*_seurat_R2.rds  (22 个)
#
# 输出 (under RESULTS_DIR/R9_Results):
#   - per_sample_signal_profile.csv, signature_auc_per_sample.csv
#   - reproducibility_summary.csv
#   - TIER_DECISION_REPORT.txt / TIER_DECISION_REPORT.json
#
# Run: Rscript R09_tier_decision.R   (paths via env vars; see config/paths.R)
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# null-coalesce 运算符 (R 没有原生的, 自定义一个)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ============================================================================
# 配置
# ============================================================================

DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
REGISTRY_PATH <- file.path(DATA_DIR, "sample_registry.json")
R2_DIR  <- file.path(RESULTS_DIR, "R2_Results")
R3_DIR  <- file.path(RESULTS_DIR, "R3_Results")
R4_DIR  <- file.path(RESULTS_DIR, "R4_Results")
R7_DIR  <- file.path(RESULTS_DIR, "R7_Results")
R8_DIR  <- file.path(RESULTS_DIR, "R8_Results")
OUT_DIR <- file.path(RESULTS_DIR, "R9_Results")

DENSITY_COL <- "density_knn_main_piecewise"

# 密度高/低组定义: 上下 20%
DENSITY_PCTL_LOW  <- 0.20
DENSITY_PCTL_HIGH <- 0.80

# AUC 判定阈值
AUC_STRONG <- 0.75
AUC_MODERATE <- 0.65

# Tier I 多个标准的阈值
T1_SINGLE_ABS_RHO <- 0.30         # 单样本 max|ρ| 达标线
T1_MEDIAN_CROSS <- 0.30           # 跨数据集 median ρ 达标线
T1_CONSISTENCY_WITHIN <- 0.40     # 数据集内一致率达标线

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat("R9: Tier 判断报告生成\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 读取 registry
# ============================================================================

registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
n_samples <- length(sample_names)
cat("样本数:", n_samples, "\n\n")

parse_sample <- function(sname) {
  parts <- strsplit(sname, "/")[[1]]
  list(dataset = parts[1], subname = parts[2])
}

# ============================================================================
# Part A: 信号强度体检 (从 R3 结果)
# ============================================================================

cat("============================================\n")
cat("Part A: 每样本信号强度体检 (从 R3 相关系数表)\n")
cat("============================================\n")

partA_rows <- list()

for (i in seq_along(sample_names)) {
  sname <- sample_names[i]
  ps <- parse_sample(sname)
  
  r3_path <- file.path(R3_DIR, ps$dataset, ps$subname, "density_gene_correlations.csv")
  if (!file.exists(r3_path)) {
    cat(sprintf("  [%d/%d] %s : R3 文件不存在, 跳过\n", i, n_samples, sname))
    next
  }
  
  r3 <- fread(r3_path)
  # 关键列: gene, rho_knn_main (或类似), qvalue_knn_main, convergence_label
  # 为了兼容不同的列命名, 自动探测
  rho_col <- grep("^rho.*knn.*main|^rho_knn_main", names(r3), value = TRUE)[1]
  q_col   <- grep("^q.*knn.*main|^qvalue_knn_main|^qval_knn_main", names(r3), value = TRUE)[1]
  conv_col <- grep("convergence|conv_label|^label", names(r3), value = TRUE)[1]
  
  if (is.na(rho_col)) {
    # 退而求其次: 用第一个 rho 列
    rho_col <- grep("^rho", names(r3), value = TRUE)[1]
    q_col   <- grep("^q", names(r3), value = TRUE)[1]
  }
  
  rhos <- r3[[rho_col]]
  qs <- if (!is.na(q_col)) r3[[q_col]] else rep(NA_real_, length(rhos))
  
  abs_rhos <- abs(rhos)
  n_total <- sum(!is.na(rhos))
  
  row <- data.frame(
    sample = sname,
    dataset = ps$dataset,
    n_genes = n_total,
    median_abs_rho = round(median(abs_rhos, na.rm = TRUE), 4),
    p90_abs_rho = round(quantile(abs_rhos, 0.90, na.rm = TRUE), 4),
    max_abs_rho = round(max(abs_rhos, na.rm = TRUE), 4),
    n_abs_rho_gt_0_1 = sum(abs_rhos >= 0.10, na.rm = TRUE),
    n_abs_rho_gt_0_2 = sum(abs_rhos >= 0.20, na.rm = TRUE),
    n_abs_rho_gt_0_3 = sum(abs_rhos >= 0.30, na.rm = TRUE),
    pct_abs_rho_gt_0_2 = round(100 * sum(abs_rhos >= 0.20, na.rm = TRUE) / n_total, 2),
    n_q_sig = if (!all(is.na(qs))) sum(qs < 0.05, na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
  partA_rows[[sname]] <- row
  
  if (i %% 5 == 0 || i == n_samples) {
    cat(sprintf("  [%d/%d] 完成\n", i, n_samples))
  }
}

partA_df <- rbindlist(partA_rows, fill = TRUE)
fwrite(partA_df, file.path(OUT_DIR, "per_sample_signal_profile.csv"))
cat(sprintf("\n  per_sample_signal_profile.csv 已写入 (%d 行)\n\n", nrow(partA_df)))

# ============================================================================
# Part B: 稳健性 (重新汇总 R7/R8)
# ============================================================================

cat("============================================\n")
cat("Part B: 跨样本/跨数据集稳健性\n")
cat("============================================\n")

# B1 - 数据集内一致性 (R7)
r7_consist_path <- file.path(R7_DIR, "ALL_DATASETS_R7_CONSISTENCY.csv")
if (file.exists(r7_consist_path)) {
  r7_consist <- fread(r7_consist_path)
  cat("  R7 数据集内一致性表已读取\n")
  print(r7_consist)
} else {
  cat("  ⚠ R7 一致性表不存在\n")
  r7_consist <- data.table()
}

# B2 - 跨数据集可重复性 (R8)
r8_summary_path <- file.path(R8_DIR, "ALL_COMPARISONS_R8_SUMMARY.csv")
if (file.exists(r8_summary_path)) {
  r8_summary <- fread(r8_summary_path)
  cat("\n  R8 跨数据集汇总表已读取\n")
} else {
  cat("\n  ⚠ R8 汇总表不存在\n")
  r8_summary <- data.table()
}

# 把 R8 分成同物种 vs 跨物种两组
if (nrow(r8_summary) > 0) {
  # 兼容列名
  type_col <- grep("type", names(r8_summary), value = TRUE)[1]
  rho_col <- grep("^rank.*rho|^rho|^spearman", names(r8_summary), value = TRUE, ignore.case = TRUE)[1]
  dir_col <- grep("^dir|direction|agree", names(r8_summary), value = TRUE, ignore.case = TRUE)[1]
  
  if (!is.na(type_col) && !is.na(rho_col)) {
    within_mask <- r8_summary[[type_col]] %in% c("within_species", "same", "same_species")
    cross_mask  <- r8_summary[[type_col]] %in% c("cross_species", "cross")
    
    b2_summary <- data.frame(
      comparison_type = c("within_species", "cross_species"),
      n_pairs = c(sum(within_mask), sum(cross_mask)),
      median_rho = c(
        round(median(r8_summary[[rho_col]][within_mask], na.rm = TRUE), 4),
        round(median(r8_summary[[rho_col]][cross_mask], na.rm = TRUE), 4)
      ),
      max_rho = c(
        round(max(r8_summary[[rho_col]][within_mask], na.rm = TRUE), 4),
        round(max(r8_summary[[rho_col]][cross_mask], na.rm = TRUE), 4)
      ),
      n_pairs_rho_gt_0_3 = c(
        sum(r8_summary[[rho_col]][within_mask] >= 0.30, na.rm = TRUE),
        sum(r8_summary[[rho_col]][cross_mask] >= 0.30, na.rm = TRUE)
      ),
      median_direction_pct = if (!is.na(dir_col)) c(
        round(median(r8_summary[[dir_col]][within_mask], na.rm = TRUE), 1),
        round(median(r8_summary[[dir_col]][cross_mask], na.rm = TRUE), 1)
      ) else c(NA, NA),
      stringsAsFactors = FALSE
    )
    fwrite(b2_summary, file.path(OUT_DIR, "reproducibility_summary.csv"))
    cat("\n  reproducibility_summary.csv:\n")
    print(b2_summary)
  }
}

cat("\n")

# ============================================================================
# Part C: Signature AUC 扫描 (全量 22 样本)
# ============================================================================

cat("============================================\n")
cat("Part C: Signature AUC 扫描 (22 个 .rds)\n")
cat("============================================\n\n")

# C0: 定义 3 个候选 signature

# C0.1 - Mouse maturation signature (Alzheimer_Mouse ∩ Brain_Mouse tier1)
MOUSE_MATURATION_SIG <- c(
  "Olig2", "Cnp", "Sox10", "Bdnf", "Epha4",
  "Nrn1", "Nrep", "Neurod6", "Trpc4"
)

# C0.2 - Cross-species conserved signature (R8 landscape >= 2 datasets)
# 用大写存, 实际读数据时会两种 case 都试
CROSS_SPECIES_SIG <- c(
  "APOE", "AQP4", "OLIG2", "CNP", "SOX10", "BDNF",
  "C1QB", "EPHA4", "FN1", "TRPC4", "NRN1", "NREP",
  "NEUROD6", "CRYM", "SLC17A7", "RBFOX3", "SOX9"
)

# C0.3 - MB rhombic lip candidate signature
# 从 R4 的 4 个 MB 样本里自动筛 "至少 2 个样本 tier1" 的基因
mb_samples <- sample_names[startsWith(sample_names, "Medulloblastoma_Human")]
cat("  从 R4 提取 MB rhombic lip candidate signature...\n")
mb_tier1_lists <- list()
for (mb_s in mb_samples) {
  ps <- parse_sample(mb_s)
  r4_path <- file.path(R4_DIR, ps$dataset, ps$subname, "filtered_density_genes.csv")
  if (file.exists(r4_path)) {
    r4 <- fread(r4_path)
    # tier 列可能叫 tier / tier_label
    tier_col <- grep("^tier$|tier_label", names(r4), value = TRUE)[1]
    gene_col <- grep("^gene$", names(r4), value = TRUE)[1]
    if (!is.na(tier_col) && !is.na(gene_col)) {
      t1_genes <- r4[[gene_col]][r4[[tier_col]] == "tier1" | r4[[tier_col]] == "1"]
      mb_tier1_lists[[mb_s]] <- unique(t1_genes)
    }
  }
}
if (length(mb_tier1_lists) > 0) {
  all_mb_t1 <- unlist(mb_tier1_lists)
  tbl <- table(all_mb_t1)
  MB_RL_SIG <- names(tbl)[tbl >= 2]
  # 附加先验 rhombic lip / SHH / Group3-4 标志
  MB_RL_PRIOR <- c("HHIP", "CXCR4", "OTX2", "GLI2", "PTCH1", "MYC",
                   "EOMES", "SLC17A7", "TUBB4A", "RBFOX3", "NRXN3",
                   "SV2B", "DCN", "CCR7", "BOC")
  MB_RL_SIG <- unique(c(MB_RL_SIG, MB_RL_PRIOR))
  cat(sprintf("  MB signature 基因数: %d (自动筛 %d + 先验 %d, 去重)\n\n",
              length(MB_RL_SIG), sum(tbl >= 2), length(MB_RL_PRIOR)))
} else {
  MB_RL_SIG <- c("HHIP", "CXCR4", "OTX2", "GLI2", "PTCH1", "MYC",
                 "EOMES", "SLC17A7", "TUBB4A", "RBFOX3", "NRXN3")
  cat("  ⚠ 未能从 R4 读取 MB tier1 基因, 使用先验基因集\n\n")
}

# 打印 3 个 signature
cat("  Signature 1 (Mouse maturation):", length(MOUSE_MATURATION_SIG), "genes\n")
cat("    ", paste(MOUSE_MATURATION_SIG, collapse = ", "), "\n\n")
cat("  Signature 2 (Cross-species conserved):", length(CROSS_SPECIES_SIG), "genes\n")
cat("    ", paste(CROSS_SPECIES_SIG, collapse = ", "), "\n\n")
cat("  Signature 3 (MB rhombic lip):", length(MB_RL_SIG), "genes\n")
cat("    ", paste(MB_RL_SIG, collapse = ", "), "\n\n")

# C1: 定义辅助函数

# 简单 AUC (Mann-Whitney U, 不依赖 pROC)
fast_auc <- function(scores, labels) {
  # labels: 0/1, 1 = "high density group"
  pos <- scores[labels == 1]
  neg <- scores[labels == 0]
  n_pos <- length(pos)
  n_neg <- length(neg)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  # 合并排名
  all_scores <- c(pos, neg)
  r <- rank(all_scores)
  sum_rank_pos <- sum(r[1:n_pos])
  auc <- (sum_rank_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  return(auc)
}

# 给定一个 Seurat 对象和一个基因集, 计算 signature score = 该基因集 SCT 数据的按细胞平均
compute_signature_score <- function(obj, gene_set) {
  # 尝试 SCT, 如果没有就用 RNA
  assay_use <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"
  
  # SCT layer 可能叫 data / scale.data. 为了不被负残差干扰, 用 "data" (log-normalized)
  mat <- tryCatch(
    GetAssayData(obj, assay = assay_use, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay_use, slot = "data")
  )
  
  # 基因名大小写兼容: mouse 用驼峰, human 用大写
  species <- obj@meta.data$species[1]
  if (!is.null(species) && species == "mouse") {
    # 小鼠 panel: 驼峰
    genes_try <- unique(c(
      gene_set,
      paste0(substr(gene_set, 1, 1), tolower(substr(gene_set, 2, nchar(gene_set))))
    ))
  } else {
    # 人 panel: 大写
    genes_try <- unique(c(gene_set, toupper(gene_set)))
  }
  
  found <- intersect(genes_try, rownames(mat))
  if (length(found) < 3) {
    return(list(score = NULL, n_found = length(found), found = found))
  }
  
  sub <- mat[found, , drop = FALSE]
  # 按基因 z-score 再取平均 (避免高表达基因主导)
  sub_z <- t(scale(t(as.matrix(sub))))
  sub_z[is.na(sub_z)] <- 0
  score <- colMeans(sub_z)
  
  return(list(score = score, n_found = length(found), found = found))
}

# 评估一个 signature 在一个样本上的 AUC
evaluate_signature <- function(obj, gene_set, density_col,
                               pctl_low = DENSITY_PCTL_LOW,
                               pctl_high = DENSITY_PCTL_HIGH) {
  meta <- obj@meta.data
  if (!(density_col %in% names(meta))) {
    return(list(auc = NA, n_found = NA, direction = NA, n_low = NA, n_high = NA))
  }
  
  d <- meta[[density_col]]
  valid <- !is.na(d)
  
  sig_res <- compute_signature_score(obj, gene_set)
  if (is.null(sig_res$score)) {
    return(list(auc = NA, n_found = sig_res$n_found, direction = NA,
                n_low = NA, n_high = NA))
  }
  
  score <- sig_res$score
  # 对齐 valid 细胞
  score <- score[valid]
  d <- d[valid]
  
  q_low <- quantile(d, pctl_low, na.rm = TRUE)
  q_high <- quantile(d, pctl_high, na.rm = TRUE)
  
  low_mask <- d <= q_low
  high_mask <- d >= q_high
  
  if (sum(low_mask) < 50 || sum(high_mask) < 50) {
    return(list(auc = NA, n_found = sig_res$n_found, direction = NA,
                n_low = sum(low_mask), n_high = sum(high_mask)))
  }
  
  # AUC: "signature 区分 high vs low density 的能力"
  # 1 = high-density group (正类)
  labels <- c(rep(1, sum(high_mask)), rep(0, sum(low_mask)))
  scores <- c(score[high_mask], score[low_mask])
  
  auc <- fast_auc(scores, labels)
  # 方向: AUC > 0.5 说明 signature 在 high-density 中评分更高
  direction <- if (is.na(auc)) NA else if (auc >= 0.5) "positive" else "negative"
  # 统一到 "单调强度" (0.5-1 区间)
  auc_adj <- if (is.na(auc)) NA else max(auc, 1 - auc)
  
  return(list(
    auc_raw = auc,
    auc = auc_adj,
    direction = direction,
    n_found = sig_res$n_found,
    genes_found = paste(sig_res$found, collapse = ";"),
    n_low = sum(low_mask),
    n_high = sum(high_mask)
  ))
}

# C2: 逐样本扫描

partC_rows <- list()

# 断点续跑
partC_cache <- file.path(OUT_DIR, "signature_auc_per_sample.csv")
done_samples <- c()
if (file.exists(partC_cache)) {
  cached <- fread(partC_cache)
  done_samples <- unique(cached$sample)
  partC_rows <- split(cached, cached$sample)
  cat(sprintf("  发现缓存: %d 样本已完成, 将跳过\n", length(done_samples)))
}

for (i in seq_along(sample_names)) {
  sname <- sample_names[i]
  ps <- parse_sample(sname)
  
  if (sname %in% done_samples) {
    cat(sprintf("  [%d/%d] %s : 已缓存, 跳过\n", i, n_samples, sname))
    next
  }
  
  rds_path <- file.path(R2_DIR, ps$dataset, ps$subname,
                        paste0(ps$subname, "_seurat_R2.rds"))
  if (!file.exists(rds_path)) {
    cat(sprintf("  [%d/%d] %s : .rds 不存在, 跳过\n", i, n_samples, sname))
    next
  }
  
  cat(sprintf("  [%d/%d] %s ", i, n_samples, sname))
  t0 <- Sys.time()
  
  obj <- tryCatch(
    readRDS(rds_path),
    error = function(e) { cat(" 读取失败\n"); NULL }
  )
  if (is.null(obj)) next
  
  cat(sprintf("(%d cells, %.1fs 读取) ",
              ncol(obj), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  # 3 个 signature 逐个测
  res_mouse <- evaluate_signature(obj, MOUSE_MATURATION_SIG, DENSITY_COL)
  res_cross <- evaluate_signature(obj, CROSS_SPECIES_SIG, DENSITY_COL)
  res_mb    <- evaluate_signature(obj, MB_RL_SIG, DENSITY_COL)
  
  row <- data.frame(
    sample = sname,
    dataset = ps$dataset,
    species = registry[[sname]]$species,
    n_cells = ncol(obj),
    
    mouse_sig_n_found = res_mouse$n_found,
    mouse_sig_auc = round(res_mouse$auc %||% NA, 4),
    mouse_sig_direction = res_mouse$direction %||% NA,
    
    cross_sig_n_found = res_cross$n_found,
    cross_sig_auc = round(res_cross$auc %||% NA, 4),
    cross_sig_direction = res_cross$direction %||% NA,
    
    mb_sig_n_found = res_mb$n_found,
    mb_sig_auc = round(res_mb$auc %||% NA, 4),
    mb_sig_direction = res_mb$direction %||% NA,
    
    n_low = res_cross$n_low, n_high = res_cross$n_high,
    elapsed_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    stringsAsFactors = FALSE
  )
  partC_rows[[sname]] <- row
  
  cat(sprintf("| Mouse AUC=%.3f, Cross AUC=%.3f, MB AUC=%.3f (%.1fs)\n",
              row$mouse_sig_auc, row$cross_sig_auc, row$mb_sig_auc, row$elapsed_s))
  
  # 增量保存 (断点续跑用)
  current_df <- rbindlist(partC_rows, fill = TRUE)
  fwrite(current_df, partC_cache)
  
  rm(obj); gc(verbose = FALSE)
}

# (%||% 已在文件开头定义)

partC_df <- rbindlist(partC_rows, fill = TRUE)
fwrite(partC_df, partC_cache)
cat(sprintf("\n  signature_auc_per_sample.csv 写入 (%d 行)\n\n", nrow(partC_df)))

# ============================================================================
# Part D: 综合判断报告
# ============================================================================

cat("============================================\n")
cat("Part D: 生成 Tier 判断报告\n")
cat("============================================\n\n")

# 从 partC 汇总各 signature 的 AUC 表现
sig_summary_cross <- function(col_auc, col_dir, group_filter = NULL) {
  df <- partC_df
  if (!is.null(group_filter)) df <- df[eval(group_filter, df), ]
  vals <- df[[col_auc]]
  vals <- vals[!is.na(vals)]
  list(
    n = length(vals),
    median = round(median(vals), 4),
    max = round(max(vals), 4),
    min = round(min(vals), 4),
    n_strong = sum(vals >= AUC_STRONG),
    n_moderate = sum(vals >= AUC_MODERATE & vals < AUC_STRONG),
    n_weak = sum(vals < AUC_MODERATE),
    direction_consistent = {
      dirs <- df[[col_dir]][!is.na(df[[col_auc]])]
      pos <- sum(dirs == "positive", na.rm = TRUE)
      neg <- sum(dirs == "negative", na.rm = TRUE)
      round(max(pos, neg) / length(dirs) * 100, 1)
    }
  )
}

mouse_sum <- sig_summary_cross("mouse_sig_auc", "mouse_sig_direction",
                               quote(species == "mouse"))
cross_sum <- sig_summary_cross("cross_sig_auc", "cross_sig_direction")
mb_sum    <- sig_summary_cross("mb_sig_auc", "mb_sig_direction",
                               quote(dataset == "Medulloblastoma_Human"))

# 信号强度全局统计 (partA)
signal_summary <- list(
  n_samples_max_rho_gt_0_3 = sum(partA_df$max_abs_rho >= T1_SINGLE_ABS_RHO, na.rm = TRUE),
  n_samples_max_rho_gt_0_5 = sum(partA_df$max_abs_rho >= 0.50, na.rm = TRUE),
  median_max_abs_rho = round(median(partA_df$max_abs_rho, na.rm = TRUE), 4),
  median_median_abs_rho = round(median(partA_df$median_abs_rho, na.rm = TRUE), 4)
)

# 升舱判断
decide_upgrade <- function(sum_obj, label) {
  n_assessed <- sum_obj$n
  if (n_assessed == 0) return(list(verdict = "NO-DATA", reason = "No samples evaluated"))
  
  frac_strong <- sum_obj$n_strong / n_assessed
  frac_at_least_mod <- (sum_obj$n_strong + sum_obj$n_moderate) / n_assessed
  
  if (sum_obj$median >= 0.70 && frac_strong >= 0.5 &&
      sum_obj$direction_consistent >= 80) {
    return(list(
      verdict = "GO (strong Tier I candidate)",
      reason = sprintf("median AUC=%.3f, %d/%d samples AUC>=%.2f, direction %.1f%% consistent",
                       sum_obj$median, sum_obj$n_strong, n_assessed,
                       AUC_STRONG, sum_obj$direction_consistent)
    ))
  }
  
  if (sum_obj$median >= 0.65 && frac_at_least_mod >= 0.5 &&
      sum_obj$direction_consistent >= 70) {
    return(list(
      verdict = "MAYBE (worth pilot external validation)",
      reason = sprintf("median AUC=%.3f, %d/%d samples AUC>=%.2f, direction %.1f%% consistent",
                       sum_obj$median,
                       sum_obj$n_strong + sum_obj$n_moderate, n_assessed,
                       AUC_MODERATE, sum_obj$direction_consistent)
    ))
  }
  
  return(list(
    verdict = "NO-GO (stay in Tier II)",
    reason = sprintf("median AUC=%.3f too low, only %d/%d samples AUC>=%.2f",
                     sum_obj$median, sum_obj$n_moderate + sum_obj$n_strong,
                     n_assessed, AUC_MODERATE)
  ))
}

verdict_mouse <- decide_upgrade(mouse_sum, "Mouse maturation")
verdict_cross <- decide_upgrade(cross_sum, "Cross-species conserved")
verdict_mb    <- decide_upgrade(mb_sum, "MB rhombic lip")

# 全局 Tier 判断
overall_tier <- if (
  verdict_mouse$verdict == "GO (strong Tier I candidate)" ||
  verdict_cross$verdict == "GO (strong Tier I candidate)" ||
  verdict_mb$verdict == "GO (strong Tier I candidate)"
) {
  "TIER I local (overall Tier II)"
} else if (
  verdict_mouse$verdict == "MAYBE (worth pilot external validation)" ||
  verdict_cross$verdict == "MAYBE (worth pilot external validation)" ||
  verdict_mb$verdict == "MAYBE (worth pilot external validation)"
) {
  "TIER II (with upgrade potential pending external data)"
} else {
  "TIER II (stable) or TIER III"
}

# 写纯文本报告
report_path <- file.path(OUT_DIR, "TIER_DECISION_REPORT.txt")
sink(report_path)

cat("================================================================\n")
cat("  TIER DECISION REPORT — Master Thesis: Spatial Density Genes\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("OVERALL VERDICT:", overall_tier, "\n\n")
cat("----------------------------------------------------------------\n")
cat("Part A: 信号强度体检 (来自 R3)\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  样本总数: %d\n", nrow(partA_df)))
cat(sprintf("  samples with max|ρ| >= 0.30: %d\n", signal_summary$n_samples_max_rho_gt_0_3))
cat(sprintf("  samples with max|ρ| >= 0.50: %d\n", signal_summary$n_samples_max_rho_gt_0_5))
cat(sprintf("  median(max|ρ|) across samples: %.4f\n", signal_summary$median_max_abs_rho))
cat(sprintf("  median(median|ρ|) across samples: %.4f\n", signal_summary$median_median_abs_rho))

cat("\n  Per-sample breakdown:\n")
printf <- function(...) cat(sprintf(...))
printf("  %-45s %8s %8s %8s %8s\n",
       "sample", "med|ρ|", "p90|ρ|", "max|ρ|", "n|ρ|>.2")
printf("  %s\n", paste(rep("-", 80), collapse = ""))
for (k in seq_len(nrow(partA_df))) {
  printf("  %-45s %8.4f %8.4f %8.4f %8d\n",
         partA_df$sample[k],
         partA_df$median_abs_rho[k],
         partA_df$p90_abs_rho[k],
         partA_df$max_abs_rho[k],
         partA_df$n_abs_rho_gt_0_2[k])
}

cat("\n----------------------------------------------------------------\n")
cat("Part B: 稳健性 (R7 + R8)\n")
cat("----------------------------------------------------------------\n")

if (nrow(r7_consist) > 0) {
  cat("\n  数据集内 tier1 一致性 (R7):\n")
  print(r7_consist)
}

if (exists("b2_summary")) {
  cat("\n  跨数据集可重复性 (R8, 按物种类型):\n")
  print(b2_summary)
}

cat("\n----------------------------------------------------------------\n")
cat("Part C: Signature AUC 扫描 (3 个候选方向)\n")
cat("----------------------------------------------------------------\n")

print_sig <- function(label, sum_obj, verdict, gene_set) {
  cat("\n  ==", label, "==\n")
  cat(sprintf("  Gene set size: %d\n", length(gene_set)))
  cat(sprintf("  Samples evaluated: %d\n", sum_obj$n))
  cat(sprintf("  AUC: median=%.4f, range=[%.4f, %.4f]\n",
              sum_obj$median, sum_obj$min, sum_obj$max))
  cat(sprintf("  Strong (AUC>=%.2f): %d | Moderate (%.2f-%.2f): %d | Weak (<%.2f): %d\n",
              AUC_STRONG, sum_obj$n_strong,
              AUC_MODERATE, AUC_STRONG, sum_obj$n_moderate,
              AUC_MODERATE, sum_obj$n_weak))
  cat(sprintf("  Direction consistency: %.1f%%\n", sum_obj$direction_consistent))
  cat(sprintf("  VERDICT: %s\n", verdict$verdict))
  cat(sprintf("    reason: %s\n", verdict$reason))
}

print_sig("Signature 1 — Mouse maturation (evaluated on mouse samples)",
          mouse_sum, verdict_mouse, MOUSE_MATURATION_SIG)
print_sig("Signature 2 — Cross-species conserved (evaluated on all samples)",
          cross_sum, verdict_cross, CROSS_SPECIES_SIG)
print_sig("Signature 3 — MB rhombic lip (evaluated on MB samples)",
          mb_sum, verdict_mb, MB_RL_SIG)

cat("\n----------------------------------------------------------------\n")
cat("Part D: Per-sample signature AUC table\n")
cat("----------------------------------------------------------------\n\n")

printf("  %-45s %4s %6s %4s %6s %4s %6s\n",
       "sample", "MsN", "MsAUC", "CsN", "CsAUC", "MBN", "MBAUC")
printf("  %s\n", paste(rep("-", 90), collapse = ""))
for (k in seq_len(nrow(partC_df))) {
  printf("  %-45s %4d %6.3f %4d %6.3f %4d %6.3f\n",
         partC_df$sample[k],
         partC_df$mouse_sig_n_found[k] %||% 0,
         partC_df$mouse_sig_auc[k] %||% NA,
         partC_df$cross_sig_n_found[k] %||% 0,
         partC_df$cross_sig_auc[k] %||% NA,
         partC_df$mb_sig_n_found[k] %||% 0,
         partC_df$mb_sig_auc[k] %||% NA)
}

cat("\n================================================================\n")
cat("INTERPRETATION GUIDE (paste to Claude):\n")
cat("================================================================\n")
cat("1. If ANY signature got 'GO' verdict  -> project has Tier I potential\n")
cat("   in that specific axis, should build predictive model against\n")
cat("   external data (TCGA / scRNA-seq development atlas).\n\n")
cat("2. If 'MAYBE' is the best verdict     -> Tier II main conclusion,\n")
cat("   with one targeted external validation effort worth trying.\n\n")
cat("3. If all are 'NO-GO'                  -> Lock in Tier II (gene-set +\n")
cat("   pathway description) or Tier III (case-study approach).\n")
cat("================================================================\n")

sink()

# JSON 版本
json_out <- list(
  generated_at = format(Sys.time()),
  overall_tier = overall_tier,
  partA_signal_summary = signal_summary,
  partC_signatures = list(
    mouse_maturation = list(stats = mouse_sum, verdict = verdict_mouse,
                            genes = MOUSE_MATURATION_SIG),
    cross_species = list(stats = cross_sum, verdict = verdict_cross,
                         genes = CROSS_SPECIES_SIG),
    mb_rhombic_lip = list(stats = mb_sum, verdict = verdict_mb,
                          genes = MB_RL_SIG)
  )
)
write_json(json_out, file.path(OUT_DIR, "TIER_DECISION_REPORT.json"),
           pretty = TRUE, auto_unbox = TRUE)

cat("============================================\n")
cat("R9 全部完成!\n")
cat("结束时间:", format(Sys.time()), "\n")
cat("============================================\n")
cat("\n关键输出 (请把这个文件的内容贴给 Claude):\n")
cat("  ", report_path, "\n")
