#!/usr/bin/env Rscript
# ============================================================
# R7: 样本内综合 + 数据集内跨样本一致性
# (EN) R7_sample_integration: per-sample integration + within-dataset
#      cross-sample consistency of density genes. Paths via env vars
#      (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================
# 输入: R3/R4/R6 的输出 csv
# 输出: 每样本完整档案, 数据集一致性评估, 全局汇总
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

options(future.globals.maxSize = Inf)

# --- 路径配置 ---
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
REGISTRY_PATH <- file.path(DATA_DIR, "sample_registry.json")
R3_DIR <- file.path(RESULTS_DIR, "R3_Results")
R4_DIR <- file.path(RESULTS_DIR, "R4_Results")
R6_DIR <- file.path(RESULTS_DIR, "R6_Results")
OUTPUT_DIR <- file.path(RESULTS_DIR, "R7_Results")
CONSISTENCY_THRESHOLD <- 0.5  # 在≥50%样本中为tier1才算一致

cat("============================================\n")
cat("R7: 样本内综合 + 跨样本一致性\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time()), "\n")
cat("一致性阈值:", CONSISTENCY_THRESHOLD, "\n\n")

# --- 读取 registry ---
registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
cat("样本数:", length(sample_names), "\n\n")

# --- 按数据集分组 ---
dataset_map <- sapply(sample_names, function(s) strsplit(s, "/")[[1]][1])
datasets <- unique(dataset_map)

# ============================================================
# 阶段 1: 逐样本整合
# ============================================================

profile_list <- list()

for (i in seq_along(sample_names)) {
  sname <- sample_names[i]
  parts <- strsplit(sname, "/")[[1]]
  dataset_name <- parts[1]
  sample_subname <- parts[2]
  
  cat("========================================\n")
  cat(sprintf("[%d/%d] %s\n", i, length(sample_names), sname))
  cat("========================================\n")
  
  # --- 读取 R3 结果 ---
  r3_path <- file.path(R3_DIR, dataset_name, sample_subname, "density_gene_correlations.csv")
  if (!file.exists(r3_path)) {
    cat("  ⚠ R3 文件不存在, 跳过\n\n")
    next
  }
  r3 <- fread(r3_path)
  
  # --- 读取 R4 结果 ---
  r4_path <- file.path(R4_DIR, dataset_name, sample_subname, "filtered_density_genes.csv")
  if (!file.exists(r4_path)) {
    cat("  ⚠ R4 文件不存在, 跳过\n\n")
    next
  }
  r4 <- fread(r4_path)
  
  # --- 读取 R6 结果 ---
  r6_path <- file.path(R6_DIR, dataset_name, sample_subname, "cell_state_coupling.csv")
  r6_exists <- file.exists(r6_path)
  if (r6_exists) {
    r6 <- fread(r6_path)
  }
  
  # --- 合并 ---
  # R4 包含 R3 的核心信息 (rho, q, tier) + 筛选标签
  # R6 包含 effect decomposition
  merged <- copy(r4)
  
  if (r6_exists && nrow(r6) > 0) {
    # R6 只有 tier1 基因, 左连接
    r6_cols <- intersect(names(r6), c("gene", "effect_class", "median_within_rho", 
                                       "n_clusters_analyzed", "n_reg_clusters"))
    if (length(r6_cols) > 1) {
      r6_sub <- r6[, ..r6_cols]
      merged <- merge(merged, r6_sub, by = "gene", all.x = TRUE)
    }
  }
  
  # --- 输出样本完整档案 ---
  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(merged, file.path(out_dir, "sample_density_profile.csv"))
  
  # --- 样本级摘要 ---
  n_total <- nrow(merged)
  n_tier1 <- sum(merged$tier == "tier1_strict", na.rm = TRUE)
  n_tier2 <- sum(merged$tier == "tier2_moderate", na.rm = TRUE)
  n_tier3 <- sum(merged$tier == "tier3_lenient", na.rm = TRUE)
  n_ns <- sum(merged$tier == "not_significant", na.rm = TRUE)
  
  # 用 knn_main 的 rho
  rho_col <- grep("rho_knn_main", names(merged), value = TRUE)[1]
  if (is.null(rho_col) || is.na(rho_col)) rho_col <- "rho_knn_main"
  
  tier1_rows <- merged[tier == "tier1_strict"]
  if (nrow(tier1_rows) > 0) {
    n_pos <- sum(tier1_rows[[rho_col]] > 0, na.rm = TRUE)
    n_neg <- sum(tier1_rows[[rho_col]] < 0, na.rm = TRUE)
    med_rho <- round(median(abs(tier1_rows[[rho_col]]), na.rm = TRUE), 4)
    max_rho <- round(max(abs(tier1_rows[[rho_col]]), na.rm = TRUE), 4)
    
    # top5 基因
    top5 <- tier1_rows[order(-abs(get(rho_col)))][1:min(5, nrow(tier1_rows))]
    top5_str <- paste(sprintf("%s(%.3f)", top5$gene, top5[[rho_col]]), collapse = ", ")
    
    # R6 标签统计 (如果有)
    if ("effect_class" %in% names(tier1_rows)) {
      n_comp <- sum(tier1_rows$effect_class == "composition_driven", na.rm = TRUE)
      n_reg <- sum(tier1_rows$effect_class == "regulation_present", na.rm = TRUE)
      n_het <- sum(tier1_rows$effect_class == "heterogeneous", na.rm = TRUE)
      reg_pct <- round(n_reg / n_tier1 * 100, 1)
    } else {
      n_comp <- NA; n_reg <- NA; n_het <- NA; reg_pct <- NA
    }
  } else {
    n_pos <- 0; n_neg <- 0; med_rho <- NA; max_rho <- NA
    top5_str <- ""; n_comp <- NA; n_reg <- NA; n_het <- NA; reg_pct <- NA
  }
  
  profile_list[[sname]] <- data.frame(
    sample_name = sname,
    dataset = dataset_name,
    n_genes_total = n_total,
    n_tier1 = n_tier1,
    n_tier2 = n_tier2,
    n_tier3 = n_tier3,
    n_not_sig = n_ns,
    tier1_pos = n_pos,
    tier1_neg = n_neg,
    tier1_neg_pct = ifelse(n_tier1 > 0, round(n_neg / n_tier1 * 100, 1), NA),
    tier1_median_abs_rho = med_rho,
    tier1_max_abs_rho = max_rho,
    tier1_regulation_pct = reg_pct,
    tier1_composition = ifelse(!is.na(n_comp), n_comp, NA),
    top5_genes = top5_str,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  总基因=%d, tier1=%d(+%d/-%d), tier2=%d, med|ρ|=%.4f, max|ρ|=%.4f\n",
              n_total, n_tier1, n_pos, n_neg, n_tier2, 
              ifelse(is.na(med_rho), 0, med_rho),
              ifelse(is.na(max_rho), 0, max_rho)))
  cat(sprintf("  regulation=%s%%, top5: %s\n", 
              ifelse(is.na(reg_pct), "N/A", reg_pct), top5_str))
  cat("\n")
}

# --- 输出全局样本摘要 ---
profile_df <- do.call(rbind, profile_list)
fwrite(profile_df, file.path(OUTPUT_DIR, "ALL_SAMPLES_R7_PROFILE.csv"))
cat("✓ 样本摘要已保存:", file.path(OUTPUT_DIR, "ALL_SAMPLES_R7_PROFILE.csv"), "\n\n")

# ============================================================
# 阶段 2: 数据集内跨样本一致性
# ============================================================

cat("============================================\n")
cat("阶段 2: 数据集内跨样本一致性\n")
cat("============================================\n\n")

consistency_list <- list()

for (ds in datasets) {
  ds_samples <- sample_names[dataset_map == ds]
  n_samples <- length(ds_samples)
  
  cat(sprintf("--- %s (%d 样本) ---\n", ds, n_samples))
  
  ds_out_dir <- file.path(OUTPUT_DIR, ds)
  dir.create(ds_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (n_samples == 1) {
    cat("  单样本数据集, 一致性分析 N/A\n")
    
    # 读取该样本的 tier1 列表
    parts <- strsplit(ds_samples[1], "/")[[1]]
    r4_path <- file.path(R4_DIR, parts[1], parts[2], "filtered_density_genes.csv")
    if (file.exists(r4_path)) {
      r4 <- fread(r4_path)
      tier1_genes <- r4[tier == "tier1_strict"]$gene
      
      consistency_list[[ds]] <- data.frame(
        dataset = ds,
        n_samples = 1,
        n_tier1_union = length(tier1_genes),
        n_consistent = NA,
        consistency_pct = NA,
        consistent_genes = paste(head(tier1_genes, 10), collapse = ", "),
        note = "single_sample",
        stringsAsFactors = FALSE
      )
    }
    
    # 输出
    fwrite(data.frame(gene = tier1_genes, n_samples_tier1 = 1, pct = 100,
                      consistent = "single_sample", stringsAsFactors = FALSE),
           file.path(ds_out_dir, "dataset_consistency.csv"))
    
    cat("\n")
    next
  }
  
  # --- 多样本数据集: 收集每个样本的 tier1 基因列表 ---
  rho_col <- "rho_knn_main"
  
  tier1_by_sample <- list()
  rho_by_sample <- list()
  
  for (sname in ds_samples) {
    parts <- strsplit(sname, "/")[[1]]
    r4_path <- file.path(R4_DIR, parts[1], parts[2], "filtered_density_genes.csv")
    if (!file.exists(r4_path)) next
    
    r4 <- fread(r4_path)
    tier1 <- r4[tier == "tier1_strict"]
    tier1_by_sample[[sname]] <- tier1$gene
    
    # 记录每个基因的 rho
    rho_vals <- setNames(tier1[[rho_col]], tier1$gene)
    rho_by_sample[[sname]] <- rho_vals
  }
  
  # --- tier1 并集 ---
  all_tier1_genes <- unique(unlist(tier1_by_sample))
  n_union <- length(all_tier1_genes)
  
  if (n_union == 0) {
    cat("  无 tier1 基因\n\n")
    consistency_list[[ds]] <- data.frame(
      dataset = ds, n_samples = n_samples, n_tier1_union = 0,
      n_consistent = 0, consistency_pct = NA, consistent_genes = "",
      note = "no_tier1", stringsAsFactors = FALSE)
    next
  }
  
  # --- 统计每个基因在几个样本中是 tier1 ---
  gene_counts <- data.frame(gene = all_tier1_genes, stringsAsFactors = FALSE)
  gene_counts$n_samples_tier1 <- sapply(gene_counts$gene, function(g) {
    sum(sapply(tier1_by_sample, function(x) g %in% x))
  })
  gene_counts$pct <- round(gene_counts$n_samples_tier1 / n_samples * 100, 1)
  
  # --- 每个基因在各样本中的 rho ---
  for (sname in ds_samples) {
    col_name <- paste0("rho_", gsub("/", "_", sname))
    gene_counts[[col_name]] <- sapply(gene_counts$gene, function(g) {
      rv <- rho_by_sample[[sname]]
      if (!is.null(rv) && g %in% names(rv)) round(rv[g], 4) else NA
    })
  }
  
  # --- 标记一致性 ---
  min_samples <- ceiling(n_samples * CONSISTENCY_THRESHOLD)
  gene_counts$consistent <- ifelse(gene_counts$n_samples_tier1 >= min_samples, "consistent", "sporadic")
  
  # 按出现频率排序
  gene_counts <- gene_counts[order(-gene_counts$n_samples_tier1, -gene_counts$pct), ]
  
  n_consistent <- sum(gene_counts$consistent == "consistent")
  consistent_genes <- gene_counts$gene[gene_counts$consistent == "consistent"]
  
  cat(sprintf("  tier1并集=%d, 一致(≥%d/%d样本)=%d (%.1f%%)\n",
              n_union, min_samples, n_samples, n_consistent,
              ifelse(n_union > 0, n_consistent / n_union * 100, 0)))
  cat(sprintf("  一致基因: %s\n", 
              paste(head(consistent_genes, 10), collapse = ", ")))
  
  # --- 输出 ---
  fwrite(gene_counts, file.path(ds_out_dir, "dataset_consistency.csv"))
  
  consistency_list[[ds]] <- data.frame(
    dataset = ds,
    n_samples = n_samples,
    n_tier1_union = n_union,
    n_consistent = n_consistent,
    consistency_pct = round(n_consistent / n_union * 100, 1),
    consistent_genes = paste(head(consistent_genes, 15), collapse = ", "),
    note = "",
    stringsAsFactors = FALSE
  )
  
  cat("\n")
}

# --- 输出全局一致性汇总 ---
consistency_df <- do.call(rbind, consistency_list)
fwrite(consistency_df, file.path(OUTPUT_DIR, "ALL_DATASETS_R7_CONSISTENCY.csv"))
cat("✓ 一致性汇总已保存:", file.path(OUTPUT_DIR, "ALL_DATASETS_R7_CONSISTENCY.csv"), "\n\n")

# ============================================================
# 阶段 3: 打印汇总
# ============================================================

cat("============================================\n")
cat("全局汇总\n")
cat("============================================\n\n")

cat("--- 样本摘要 ---\n")
cat(sprintf("%-45s %5s %5s %5s %6s %6s %6s %s\n",
            "sample", "tier1_strict", "+", "-", "med|ρ|", "max|ρ|", "reg%", "top3"))
cat(paste(rep("-", 130), collapse = ""), "\n")
for (i in 1:nrow(profile_df)) {
  row <- profile_df[i, ]
  top3 <- paste(head(strsplit(row$top5_genes, ", ")[[1]], 3), collapse = ", ")
  cat(sprintf("%-45s %5d %5d %5d %6s %6s %6s %s\n",
              row$sample_name, row$n_tier1, row$tier1_pos, row$tier1_neg,
              ifelse(is.na(row$tier1_median_abs_rho), "N/A", sprintf("%.4f", row$tier1_median_abs_rho)),
              ifelse(is.na(row$tier1_max_abs_rho), "N/A", sprintf("%.4f", row$tier1_max_abs_rho)),
              ifelse(is.na(row$tier1_regulation_pct), "N/A", sprintf("%.1f", row$tier1_regulation_pct)),
              top3))
}

cat("\n--- 数据集一致性 ---\n")
for (i in 1:nrow(consistency_df)) {
  row <- consistency_df[i, ]
  cat(sprintf("  %-30s %d样本, tier1并集=%s, 一致=%s (%s%%)\n",
              row$dataset, row$n_samples,
              ifelse(is.na(row$n_tier1_union), "N/A", row$n_tier1_union),
              ifelse(is.na(row$n_consistent), "N/A", row$n_consistent),
              ifelse(is.na(row$consistency_pct), "N/A", row$consistency_pct)))
  if (!is.na(row$consistent_genes) && nchar(row$consistent_genes) > 0 && row$note != "single_sample") {
    cat(sprintf("    一致基因: %s\n", row$consistent_genes))
  }
}

cat("\n============================================\n")
cat("R7 全部完成!\n")
cat("结束时间:", format(Sys.time()), "\n")
cat("============================================\n")
