#!/usr/bin/env Rscript
# ============================================================
# R8: 跨数据集 / 跨物种 density gene 可重复性分析
# (EN) R08_cross_dataset_comparison: cross-dataset / cross-species
#      reproducibility of density genes. Paths via env vars
#      (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================
# 输入: R3/R4/R7 的输出, P1 基因交集分析结果
# 输出: 两两比较结果, 全局汇总, 可视化
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

options(future.globals.maxSize = Inf)

# --- 路径配置 ---
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
REGISTRY_PATH <- file.path(DATA_DIR, "sample_registry.json")
R3_DIR <- file.path(RESULTS_DIR, "R3_Results")
R4_DIR <- file.path(RESULTS_DIR, "R4_Results")
R7_DIR <- file.path(RESULTS_DIR, "R7_Results")
P1_INTERSECTION_DIR <- file.path(RESULTS_DIR, "P1_Results", "Gene_Intersection")
OUTPUT_DIR <- file.path(RESULTS_DIR, "R8_Results")

MIN_SHARED_GENES <- 50  # 共有基因≥50才做比较

cat("============================================\n")
cat("R8: 跨数据集 / 跨物种可重复性分析\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time()), "\n")
cat("最小共有基因数:", MIN_SHARED_GENES, "\n\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "pairwise_comparisons"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cross_species_comparisons"), showWarnings = FALSE)

# --- 读取 registry ---
registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
dataset_map <- sapply(sample_names, function(s) strsplit(s, "/")[[1]][1])
datasets <- unique(dataset_map)
species_map <- sapply(sample_names, function(s) registry[[s]]$species)

cat("数据集:", paste(datasets, collapse = ", "), "\n\n")

# ============================================================
# 阶段 0: 收集每个数据集的基因列表和density gene信息
# ============================================================

cat("============================================\n")
cat("阶段 0: 收集数据集信息\n")
cat("============================================\n\n")

dataset_info <- list()

for (ds in datasets) {
  ds_samples <- sample_names[dataset_map == ds]
  ds_species <- species_map[ds_samples[1]]
  
  # --- 收集该数据集所有样本的基因列表和R3/R4结果 ---
  all_genes <- c()
  sample_r3_list <- list()
  sample_r4_list <- list()
  
  for (sname in ds_samples) {
    parts <- strsplit(sname, "/")[[1]]
    
    # R3 结果
    r3_path <- file.path(R3_DIR, parts[1], parts[2], "density_gene_correlations.csv")
    if (file.exists(r3_path)) {
      r3 <- fread(r3_path)
      sample_r3_list[[sname]] <- r3
      all_genes <- union(all_genes, r3$gene)
    }
    
    # R4 结果
    r4_path <- file.path(R4_DIR, parts[1], parts[2], "filtered_density_genes.csv")
    if (file.exists(r4_path)) {
      r4 <- fread(r4_path)
      sample_r4_list[[sname]] <- r4
    }
  }
  
  # --- 读取R7一致性结果 ---
  consistency_path <- file.path(R7_DIR, ds, "dataset_consistency.csv")
  consistent_genes_tier1 <- c()
  consistent_genes_tier2 <- c()
  
  if (file.exists(consistency_path)) {
    cons <- fread(consistency_path)
    if ("consistent" %in% names(cons)) {
      consistent_genes_tier1 <- cons[consistent == "consistent" | consistent == "single_sample"]$gene
    }
  }
  
  # --- 计算数据集级别的 median rho (跨样本取中位数) ---
  # 用 rho_knn_main
  rho_col <- "rho_knn_main"
  
  gene_median_rho <- data.table(gene = all_genes)
  
  for (sname in names(sample_r3_list)) {
    r3 <- sample_r3_list[[sname]]
    if (rho_col %in% names(r3)) {
      col_name <- paste0("rho_", gsub("/", "_", sname))
      gene_median_rho <- merge(gene_median_rho, 
                                r3[, .(gene, rho_val = get(rho_col))],
                                by = "gene", all.x = TRUE)
      setnames(gene_median_rho, "rho_val", col_name)
    }
  }
  
  # 计算跨样本 median rho
  rho_cols <- grep("^rho_", names(gene_median_rho), value = TRUE)
  if (length(rho_cols) > 0) {
    gene_median_rho$median_rho <- apply(gene_median_rho[, ..rho_cols], 1, 
                                         function(x) median(x, na.rm = TRUE))
  } else {
    gene_median_rho$median_rho <- NA
  }
  
  # --- 确定该数据集的 density genes ---
  # tier1: 一致基因 (R7) 或 单样本的 tier1
  # tier2: 在任意样本中是 tier2 的基因 (更宽松)
  all_tier1 <- c()
  all_tier2 <- c()
  for (sname in names(sample_r4_list)) {
    r4 <- sample_r4_list[[sname]]
    all_tier1 <- union(all_tier1, r4[tier == "tier1_strict"]$gene)
    all_tier2 <- union(all_tier2, r4[tier %in% c("tier1_strict", "tier2_moderate")]$gene)
  }
  
  # 用R7一致基因作为数据集代表 (如果有多样本)
  if (length(ds_samples) > 1 && length(consistent_genes_tier1) > 0) {
    representative_genes_tier1 <- consistent_genes_tier1
  } else {
    representative_genes_tier1 <- all_tier1
  }
  
  dataset_info[[ds]] <- list(
    name = ds,
    species = ds_species,
    n_samples = length(ds_samples),
    all_genes = all_genes,
    all_genes_upper = toupper(all_genes),  # 跨物种用
    gene_name_map = setNames(all_genes, toupper(all_genes)),  # 大写→原始名映射
    representative_tier1 = representative_genes_tier1,
    representative_tier1_upper = toupper(representative_genes_tier1),
    all_tier1 = all_tier1,
    all_tier2 = all_tier2,
    gene_median_rho = gene_median_rho,
    sample_r3 = sample_r3_list,
    sample_r4 = sample_r4_list
  )
  
  cat(sprintf("  %-30s species=%-5s samples=%d genes=%d tier1_rep=%d tier1_all=%d tier2_all=%d\n",
              ds, ds_species, length(ds_samples), length(all_genes),
              length(representative_genes_tier1), length(all_tier1), length(all_tier2)))
}

cat("\n")

# ============================================================
# 阶段 1: 生成所有可行的比较对
# ============================================================

cat("============================================\n")
cat("阶段 1: 确定比较对\n")
cat("============================================\n\n")

comparisons <- list()
comp_idx <- 0

for (i in 1:(length(datasets) - 1)) {
  for (j in (i + 1):length(datasets)) {
    ds_a <- datasets[i]
    ds_b <- datasets[j]
    info_a <- dataset_info[[ds_a]]
    info_b <- dataset_info[[ds_b]]
    
    # 确定是否跨物种
    same_species <- (info_a$species == info_b$species)
    
    if (same_species) {
      # 同物种: 直接比较基因名
      shared_genes <- intersect(info_a$all_genes, info_b$all_genes)
      comparison_type <- "within_species"
    } else {
      # 跨物种: 大写统一后比较
      shared_upper <- intersect(info_a$all_genes_upper, info_b$all_genes_upper)
      shared_genes <- shared_upper  # 用大写版本作为共有基因标识
      comparison_type <- "cross_species"
    }
    
    n_shared <- length(shared_genes)
    
    if (n_shared >= MIN_SHARED_GENES) {
      comp_idx <- comp_idx + 1
      comparisons[[comp_idx]] <- list(
        ds_a = ds_a,
        ds_b = ds_b,
        shared_genes = shared_genes,
        n_shared = n_shared,
        same_species = same_species,
        comparison_type = comparison_type
      )
      cat(sprintf("  [%d] %-25s vs %-25s shared=%d type=%s\n",
                  comp_idx, ds_a, ds_b, n_shared, comparison_type))
    } else {
      cat(sprintf("  [x] %-25s vs %-25s shared=%d < %d, 跳过\n",
                  ds_a, ds_b, n_shared, MIN_SHARED_GENES))
    }
  }
}

cat(sprintf("\n共 %d 个可行比较对\n\n", length(comparisons)))

# ============================================================
# 阶段 2: 逐对比较
# ============================================================

cat("============================================\n")
cat("阶段 2: 逐对比较\n")
cat("============================================\n\n")

comparison_results <- list()

for (ci in seq_along(comparisons)) {
  comp <- comparisons[[ci]]
  ds_a <- comp$ds_a
  ds_b <- comp$ds_b
  info_a <- dataset_info[[ds_a]]
  info_b <- dataset_info[[ds_b]]
  shared <- comp$shared_genes
  
  cat("========================================\n")
  cat(sprintf("[%d/%d] %s vs %s (%d shared genes, %s)\n",
              ci, length(comparisons), ds_a, ds_b, comp$n_shared, comp$comparison_type))
  cat("========================================\n")
  
  # --- 获取两侧在共有基因上的 median rho ---
  if (comp$same_species) {
    # 同物种: 直接匹配
    rho_a <- info_a$gene_median_rho[gene %in% shared, .(gene, rho_a = median_rho)]
    rho_b <- info_b$gene_median_rho[gene %in% shared, .(gene, rho_b = median_rho)]
    merged_rho <- merge(rho_a, rho_b, by = "gene")
  } else {
    # 跨物种: 用大写匹配
    rho_a_dt <- copy(info_a$gene_median_rho)
    rho_a_dt[, gene_upper := toupper(gene)]
    rho_a_dt <- rho_a_dt[gene_upper %in% shared, .(gene_upper, gene_a = gene, rho_a = median_rho)]
    
    rho_b_dt <- copy(info_b$gene_median_rho)
    rho_b_dt[, gene_upper := toupper(gene)]
    rho_b_dt <- rho_b_dt[gene_upper %in% shared, .(gene_upper, gene_b = gene, rho_b = median_rho)]
    
    merged_rho <- merge(rho_a_dt, rho_b_dt, by = "gene_upper")
    merged_rho[, gene := gene_upper]
  }
  
  # 去除 NA
  merged_rho <- merged_rho[!is.na(rho_a) & !is.na(rho_b)]
  n_valid <- nrow(merged_rho)
  
  if (n_valid < 10) {
    cat("  ⚠ 有效基因数不足 (<10), 跳过\n\n")
    next
  }
  
  # --- 方法1: 排名相关 (两侧 median rho 的 Spearman) ---
  rank_cor <- cor.test(merged_rho$rho_a, merged_rho$rho_b, method = "spearman")
  rank_rho <- round(rank_cor$estimate, 4)
  rank_p <- rank_cor$p.value
  
  # --- 方法2: 列表重叠 (tier1) ---
  if (comp$same_species) {
    tier1_a <- intersect(info_a$representative_tier1, shared)
    tier1_b <- intersect(info_b$representative_tier1, shared)
  } else {
    tier1_a <- intersect(info_a$representative_tier1_upper, shared)
    tier1_b <- intersect(info_b$representative_tier1_upper, shared)
  }
  
  overlap_tier1 <- intersect(tier1_a, tier1_b)
  n_overlap_tier1 <- length(overlap_tier1)
  
  # Fisher exact test
  n_total <- comp$n_shared
  n_a <- length(tier1_a)
  n_b <- length(tier1_b)
  
  if (n_a > 0 && n_b > 0 && n_total > 0) {
    # 构建 2x2 列联表
    a_and_b <- n_overlap_tier1
    a_not_b <- n_a - a_and_b
    b_not_a <- n_b - a_and_b
    neither <- n_total - n_a - n_b + a_and_b
    
    # 确保非负
    neither <- max(0, neither)
    
    contingency <- matrix(c(a_and_b, a_not_b, b_not_a, neither), nrow = 2)
    fisher_result <- tryCatch(fisher.test(contingency, alternative = "greater"),
                               error = function(e) list(p.value = NA, estimate = NA))
    fisher_p <- fisher_result$p.value
    fisher_or <- ifelse(is.null(fisher_result$estimate), NA, round(fisher_result$estimate, 2))
    expected_overlap <- round(n_a * n_b / n_total, 1)
  } else {
    fisher_p <- NA
    fisher_or <- NA
    expected_overlap <- 0
  }
  
  # --- 方法2b: 列表重叠 (tier2, 敏感性分析) ---
  if (comp$same_species) {
    tier2_a <- intersect(info_a$all_tier2, shared)
    tier2_b <- intersect(info_b$all_tier2, shared)
  } else {
    tier2_a <- intersect(toupper(info_a$all_tier2), shared)
    tier2_b <- intersect(toupper(info_b$all_tier2), shared)
  }
  overlap_tier2 <- intersect(tier2_a, tier2_b)
  n_overlap_tier2 <- length(overlap_tier2)
  
  # --- 方法3: 方向一致性 ---
  # 在两侧都是 tier1 或 tier2 的基因中, 检查 rho 符号是否一致
  sig_both <- merged_rho[gene %in% union(overlap_tier1, overlap_tier2)]
  if (nrow(sig_both) > 0) {
    n_same_dir <- sum(sign(sig_both$rho_a) == sign(sig_both$rho_b))
    direction_consistency <- round(n_same_dir / nrow(sig_both) * 100, 1)
  } else {
    n_same_dir <- 0
    direction_consistency <- NA
  }
  
  # --- 输出结果 ---
  cat(sprintf("  排名相关: ρ=%.4f, p=%.2e\n", rank_rho, rank_p))
  cat(sprintf("  Tier1重叠: %d/%d(A) ∩ %d/%d(B) = %d (期望=%.1f, Fisher p=%.2e, OR=%.2f)\n",
              n_a, n_total, n_b, n_total, n_overlap_tier1, expected_overlap,
              ifelse(is.na(fisher_p), 0, fisher_p),
              ifelse(is.na(fisher_or), 0, fisher_or)))
  cat(sprintf("  Tier2重叠: %d ∩ %d = %d\n",
              length(tier2_a), length(tier2_b), n_overlap_tier2))
  cat(sprintf("  方向一致性: %d/%d = %s%%\n",
              n_same_dir, nrow(sig_both),
              ifelse(is.na(direction_consistency), "N/A", direction_consistency)))
  
  if (n_overlap_tier1 > 0) {
    cat(sprintf("  重叠tier1基因: %s\n", paste(head(overlap_tier1, 15), collapse = ", ")))
  }
  
  # --- 保存比较对级别的结果 ---
  comp_dir <- file.path(OUTPUT_DIR, "pairwise_comparisons", 
                         paste0(ds_a, "_vs_", ds_b))
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 基因级别结果
  gene_level <- copy(merged_rho)
  if (comp$same_species) {
    gene_level[, tier1_in_A := gene %in% tier1_a]
    gene_level[, tier1_in_B := gene %in% tier1_b]
    gene_level[, tier2_in_A := gene %in% tier2_a]
    gene_level[, tier2_in_B := gene %in% tier2_b]
  } else {
    gene_level[, tier1_in_A := gene_upper %in% tier1_a]
    gene_level[, tier1_in_B := gene_upper %in% tier1_b]
    gene_level[, tier2_in_A := gene_upper %in% tier2_a]
    gene_level[, tier2_in_B := gene_upper %in% tier2_b]
  }
  gene_level[, same_direction := sign(rho_a) == sign(rho_b)]
  gene_level[, both_tier1 := tier1_in_A & tier1_in_B]
  
  fwrite(gene_level, file.path(comp_dir, "gene_level_comparison.csv"))
  
  # 汇总结果
  summary_row <- data.frame(
    dataset_A = ds_a,
    dataset_B = ds_b,
    species_A = info_a$species,
    species_B = info_b$species,
    comparison_type = comp$comparison_type,
    n_shared_genes = comp$n_shared,
    n_valid_genes = n_valid,
    rank_correlation_rho = rank_rho,
    rank_correlation_p = rank_p,
    tier1_A = n_a,
    tier1_B = n_b,
    tier1_overlap = n_overlap_tier1,
    tier1_expected = expected_overlap,
    tier1_fisher_p = fisher_p,
    tier1_odds_ratio = fisher_or,
    tier2_A = length(tier2_a),
    tier2_B = length(tier2_b),
    tier2_overlap = n_overlap_tier2,
    direction_consistency_pct = direction_consistency,
    overlap_genes_tier1 = paste(head(overlap_tier1, 20), collapse = "; "),
    stringsAsFactors = FALSE
  )
  
  fwrite(summary_row, file.path(comp_dir, "comparison_summary.csv"))
  comparison_results[[ci]] <- summary_row
  
  # --- 散点图: A侧 median rho vs B侧 median rho ---
  tryCatch({
    p <- ggplot(merged_rho, aes(x = rho_a, y = rho_b)) +
      geom_point(alpha = 0.4, size = 1.5, color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
      labs(
        title = sprintf("%s vs %s (n=%d genes)", ds_a, ds_b, n_valid),
        subtitle = sprintf("Rank ρ=%.3f, p=%.2e | Tier1 overlap=%d (expected=%.1f, p=%.2e)",
                           rank_rho, rank_p, n_overlap_tier1, expected_overlap,
                           ifelse(is.na(fisher_p), 1, fisher_p)),
        x = paste0("Median ρ (", ds_a, ")"),
        y = paste0("Median ρ (", ds_b, ")")
      ) +
      theme_minimal(base_size = 11) +
      coord_fixed()
    
    # 标注重叠的 tier1 基因
    if (n_overlap_tier1 > 0) {
      label_genes <- head(overlap_tier1, 10)
      label_dt <- merged_rho[gene %in% label_genes]
      p <- p + geom_point(data = label_dt, color = "red", size = 2.5) +
        geom_text(data = label_dt, aes(label = gene), 
                  hjust = -0.1, vjust = -0.3, size = 2.8, color = "red")
    }
    
    ggsave(file.path(comp_dir, "rho_scatter.png"), p, width = 8, height = 8, dpi = 150)
    cat("  ✓ 散点图已保存\n")
  }, error = function(e) {
    cat("  ⚠ 散点图生成失败:", conditionMessage(e), "\n")
  })
  
  cat("\n")
}

# ============================================================
# 阶段 3: 全局汇总
# ============================================================

cat("============================================\n")
cat("阶段 3: 全局汇总\n")
cat("============================================\n\n")

if (length(comparison_results) > 0) {
  all_comparisons <- rbindlist(comparison_results, fill = TRUE)
  fwrite(all_comparisons, file.path(OUTPUT_DIR, "ALL_COMPARISONS_R8_SUMMARY.csv"))
  cat("✓ 汇总已保存:", file.path(OUTPUT_DIR, "ALL_COMPARISONS_R8_SUMMARY.csv"), "\n\n")
  
  # --- 打印汇总表 ---
  cat(sprintf("%-25s %-25s %5s %5s %6s %6s %5s %5s %5s %8s %6s\n",
              "Dataset_A", "Dataset_B", "type", "shared", "rank_ρ", "rank_p",
              "t1_A", "t1_B", "t1_ov", "fisher_p", "dir%"))
  cat(paste(rep("-", 140), collapse = ""), "\n")
  
  for (i in 1:nrow(all_comparisons)) {
    row <- all_comparisons[i, ]
    type_short <- ifelse(row$comparison_type == "within_species", "same", "cross")
    cat(sprintf("%-25s %-25s %5s %5d %6.3f %6.1e %5d %5d %5d %8.1e %6s\n",
                row$dataset_A, row$dataset_B, type_short, row$n_shared_genes,
                row$rank_correlation_rho, row$rank_correlation_p,
                row$tier1_A, row$tier1_B, row$tier1_overlap,
                ifelse(is.na(row$tier1_fisher_p), 1, row$tier1_fisher_p),
                ifelse(is.na(row$direction_consistency_pct), "N/A", 
                       sprintf("%.1f", row$direction_consistency_pct))))
  }
  
  # --- 按比较类型分组统计 ---
  cat("\n--- 按类型分组 ---\n")
  within <- all_comparisons[comparison_type == "within_species"]
  cross <- all_comparisons[comparison_type == "cross_species"]
  
  if (nrow(within) > 0) {
    cat(sprintf("  同物种比较: %d 对, 排名相关ρ: %.3f ~ %.3f (median=%.3f)\n",
                nrow(within), min(within$rank_correlation_rho), 
                max(within$rank_correlation_rho),
                median(within$rank_correlation_rho)))
    sig_within <- sum(within$tier1_fisher_p < 0.05, na.rm = TRUE)
    cat(sprintf("  tier1重叠显著(p<0.05): %d/%d 对\n", sig_within, nrow(within)))
  }
  
  if (nrow(cross) > 0) {
    cat(sprintf("  跨物种比较: %d 对, 排名相关ρ: %.3f ~ %.3f (median=%.3f)\n",
                nrow(cross), min(cross$rank_correlation_rho),
                max(cross$rank_correlation_rho),
                median(cross$rank_correlation_rho)))
    sig_cross <- sum(cross$tier1_fisher_p < 0.05, na.rm = TRUE)
    cat(sprintf("  tier1重叠显著(p<0.05): %d/%d 对\n", sig_cross, nrow(cross)))
  }
  
} else {
  cat("⚠ 无有效比较对\n")
}

# ============================================================
# 阶段 4: 全局 density gene 景观
# ============================================================

cat("\n============================================\n")
cat("阶段 4: 全局 density gene 景观\n")
cat("============================================\n\n")

# 收集所有数据集的一致 density genes, 构建全局视图
landscape <- list()

for (ds in datasets) {
  info <- dataset_info[[ds]]
  rep_genes <- info$representative_tier1
  
  for (g in rep_genes) {
    g_upper <- toupper(g)
    rho_row <- info$gene_median_rho[gene == g]
    med_rho <- ifelse(nrow(rho_row) > 0, round(rho_row$median_rho[1], 4), NA)
    
    landscape[[paste0(ds, "_", g)]] <- data.frame(
      gene_original = g,
      gene_upper = g_upper,
      dataset = ds,
      species = info$species,
      median_rho = med_rho,
      direction = ifelse(!is.na(med_rho) && med_rho > 0, "positive", 
                          ifelse(!is.na(med_rho) && med_rho < 0, "negative", "NA")),
      stringsAsFactors = FALSE
    )
  }
}

if (length(landscape) > 0) {
  landscape_df <- rbindlist(landscape, fill = TRUE)
  
  # 统计每个基因（大写统一）在几个数据集中是一致 density gene
  gene_summary <- landscape_df[, .(
    n_datasets = .N,
    datasets = paste(dataset, collapse = "; "),
    species_list = paste(unique(species), collapse = "; "),
    directions = paste(direction, collapse = "; "),
    median_rhos = paste(round(median_rho, 4), collapse = "; ")
  ), by = gene_upper]
  
  gene_summary <- gene_summary[order(-n_datasets)]
  
  fwrite(landscape_df, file.path(OUTPUT_DIR, "global_density_gene_landscape.csv"))
  fwrite(gene_summary, file.path(OUTPUT_DIR, "global_gene_summary.csv"))
  
  cat("全局 density gene 景观:\n")
  cat(sprintf("  总基因数(大写统一): %d\n", nrow(gene_summary)))
  cat(sprintf("  在≥2个数据集中出现: %d\n", sum(gene_summary$n_datasets >= 2)))
  cat(sprintf("  在≥3个数据集中出现: %d\n", sum(gene_summary$n_datasets >= 3)))
  
  # 打印在多个数据集中出现的基因
  multi_ds <- gene_summary[n_datasets >= 2]
  if (nrow(multi_ds) > 0) {
    cat("\n  跨数据集重复的 density genes:\n")
    for (i in 1:min(30, nrow(multi_ds))) {
      row <- multi_ds[i, ]
      cat(sprintf("    %-15s %d个数据集: %s | ρ: %s | 方向: %s\n",
                  row$gene_upper, row$n_datasets, row$datasets,
                  row$median_rhos, row$directions))
    }
  }
} else {
  cat("⚠ 无 density gene 数据\n")
}

cat("\n============================================\n")
cat("R8 全部完成!\n")
cat("结束时间:", format(Sys.time()), "\n")
cat("============================================\n")
