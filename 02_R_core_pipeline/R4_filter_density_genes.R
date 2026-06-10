#!/usr/bin/env Rscript
# ============================================================================
# R4_filter_density_genes.R
# 功能：从R3的相关系数表中筛选density genes，分级标记，生成诊断图
# 输入：R3输出的 density_gene_correlations.csv（22个样本）
# 输出：filtered_density_genes.csv + 诊断图 + 全局汇总
# Run (EN): Rscript R4_filter_density_genes.R
#   Purpose: filter and tier density genes from the R3 correlation tables; diagnostics.
#   Paths configured via env vars (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================================

suppressPackageStartupMessages({
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
R3_DIR         <- file.path(RESULTS_DIR, "R3_Results")
OUTPUT_DIR     <- file.path(RESULTS_DIR, "R4_Results")

# 效应量阈值（初始值，看结果后可调整重跑）
RHO_TIER1 <- 0.10   # strict
RHO_TIER2 <- 0.05   # moderate

# 显著性阈值
Q_TIER1 <- 0.01
Q_TIER2 <- 0.05

cat("============================================\n")
cat("R4: Density Gene 筛选与分级\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("阈值: tier1 = q<%.2f & |rho|>=%.2f & (method_robust|high_confidence)\n", Q_TIER1, RHO_TIER1))
cat(sprintf("       tier2 = q<%.2f & |rho|>=%.2f\n", Q_TIER2, RHO_TIER2))
cat(sprintf("       tier3 = q<%.2f\n", Q_TIER2))
cat("\n")

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

  # --- 确定文件路径 ---
  parts <- strsplit(sample_name, "/")[[1]]
  dataset_name <- parts[1]
  sample_subname <- parts[2]

  r3_path <- file.path(R3_DIR, dataset_name, sample_subname, "density_gene_correlations.csv")

  if (!file.exists(r3_path)) {
    cat("  ✗ R3 结果不存在:", r3_path, "\n\n")
    next
  }

  # --- 创建输出目录 ---
  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------
  # 第一步：读取 R3 结果
  # ---------------------------------------------------------------

  cat("  [1/5] 读取 R3 结果 ...\n")
  df <- fread(r3_path)
  n_genes <- nrow(df)
  cat(sprintf("        %d 基因\n", n_genes))

  # ---------------------------------------------------------------
  # 第二步：描述统计
  # ---------------------------------------------------------------

  cat("  [2/5] |ρ| 分布描述统计 ...\n")

  abs_rho <- abs(df$rho_knn_main)

  # 基本统计
  rho_stats <- list(
    min    = round(min(abs_rho, na.rm = TRUE), 4),
    q25    = round(quantile(abs_rho, 0.25, na.rm = TRUE), 4),
    median = round(median(abs_rho, na.rm = TRUE), 4),
    mean   = round(mean(abs_rho, na.rm = TRUE), 4),
    q75    = round(quantile(abs_rho, 0.75, na.rm = TRUE), 4),
    max    = round(max(abs_rho, na.rm = TRUE), 4)
  )

  # 百分位阈值（从大到小排列时，top X% 的阈值）
  pct_thresholds <- list(
    top_5  = round(quantile(abs_rho, 0.95, na.rm = TRUE), 4),
    top_10 = round(quantile(abs_rho, 0.90, na.rm = TRUE), 4),
    top_20 = round(quantile(abs_rho, 0.80, na.rm = TRUE), 4),
    top_25 = round(quantile(abs_rho, 0.75, na.rm = TRUE), 4)
  )

  cat(sprintf("        min=%.4f  Q1=%.4f  median=%.4f  mean=%.4f  Q3=%.4f  max=%.4f\n",
              rho_stats$min, rho_stats$q25, rho_stats$median,
              rho_stats$mean, rho_stats$q75, rho_stats$max))
  cat(sprintf("        top5%%≥%.4f  top10%%≥%.4f  top20%%≥%.4f\n",
              pct_thresholds$top_5, pct_thresholds$top_10, pct_thresholds$top_20))

  # 正 vs 负相关的 |ρ| 统计
  pos_rho <- abs_rho[df$rho_knn_main > 0]
  neg_rho <- abs_rho[df$rho_knn_main < 0]
  cat(sprintf("        正相关基因: n=%d, median|ρ|=%.4f\n",
              length(pos_rho), round(median(pos_rho, na.rm = TRUE), 4)))
  cat(sprintf("        负相关基因: n=%d, median|ρ|=%.4f\n",
              length(neg_rho), round(median(neg_rho, na.rm = TRUE), 4)))

  # ---------------------------------------------------------------
  # 第三步：分级筛选
  # ---------------------------------------------------------------

  cat("  [3/5] 分级筛选 ...\n")

  # 绝对值列
  df[, abs_rho_main := abs(rho_knn_main)]

  # tier 标记
  df[, tier := "not_significant"]

  # tier3: q<0.05
  df[q_knn_main < Q_TIER2, tier := "tier3_lenient"]

  # tier2: q<0.05 & |rho|>=0.05
  df[q_knn_main < Q_TIER2 & abs_rho_main >= RHO_TIER2, tier := "tier2_moderate"]

  # tier1: q<0.01 & |rho|>=0.10 & (method_robust | high_confidence)
  df[q_knn_main < Q_TIER1 & abs_rho_main >= RHO_TIER1 &
     convergence %in% c("method_robust", "high_confidence"),
     tier := "tier1_strict"]

  # top百分位标记
  rho_rank <- rank(-abs_rho)  # 降序排名
  df[, top_10_pct := rho_rank <= ceiling(n_genes * 0.10)]
  df[, top_20_pct := rho_rank <= ceiling(n_genes * 0.20)]

  # 方向标记
  df[, direction := ifelse(rho_knn_main > 0, "positive", ifelse(rho_knn_main < 0, "negative", "zero"))]

  # 计数
  n_tier1 <- sum(df$tier == "tier1_strict")
  n_tier2 <- sum(df$tier == "tier2_moderate")
  n_tier3 <- sum(df$tier == "tier3_lenient")
  n_not_sig <- sum(df$tier == "not_significant")
  n_top10 <- sum(df$top_10_pct)
  n_top20 <- sum(df$top_20_pct)

  # tier1 中的正/负
  n_tier1_pos <- sum(df$tier == "tier1_strict" & df$direction == "positive")
  n_tier1_neg <- sum(df$tier == "tier1_strict" & df$direction == "negative")

  cat(sprintf("        tier1_strict:   %4d 基因 (正=%d, 负=%d)\n", n_tier1, n_tier1_pos, n_tier1_neg))
  cat(sprintf("        tier2_moderate: %4d 基因\n", n_tier2))
  cat(sprintf("        tier3_lenient:  %4d 基因\n", n_tier3))
  cat(sprintf("        not_significant:%4d 基因\n", n_not_sig))
  cat(sprintf("        top_10%%:        %4d 基因\n", n_top10))
  cat(sprintf("        top_20%%:        %4d 基因\n", n_top20))

  # ---------------------------------------------------------------
  # 第四步：输出
  # ---------------------------------------------------------------

  cat("  [4/5] 输出 ...\n")

  # 按 |ρ_knn_main| 降序排列
  df <- df[order(-abs_rho_main)]

  # 输出完整标记结果
  out_path <- file.path(out_dir, "filtered_density_genes.csv")
  fwrite(df, out_path)

  # Top 10 基因列表
  top10 <- head(df, 10)
  cat("        Top 10 density genes:\n")
  for (j in 1:min(10, nrow(top10))) {
    row <- top10[j]
    cat(sprintf("          %2d. %-15s ρ=%+.4f  q=%.2e  %s  %s\n",
                j, row$gene, row$rho_knn_main, row$q_knn_main,
                row$convergence, row$tier))
  }

  # ---------------------------------------------------------------
  # 第五步：诊断图
  # ---------------------------------------------------------------

  cat("  [5/5] 诊断图 ...\n")

  # --- 图1：|ρ| 直方图 ---
  tryCatch({
    p_hist <- ggplot(df, aes(x = abs_rho_main)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
      geom_vline(xintercept = RHO_TIER1, color = "red", linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = RHO_TIER2, color = "orange", linetype = "dashed", linewidth = 0.8) +
      annotate("text", x = RHO_TIER1 + 0.005, y = Inf, label = paste0("tier1: ", RHO_TIER1),
               color = "red", vjust = 2, hjust = 0, size = 3) +
      annotate("text", x = RHO_TIER2 + 0.005, y = Inf, label = paste0("tier2: ", RHO_TIER2),
               color = "orange", vjust = 3.5, hjust = 0, size = 3) +
      labs(title = paste0(sample_name, " — |ρ| Distribution"),
           x = "|ρ| (Spearman, KNN main)", y = "Gene count") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))

    ggsave(file.path(out_dir, "rho_distribution.png"), p_hist,
           width = 8, height = 5, dpi = 150)
  }, error = function(e) {
    cat(sprintf("        ⚠ 直方图生成失败: %s\n", e$message))
  })

  # --- 图2：火山图（按tier着色） ---
  tryCatch({
    # 准备数据
    plot_df <- data.frame(
      rho = df$rho_knn_main,
      neg_log_q = -log10(pmax(df$q_knn_main, 1e-300)),  # 防止log(0)
      tier = factor(df$tier, levels = c("tier1_strict", "tier2_moderate",
                                         "tier3_lenient", "not_significant")),
      gene = df$gene
    )

    # 标记top5用于label
    plot_df$label <- ""
    top5_idx <- head(order(-abs(plot_df$rho)), 5)
    plot_df$label[top5_idx] <- plot_df$gene[top5_idx]

    tier_colors <- c("tier1_strict" = "#D62728", "tier2_moderate" = "#FF7F0E",
                     "tier3_lenient" = "#AAAAAA", "not_significant" = "#DDDDDD")

    p_volcano_tier <- ggplot(plot_df, aes(x = rho, y = neg_log_q, color = tier)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = tier_colors, name = "Tier") +
      geom_vline(xintercept = c(-RHO_TIER1, RHO_TIER1), color = "red",
                 linetype = "dashed", linewidth = 0.3) +
      geom_vline(xintercept = c(-RHO_TIER2, RHO_TIER2), color = "orange",
                 linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = -log10(Q_TIER2), color = "grey50",
                 linetype = "dotted", linewidth = 0.3) +
      labs(title = paste0(sample_name, " — Volcano Plot (Tier)"),
           x = "ρ (Spearman)", y = "-log10(q)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            legend.position = "right")

    # 添加label（如果ggrepel可用）
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
    if (has_ggrepel && sum(plot_df$label != "") > 0) {
      p_volcano_tier <- p_volcano_tier +
        ggrepel::geom_text_repel(
          data = plot_df[plot_df$label != "", ],
          aes(label = label), size = 2.5, max.overlaps = 10,
          color = "black", segment.color = "grey50")
    }

    ggsave(file.path(out_dir, "volcano_tier.png"), p_volcano_tier,
           width = 9, height = 6, dpi = 150)
  }, error = function(e) {
    cat(sprintf("        ⚠ 火山图(tier)生成失败: %s\n", e$message))
  })

  # --- 图3：火山图（按convergence着色） ---
  tryCatch({
    plot_df2 <- data.frame(
      rho = df$rho_knn_main,
      neg_log_q = -log10(pmax(df$q_knn_main, 1e-300)),
      convergence = factor(df$convergence,
                           levels = c("method_robust", "high_confidence",
                                      "K_sensitive", "not_significant", "no_data")),
      gene = df$gene
    )

    conv_colors <- c("method_robust" = "#2CA02C", "high_confidence" = "#1F77B4",
                     "K_sensitive" = "#FF7F0E", "not_significant" = "#DDDDDD",
                     "no_data" = "#999999")

    p_volcano_conv <- ggplot(plot_df2, aes(x = rho, y = neg_log_q, color = convergence)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = conv_colors, name = "Convergence") +
      geom_vline(xintercept = c(-RHO_TIER1, RHO_TIER1), color = "red",
                 linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = -log10(Q_TIER2), color = "grey50",
                 linetype = "dotted", linewidth = 0.3) +
      labs(title = paste0(sample_name, " — Volcano Plot (Convergence)"),
           x = "ρ (Spearman)", y = "-log10(q)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            legend.position = "right")

    ggsave(file.path(out_dir, "volcano_convergence.png"), p_volcano_conv,
           width = 9, height = 6, dpi = 150)
  }, error = function(e) {
    cat(sprintf("        ⚠ 火山图(convergence)生成失败: %s\n", e$message))
  })

  # --- 图4：正 vs 负相关 |ρ| 箱线图 ---
  tryCatch({
    box_df <- data.frame(
      abs_rho = df$abs_rho_main,
      direction = df$direction
    )
    box_df <- box_df[box_df$direction != "zero", ]

    p_box <- ggplot(box_df, aes(x = direction, y = abs_rho, fill = direction)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      scale_fill_manual(values = c("positive" = "#D62728", "negative" = "#1F77B4")) +
      labs(title = paste0(sample_name, " — |ρ| by Direction"),
           x = "Direction", y = "|ρ|") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            legend.position = "none")

    ggsave(file.path(out_dir, "rho_by_direction.png"), p_box,
           width = 6, height = 5, dpi = 150)
  }, error = function(e) {
    cat(sprintf("        ⚠ 箱线图生成失败: %s\n", e$message))
  })

  cat("        ✓ 诊断图完成\n")

  cat("\n")

  # ---------------------------------------------------------------
  # 收集汇总
  # ---------------------------------------------------------------

  # tier1的top基因
  tier1_genes <- df[df$tier == "tier1_strict", ]
  top1_gene <- ifelse(nrow(tier1_genes) > 0, tier1_genes$gene[1], NA)
  top1_rho <- ifelse(nrow(tier1_genes) > 0, round(tier1_genes$rho_knn_main[1], 4), NA)

  summary_list[[sample_name]] <- data.frame(
    sample_name = sample_name,
    dataset = dataset_name,
    species = ifelse(!is.null(sample_info$species), sample_info$species, NA),
    condition = ifelse(!is.null(sample_info$condition), sample_info$condition, NA),
    n_genes = n_genes,
    n_tier1 = n_tier1,
    n_tier1_pos = n_tier1_pos,
    n_tier1_neg = n_tier1_neg,
    n_tier2 = n_tier2,
    n_tier3 = n_tier3,
    n_not_sig = n_not_sig,
    n_top10pct = n_top10,
    n_top20pct = n_top20,
    median_abs_rho = as.numeric(rho_stats$median),
    max_abs_rho = as.numeric(rho_stats$max),
    top10pct_threshold = as.numeric(pct_thresholds$top_10),
    top20pct_threshold = as.numeric(pct_thresholds$top_20),
    top1_tier1_gene = top1_gene,
    top1_tier1_rho = top1_rho,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# 全局汇总
# ============================================================================

cat("============================================\n")
cat("全局汇总\n")
cat("============================================\n")

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

summary_path <- file.path(OUTPUT_DIR, "ALL_SAMPLES_R4_SUMMARY.csv")
fwrite(summary_df, summary_path)
cat("汇总已保存:", summary_path, "\n\n")

# 打印汇总表
cat("--- 22 样本 R4 汇总 ---\n")
cat(sprintf("%-45s %5s %5s %5s %5s %5s %5s %7s %7s %8s\n",
            "sample", "genes", "tier1", "t1(+)", "t1(-)", "tier2", "tier3",
            "med|ρ|", "max|ρ|", "top1"))
cat(paste(rep("-", 130), collapse = ""), "\n")

for (j in 1:nrow(summary_df)) {
  row <- summary_df[j, ]
  top1_str <- ifelse(is.na(row$top1_tier1_gene), "-",
                     sprintf("%s(%.3f)", row$top1_tier1_gene, row$top1_tier1_rho))
  cat(sprintf("%-45s %5d %5d %5d %5d %5d %5d %7.4f %7.4f %8s\n",
              row$sample_name,
              row$n_genes,
              row$n_tier1,
              row$n_tier1_pos,
              row$n_tier1_neg,
              row$n_tier2,
              row$n_tier3,
              row$median_abs_rho,
              row$max_abs_rho,
              top1_str))
}

# 按数据集分组
cat("\n--- 按数据集分组 ---\n")
for (ds in unique(summary_df$dataset)) {
  sub <- summary_df[summary_df$dataset == ds, ]
  cat(sprintf("  %-30s %2d 样本, tier1: %d-%d, tier2: %d-%d, median|ρ|: %.4f-%.4f\n",
              ds,
              nrow(sub),
              min(sub$n_tier1), max(sub$n_tier1),
              min(sub$n_tier2), max(sub$n_tier2),
              min(sub$median_abs_rho), max(sub$median_abs_rho)))
}

# 阈值合理性检查
cat("\n--- 阈值合理性检查 ---\n")
for (j in 1:nrow(summary_df)) {
  row <- summary_df[j, ]
  if (row$n_tier1 == 0) {
    cat(sprintf("  ⚠ %s: tier1=0 基因! 考虑降低 RHO_TIER1 (当前=%.2f, top10%%≥%.4f)\n",
                row$sample_name, RHO_TIER1, row$top10pct_threshold))
  } else if (row$n_tier1 > row$n_genes * 0.5) {
    cat(sprintf("  ⚠ %s: tier1=%d (>50%%基因)! 考虑提高 RHO_TIER1\n",
                row$sample_name, row$n_tier1))
  }
}

cat("\n============================================\n")
cat("R4 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
