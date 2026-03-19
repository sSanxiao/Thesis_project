############################################################
# R4_threshold_sensitivity.R
# 功能: 对 R3 选定的 KNN 方法, 测试不同 Spearman ρ 阈值
#       的敏感性, 确定最终效应量阈值
# 输入: density_results_KNN.csv (R3 输出, 每样本一份)
# 输出: 阈值决策表 + 诊断图
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R4_threshold_sensitivity.R > R4_run.log 2>&1 &
############################################################

library(dplyr)
library(ggplot2)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"
FDR_THRESHOLD <- 0.05
SELECTED_METHOD <- "KNN"  # R3 选定

# 候选阈值
COR_CANDIDATES <- c(0.00, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20)

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
# Jaccard 函数
# ===========================================================

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# ===========================================================
# 读取所有样本的 KNN 结果
# ===========================================================

cat("============================================================\n")
cat("  R4: SPEARMAN THRESHOLD SENSITIVITY ANALYSIS\n")
cat("  Method:", SELECTED_METHOD, "\n")
cat("  Thresholds:", paste(COR_CANDIDATES, collapse = ", "), "\n")
cat("============================================================\n\n")

all_results <- list()

for (s in sample_list) {
  label    <- paste0(s$project, "/", s$sample)
  res_file <- file.path(RESULTS_ROOT, s$project, s$sample,
                        paste0("density_results_", SELECTED_METHOD, ".csv"))

  if (!file.exists(res_file)) {
    cat("  [SKIP]", label, "- file not found\n")
    next
  }

  df <- read.csv(res_file, stringsAsFactors = FALSE)
  all_results[[label]] <- list(
    data      = df,
    species   = s$species,
    condition = s$condition
  )
}

cat("  Loaded", length(all_results), "samples\n\n")

# ===========================================================
# 对每个阈值, 在每个样本上筛选基因
# ===========================================================

# 存储: threshold × sample × gene list
threshold_data <- list()

for (thr in COR_CANDIDATES) {
  thr_label <- as.character(thr)
  threshold_data[[thr_label]] <- list()

  for (label in names(all_results)) {
    df <- all_results[[label]]$data
    sig <- df$gene[df$FDR < FDR_THRESHOLD & abs(df$spearman_cor) > thr]
    threshold_data[[thr_label]][[label]] <- sig
  }
}

# ===========================================================
# 分析 1: 每个阈值下的基因数量统计
# ===========================================================

gene_count_rows <- list()

for (thr in COR_CANDIDATES) {
  thr_label <- as.character(thr)
  for (label in names(all_results)) {
    n_sig <- length(threshold_data[[thr_label]][[label]])
    gene_count_rows[[length(gene_count_rows) + 1]] <- data.frame(
      threshold = thr,
      sample    = label,
      species   = all_results[[label]]$species,
      condition = all_results[[label]]$condition,
      n_sig     = n_sig,
      stringsAsFactors = FALSE
    )
  }
}

gc_df <- do.call(rbind, gene_count_rows)

# 汇总: 每阈值的平均/中位基因数
thr_summary <- gc_df %>%
  group_by(threshold) %>%
  summarise(
    mean_n_sig   = round(mean(n_sig), 1),
    median_n_sig = median(n_sig),
    sd_n_sig     = round(sd(n_sig), 1),
    min_n_sig    = min(n_sig),
    max_n_sig    = max(n_sig),
    .groups = "drop"
  )

cat("  Gene count by threshold:\n")
print(as.data.frame(thr_summary))

# ===========================================================
# 分析 2: 相邻阈值之间的 Jaccard 稳定性
# ===========================================================

# 对每对相邻阈值, 跨样本平均 Jaccard
adjacent_jac <- list()

for (i in seq_along(COR_CANDIDATES)[-1]) {
  thr_a <- as.character(COR_CANDIDATES[i - 1])
  thr_b <- as.character(COR_CANDIDATES[i])

  jac_vals <- c()
  for (label in names(all_results)) {
    genes_a <- threshold_data[[thr_a]][[label]]
    genes_b <- threshold_data[[thr_b]][[label]]
    jac_vals <- c(jac_vals, jaccard(genes_a, genes_b))
  }

  adjacent_jac[[length(adjacent_jac) + 1]] <- data.frame(
    pair         = paste0(thr_a, " → ", thr_b),
    thr_from     = COR_CANDIDATES[i - 1],
    thr_to       = COR_CANDIDATES[i],
    mean_jaccard = round(mean(jac_vals, na.rm = TRUE), 4),
    sd_jaccard   = round(sd(jac_vals, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}

jac_df <- do.call(rbind, adjacent_jac)

cat("\n  Adjacent threshold Jaccard stability:\n")
print(jac_df)

# ===========================================================
# 分析 3: top 50 基因在不同阈值下的保留率
# ===========================================================

# 用 threshold=0 (只有 FDR 筛选) 的 top 50 作为参考
top50_retention <- list()

for (label in names(all_results)) {
  df <- all_results[[label]]$data
  # 按 |spearman_cor| 降序取 top 50
  df$abs_cor <- abs(df$spearman_cor)
  top50 <- head(df$gene[order(-df$abs_cor)], 50)

  for (thr in COR_CANDIDATES) {
    thr_label <- as.character(thr)
    sig_at_thr <- threshold_data[[thr_label]][[label]]
    retained <- sum(top50 %in% sig_at_thr)

    top50_retention[[length(top50_retention) + 1]] <- data.frame(
      threshold   = thr,
      sample      = label,
      n_retained  = retained,
      pct_retained = round(retained / 50 * 100, 1),
      stringsAsFactors = FALSE
    )
  }
}

ret_df <- do.call(rbind, top50_retention)

ret_summary <- ret_df %>%
  group_by(threshold) %>%
  summarise(
    mean_pct_retained = round(mean(pct_retained), 1),
    min_pct_retained  = min(pct_retained),
    .groups = "drop"
  )

cat("\n  Top-50 gene retention by threshold:\n")
print(as.data.frame(ret_summary))

# ===========================================================
# 分析 4: 跨样本可重复性 (每阈值下, 在 ≥30% 样本中显著的基因数)
# ===========================================================

repro_rows <- list()
n_samples <- length(all_results)

for (thr in COR_CANDIDATES) {
  thr_label <- as.character(thr)

  # 统计每个基因在多少样本中显著
  all_genes <- unique(unlist(threshold_data[[thr_label]]))
  gene_counts <- sapply(all_genes, function(g) {
    sum(sapply(threshold_data[[thr_label]], function(gl) g %in% gl))
  })

  n_repro_30pct <- sum(gene_counts >= ceiling(n_samples * 0.3))
  n_repro_50pct <- sum(gene_counts >= ceiling(n_samples * 0.5))

  repro_rows[[length(repro_rows) + 1]] <- data.frame(
    threshold     = thr,
    n_total_union = length(all_genes),
    n_repro_30pct = n_repro_30pct,
    n_repro_50pct = n_repro_50pct,
    stringsAsFactors = FALSE
  )
}

repro_df <- do.call(rbind, repro_rows)

cat("\n  Cross-sample reproducibility by threshold:\n")
print(repro_df)

# ===========================================================
# 决策逻辑
# ===========================================================

# 选择标准:
# 1. 排除 top50 保留率 < 90% 的阈值 (太严格, 丢了核心基因)
# 2. 在剩余阈值中, 选 Jaccard 稳定性开始下降的拐点
#    (即相邻 Jaccard 首次 < 0.85 之前的阈值)
# 3. 如果多个阈值满足, 选 n_repro_30pct 最大的

valid_thr <- ret_summary$threshold[ret_summary$mean_pct_retained >= 90]
cat("\n  Thresholds retaining ≥90% of top-50:", paste(valid_thr, collapse = ", "), "\n")

if (length(valid_thr) == 0) {
  selected_thr <- 0.05
  reason <- "Fallback: no threshold retains >=90% top-50, using 0.05"
} else {
  # 找 Jaccard 稳定区间的右边界
  stable_thr <- valid_thr
  for (i in seq_len(nrow(jac_df))) {
    if (jac_df$mean_jaccard[i] < 0.85) {
      # 这个阈值对之后就不稳定了
      cutoff <- jac_df$thr_from[i]
      stable_thr <- valid_thr[valid_thr <= cutoff]
      break
    }
  }

  if (length(stable_thr) == 0) stable_thr <- valid_thr

  # 在稳定区间内选 repro 最大的
  repro_sub <- repro_df[repro_df$threshold %in% stable_thr, ]
  best_idx  <- which.max(repro_sub$n_repro_30pct)
  selected_thr <- repro_sub$threshold[best_idx]
  reason <- paste0(
    "In stable Jaccard zone (", paste(stable_thr, collapse=","),
    "), threshold=", selected_thr,
    " maximizes cross-sample reproducibility (",
    repro_sub$n_repro_30pct[best_idx], " genes in >=30% samples)."
  )
}

cat("\n  >>> SELECTED THRESHOLD:", selected_thr, "\n")
cat("  >>> Reason:", reason, "\n")

# ===========================================================
# 保存结果
# ===========================================================

qc_dir <- file.path(RESULTS_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(gc_df,
          file.path(qc_dir, "R4_gene_count_by_threshold.csv"),
          row.names = FALSE)
write.csv(as.data.frame(thr_summary),
          file.path(qc_dir, "R4_threshold_summary.csv"),
          row.names = FALSE)
write.csv(jac_df,
          file.path(qc_dir, "R4_adjacent_jaccard.csv"),
          row.names = FALSE)
write.csv(as.data.frame(ret_summary),
          file.path(qc_dir, "R4_top50_retention.csv"),
          row.names = FALSE)
write.csv(repro_df,
          file.path(qc_dir, "R4_reproducibility_by_threshold.csv"),
          row.names = FALSE)

decision <- data.frame(
  selected_threshold = selected_thr,
  reason             = reason,
  fdr_threshold      = FDR_THRESHOLD,
  method             = SELECTED_METHOD,
  stringsAsFactors   = FALSE
)
write.csv(decision,
          file.path(qc_dir, "R4_THRESHOLD_SELECTION.csv"),
          row.names = FALSE)

# ===========================================================
# 诊断图
# ===========================================================

# 图1: 阈值 vs 平均基因数 + 跨样本 SD
p1 <- ggplot(thr_summary, aes(x = threshold, y = mean_n_sig)) +
  geom_ribbon(aes(ymin = mean_n_sig - sd_n_sig,
                  ymax = mean_n_sig + sd_n_sig), alpha = 0.2, fill = "#2c7bb6") +
  geom_line(color = "#2c7bb6", linewidth = 1) +
  geom_point(color = "#2c7bb6", size = 3) +
  geom_vline(xintercept = selected_thr, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = selected_thr + 0.01, y = max(thr_summary$mean_n_sig) * 0.95,
           label = paste0("Selected: ", selected_thr), color = "red", hjust = 0, size = 4) +
  labs(title = "Mean significant genes vs Spearman threshold",
       subtitle = "Blue band = ±1 SD across 21 samples",
       x = "|Spearman ρ| threshold", y = "Mean sig. genes") +
  theme_minimal(base_size = 13) +
  scale_x_continuous(breaks = COR_CANDIDATES)

ggsave(file.path(qc_dir, "R4_gene_count_vs_threshold.png"),
       p1, width = 8, height = 5, dpi = 200)

# 图2: 相邻阈值 Jaccard
p2 <- ggplot(jac_df, aes(x = thr_to, y = mean_jaccard)) +
  geom_line(color = "#d7191c", linewidth = 1) +
  geom_point(color = "#d7191c", size = 3) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(COR_CANDIDATES) - 0.02, y = 0.86,
           label = "Stability threshold (0.85)", color = "gray50", size = 3.5) +
  geom_vline(xintercept = selected_thr, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Gene list stability between adjacent thresholds",
       x = "Threshold (to)", y = "Mean Jaccard with previous threshold") +
  theme_minimal(base_size = 13) +
  scale_x_continuous(breaks = COR_CANDIDATES) +
  ylim(0, 1)

ggsave(file.path(qc_dir, "R4_jaccard_stability.png"),
       p2, width = 8, height = 5, dpi = 200)

# 图3: Top-50 保留率
p3 <- ggplot(ret_summary, aes(x = threshold, y = mean_pct_retained)) +
  geom_line(color = "#1a9850", linewidth = 1) +
  geom_point(color = "#1a9850", size = 3) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = selected_thr, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Top-50 gene retention rate by threshold",
       x = "|Spearman ρ| threshold", y = "Mean % retained") +
  theme_minimal(base_size = 13) +
  scale_x_continuous(breaks = COR_CANDIDATES) +
  ylim(0, 105)

ggsave(file.path(qc_dir, "R4_top50_retention.png"),
       p3, width = 8, height = 5, dpi = 200)

# 图4: 跨样本可重复基因数
repro_long <- repro_df %>%
  tidyr::pivot_longer(cols = c(n_repro_30pct, n_repro_50pct),
                      names_to = "level", values_to = "n_genes") %>%
  mutate(level = ifelse(level == "n_repro_30pct", "≥30% samples", "≥50% samples"))

p4 <- ggplot(repro_long, aes(x = threshold, y = n_genes, color = level)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = selected_thr, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Reproducible gene count by threshold",
       x = "|Spearman ρ| threshold", y = "Number of genes", color = "Reproducibility") +
  theme_minimal(base_size = 13) +
  scale_x_continuous(breaks = COR_CANDIDATES) +
  scale_color_manual(values = c("≥30% samples" = "#2c7bb6", "≥50% samples" = "#d7191c"))

ggsave(file.path(qc_dir, "R4_reproducibility_vs_threshold.png"),
       p4, width = 8, height = 5, dpi = 200)

cat("\n============================================================\n")
cat("  R4 COMPLETE\n")
cat("  Selected threshold:", selected_thr, "\n")
cat("  Plots saved to:", qc_dir, "\n")
cat("============================================================\n")
