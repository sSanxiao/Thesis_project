############################################################
# R8_cross_sample_summary.R
# 功能: 跨样本汇总 (Stage 2.5 Task 3)
#   1. Gene × sample Spearman 相关系数矩阵 + FDR 矩阵
#   2. 可重复 density-associated genes 识别
#   3. Human vs Mouse 比较
#   4. Per-condition 汇总
#   5. 层次聚类热力图 (带 species/condition 注释)
#   6. Gene panel 覆盖矩阵 (跨数据集基因交集)
#   7. 整合 R5 tier 信息到可重复基因表
#
# 输入: density_results_KNN.csv (R3, 每样本)
#       gene_tier_classification.csv (R5, 每样本)
#       R6_all_coupling.csv
# 输出: 跨样本矩阵, 可重复基因表, 比较报告
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R8_cross_sample_summary.R > R8_run.log 2>&1 &
############################################################

library(dplyr)
library(pheatmap)
library(ggplot2)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT  <- "/home/disk/wangqilu/Stage2_new/Results"
FDR_THRESHOLD <- 0.05
COR_THRESHOLD <- 0.05

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

# 构建 label, species, condition 向量
labels <- sapply(sample_list, function(s) paste0(s$project, "/", s$sample))
species_vec <- sapply(sample_list, function(s) s$species)
cond_vec    <- sapply(sample_list, function(s) s$condition)
n_samples   <- length(sample_list)

cat("============================================================\n")
cat("  R8: CROSS-SAMPLE SUMMARY (Task 3)\n")
cat("  Samples:", n_samples, "\n")
cat("============================================================\n\n")

out_dir <- file.path(RESULTS_ROOT, "cross_sample_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================
# 1. 读取所有样本的 KNN 结果
# ===========================================================

cat("  Loading all density_results_KNN.csv...\n")
all_res <- list()
for (i in seq_along(sample_list)) {
  s <- sample_list[[i]]
  f <- file.path(RESULTS_ROOT, s$project, s$sample, "density_results_KNN.csv")
  if (file.exists(f)) {
    all_res[[labels[i]]] <- read.csv(f, stringsAsFactors = FALSE)
  } else {
    cat("    [SKIP]", labels[i], "- file not found\n")
  }
}
cat("  Loaded:", length(all_res), "samples\n\n")

# ===========================================================
# 2. 基因名标准化 (大写, 跨物种兼容)
# ===========================================================

cat("  Standardizing gene names to uppercase...\n")
for (lab in names(all_res)) {
  all_res[[lab]]$gene_original <- all_res[[lab]]$gene
  all_res[[lab]]$gene <- toupper(all_res[[lab]]$gene)
}

all_genes <- sort(unique(unlist(lapply(all_res, function(x) x$gene))))
cat("  Total unique genes (uppercase):", length(all_genes), "\n\n")

# ===========================================================
# 3. Gene × Sample 矩阵 (cor + FDR)
# ===========================================================

cat("  Building gene x sample matrices...\n")
cor_mat <- matrix(NA_real_, nrow = length(all_genes), ncol = n_samples,
                  dimnames = list(all_genes, names(all_res)))
fdr_mat <- cor_mat

for (j in seq_along(all_res)) {
  lab <- names(all_res)[j]
  df  <- all_res[[lab]]
  idx <- match(df$gene, all_genes)
  cor_mat[idx, j] <- df$spearman_cor
  fdr_mat[idx, j] <- df$FDR
}

write.csv(cor_mat, file.path(out_dir, "gene_sample_correlation_matrix.csv"))
write.csv(fdr_mat, file.path(out_dir, "gene_sample_fdr_matrix.csv"))

# ===========================================================
# 4. Gene panel 覆盖矩阵
# ===========================================================

cat("  Building gene panel presence matrix...\n")
presence_mat <- matrix(FALSE, nrow = length(all_genes), ncol = n_samples,
                       dimnames = list(all_genes, names(all_res)))
for (j in seq_along(all_res)) {
  presence_mat[match(all_res[[j]]$gene, all_genes), j] <- TRUE
}
write.csv(presence_mat, file.path(out_dir, "gene_panel_presence_matrix.csv"))

# 基因在多少样本的 panel 中
gene_coverage <- rowSums(presence_mat)
cat("  Genes in all 21 samples:", sum(gene_coverage == n_samples), "\n")
cat("  Genes in >=10 samples:", sum(gene_coverage >= 10), "\n")
cat("  Genes in only 1 sample:", sum(gene_coverage == 1), "\n\n")

# ===========================================================
# 5. 可重复 density genes
# ===========================================================

cat("  Identifying reproducible density genes...\n")
sig_mat <- (!is.na(fdr_mat)) & (fdr_mat < FDR_THRESHOLD) &
           (abs(cor_mat) > COR_THRESHOLD)
pos_mat <- sig_mat & (cor_mat > 0)
neg_mat <- sig_mat & (cor_mat < 0)

repro <- data.frame(
  gene         = all_genes,
  n_tested     = rowSums(!is.na(cor_mat)),
  n_sig        = rowSums(sig_mat, na.rm = TRUE),
  n_pos        = rowSums(pos_mat, na.rm = TRUE),
  n_neg        = rowSums(neg_mat, na.rm = TRUE),
  mean_cor     = rowMeans(cor_mat, na.rm = TRUE),
  median_cor   = apply(cor_mat, 1, median, na.rm = TRUE),
  sd_cor       = apply(cor_mat, 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# 方向一致性: 在显著样本中, 方向最多的那个占比
repro$dir_consistency <- ifelse(
  repro$n_sig > 0,
  pmax(repro$n_pos, repro$n_neg) / repro$n_sig,
  NA
)

repro$pct_sig <- round(repro$n_sig / repro$n_tested * 100, 1)

# ===========================================================
# 6. 整合 R5 tier 信息
# ===========================================================

cat("  Integrating R5 tier information...\n")
tier_counts <- data.frame(gene = all_genes,
                          tier_strong_count = 0L,
                          tier_moderate_count = 0L,
                          tier_weak_count = 0L,
                          stringsAsFactors = FALSE)

for (s in sample_list) {
  tf <- file.path(RESULTS_ROOT, s$project, s$sample, "gene_tier_classification.csv")
  if (!file.exists(tf)) next
  td <- read.csv(tf, stringsAsFactors = FALSE)
  td$gene_upper <- toupper(td$gene)
  for (tier in c("strong", "moderate", "weak")) {
    genes_in_tier <- td$gene_upper[td$tier == tier]
    idx <- match(genes_in_tier, all_genes)
    idx <- idx[!is.na(idx)]
    col_name <- paste0("tier_", tier, "_count")
    tier_counts[[col_name]][idx] <- tier_counts[[col_name]][idx] + 1L
  }
}

repro <- merge(repro, tier_counts, by = "gene", all.x = TRUE)
repro <- repro[order(-repro$n_sig, -abs(repro$mean_cor)), ]

write.csv(repro, file.path(out_dir, "reproducible_density_genes.csv"),
          row.names = FALSE)

cat("  Genes sig in >=3 samples:", sum(repro$n_sig >= 3), "\n")
cat("  Genes sig in >=5 samples:", sum(repro$n_sig >= 5), "\n")
cat("  Genes sig in >=10 samples:", sum(repro$n_sig >= 10), "\n")
cat("  Genes sig in >=50% tested:", sum(repro$pct_sig >= 50, na.rm = TRUE), "\n\n")

# ===========================================================
# 7. Human vs Mouse 比较
# ===========================================================

cat("  Human vs Mouse comparison...\n")
h_idx <- which(species_vec == "human")
m_idx <- which(species_vec == "mouse")
h_labels <- names(all_res)[h_idx]
m_labels <- names(all_res)[m_idx]

hmc <- data.frame(
  gene           = all_genes,
  human_mean_cor = rowMeans(cor_mat[, h_labels, drop = FALSE], na.rm = TRUE),
  mouse_mean_cor = rowMeans(cor_mat[, m_labels, drop = FALSE], na.rm = TRUE),
  human_n_sig    = rowSums(sig_mat[, h_labels, drop = FALSE], na.rm = TRUE),
  mouse_n_sig    = rowSums(sig_mat[, m_labels, drop = FALSE], na.rm = TRUE),
  human_n_tested = rowSums(!is.na(cor_mat[, h_labels, drop = FALSE])),
  mouse_n_tested = rowSums(!is.na(cor_mat[, m_labels, drop = FALSE])),
  stringsAsFactors = FALSE
)

hmc$both_species   <- (hmc$human_n_sig > 0) & (hmc$mouse_n_sig > 0)
hmc$human_only     <- (hmc$human_n_sig > 0) & (hmc$mouse_n_sig == 0)
hmc$mouse_only     <- (hmc$human_n_sig == 0) & (hmc$mouse_n_sig > 0)
hmc$neither        <- (hmc$human_n_sig == 0) & (hmc$mouse_n_sig == 0)

# 方向一致性检查 (跨物种)
hmc$human_direction <- ifelse(hmc$human_mean_cor > 0, "positive",
                       ifelse(hmc$human_mean_cor < 0, "negative", "none"))
hmc$mouse_direction <- ifelse(hmc$mouse_mean_cor > 0, "positive",
                       ifelse(hmc$mouse_mean_cor < 0, "negative", "none"))
hmc$cross_species_consistent <- (hmc$both_species) &
                                (hmc$human_direction == hmc$mouse_direction)

hmc <- hmc[order(-hmc$both_species,
                 -(hmc$human_n_sig + hmc$mouse_n_sig),
                 -abs(hmc$human_mean_cor + hmc$mouse_mean_cor)), ]

write.csv(hmc, file.path(out_dir, "human_mouse_comparison.csv"),
          row.names = FALSE)

cat("  Genes sig in both species:", sum(hmc$both_species), "\n")
cat("    Direction consistent:", sum(hmc$cross_species_consistent, na.rm = TRUE), "\n")
cat("  Human-only:", sum(hmc$human_only), "\n")
cat("  Mouse-only:", sum(hmc$mouse_only), "\n\n")

# ===========================================================
# 8. Per-condition 汇总
# ===========================================================

cat("  Per-condition summary...\n")
unique_conds <- unique(cond_vec)
cond_summary_list <- list()

for (cd in unique_conds) {
  cd_idx    <- which(cond_vec == cd)
  cd_labels <- names(all_res)[cd_idx]
  cd_labels <- intersect(cd_labels, names(all_res))
  if (length(cd_labels) == 0) next

  cond_summary_list[[cd]] <- data.frame(
    gene      = all_genes,
    condition = cd,
    n_samples = length(cd_labels),
    mean_cor  = rowMeans(cor_mat[, cd_labels, drop = FALSE], na.rm = TRUE),
    sd_cor    = apply(cor_mat[, cd_labels, drop = FALSE], 1, sd, na.rm = TRUE),
    n_sig     = rowSums(sig_mat[, cd_labels, drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cond_df <- do.call(rbind, cond_summary_list)
write.csv(cond_df, file.path(out_dir, "per_condition_density_summary.csv"),
          row.names = FALSE)

# Per-condition top 5 genes
cat("\n  Top 5 genes per condition:\n")
for (cd in unique_conds) {
  sub <- cond_df[cond_df$condition == cd, ]
  sub <- sub[order(-abs(sub$mean_cor)), ]
  cat("  [", cd, "] ", paste(head(sub$gene, 5), collapse = ", "), "\n")
}

# ===========================================================
# 9. Sample registry
# ===========================================================

reg <- data.frame(
  sample_id = labels,
  project   = sapply(sample_list, function(s) s$project),
  sample    = sapply(sample_list, function(s) s$sample),
  species   = species_vec,
  condition = cond_vec,
  stringsAsFactors = FALSE
)
write.csv(reg, file.path(out_dir, "sample_registry.csv"), row.names = FALSE)

# ===========================================================
# 10. 层次聚类热力图
# ===========================================================

cat("\n  Generating heatmaps...\n")

# 取可重复基因子集 (sig in >=3 samples, 且 direction consistency > 0.8)
repro_genes <- repro$gene[repro$n_sig >= 3 &
                          !is.na(repro$dir_consistency) &
                          repro$dir_consistency >= 0.8]
cat("  Heatmap genes (sig >=3 & dir consistency >=0.8):", length(repro_genes), "\n")

if (length(repro_genes) >= 5 && length(repro_genes) <= 500) {
  hm_mat <- cor_mat[repro_genes, , drop = FALSE]

  # 替换 NA 为 0 (基因不在该样本 panel 中)
  hm_mat[is.na(hm_mat)] <- 0

  # 样本注释
  annot_col <- data.frame(
    species   = species_vec,
    condition = cond_vec,
    row.names = names(all_res)
  )

  # tier 注释 (行注释)
  tier_info <- repro[match(repro_genes, repro$gene), ]
  annot_row <- data.frame(
    n_sig        = tier_info$n_sig,
    n_strong     = tier_info$tier_strong_count,
    row.names    = repro_genes
  )

  # 颜色
  ann_colors <- list(
    species   = c(human = "#2c7bb6", mouse = "#d7191c"),
    condition = c(wild_type = "#1a9850", TgCRND8_AD = "#d7191c",
                  alzheimer = "#fdae61", glioblastoma = "#542788",
                  healthy_brain = "#66c2a5", normal_brain = "#3288bd",
                  ATRT = "#f46d43", medulloblastoma = "#e08214")
  )

  # 限制热力图大小
  n_show <- min(length(repro_genes), 100)
  hm_show <- hm_mat[seq_len(n_show), ]

  png(file.path(out_dir, "reproducible_genes_heatmap.png"),
      width = 1200, height = max(600, n_show * 12), res = 150)
  pheatmap(
    hm_show,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-0.3, 0.3, length.out = 101),
    annotation_col = annot_col,
    annotation_colors = ann_colors,
    show_rownames = (n_show <= 60),
    show_colnames = TRUE,
    fontsize = 7, fontsize_col = 6, fontsize_row = 6,
    main = paste0("Reproducible Density Genes (n=", n_show,
                  ", sig in >=3 samples, dir consistency >=0.8)"),
    na_col = "gray90"
  )
  dev.off()
  cat("  Heatmap saved.\n")

} else if (length(repro_genes) > 500) {
  cat("  [NOTE] Too many genes (", length(repro_genes),
      ") for heatmap, taking top 100 by n_sig\n")
  top100 <- head(repro$gene[repro$gene %in% repro_genes], 100)
  hm_mat <- cor_mat[top100, , drop = FALSE]
  hm_mat[is.na(hm_mat)] <- 0

  annot_col <- data.frame(
    species = species_vec, condition = cond_vec,
    row.names = names(all_res)
  )

  png(file.path(out_dir, "reproducible_genes_heatmap.png"),
      width = 1200, height = 1000, res = 150)
  pheatmap(
    hm_mat, cluster_rows = TRUE, cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-0.3, 0.3, length.out = 101),
    annotation_col = annot_col,
    show_rownames = TRUE, show_colnames = TRUE,
    fontsize = 7, fontsize_col = 6, fontsize_row = 5,
    main = "Top 100 Reproducible Density Genes",
    na_col = "gray90"
  )
  dev.off()
  cat("  Heatmap saved (top 100).\n")

} else {
  cat("  [NOTE] <5 reproducible genes, skipping heatmap\n")
}

# ===========================================================
# 11. Human vs Mouse 散点图
# ===========================================================

cat("  Human vs Mouse scatter plot...\n")
hmc_plot <- hmc[hmc$human_n_tested > 0 & hmc$mouse_n_tested > 0, ]
hmc_plot$category <- "Neither"
hmc_plot$category[hmc_plot$both_species] <- "Both species"
hmc_plot$category[hmc_plot$human_only]   <- "Human only"
hmc_plot$category[hmc_plot$mouse_only]   <- "Mouse only"
hmc_plot$category <- factor(hmc_plot$category,
                            levels = c("Both species", "Human only", "Mouse only", "Neither"))

# Label top genes
hmc_plot$label <- ""
top_both <- head(hmc_plot$gene[hmc_plot$category == "Both species"], 10)
hmc_plot$label[hmc_plot$gene %in% top_both] <- hmc_plot$gene[hmc_plot$gene %in% top_both]

p_hm <- ggplot(hmc_plot, aes(x = human_mean_cor, y = mouse_mean_cor, color = category)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_text(aes(label = label), size = 2.5, nudge_y = 0.005, color = "black",
            check_overlap = TRUE) +
  scale_color_manual(values = c("Both species" = "#d7191c",
                                "Human only" = "#2c7bb6",
                                "Mouse only" = "#1a9850",
                                "Neither" = "gray80")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray40") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold")) +
  labs(title = "Density-Gene Correlation: Human vs Mouse",
       x = "Mean Spearman ρ (Human samples)",
       y = "Mean Spearman ρ (Mouse samples)",
       color = "Significance")
ggsave(file.path(out_dir, "human_vs_mouse_scatter.png"),
       p_hm, width = 8, height = 7, dpi = 200)

# ===========================================================
# 12. 汇总统计
# ===========================================================

cat("\n============================================================\n")
cat("  CROSS-SAMPLE SUMMARY STATISTICS\n")
cat("============================================================\n\n")

cat("  Total unique genes:", length(all_genes), "\n")
cat("  Genes in all samples:", sum(gene_coverage == n_samples), "\n\n")

cat("  Reproducible density genes (global Spearman):\n")
cat("    Sig in >=3 samples:", sum(repro$n_sig >= 3), "\n")
cat("    Sig in >=5 samples:", sum(repro$n_sig >= 5), "\n")
cat("    Sig in >=10 samples:", sum(repro$n_sig >= 10), "\n\n")

cat("  Tier-based (from R5 multilevel analysis):\n")
cat("    Strong in >=3 samples:", sum(repro$tier_strong_count >= 3, na.rm = TRUE), "\n")
cat("    Moderate in >=3 samples:", sum(repro$tier_moderate_count >= 3, na.rm = TRUE), "\n\n")

cat("  Cross-species:\n")
cat("    Both human & mouse:", sum(hmc$both_species), "\n")
cat("    Direction consistent:", sum(hmc$cross_species_consistent, na.rm = TRUE), "\n")
cat("    Human-only:", sum(hmc$human_only), "\n")
cat("    Mouse-only:", sum(hmc$mouse_only), "\n\n")

cat("  Top 10 most reproducible genes:\n")
top10 <- head(repro, 10)
for (i in seq_len(nrow(top10))) {
  cat("    ", top10$gene[i],
      " (sig:", top10$n_sig[i], "/", top10$n_tested[i],
      ", mean_ρ=", round(top10$mean_cor[i], 3),
      ", strong:", top10$tier_strong_count[i],
      ", dir:", round(top10$dir_consistency[i], 2), ")\n")
}

cat("\n============================================================\n")
cat("  R8 COMPLETE\n")
cat("  Output:", out_dir, "\n")
cat("============================================================\n")
