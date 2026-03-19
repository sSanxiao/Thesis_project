############################################################
# R7_visualization.R
# 功能: 批量生成 Stage 2 的所有可视化
#   1. Density 空间图 (KNN)
#   2. Density Signature 空间图
#   3. Gene tier 空间图 (strong 基因示例)
#   4. Cluster 级别 density 小提琴图
#   5. Top 10 density 基因空间表达
#   6. Binned DE 火山图
#   7. 全局汇总: coupling 热力图, tier 分布柱状图
#
# 输入: seurat_final.rds (R6), gene_tier_classification.csv (R5),
#       density_results_KNN.csv (R3), binned_DE_results.csv (R5)
# 输出: PDF/PNG per sample + summary plots
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R7_visualization.R > R7_run.log 2>&1 &
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"

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
# 通用空间绘图函数
# ===========================================================

plot_spatial <- function(x, y, values, title, color_option = "magma",
                         point_size = 0.2, clip_pct = c(1, 99)) {
  valid <- !is.na(values)
  df <- data.frame(x = x[valid], y = y[valid], v = values[valid])

  if (nrow(df) == 0) return(ggplot() + labs(title = paste(title, "(no data)")))

  vmin <- quantile(df$v, clip_pct[1] / 100)
  vmax <- quantile(df$v, clip_pct[2] / 100)
  df$v_clip <- pmax(pmin(df$v, vmax), vmin)

  ggplot(df, aes(x = x, y = y, color = v_clip)) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_viridis_c(option = color_option) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = title, x = "X (μm)", y = "Y (μm)", color = "")
}

# ===========================================================
# 单样本可视化
# ===========================================================

visualize_one_sample <- function(s) {

  project <- s$project
  sname   <- s$sample
  r_label <- paste0(project, "/", sname)
  sample_dir <- file.path(RESULTS_ROOT, project, sname)
  plot_dir   <- file.path(sample_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  rds_path <- file.path(sample_dir, "seurat_final.rds")
  if (!file.exists(rds_path)) { cat("  [SKIP] no seurat_final.rds\n"); return() }

  cat("  Loading...\n")
  so <- readRDS(rds_path)

  x <- so$x_centroid
  y <- so$y_centroid
  expr_mat <- GetAssayData(so, assay = "SCT", layer = "data")

  # --- 1. Density 空间图 ---
  cat("  Plot 1: density spatial...\n")
  p1 <- plot_spatial(x, y, so$density_knn,
                     paste0(r_label, " — KNN Density"), "inferno")
  ggsave(file.path(plot_dir, "density_spatial.png"), p1, width = 8, height = 6, dpi = 150)

  # --- 2. Density Signature 空间图 ---
  if (!all(is.na(so$Density_Signature))) {
    cat("  Plot 2: signature spatial...\n")
    p2 <- plot_spatial(x, y, so$Density_Signature,
                       paste0(r_label, " — Density Signature"), "magma")
    ggsave(file.path(plot_dir, "signature_spatial.png"), p2, width = 8, height = 6, dpi = 150)
  }

  # --- 3. Cluster 级别 density 小提琴图 ---
  cat("  Plot 3: cluster violin...\n")
  vln_df <- data.frame(
    cluster = as.character(Idents(so)),
    density = so$density_knn
  )
  # 按中位密度排序 cluster
  cl_order <- vln_df %>% group_by(cluster) %>%
    summarise(med = median(density, na.rm = TRUE)) %>%
    arrange(med) %>% pull(cluster)
  vln_df$cluster <- factor(vln_df$cluster, levels = cl_order)

  p3 <- ggplot(vln_df, aes(x = cluster, y = density, fill = cluster)) +
    geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = paste0(r_label, " — Density by Cluster"),
         x = "Cluster", y = "KNN Density")
  ggsave(file.path(plot_dir, "cluster_density_violin.png"), p3, width = 8, height = 5, dpi = 150)

  # --- 4. Density 分布直方图 ---
  cat("  Plot 4: density histogram...\n")
  hist_df <- data.frame(d = so$density_knn[!is.na(so$density_knn)])
  clip_max <- quantile(hist_df$d, 0.99)
  p4 <- ggplot(hist_df[hist_df$d <= clip_max, ], aes(x = d)) +
    geom_histogram(bins = 60, fill = "#2c7bb6", color = "white", linewidth = 0.2) +
    theme_minimal(base_size = 11) +
    labs(title = paste0(r_label, " — Density Distribution"),
         x = "KNN Density", y = "Cell Count")
  ggsave(file.path(plot_dir, "density_histogram.png"), p4, width = 6, height = 4, dpi = 150)

  # --- 5. Top 10 density 基因空间表达 ---
  cat("  Plot 5: top 10 genes...\n")
  res_file <- file.path(sample_dir, "density_results_KNN.csv")
  if (file.exists(res_file)) {
    res <- read.csv(res_file, stringsAsFactors = FALSE)
    # Top 5 正相关 + Top 5 负相关
    top_pos <- head(res$gene[res$spearman_cor > 0], 5)
    top_neg <- head(res$gene[order(res$spearman_cor)], 5)
    top_genes <- unique(c(top_pos, top_neg))
    top_genes <- top_genes[top_genes %in% rownames(expr_mat)]

    if (length(top_genes) > 0) {
      pdf(file.path(plot_dir, "top_density_genes_spatial.pdf"), width = 8, height = 6)
      for (g in top_genes) {
        rho_val <- round(res$spearman_cor[res$gene == g], 3)
        p <- plot_spatial(x, y, as.numeric(expr_mat[g, ]),
                          paste0(g, " (ρ=", rho_val, ")"), "magma")
        print(p)
      }
      dev.off()
    }
  }

  # --- 6. Binned DE 火山图 ---
  cat("  Plot 6: binned DE volcano...\n")
  de_file <- file.path(sample_dir, "binned_DE_results.csv")
  if (file.exists(de_file)) {
    de <- read.csv(de_file, stringsAsFactors = FALSE)
    de$neg_log10_fdr <- -log10(pmax(de$FDR, 1e-300))
    de$sig_label <- "NS"
    de$sig_label[de$FDR < 0.05 & abs(de$log2FC) > 0.25] <- "Moderate"
    de$sig_label[de$FDR < 0.05 & abs(de$log2FC) > 1.0]  <- "Strong"
    de$sig_label <- factor(de$sig_label, levels = c("NS", "Moderate", "Strong"))

    # 标注 top 基因名
    de$label <- ""
    top_show <- head(de$gene[de$sig_label == "Strong"][order(-abs(de$log2FC[de$sig_label == "Strong"]))], 8)
    de$label[de$gene %in% top_show] <- de$gene[de$gene %in% top_show]

    p6 <- ggplot(de, aes(x = log2FC, y = neg_log10_fdr, color = sig_label)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("NS" = "gray70", "Moderate" = "#2c7bb6", "Strong" = "#d7191c")) +
      geom_text(aes(label = label), size = 2.5, nudge_y = 1, color = "black", check_overlap = TRUE) +
      geom_vline(xintercept = c(-1, -0.25, 0.25, 1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(size = 12, face = "bold")) +
      labs(title = paste0(r_label, " — High vs Low Density DE"),
           x = "log2 Fold Change (Q4 vs Q1)", y = "-log10(FDR)", color = "Tier")
    ggsave(file.path(plot_dir, "binned_DE_volcano.png"), p6, width = 8, height = 6, dpi = 150)
  }

  # --- 7. Gene tier 空间图 (top 3 strong 基因) ---
  cat("  Plot 7: strong gene spatial...\n")
  tier_file <- file.path(sample_dir, "gene_tier_classification.csv")
  if (file.exists(tier_file)) {
    tier <- read.csv(tier_file, stringsAsFactors = FALSE)
    strong_genes <- tier$gene[tier$tier == "strong"]
    strong_genes <- strong_genes[strong_genes %in% rownames(expr_mat)]

    if (length(strong_genes) >= 1) {
      show_genes <- head(strong_genes, 3)
      pdf(file.path(plot_dir, "strong_tier_genes_spatial.pdf"), width = 8, height = 6)
      for (g in show_genes) {
        fc_val <- round(tier$binned_log2FC[tier$gene == g], 2)
        rho_val <- round(tier$global_rho[tier$gene == g], 3)
        p <- plot_spatial(x, y, as.numeric(expr_mat[g, ]),
                          paste0(g, " [STRONG] (ρ=", rho_val, ", FC=", fc_val, ")"), "magma")
        print(p)
      }
      dev.off()
    }
  }

  cat("  Plots saved to:", plot_dir, "\n")
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R7: VISUALIZATION\n")
cat("  Samples:", length(sample_list), "\n")
cat("============================================================\n\n")

t_total <- Sys.time()

for (i in seq_along(sample_list)) {
  s     <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("\n[%d/%d] %s\n", i, length(sample_list), label))
  visualize_one_sample(s)
}

# ===========================================================
# 全局汇总图
# ===========================================================

cat("\n============================================================\n")
cat("  GLOBAL SUMMARY PLOTS\n")
cat("============================================================\n")

summary_dir <- file.path(RESULTS_ROOT, "QC", "summary_plots")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# --- A. Coupling 热力图 ---
cat("  Coupling heatmap...\n")
coupling_file <- file.path(RESULTS_ROOT, "QC", "R6_all_coupling.csv")
if (file.exists(coupling_file)) {
  cp <- read.csv(coupling_file, stringsAsFactors = FALSE)
  rownames(cp) <- paste0(cp$condition, ": ", cp$sample)

  # 只取数值列
  num_cols <- c("rho_PC1", "rho_PC2", "rho_cluster_mean", "rho_signature")
  # 加入 cell cycle / inflammation (如果有非 NaN 值)
  if (any(!is.na(cp$rho_S_Score))) num_cols <- c("rho_S_Score", "rho_G2M_Score", num_cols)
  if (any(!is.na(cp$rho_inflammation))) num_cols <- c("rho_inflammation", num_cols)

  cp_mat <- as.matrix(cp[, num_cols])
  rownames(cp_mat) <- rownames(cp)

  # 注释条: species + condition
  annot <- data.frame(
    species   = cp$species,
    condition = cp$condition,
    row.names = rownames(cp)
  )

  png(file.path(summary_dir, "coupling_heatmap.png"), width = 1000, height = 800, res = 150)
  pheatmap(
    cp_mat,
    cluster_rows = TRUE, cluster_cols = FALSE,
    display_numbers = TRUE, number_format = "%.2f",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-0.6, 0.6, length.out = 101),
    annotation_row = annot,
    main = "Density–State Coupling (Spearman ρ)",
    fontsize = 8, fontsize_number = 7,
    na_col = "gray90"
  )
  dev.off()
}

# --- B. Tier 分布柱状图 (按 condition) ---
cat("  Tier distribution bar chart...\n")
tier_qc_file <- file.path(RESULTS_ROOT, "QC", "R5_multilevel_qc.csv")
if (file.exists(tier_qc_file)) {
  tqc <- read.csv(tier_qc_file, stringsAsFactors = FALSE)

  tier_long <- tqc %>%
    select(sample, condition, tier_strong, tier_moderate, tier_weak) %>%
    tidyr::pivot_longer(cols = starts_with("tier_"),
                        names_to = "tier", values_to = "n_genes") %>%
    mutate(tier = gsub("tier_", "", tier),
           tier = factor(tier, levels = c("strong", "moderate", "weak")))

  pB <- ggplot(tier_long, aes(x = reorder(sample, -n_genes), y = n_genes, fill = tier)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("strong" = "#d7191c", "moderate" = "#2c7bb6", "weak" = "#cccccc")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(size = 13, face = "bold")) +
    labs(title = "Density-Gene Association Tiers per Sample",
         x = "", y = "Number of Genes", fill = "Tier")
  ggsave(file.path(summary_dir, "tier_distribution_by_sample.png"),
         pB, width = 12, height = 6, dpi = 200)
}

# --- C. Cross-sample strong genes 条形图 ---
cat("  Cross-sample strong genes...\n")
cross_file <- file.path(RESULTS_ROOT, "QC", "R5_cross_sample_gene_tiers.csv")
if (file.exists(cross_file)) {
  cg <- read.csv(cross_file, stringsAsFactors = FALSE)
  top_cross <- head(cg[order(-cg$n_strong, -cg$n_moderate), ], 30)
  top_cross$gene <- factor(top_cross$gene, levels = rev(top_cross$gene))

  pC <- ggplot(top_cross, aes(x = gene, y = n_strong)) +
    geom_bar(stat = "identity", fill = "#d7191c", width = 0.6) +
    geom_bar(aes(y = n_moderate), stat = "identity", fill = "#2c7bb6", width = 0.6, alpha = 0.5) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 13, face = "bold")) +
    labs(title = "Top 30 Cross-Sample Density Genes",
         subtitle = "Red = strong tier count, Blue overlay = moderate tier count",
         x = "", y = "Number of Samples")
  ggsave(file.path(summary_dir, "top30_cross_sample_genes.png"),
         pC, width = 8, height = 8, dpi = 200)
}

total_time <- round(difftime(Sys.time(), t_total, units = "mins"), 1)

cat("\n============================================================\n")
cat("  R7 COMPLETE\n")
cat("  Total time:", total_time, "min\n")
cat("  Per-sample plots: Results/{project}/{sample}/plots/\n")
cat("  Summary plots:", summary_dir, "\n")
cat("============================================================\n")
