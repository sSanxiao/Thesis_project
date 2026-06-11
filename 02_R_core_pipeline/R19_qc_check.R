# ============================================================
# R19: 4-sample MBEN QC and batch effect check (Item 3)
# ------------------------------------------------------------
# Goal: rule out the criticism that "MB266 has stronger signal
# because of better technical QC, not biology"
#
# Tests:
#   3a: Basic QC metrics across 4 samples (boxplot)
#   3b: sig_94 score vs nCount_RNA correlation (within-sample)
#   3c: transcript_counts distribution comparison
#
# Key sanity checks:
#   - If MB266 has 2x the median nCount_RNA → technical confound
#   - If sig_94 strongly correlates with nCount_RNA in MB266 →
#     "high-expressing cells" inflates sig_94, not biology
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# 4 sample paths
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
R2_MB <- file.path(RESULTS_DIR, "R2_Results", "Medulloblastoma_Human")
MB_SAMPLES <- list(
  MB263 = file.path(R2_MB, "GSM8840046", "GSM8840046_seurat_R2.rds"),
  MB266 = file.path(R2_MB, "GSM8840047", "GSM8840047_seurat_R2.rds"),
  MB295 = file.path(R2_MB, "GSM8840048", "GSM8840048_seurat_R2.rds"),
  MB299 = file.path(R2_MB, "GSM8840049", "GSM8840049_seurat_R2.rds")
)

OUT_DIR <- file.path(RESULTS_DIR, "R19_QC_check")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("R19: MB sample QC and batch effect check\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Read sig_94 score data from R18 outputs (already computed)
# ============================================================
cat("[1] Loading existing sig_94 scores from R18 outputs...\n")

R18_ROOT <- file.path(RESULTS_DIR, "R18_SpatialSignature")
R18_DATA <- list(
  MB263 = file.path(R18_ROOT, "MB263", "cell_data.csv"),
  MB266 = file.path(R18_ROOT, "MB266", "cell_data.csv"),
  MB295 = file.path(R18_ROOT, "MB295", "cell_data.csv"),
  MB299 = file.path(R18_ROOT, "MB299", "cell_data.csv")
)

# ============================================================
# 2. Loop through samples, extract QC metrics + sig_94
# ============================================================
all_data <- list()

for (sn in names(MB_SAMPLES)) {
  cat(sprintf("\n  Processing %s...\n", sn))
  
  # Load Seurat for QC metrics
  obj <- readRDS(MB_SAMPLES[[sn]])
  md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  rm(obj); gc(verbose = FALSE)
  
  # Add sample label
  md$sample <- sn
  
  # Load sig_94 score from R18
  sig_data <- fread(R18_DATA[[sn]])
  
  # Merge sig_94 / sig_core into QC table
  md_keep <- md[, .(cell_id, sample, nCount_RNA, nFeature_RNA,
                     transcript_counts, cell_area, nucleus_area,
                     density_knn_main_piecewise,
                     data_quality_tier, seurat_clusters)]
  
  merged <- merge(md_keep, sig_data[, .(cell_id, sig_94, sig_core)],
                   by = "cell_id", all.x = TRUE)
  
  cat(sprintf("    %d cells, %d with sig_94\n",
              nrow(merged), sum(!is.na(merged$sig_94))))
  
  all_data[[sn]] <- merged
}

# Combine all samples
big_dt <- rbindlist(all_data)
big_dt[, sample := factor(sample, levels = c("MB263", "MB266", "MB295", "MB299"))]

cat(sprintf("\n  Combined dataset: %d cells\n", nrow(big_dt)))

# ============================================================
# 3a: Per-sample QC summary table
# ============================================================
cat("\n[3a] Computing per-sample QC summary...\n")

qc_summary <- big_dt[, .(
  n_cells = .N,
  median_nCount_RNA = round(median(nCount_RNA, na.rm = TRUE), 1),
  mean_nCount_RNA = round(mean(nCount_RNA, na.rm = TRUE), 1),
  median_nFeature_RNA = round(median(nFeature_RNA, na.rm = TRUE), 1),
  mean_nFeature_RNA = round(mean(nFeature_RNA, na.rm = TRUE), 1),
  median_cell_area = round(median(cell_area, na.rm = TRUE), 1),
  median_nucleus_area = round(median(nucleus_area, na.rm = TRUE), 1),
  median_transcript_counts = round(median(transcript_counts, na.rm = TRUE), 1),
  median_sig_94 = round(median(sig_94, na.rm = TRUE), 4)
), by = sample]

cat("\nPer-sample QC summary:\n")
print(qc_summary)

fwrite(qc_summary, file.path(OUT_DIR, "qc_summary_per_sample.csv"))

# Compute MB266 vs others ratio for nCount_RNA
mb266_nCount <- qc_summary[sample == "MB266", median_nCount_RNA]
others_nCount <- qc_summary[sample != "MB266", median_nCount_RNA]

cat("\n[Sanity check]\n")
cat(sprintf("  MB266 median nCount_RNA: %.1f\n", mb266_nCount))
cat(sprintf("  Other samples median nCount_RNA: %s\n",
            paste(round(others_nCount, 1), collapse = ", ")))
cat(sprintf("  MB266 vs others (max ratio): %.2f×\n",
            mb266_nCount / min(others_nCount)))
if (mb266_nCount / min(others_nCount) > 2.0) {
  cat("  ⚠️ WARNING: MB266 nCount_RNA > 2x others. Possible technical confound.\n")
} else {
  cat("  ✓ MB266 not dramatically different in nCount. Technical confound unlikely.\n")
}

# ============================================================
# 3b: 4-panel QC boxplot
# ============================================================
cat("\n[3b] Rendering 4-panel QC boxplot...\n")

# Pick 4 key QC metrics
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "cell_area", "transcript_counts")
plots <- list()

for (m in qc_metrics) {
  p <- ggplot(big_dt, aes(x = sample, y = .data[[m]], fill = sample)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = m, x = NULL, y = m) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
  
  # log scale if appropriate
  if (m %in% c("nCount_RNA", "nFeature_RNA", "transcript_counts")) {
    p <- p + scale_y_continuous(trans = "log10")
  }
  
  plots[[m]] <- p
}

p_4panel <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "4-sample QC comparison",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "qc_4panel_boxplot.png"),
       p_4panel, width = 12, height = 10, dpi = 150, bg = "white")

# ============================================================
# 3c: sig_94 vs nCount_RNA per-sample correlation
# ============================================================
cat("\n[3c] sig_94 vs nCount_RNA within-sample correlation...\n")

cor_per_sample <- big_dt[!is.na(sig_94), .(
  rho_sig94_nCount = cor(sig_94, nCount_RNA,
                          method = "spearman", use = "complete.obs"),
  rho_sig94_nFeature = cor(sig_94, nFeature_RNA,
                            method = "spearman", use = "complete.obs"),
  rho_sig94_cell_area = cor(sig_94, cell_area,
                             method = "spearman", use = "complete.obs"),
  n = .N
), by = sample]

cat("\nWithin-sample sig_94 vs technical metrics correlation:\n")
print(cor_per_sample)

fwrite(cor_per_sample, file.path(OUT_DIR, "sig94_vs_QC_correlation.csv"))

# Interpretation
cat("\n[Interpretation]\n")
cat("  If sig_94 strongly correlates with nCount_RNA (|ρ| > 0.5),\n")
cat("  it suggests sig_94 score is partially driven by total transcript count.\n")
cat("  Each row above shows the within-sample ρ:\n")
for (i in 1:nrow(cor_per_sample)) {
  rho <- cor_per_sample$rho_sig94_nCount[i]
  status <- if (abs(rho) > 0.5) "⚠ HIGH" else if (abs(rho) > 0.3) "moderate" else "OK"
  cat(sprintf("    %s: ρ(sig_94, nCount_RNA) = %+.3f  [%s]\n",
              cor_per_sample$sample[i], rho, status))
}

# ============================================================
# 3d: Scatter plot — sig_94 vs nCount_RNA per sample
# ============================================================
cat("\n[3d] Rendering scatter plot...\n")

# Subsample for scatter clarity
sub_dt <- big_dt[!is.na(sig_94), .SD[sample(.N, min(5000, .N))], by = sample]

p_scatter <- ggplot(sub_dt, aes(x = nCount_RNA, y = sig_94)) +
  geom_point(size = 0.2, alpha = 0.3, color = "#377eb8") +
  geom_smooth(method = "loess", color = "red", linewidth = 0.5, se = TRUE) +
  facet_wrap(~ sample, ncol = 2, scales = "free") +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "sig_94 vs nCount_RNA (technical confound check)",
       subtitle = "If sig_94 strongly tracks nCount_RNA, signature partially driven by sequencing depth",
       x = "nCount_RNA (log10)", y = "sig_94 score") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "sig94_vs_nCount_scatter.png"),
       p_scatter, width = 12, height = 8, dpi = 150, bg = "white")

# ============================================================
# 3e: SUMMARY
# ============================================================
cat("\n[Summary] Writing R19_SUMMARY.txt...\n")

summary_lines <- c(
  "================================================================",
  "R19 — MB Sample QC + Batch Effect Check",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "PURPOSE",
  "  Rule out technical (not biological) explanation for MB266's",
  "  strong signature signal in R9 AUC and R18 spatial coupling.",
  "",
  "================================================================",
  "PER-SAMPLE QC SUMMARY",
  "================================================================",
  ""
)

for (i in 1:nrow(qc_summary)) {
  summary_lines <- c(summary_lines,
    sprintf("[%s]  n=%d", qc_summary$sample[i], qc_summary$n_cells[i]),
    sprintf("  median nCount_RNA: %.1f  (mean: %.1f)",
            qc_summary$median_nCount_RNA[i], qc_summary$mean_nCount_RNA[i]),
    sprintf("  median nFeature_RNA: %.1f", qc_summary$median_nFeature_RNA[i]),
    sprintf("  median cell_area: %.1f, nucleus_area: %.1f",
            qc_summary$median_cell_area[i], qc_summary$median_nucleus_area[i]),
    sprintf("  median transcript_counts: %.1f", qc_summary$median_transcript_counts[i]),
    sprintf("  median sig_94: %+.4f", qc_summary$median_sig_94[i]),
    "")
}

# Verdict
mb266_factor <- mb266_nCount / min(others_nCount)
summary_lines <- c(summary_lines,
  "================================================================",
  "VERDICT 1: nCount_RNA comparison",
  "================================================================",
  sprintf("  MB266 median nCount_RNA = %.1f", mb266_nCount),
  sprintf("  Min of others        = %.1f", min(others_nCount)),
  sprintf("  Ratio (MB266/min)    = %.2f×", mb266_factor),
  ""
)

if (mb266_factor > 2.0) {
  summary_lines <- c(summary_lines,
    "  ⚠️ WARNING: MB266 has > 2× nCount_RNA vs other samples.",
    "  This may partially explain stronger signature signal.",
    "  Discussion must acknowledge technical confound.")
} else if (mb266_factor > 1.5) {
  summary_lines <- c(summary_lines,
    "  ⚠ MB266 modestly higher nCount_RNA (1.5-2× others).",
    "  Some technical contribution possible but unlikely dominant.")
} else {
  summary_lines <- c(summary_lines,
    "  ✓ MB266 nCount_RNA NOT dramatically different from others.",
    "  Technical confound is unlikely to explain MB266's signal strength.")
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "VERDICT 2: sig_94 vs nCount_RNA correlation",
  "================================================================",
  ""
)

for (i in 1:nrow(cor_per_sample)) {
  rho <- cor_per_sample$rho_sig94_nCount[i]
  status <- if (abs(rho) > 0.5) "⚠ HIGH (technical confound likely)" 
            else if (abs(rho) > 0.3) "moderate (some confound possible)"
            else "OK (independent of seq depth)"
  summary_lines <- c(summary_lines,
    sprintf("  %s: ρ(sig_94, nCount_RNA) = %+.3f  [%s]",
            cor_per_sample$sample[i], rho, status))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  sprintf("  qc_summary_per_sample.csv"),
  sprintf("  sig94_vs_QC_correlation.csv"),
  sprintf("  qc_4panel_boxplot.png"),
  sprintf("  sig94_vs_nCount_scatter.png"),
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R19_SUMMARY.txt"))

cat("\n=== SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R19 DONE\n")
cat("================================================================\n")
