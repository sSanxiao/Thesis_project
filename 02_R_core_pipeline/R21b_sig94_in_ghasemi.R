# ============================================================
# R21b: sig_94 cross-validation in Ghasemi 2024 snRNA-seq
# ------------------------------------------------------------
# Builds on R21a (ghasemi_seurat_R21a.rds with 13 clusters).
#
# Main analyses:
#   1. Compute sig_94 + sig_core scores
#   2. Compute 5 Ghasemi cell-state module scores
#   3. Boxplot sig_94 across 13 clusters + ANOVA + Tukey HSD
#   4. sig_94 vs 5 Ghasemi modules (cell-level Spearman)
#   5. UMAP colored by sig_94 + sig_core (FeaturePlot)
#   6. Cluster-level summary table
#
# Outputs:
#   <RESULTS_DIR>/R21_Ghasemi/
#     R21b_SUMMARY.txt
#     sig94_by_cluster_boxplot.png
#     sigcore_by_cluster_boxplot.png
#     umap_sig94.png / umap_sigcore.png
#     sig94_vs_ghasemi_modules.png
#     cluster_summary_table.csv
#     tukey_HSD_results.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
})

set.seed(42)

RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(RESULTS_DIR, "R21_Ghasemi")
SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")

cat("================================================================\n")
cat("R21b: sig_94 cross-validation in Ghasemi 2024\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Load R21a Seurat object
# ============================================================
cat("[1] Loading R21a Seurat object...\n")
combined <- readRDS(file.path(OUT_DIR, "ghasemi_seurat_R21a.rds"))
cat(sprintf("  Loaded: %d cells × %d genes, %d clusters\n",
            ncol(combined), nrow(combined),
            length(unique(combined$seurat_clusters))))

# ============================================================
# 2. Load sig_94 provenance and define gene sets
# ============================================================
cat("\n[2] Loading sig_94 provenance...\n")

CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
SIG_CORE_GENES <- c("TUBB4A", "APOE", "CCR7", "EOMES", "ST18", "NES", "AQP4", "QKI")

sig_prov <- fread(SIG_PROV)
sig_94_full <- setdiff(sig_prov$gene, CONFLICT_GENES)

cat(sprintf("  sig_94 total: %d genes\n", length(sig_94_full)))

# Match to Ghasemi data
match_94 <- intersect(sig_94_full, rownames(combined))
match_core <- intersect(SIG_CORE_GENES, rownames(combined))
cat(sprintf("  sig_94 matched in Ghasemi: %d / %d (%.1f%%)\n",
            length(match_94), length(sig_94_full),
            100 * length(match_94) / length(sig_94_full)))
cat(sprintf("  sig_core matched: %d / %d\n",
            length(match_core), length(SIG_CORE_GENES)))

# Direction map
dir_table <- sig_prov[gene %in% match_94, .(gene, direction_final)]
setkey(dir_table, gene)
dir_vec <- dir_table[J(match_94), direction_final]
sig_94_pos <- match_94[dir_vec == "positive"]
sig_94_neg <- match_94[dir_vec == "negative"]

core_dir <- sig_prov[gene %in% match_core, .(gene, direction_final)]
setkey(core_dir, gene)
core_dir_vec <- core_dir[J(match_core), direction_final]
sig_core_pos <- match_core[core_dir_vec == "positive"]
sig_core_neg <- match_core[core_dir_vec == "negative"]

cat(sprintf("  sig_94: %d pos + %d neg directions\n",
            length(sig_94_pos), length(sig_94_neg)))
cat(sprintf("  sig_core: %d pos + %d neg directions\n",
            length(sig_core_pos), length(sig_core_neg)))

# ============================================================
# 3. Compute sig_94 / sig_core via AddModuleScore (paired pos/neg)
# ============================================================
cat("\n[3] Computing signature scores...\n")

DefaultAssay(combined) <- "RNA"

# AddModuleScore for sig_94 pos and neg separately
combined <- AddModuleScore(combined,
                            features = list(sig_94_pos),
                            name = "sig_94_pos",
                            ctrl = 50,
                            assay = "RNA")
combined <- AddModuleScore(combined,
                            features = list(sig_94_neg),
                            name = "sig_94_neg",
                            ctrl = 50,
                            assay = "RNA")
combined$sig_94 <- combined$sig_94_pos1 - combined$sig_94_neg1

# sig_core
combined <- AddModuleScore(combined,
                            features = list(sig_core_neg),
                            name = "sig_core_neg",
                            ctrl = 50,
                            assay = "RNA")
combined$sig_core <- -combined$sig_core_neg1  # all neg direction → negate

cat(sprintf("  sig_94: range [%.3f, %.3f], median=%.3f\n",
            min(combined$sig_94), max(combined$sig_94), median(combined$sig_94)))
cat(sprintf("  sig_core: range [%.3f, %.3f], median=%.3f\n",
            min(combined$sig_core), max(combined$sig_core), median(combined$sig_core)))

# ============================================================
# 4. Compute 5 Ghasemi cell-state module scores
# ============================================================
cat("\n[4] Computing Ghasemi cell-state module scores...\n")

ghasemi_markers <- list(
  early_CGNP_proliferating = c("MKI67", "TOP2A", "ATOH1", "BARHL1", "ZIC1", "ZIC3"),
  early_CGNP_quiescent = c("PTCH1", "SMO", "HHIP", "GLI1", "GLI2", "PTPRK"),
  migrating = c("GRIN2B", "CNTN2", "ASTN1", "SEMA6A"),
  postmitotic_differentiated = c("GABRA1", "GABRA6", "GRIN2C", "RBFOX3"),
  astrocytic_like = c("LAMA2", "SOX2", "SOX9", "GFAP", "AQP4")
)

# Filter to present
ghasemi_present <- lapply(ghasemi_markers, function(g) intersect(g, rownames(combined)))

for (state_name in names(ghasemi_present)) {
  genes <- ghasemi_present[[state_name]]
  if (length(genes) == 0) next
  cat(sprintf("  Module %s: %d / %d genes present\n",
              state_name, length(genes), length(ghasemi_markers[[state_name]])))
  combined <- AddModuleScore(combined,
                              features = list(genes),
                              name = state_name,
                              ctrl = 50,
                              assay = "RNA")
}

# Rename module score columns (AddModuleScore appends "1" suffix)
for (state_name in names(ghasemi_present)) {
  old_col <- paste0(state_name, "1")
  new_col <- paste0("module_", state_name)
  combined[[new_col]] <- combined[[old_col]]
  combined[[old_col]] <- NULL
}

# ============================================================
# 5. Build cell-level data table
# ============================================================
cat("\n[5] Building cell-level data table...\n")

cell_dt <- data.table(
  cell_id = colnames(combined),
  sample = combined$sample,
  cluster = factor(combined$seurat_clusters),
  sig_94 = combined$sig_94,
  sig_core = combined$sig_core,
  module_early_CGNP_proliferating = combined$module_early_CGNP_proliferating,
  module_early_CGNP_quiescent = combined$module_early_CGNP_quiescent,
  module_migrating = combined$module_migrating,
  module_postmitotic_differentiated = combined$module_postmitotic_differentiated,
  module_astrocytic_like = combined$module_astrocytic_like
)

fwrite(cell_dt, file.path(OUT_DIR, "cell_data_with_scores.csv"))
cat(sprintf("  Saved cell-level data: %d cells × %d columns\n",
            nrow(cell_dt), ncol(cell_dt)))

# ============================================================
# 6. Cluster-level summary
# ============================================================
cat("\n[6] Cluster-level summary...\n")

cluster_summary <- cell_dt[, .(
  n = .N,
  sig_94_mean = round(mean(sig_94), 4),
  sig_94_median = round(median(sig_94), 4),
  sig_core_mean = round(mean(sig_core), 4),
  prolif_mean = round(mean(module_early_CGNP_proliferating), 3),
  quiescent_mean = round(mean(module_early_CGNP_quiescent), 3),
  migrating_mean = round(mean(module_migrating), 3),
  diff_mean = round(mean(module_postmitotic_differentiated), 3),
  astro_mean = round(mean(module_astrocytic_like), 3)
), by = cluster]
setorder(cluster_summary, -sig_94_mean)

cat("\n=== Cluster summary (sorted by sig_94 mean, descending) ===\n")
print(cluster_summary)

fwrite(cluster_summary, file.path(OUT_DIR, "cluster_summary_table.csv"))

# Identify dominant Ghasemi state per cluster
state_cols <- c("prolif_mean", "quiescent_mean", "migrating_mean", "diff_mean", "astro_mean")
state_labels <- c("proliferating", "quiescent", "migrating", "differentiated", "astrocytic")

cluster_summary[, dominant_state := state_labels[apply(.SD, 1, which.max)],
                 .SDcols = state_cols]

cat("\n=== Cluster x Dominant Ghasemi state ===\n")
print(cluster_summary[, .(cluster, n, sig_94_mean, dominant_state)])

# ============================================================
# 7. ANOVA + Tukey HSD: sig_94 across 13 clusters
# ============================================================
cat("\n[7] ANOVA + Tukey HSD on sig_94...\n")

aov_sig94 <- aov(sig_94 ~ cluster, data = cell_dt)
anova_sig94 <- summary(aov_sig94)
cat(sprintf("  sig_94 cluster ANOVA: F=%.1f, p=%.3e\n",
            anova_sig94[[1]]$`F value`[1],
            anova_sig94[[1]]$`Pr(>F)`[1]))

# Tukey HSD
tukey_sig94 <- TukeyHSD(aov_sig94)
tukey_dt <- as.data.table(tukey_sig94$cluster, keep.rownames = "comparison")
setnames(tukey_dt, c("comparison", "diff", "lwr", "upr", "p_adj"))
tukey_dt[, sig := fcase(
  p_adj < 0.001, "***",
  p_adj < 0.01, "**",
  p_adj < 0.05, "*",
  default = "ns"
)]
fwrite(tukey_dt, file.path(OUT_DIR, "tukey_HSD_results.csv"))
cat(sprintf("  Tukey HSD: %d / %d pairwise comparisons significant (p<0.05)\n",
            sum(tukey_dt$p_adj < 0.05), nrow(tukey_dt)))

# Same for sig_core
aov_sigcore <- aov(sig_core ~ cluster, data = cell_dt)
anova_sigcore <- summary(aov_sigcore)
cat(sprintf("  sig_core cluster ANOVA: F=%.1f, p=%.3e\n",
            anova_sigcore[[1]]$`F value`[1],
            anova_sigcore[[1]]$`Pr(>F)`[1]))

# ============================================================
# 8. sig_94 vs Ghasemi modules: Spearman correlation
# ============================================================
cat("\n[8] sig_94 vs Ghasemi modules (cell-level Spearman)...\n")

sp_table <- data.table(
  module = c("early_CGNP_proliferating", "early_CGNP_quiescent",
             "migrating", "postmitotic_differentiated", "astrocytic_like"),
  rho_sig94 = NA_real_,
  rho_sigcore = NA_real_
)

for (i in 1:nrow(sp_table)) {
  module_col <- paste0("module_", sp_table$module[i])
  sp_table$rho_sig94[i] <- cor(cell_dt$sig_94, cell_dt[[module_col]],
                                 method = "spearman", use = "complete.obs")
  sp_table$rho_sigcore[i] <- cor(cell_dt$sig_core, cell_dt[[module_col]],
                                   method = "spearman", use = "complete.obs")
}
sp_table[, rho_sig94 := round(rho_sig94, 3)]
sp_table[, rho_sigcore := round(rho_sigcore, 3)]

cat("\n  sig_94 / sig_core vs Ghasemi modules:\n")
print(sp_table)

fwrite(sp_table, file.path(OUT_DIR, "sig94_vs_ghasemi_modules_spearman.csv"))

# ============================================================
# 9. Visualizations
# ============================================================
cat("\n[9] Rendering plots...\n")

# Sort cluster levels by sig_94 mean for boxplot
cluster_order <- as.character(cluster_summary$cluster)
cell_dt[, cluster_ordered := factor(cluster, levels = cluster_order)]

# Annotation labels for x axis
cluster_labels <- sapply(cluster_order, function(cl) {
  ds <- cluster_summary[cluster == cl, dominant_state]
  sprintf("C%s\n(%s)", cl, ds)
})

# Plot 1: sig_94 boxplot across clusters
p1 <- ggplot(cell_dt, aes(x = cluster_ordered, y = sig_94, fill = cluster_ordered)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_x_discrete(labels = cluster_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = sprintf("R21b: sig_94 across Ghasemi clusters (n=%d cells)",
                        nrow(cell_dt)),
       subtitle = sprintf("ANOVA F=%.1f, p=%.2e | Sorted by sig_94 mean (high → low)",
                           anova_sig94[[1]]$`F value`[1],
                           anova_sig94[[1]]$`Pr(>F)`[1]),
       x = "Cluster (annotated by dominant Ghasemi state)",
       y = "sig_94 score") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(size = 9))

ggsave(file.path(OUT_DIR, "sig94_by_cluster_boxplot.png"),
       p1, width = 14, height = 7, dpi = 150, bg = "white")

# Plot 2: sig_core boxplot
p2 <- ggplot(cell_dt, aes(x = cluster_ordered, y = sig_core, fill = cluster_ordered)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_x_discrete(labels = cluster_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = sprintf("R21b: sig_core across Ghasemi clusters"),
       subtitle = sprintf("ANOVA F=%.1f, p=%.2e",
                           anova_sigcore[[1]]$`F value`[1],
                           anova_sigcore[[1]]$`Pr(>F)`[1]),
       x = "Cluster",
       y = "sig_core score") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(size = 9))

ggsave(file.path(OUT_DIR, "sigcore_by_cluster_boxplot.png"),
       p2, width = 14, height = 7, dpi = 150, bg = "white")

# Plot 3: UMAP colored by sig_94
DefaultAssay(combined) <- "RNA"
p3a <- FeaturePlot(combined, features = "sig_94", reduction = "umap",
                    pt.size = 0.3) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                          midpoint = 0,
                          limits = c(-max(abs(range(combined$sig_94))),
                                     max(abs(range(combined$sig_94)))),
                          oob = scales::squish) +
  labs(title = "R21b: UMAP colored by sig_94",
       subtitle = sprintf("Range [%.2f, %.2f], median=%.3f",
                           min(combined$sig_94), max(combined$sig_94),
                           median(combined$sig_94))) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "umap_sig94.png"),
       p3a, width = 10, height = 8, dpi = 150, bg = "white")

# Plot 4: UMAP colored by sig_core
p3b <- FeaturePlot(combined, features = "sig_core", reduction = "umap",
                    pt.size = 0.3) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                          midpoint = 0,
                          limits = c(-max(abs(range(combined$sig_core))),
                                     max(abs(range(combined$sig_core)))),
                          oob = scales::squish) +
  labs(title = "R21b: UMAP colored by sig_core",
       subtitle = sprintf("Range [%.2f, %.2f]",
                           min(combined$sig_core), max(combined$sig_core))) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "umap_sigcore.png"),
       p3b, width = 10, height = 8, dpi = 150, bg = "white")

# Plot 5: sig_94 vs Ghasemi modules (4-panel scatter)
ghasemi_module_cols <- c("module_early_CGNP_proliferating",
                          "module_migrating",
                          "module_postmitotic_differentiated",
                          "module_astrocytic_like")
ghasemi_module_labels <- c("Proliferating CGNP", "Migrating",
                            "Differentiated GN", "Astrocytic-like")

scatter_panels <- list()
for (i in seq_along(ghasemi_module_cols)) {
  module_col <- ghasemi_module_cols[i]
  module_label <- ghasemi_module_labels[i]
  rho <- sp_table[module == sub("module_", "", module_col), rho_sig94]
  
  # Subsample for plot clarity
  n_sub <- min(8000, nrow(cell_dt))
  sub <- cell_dt[sample(.N, n_sub)]
  
  p <- ggplot(sub, aes(x = .data[[module_col]], y = sig_94)) +
    geom_point(size = 0.3, alpha = 0.3, color = "#377eb8") +
    geom_smooth(method = "loess", color = "red", linewidth = 0.6, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("%s   ρ=%+.3f", module_label, rho),
         x = sprintf("%s module score", module_label),
         y = "sig_94 score") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  scatter_panels[[i]] <- p
}

p_scatter_all <- wrap_plots(scatter_panels, ncol = 2) +
  plot_annotation(title = "R21b: sig_94 vs Ghasemi 2024 cell-state modules",
                  subtitle = sprintf("Cell-level Spearman, n=%d (subsampled to 8000 for display)",
                                       nrow(cell_dt)),
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "sig94_vs_ghasemi_modules.png"),
       p_scatter_all, width = 14, height = 10, dpi = 150, bg = "white")

# ============================================================
# 10. SUMMARY
# ============================================================
cat("\n[10] Writing R21b_SUMMARY.txt...\n")

# Top 3 clusters by sig_94
top3 <- cluster_summary[1:3]
bot3 <- cluster_summary[(nrow(cluster_summary)-2):nrow(cluster_summary)]

summary_lines <- c(
  "================================================================",
  "R21b — sig_94 cross-validation in Ghasemi 2024 snRNA-seq",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "DATA",
  sprintf("  17,495 cells from 4 MBEN samples (Ghasemi 2024 GSE239854)"),
  sprintf("  13 clusters from R21a Harmony-integrated UMAP"),
  "",
  "================================================================",
  "SIGNATURE COVERAGE IN GHASEMI DATA",
  "================================================================",
  sprintf("  sig_94: %d / %d (%.1f%%) genes matched",
          length(match_94), length(sig_94_full),
          100 * length(match_94) / length(sig_94_full)),
  sprintf("    %d positive direction + %d negative direction",
          length(sig_94_pos), length(sig_94_neg)),
  sprintf("  sig_core: %d / %d genes matched",
          length(match_core), length(SIG_CORE_GENES)),
  "",
  "================================================================",
  "SCORE RANGES",
  "================================================================",
  sprintf("  sig_94: [%.3f, %.3f] median=%.3f",
          min(combined$sig_94), max(combined$sig_94), median(combined$sig_94)),
  sprintf("  sig_core: [%.3f, %.3f] median=%.3f",
          min(combined$sig_core), max(combined$sig_core), median(combined$sig_core)),
  "",
  "================================================================",
  "CLUSTER × sig_94 SUMMARY (sorted by sig_94 mean, descending)",
  "================================================================",
  ""
)

# Cluster table formatted
header <- sprintf("  %-10s %7s %12s %12s %20s",
                   "Cluster", "n", "sig_94_mean", "sig_core_mean", "dominant_state")
summary_lines <- c(summary_lines, header,
                    paste(rep("-", nchar(header)), collapse=""))

for (i in 1:nrow(cluster_summary)) {
  row_str <- sprintf("  %-10s %7d %+12.4f %+12.4f %20s",
                      paste0("Cluster_", cluster_summary$cluster[i]),
                      cluster_summary$n[i],
                      cluster_summary$sig_94_mean[i],
                      cluster_summary$sig_core_mean[i],
                      cluster_summary$dominant_state[i])
  summary_lines <- c(summary_lines, row_str)
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "ANOVA RESULTS",
  "================================================================",
  sprintf("  sig_94 ~ cluster: F=%.1f, p=%.3e",
          anova_sig94[[1]]$`F value`[1],
          anova_sig94[[1]]$`Pr(>F)`[1]),
  sprintf("  sig_core ~ cluster: F=%.1f, p=%.3e",
          anova_sigcore[[1]]$`F value`[1],
          anova_sigcore[[1]]$`Pr(>F)`[1]),
  sprintf("  Tukey HSD: %d / %d pairwise comparisons significant (p<0.05)",
          sum(tukey_dt$p_adj < 0.05), nrow(tukey_dt)),
  "",
  "================================================================",
  "sig_94 vs GHASEMI 2024 MODULES (cell-level Spearman)",
  "================================================================",
  ""
)

for (i in 1:nrow(sp_table)) {
  status <- fcase(
    abs(sp_table$rho_sig94[i]) > 0.5, "★★★ very strong",
    abs(sp_table$rho_sig94[i]) > 0.3, "★★ strong",
    abs(sp_table$rho_sig94[i]) > 0.1, "★ moderate",
    default = "weak"
  )
  summary_lines <- c(summary_lines,
    sprintf("  %-30s  ρ_sig94=%+.3f  ρ_sigcore=%+.3f  [%s]",
            sp_table$module[i],
            sp_table$rho_sig94[i],
            sp_table$rho_sigcore[i],
            status))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "TOP 3 CLUSTERS BY sig_94 (highest)",
  "================================================================",
  ""
)
for (i in 1:nrow(top3)) {
  summary_lines <- c(summary_lines,
    sprintf("  %d. Cluster %s (n=%d): sig_94=%+.4f, dominant=%s",
            i, top3$cluster[i], top3$n[i],
            top3$sig_94_mean[i], top3$dominant_state[i]))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "BOTTOM 3 CLUSTERS BY sig_94 (lowest)",
  "================================================================",
  ""
)
for (i in 1:nrow(bot3)) {
  summary_lines <- c(summary_lines,
    sprintf("  %d. Cluster %s (n=%d): sig_94=%+.4f, dominant=%s",
            i, bot3$cluster[i], bot3$n[i],
            bot3$sig_94_mean[i], bot3$dominant_state[i]))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "INTERPRETATION GUIDE",
  "================================================================",
  "",
  "  Expected if sig_94 captures 'differentiation/RL-derived' program:",
  "    - Highest sig_94 in clusters with 'differentiated' or 'astrocytic' dominant state",
  "    - Lowest sig_94 in clusters with 'proliferating' dominant state",
  "    - Negative ρ between sig_94 and proliferating module",
  "    - Positive ρ between sig_94 and differentiated/astrocytic modules",
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  R21b_SUMMARY.txt                            - this file",
  "  cell_data_with_scores.csv                    - per-cell scores",
  "  cluster_summary_table.csv                    - per-cluster summary",
  "  tukey_HSD_results.csv                        - all pairwise tests",
  "  sig94_vs_ghasemi_modules_spearman.csv        - module correlations",
  "  sig94_by_cluster_boxplot.png                 - main boxplot ★",
  "  sigcore_by_cluster_boxplot.png",
  "  umap_sig94.png                               - UMAP gradient ★",
  "  umap_sigcore.png",
  "  sig94_vs_ghasemi_modules.png                 - 4-panel scatter ★",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21b_SUMMARY.txt"))

cat("\n=== R21b SUMMARY (preview) ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21b DONE\n")
cat("================================================================\n")
