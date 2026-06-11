# ============================================================
# R21d v2: sig_94 + SHH pathway co-localization in Aldinger
# Fix: Aldinger has only RNA assay (no SCT), use sig_94_zscore
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(RESULTS_DIR, "R21d_Aldinger_SHH")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("R21d v2: sig_94 + SHH pathway co-localization in Aldinger\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Load
# ============================================================
cat("[1] Loading Aldinger Seurat object...\n")
aldinger <- readRDS(file.path(RESULTS_DIR, "R17_Aldinger2021", "aldinger2021_scored.rds"))

cat(sprintf("  Cells: %d, Genes: %d\n", ncol(aldinger), nrow(aldinger)))
cat(sprintf("  Assays: %s\n", paste(Assays(aldinger), collapse=", ")))
cat(sprintf("  Default assay: %s\n", DefaultAssay(aldinger)))

# Use sig_94_zscore (the zscore-weighted version, consistent with R18/R21)
aldinger$sig_94 <- aldinger$sig_94_zscore
aldinger$sig_core <- aldinger$sig_core_zscore
aldinger$cluster <- as.character(aldinger$figure_clusters)

cat(sprintf("\n  Using sig_94_zscore: range [%.3f, %.3f]\n",
            min(aldinger$sig_94, na.rm=TRUE),
            max(aldinger$sig_94, na.rm=TRUE)))
cat(sprintf("  Using sig_core_zscore: range [%.3f, %.3f]\n",
            min(aldinger$sig_core, na.rm=TRUE),
            max(aldinger$sig_core, na.rm=TRUE)))

# Drop cells with NA scores
keep_cells <- !is.na(aldinger$sig_94)
aldinger <- subset(aldinger, cells = colnames(aldinger)[keep_cells])
cat(sprintf("  After NA filter: %d cells\n", ncol(aldinger)))

# ============================================================
# 2. SHH pathway genes
# ============================================================
shh_genes <- c("PTCH1", "HHIP", "GLI1", "GLI2", "GLI3", "BOC", "PTPRK")
shh_present <- intersect(shh_genes, rownames(aldinger))
cat(sprintf("\n[2] SHH pathway genes available: %d (%s)\n",
            length(shh_present), paste(shh_present, collapse=", ")))

# ============================================================
# 3. Get expression data (RNA assay only available)
# ============================================================
cat("\n[3] Extracting expression data from RNA assay...\n")

# RNA assay - use data layer (log-normalized)
expr <- GetAssayData(aldinger, assay = "RNA", layer = "data")
if (is.null(expr) || nrow(expr) == 0) {
  cat("  RNA data layer empty, trying counts layer...\n")
  expr <- GetAssayData(aldinger, assay = "RNA", layer = "counts")
  cat("  Note: using raw counts (not log-normalized)\n")
}

cat(sprintf("  Expression matrix: %d genes × %d cells\n",
            nrow(expr), ncol(expr)))
cat(sprintf("  Range: [%.3f, %.3f]\n",
            min(expr, na.rm=TRUE), max(expr, na.rm=TRUE)))

# ============================================================
# 4. Cell-level data table
# ============================================================
cat("\n[4] Building cell-level data table...\n")

cell_dt <- data.table(
  cell_id = colnames(aldinger),
  cluster = aldinger$cluster,
  sig_94 = as.numeric(aldinger$sig_94),
  sig_core = as.numeric(aldinger$sig_core)
)

for (g in shh_present) {
  cell_dt[[g]] <- as.numeric(expr[g, ])
}

cat(sprintf("  Built data table: %d cells × %d cols\n",
            nrow(cell_dt), ncol(cell_dt)))

fwrite(cell_dt, file.path(OUT_DIR, "cell_data_with_SHH.csv"))

# ============================================================
# 5. Cell-level Spearman
# ============================================================
cat("\n[5] sig_94 vs SHH genes (cell-level Spearman)...\n")

cor_table <- data.table(
  gene = shh_present,
  rho_sig94 = NA_real_,
  rho_sigcore = NA_real_,
  pct_expressing = NA_real_
)

for (i in seq_along(shh_present)) {
  g <- shh_present[i]
  expr_vec <- cell_dt[[g]]
  cor_table$rho_sig94[i] <- round(cor(cell_dt$sig_94, expr_vec,
                                        method = "spearman", use = "complete.obs"), 3)
  cor_table$rho_sigcore[i] <- round(cor(cell_dt$sig_core, expr_vec,
                                          method = "spearman", use = "complete.obs"), 3)
  cor_table$pct_expressing[i] <- round(100 * sum(expr_vec > 0) / length(expr_vec), 1)
}

setorder(cor_table, -rho_sig94)
cat("\nsig_94 vs SHH gene cell-level Spearman:\n")
print(cor_table)

fwrite(cor_table, file.path(OUT_DIR, "shh_gene_correlation_table.csv"))

# ============================================================
# 6. Top 20% vs Bottom 20%
# ============================================================
cat("\n[6] Top vs bottom sig_94 cells...\n")

q20 <- quantile(cell_dt$sig_94, 0.20, na.rm = TRUE)
q80 <- quantile(cell_dt$sig_94, 0.80, na.rm = TRUE)
cell_dt[, tier := fcase(
  sig_94 >= q80, "top20",
  sig_94 <= q20, "bot20",
  default = "mid60"
)]
cell_dt[, tier := factor(tier, levels = c("top20", "mid60", "bot20"))]

cat(sprintf("  q20=%.4f, q80=%.4f\n", q20, q80))
cat(sprintf("  top20: %d cells, bot20: %d cells\n",
            sum(cell_dt$tier == "top20"), sum(cell_dt$tier == "bot20")))

tier_expr <- data.table(
  gene = shh_present,
  top20_mean = NA_real_,
  bot20_mean = NA_real_,
  fold_change = NA_real_
)

for (i in seq_along(shh_present)) {
  g <- shh_present[i]
  top_mean <- mean(cell_dt[tier == "top20", get(g)])
  bot_mean <- mean(cell_dt[tier == "bot20", get(g)])
  tier_expr$top20_mean[i] <- round(top_mean, 4)
  tier_expr$bot20_mean[i] <- round(bot_mean, 4)
  tier_expr$fold_change[i] <- round(log2((top_mean + 0.01) / (bot_mean + 0.01)), 3)
}

setorder(tier_expr, -fold_change)
cat("\nSHH gene expression: top20 vs bot20:\n")
print(tier_expr)

fwrite(tier_expr, file.path(OUT_DIR, "tier_shh_expression.csv"))

# ============================================================
# 7. Cluster-level mean
# ============================================================
cat("\n[7] Cluster-level SHH gene expression...\n")

cluster_means <- cell_dt[, lapply(.SD, mean, na.rm = TRUE),
                          .SDcols = shh_present, by = cluster]
cluster_sig94 <- cell_dt[, .(sig_94_mean = round(mean(sig_94, na.rm=TRUE), 4)),
                          by = cluster]
cluster_means <- merge(cluster_means, cluster_sig94, by = "cluster")
setorder(cluster_means, -sig_94_mean)

cat("\nCluster summary (sorted by sig_94 mean):\n")
print(cluster_means)

fwrite(cluster_means, file.path(OUT_DIR, "cluster_shh_expression.csv"))

# ============================================================
# 8. Plots
# ============================================================
cat("\n[8] Rendering plots...\n")

# --- Plot 1: 4-panel scatter ---
top_shh <- cor_table$gene[1:min(4, nrow(cor_table))]
scatter_panels <- list()

for (i in seq_along(top_shh)) {
  g <- top_shh[i]
  rho <- cor_table[gene == g, rho_sig94]
  
  n_sub <- min(8000, nrow(cell_dt))
  sub <- cell_dt[sample(.N, n_sub)]
  
  p <- ggplot(sub, aes(x = .data[[g]], y = sig_94)) +
    geom_point(size = 0.2, alpha = 0.3, color = "#377eb8") +
    geom_smooth(method = "loess", color = "red", linewidth = 0.5, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("%s   ρ=%+.3f", g, rho),
         x = sprintf("%s expression (log-norm)", g),
         y = "sig_94 zscore") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  scatter_panels[[i]] <- p
}

p_scatter <- wrap_plots(scatter_panels, ncol = 2) +
  plot_annotation(title = "R21d: sig_94 vs SHH pathway gene expression in Aldinger fetal cerebellum",
                  subtitle = sprintf("Cell-level Spearman, n=%d (subsampled to 8000 for display)",
                                      nrow(cell_dt)),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(OUT_DIR, "sig94_vs_SHH_genes_scatter.png"),
       p_scatter, width = 12, height = 10, dpi = 150, bg = "white")

# --- Plot 2: UMAP comparison ---
cat("  Building UMAP panels...\n")

umap_coords <- as.data.table(Embeddings(aldinger, reduction = "umap"))
umap_coords[, cell_id := colnames(aldinger)]
umap_coords[, sig_94 := aldinger$sig_94]
for (g in shh_present) {
  umap_coords[[g]] <- as.numeric(expr[g, ])
}

umap_lim <- max(abs(quantile(umap_coords$sig_94, c(0.01, 0.99), na.rm=TRUE)))

p_umap_sig94 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = sig_94)) +
  geom_point(size = 0.15, alpha = 0.5) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0,
                         limits = c(-umap_lim, umap_lim),
                         oob = scales::squish, name = "sig_94") +
  coord_fixed() +
  labs(title = "sig_94 zscore") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

top3_shh <- cor_table$gene[1:3]
shh_umap_panels <- list(p_umap_sig94)

for (g in top3_shh) {
  rho <- cor_table[gene == g, rho_sig94]
  p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = .data[[g]])) +
    geom_point(size = 0.15, alpha = 0.5) +
    scale_color_gradient(low = "grey90", high = "#b2182b", name = g) +
    coord_fixed() +
    labs(title = sprintf("%s expression (ρ=%+.3f)", g, rho)) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())
  shh_umap_panels[[length(shh_umap_panels) + 1]] <- p
}

p_umap_combined <- wrap_plots(shh_umap_panels, ncol = 2) +
  plot_annotation(title = "R21d: sig_94 score vs SHH pathway genes on Aldinger UMAP",
                  subtitle = "Visual co-localization check",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(OUT_DIR, "umap_sig94_vs_SHH_genes.png"),
       p_umap_combined, width = 14, height = 12, dpi = 150, bg = "white")

# --- Plot 3: Cluster heatmap ---
cat("  Cluster heatmap...\n")

heatmap_dt <- melt(cluster_means,
                    id.vars = c("cluster", "sig_94_mean"),
                    measure.vars = shh_present,
                    variable.name = "gene", value.name = "expression")

cluster_order <- cluster_means$cluster
heatmap_dt[, cluster := factor(cluster, levels = cluster_order)]

p_heatmap <- ggplot(heatmap_dt, aes(x = gene, y = cluster, fill = expression)) +
  geom_tile() +
  scale_fill_gradient(low = "#f7f7f7", high = "#b2182b") +
  labs(title = "R21d: SHH pathway gene mean expression per Aldinger cluster",
       subtitle = "Clusters sorted top-to-bottom by sig_94 mean (high to low)",
       x = "SHH pathway gene", y = "Aldinger cluster") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())

ggsave(file.path(OUT_DIR, "shh_gene_cluster_mean_heatmap.png"),
       p_heatmap, width = 10, height = max(8, length(cluster_order) * 0.4),
       dpi = 150, bg = "white")

# --- Plot 4: tier bar ---
cat("  Tier bar chart...\n")

tier_long <- melt(tier_expr,
                   id.vars = "gene",
                   measure.vars = c("top20_mean", "bot20_mean"),
                   variable.name = "tier", value.name = "mean_expr")
tier_long[, tier := gsub("_mean", "", tier)]

gene_order <- tier_expr$gene
tier_long[, gene := factor(gene, levels = gene_order)]

p_tier <- ggplot(tier_long, aes(x = gene, y = mean_expr, fill = tier)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("top20" = "#b2182b", "bot20" = "#2166ac")) +
  labs(title = "R21d: SHH gene expression — sig_94 top20% vs bot20% cells in Aldinger",
       subtitle = "Top20% sig_94 cells should have ↑ SHH gene expression if sig_94 = SHH-active program",
       x = "SHH pathway gene (sorted by fold change)",
       y = "Mean expression (log-normalized)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "sig94_top_vs_bot_SHH_expression.png"),
       p_tier, width = 10, height = 6, dpi = 150, bg = "white")

# ============================================================
# 9. SUMMARY
# ============================================================
cat("\n[9] Writing R21d_SUMMARY.txt...\n")

summary_lines <- c(
  "================================================================",
  "R21d — sig_94 + SHH pathway co-localization in Aldinger",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "PURPOSE",
  "  6th-dimension validation: direct cell-level co-localization",
  "  of sig_94 score with SHH pathway gene expression in human",
  "  fetal cerebellum (Aldinger 2021 atlas).",
  "",
  "================================================================",
  "DATA",
  "================================================================",
  sprintf("  Aldinger 2021: %d cells, %d genes (2000 var genes panel)",
          ncol(aldinger), nrow(aldinger)),
  sprintf("  SHH pathway genes available: %d (%s)",
          length(shh_present), paste(shh_present, collapse=", ")),
  sprintf("  Using sig_94_zscore (zscore-weighted, same as R18/R21)"),
  "",
  "================================================================",
  "RESULT 1: Cell-level Spearman ρ",
  "================================================================",
  ""
)

for (i in 1:nrow(cor_table)) {
  status <- fcase(
    abs(cor_table$rho_sig94[i]) > 0.5, "*** very strong",
    abs(cor_table$rho_sig94[i]) > 0.3, "** strong",
    abs(cor_table$rho_sig94[i]) > 0.1, "* moderate",
    default = "weak"
  )
  summary_lines <- c(summary_lines,
    sprintf("  %-8s  ρ_sig94=%+.3f  ρ_sigcore=%+.3f  (%.1f%% expressing)  [%s]",
            cor_table$gene[i],
            cor_table$rho_sig94[i],
            cor_table$rho_sigcore[i],
            cor_table$pct_expressing[i],
            status))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "RESULT 2: top20% vs bot20% sig_94 cells - SHH gene expression",
  "================================================================",
  ""
)

for (i in 1:nrow(tier_expr)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-8s  top20_mean=%.4f  bot20_mean=%.4f  log2FC=%+.2f",
            tier_expr$gene[i],
            tier_expr$top20_mean[i],
            tier_expr$bot20_mean[i],
            tier_expr$fold_change[i]))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "RESULT 3: Top 5 clusters by sig_94 (and their SHH gene expression)",
  "================================================================",
  ""
)

header_str <- sprintf("  %-25s %12s   %s",
                       "Cluster", "sig_94_mean",
                       paste(sprintf("%-8s", shh_present), collapse=""))
summary_lines <- c(summary_lines, header_str,
                    paste(rep("-", nchar(header_str)), collapse=""))

for (i in 1:min(5, nrow(cluster_means))) {
  expr_vals <- sprintf("%-8.3f", as.numeric(cluster_means[i, ..shh_present]))
  row_str <- sprintf("  %-25s %12.4f   %s",
                      substr(cluster_means$cluster[i], 1, 25),
                      cluster_means$sig_94_mean[i],
                      paste(expr_vals, collapse=""))
  summary_lines <- c(summary_lines, row_str)
}

summary_lines <- c(summary_lines, "",
  "  Bottom 3 clusters by sig_94:",
  ""
)

for (i in (nrow(cluster_means)-2):nrow(cluster_means)) {
  if (i < 1) next
  expr_vals <- sprintf("%-8.3f", as.numeric(cluster_means[i, ..shh_present]))
  row_str <- sprintf("  %-25s %12.4f   %s",
                      substr(cluster_means$cluster[i], 1, 25),
                      cluster_means$sig_94_mean[i],
                      paste(expr_vals, collapse=""))
  summary_lines <- c(summary_lines, row_str)
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  R21d_SUMMARY.txt                       - this file",
  "  cell_data_with_SHH.csv",
  "  shh_gene_correlation_table.csv",
  "  tier_shh_expression.csv",
  "  cluster_shh_expression.csv",
  "  sig94_vs_SHH_genes_scatter.png         scatter plot",
  "  umap_sig94_vs_SHH_genes.png            UMAP comparison",
  "  shh_gene_cluster_mean_heatmap.png      cluster heatmap",
  "  sig94_top_vs_bot_SHH_expression.png    tier bar chart",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21d_SUMMARY.txt"))

cat("\n=== R21d SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21d DONE\n")
cat("================================================================\n")
