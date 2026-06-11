# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R21d: sig_94 + SHH pathway co-localization in Aldinger fetal cerebellum
# ------------------------------------------------------------
# Closes the 6-dimension validation circle:
#   1. MBEN spatial density coupling (R18)
#   2. Bulk MB cohort prognosis (R15, R16)
#   3. Fetal cerebellar developmental cluster mapping (R17)
#   4. MBEN snRNA-seq cluster cross-validation (R21a/b)
#   5. Cell-level differential expression of SHH pathway genes (R21c)
#   6. *** Direct cell-level sig_94 ↔ SHH pathway co-localization in fetal *** ← R21d
#
# Outputs (in ./results/R21d_Aldinger_SHH/):
#   R21d_SUMMARY.txt
#   sig94_vs_SHH_genes_scatter.png       4-panel cell-level scatter
#   umap_sig94_vs_SHH_genes.png          UMAP comparison panels
#   shh_gene_cluster_mean_heatmap.png    21 cluster × 7 SHH genes
#   sig94_top_vs_bot_SHH_expression.png  bar plot: high vs low sig_94, SHH expr
#   shh_gene_correlation_table.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

OUT_DIR <- "./results/R21d_Aldinger_SHH"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("R21d: sig_94 + SHH pathway co-localization in Aldinger\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Load Aldinger Seurat (already has sig_94 / sig_core scores from R17)
# ============================================================
cat("[1] Loading Aldinger Seurat object...\n")
aldinger <- readRDS("./results/R17_Aldinger2021/aldinger2021_scored.rds")

cat(sprintf("  Cells: %d, Genes: %d\n", ncol(aldinger), nrow(aldinger)))
cat(sprintf("  Available metadata columns: %s\n",
            paste(colnames(aldinger@meta.data), collapse=", ")))

# Identify which sig_94 column to use
# R17 should have added sig_94 score - check available names
md_cols <- colnames(aldinger@meta.data)
sig94_col <- grep("^sig_94|^sig94", md_cols, value = TRUE)[1]
sigcore_col <- grep("^sig_core|^sigcore", md_cols, value = TRUE)[1]

cat(sprintf("  Using sig_94 column: %s\n", sig94_col))
cat(sprintf("  Using sig_core column: %s\n", sigcore_col))

if (is.na(sig94_col)) {
  stop("sig_94 score column not found in Aldinger metadata. Check R17 output.")
}

# Standardize names for downstream
aldinger$sig_94 <- aldinger@meta.data[[sig94_col]]
if (!is.na(sigcore_col)) {
  aldinger$sig_core <- aldinger@meta.data[[sigcore_col]]
}

# Identify cluster column
cluster_col <- grep("cluster|figure_cluster|annotation", md_cols, value = TRUE)[1]
cat(sprintf("  Using cluster column: %s\n", cluster_col))

aldinger$cluster <- aldinger@meta.data[[cluster_col]]

# ============================================================
# 2. SHH pathway genes
# ============================================================
shh_genes <- c("PTCH1", "HHIP", "GLI1", "GLI2", "GLI3", "BOC", "PTPRK")
shh_present <- intersect(shh_genes, rownames(aldinger))
cat(sprintf("\n[2] SHH pathway genes available: %d (%s)\n",
            length(shh_present), paste(shh_present, collapse=", ")))

# ============================================================
# 3. Get expression data
# ============================================================
cat("\n[3] Extracting expression data...\n")

# Try data slot first (log-normalized), fallback to scale.data
DefaultAssay(aldinger) <- "SCT"
expr <- tryCatch({
  GetAssayData(aldinger, assay = "SCT", layer = "data")
}, error = function(e) {
  GetAssayData(aldinger, assay = "RNA", layer = "data")
})

if (is.null(expr) || nrow(expr) == 0) {
  cat("  SCT data missing, trying RNA...\n")
  DefaultAssay(aldinger) <- "RNA"
  expr <- GetAssayData(aldinger, assay = "RNA", layer = "data")
}

cat(sprintf("  Expression matrix: %d genes × %d cells\n",
            nrow(expr), ncol(expr)))

# ============================================================
# 4. Cell-level data table
# ============================================================
cat("\n[4] Building cell-level data table...\n")

cell_dt <- data.table(
  cell_id = colnames(aldinger),
  cluster = as.character(aldinger$cluster),
  sig_94 = aldinger$sig_94
)

# Add SHH gene expression
for (g in shh_present) {
  cell_dt[[g]] <- as.numeric(expr[g, ])
}

if (!is.na(sigcore_col)) {
  cell_dt$sig_core <- aldinger$sig_core
}

cat(sprintf("  Built data table: %d cells × %d cols\n",
            nrow(cell_dt), ncol(cell_dt)))

fwrite(cell_dt, file.path(OUT_DIR, "cell_data_with_SHH.csv"))

# ============================================================
# 5. Cell-level Spearman correlations
# ============================================================
cat("\n[5] sig_94 vs SHH genes (cell-level Spearman)...\n")

cor_table <- data.table(
  gene = shh_present,
  rho_sig94 = NA_real_,
  pct_expressing = NA_real_
)

for (i in seq_along(shh_present)) {
  g <- shh_present[i]
  expr_vec <- cell_dt[[g]]
  cor_table$rho_sig94[i] <- round(cor(cell_dt$sig_94, expr_vec,
                                        method = "spearman", use = "complete.obs"), 3)
  cor_table$pct_expressing[i] <- round(100 * sum(expr_vec > 0) / length(expr_vec), 1)
}

setorder(cor_table, -rho_sig94)
cat("\nsig_94 vs SHH gene cell-level Spearman correlations:\n")
print(cor_table)

fwrite(cor_table, file.path(OUT_DIR, "shh_gene_correlation_table.csv"))

# ============================================================
# 6. Top 20% vs Bottom 20% sig_94 cells: SHH gene expression
# ============================================================
cat("\n[6] Top vs bottom sig_94 cell SHH gene expression...\n")

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

# Per-tier mean expression of each SHH gene
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
  # log2 fold change in expression (avoid div by 0)
  tier_expr$fold_change[i] <- round(log2((top_mean + 0.01) / (bot_mean + 0.01)), 3)
}

setorder(tier_expr, -fold_change)
cat("\nSHH gene expression: top20 vs bot20 sig_94 cells:\n")
print(tier_expr)

fwrite(tier_expr, file.path(OUT_DIR, "tier_shh_expression.csv"))

# ============================================================
# 7. Cluster-level: mean SHH gene expression per cluster
# ============================================================
cat("\n[7] Cluster-level SHH gene expression heatmap...\n")

# Compute mean expression per cluster
cluster_means <- cell_dt[, lapply(.SD, mean, na.rm = TRUE),
                          .SDcols = shh_present, by = cluster]

# Add sig_94 mean per cluster
cluster_sig94 <- cell_dt[, .(sig_94_mean = round(mean(sig_94, na.rm = TRUE), 4)),
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

# --- Plot 1: 4-panel scatter of sig_94 vs top SHH genes ---
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
         y = "sig_94 score") +
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

# --- Plot 2: UMAP comparison panels ---
cat("  Building UMAP panels...\n")

# Get UMAP coordinates
umap_coords <- as.data.table(Embeddings(aldinger, reduction = "umap"))
umap_coords[, cell_id := colnames(aldinger)]
umap_coords[, sig_94 := aldinger$sig_94]
for (g in shh_present) {
  umap_coords[[g]] <- as.numeric(expr[g, ])
}

# Pick top 3 SHH genes by ρ for UMAP comparison
top3_shh <- cor_table$gene[1:3]

# UMAP colored by sig_94
p_umap_sig94 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = sig_94)) +
  geom_point(size = 0.15, alpha = 0.5) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0,
                         limits = c(-max(abs(range(umap_coords$sig_94, na.rm=TRUE))),
                                    max(abs(range(umap_coords$sig_94, na.rm=TRUE)))),
                         oob = scales::squish, name = "sig_94") +
  coord_fixed() +
  labs(title = "sig_94 score") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

# UMAP colored by each top SHH gene
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

# Long format for heatmap
heatmap_dt <- melt(cluster_means,
                    id.vars = c("cluster", "sig_94_mean"),
                    measure.vars = shh_present,
                    variable.name = "gene", value.name = "expression")

# Cluster ordered by sig_94 mean
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

# --- Plot 4: top vs bot tier SHH expression bar ---
cat("  Tier comparison bar chart...\n")

tier_long <- melt(tier_expr,
                   id.vars = "gene",
                   measure.vars = c("top20_mean", "bot20_mean"),
                   variable.name = "tier", value.name = "mean_expr")
tier_long[, tier := gsub("_mean", "", tier)]

# Order genes by fold change
gene_order <- tier_expr$gene
tier_long[, gene := factor(gene, levels = gene_order)]

p_tier <- ggplot(tier_long, aes(x = gene, y = mean_expr, fill = tier)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("top20" = "#b2182b", "bot20" = "#2166ac")) +
  labs(title = "R21d: SHH gene expression — sig_94 top20% vs bot20% cells",
       subtitle = "If sig_94 = SHH-active program: top20 should have ↑ SHH gene expression",
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
  "  Closes the validation circle:",
  "    1. MBEN spatial density coupling (R18)",
  "    2. Bulk MB cohort prognosis (R15, R16)",
  "    3. Fetal cerebellar developmental cluster mapping (R17)",
  "    4. MBEN snRNA-seq cluster cross-validation (R21a/b)",
  "    5. Cell-level differential expression of SHH genes (R21c)",
  "    6. *** Direct cell-level sig_94 ↔ SHH co-localization ***",
  "",
  "================================================================",
  "DATA",
  "================================================================",
  sprintf("  Aldinger 2021: %d cells, %d genes (2000 var genes panel)",
          ncol(aldinger), nrow(aldinger)),
  sprintf("  SHH pathway genes available: %d (%s)",
          length(shh_present), paste(shh_present, collapse=", ")),
  "",
  "================================================================",
  "RESULT 1: Cell-level Spearman ρ (sig_94 vs SHH gene)",
  "================================================================",
  ""
)

for (i in 1:nrow(cor_table)) {
  status <- fcase(
    abs(cor_table$rho_sig94[i]) > 0.5, "★★★ very strong",
    abs(cor_table$rho_sig94[i]) > 0.3, "★★ strong",
    abs(cor_table$rho_sig94[i]) > 0.1, "★ moderate",
    default = "weak"
  )
  summary_lines <- c(summary_lines,
    sprintf("  %-8s  ρ_sig94 = %+.3f   (%.1f%% expressing)   [%s]",
            cor_table$gene[i],
            cor_table$rho_sig94[i],
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
  "RESULT 3: Cluster-level mean (top 5 by sig_94)",
  "================================================================",
  ""
)

# Top 5 clusters by sig_94
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
  "INTERPRETATION",
  "================================================================",
  "",
  "  Expected pattern if sig_94 = SHH-pathway-active CGNP program:",
  "    - Strong positive ρ between sig_94 and PTCH1 / HHIP / GLI1",
  "      (these are SHH target genes induced by active SHH signaling)",
  "    - top20 sig_94 cells should have higher SHH gene expression",
  "    - High-sig_94 clusters should be GCP / RL / eCN (SHH-active progenitors)",
  "    - Low-sig_94 clusters should be mature GN, Purkinje, glia",
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  R21d_SUMMARY.txt                       - this file",
  "  cell_data_with_SHH.csv                 - per-cell sig_94 + SHH expression",
  "  shh_gene_correlation_table.csv         - cell-level Spearman",
  "  tier_shh_expression.csv                - top20 vs bot20 expression",
  "  cluster_shh_expression.csv             - per-cluster means",
  "  sig94_vs_SHH_genes_scatter.png         ★ 4-panel scatter",
  "  umap_sig94_vs_SHH_genes.png            ★ UMAP comparison",
  "  shh_gene_cluster_mean_heatmap.png      cluster heatmap",
  "  sig94_top_vs_bot_SHH_expression.png    tier bar chart",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21d_SUMMARY.txt"))

cat("\n=== R21d SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21d DONE — 6-dimension validation circle closed\n")
cat("================================================================\n")
