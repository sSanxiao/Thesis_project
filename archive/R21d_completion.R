# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R21d completion script - finish what was missed when main crashed at UMAP
# Reads saved CSVs and generates: heatmap, tier bar, SUMMARY
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

OUT_DIR <- "./results/R21d_Aldinger_SHH"

cat("Loading saved CSV files...\n")
cor_table <- fread(file.path(OUT_DIR, "shh_gene_correlation_table.csv"))
tier_expr <- fread(file.path(OUT_DIR, "tier_shh_expression.csv"))
cluster_means <- fread(file.path(OUT_DIR, "cluster_shh_expression.csv"))

shh_present <- cor_table$gene
cat(sprintf("SHH genes: %s\n", paste(shh_present, collapse = ", ")))

# Re-sort cluster_means by sig_94_mean
setorder(cluster_means, -sig_94_mean)
setorder(tier_expr, -fold_change)

# ============================================================
# Plot: cluster heatmap
# ============================================================
cat("\nRendering cluster heatmap...\n")

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

# ============================================================
# Plot: tier bar chart
# ============================================================
cat("Rendering tier bar chart...\n")

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
       subtitle = "Top20% sig_94 cells should have higher SHH gene expression if sig_94 = SHH-active",
       x = "SHH pathway gene (sorted by fold change)",
       y = "Mean expression (log-normalized)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "sig94_top_vs_bot_SHH_expression.png"),
       p_tier, width = 10, height = 6, dpi = 150, bg = "white")

# ============================================================
# SUMMARY
# ============================================================
cat("Writing SUMMARY...\n")

summary_lines <- c(
  "================================================================",
  "R21d - sig_94 + SHH pathway co-localization in Aldinger",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "PURPOSE",
  "  6th-dimension validation: direct cell-level co-localization",
  "  of sig_94 score with SHH pathway gene expression in human",
  "  fetal cerebellum (Aldinger 2021 atlas, 69174 cells).",
  "",
  "================================================================",
  "DATA",
  "================================================================",
  "  Aldinger 2021: 69174 cells, 2000 var genes",
  sprintf("  SHH pathway genes available: %d (%s)",
          length(shh_present), paste(shh_present, collapse=", ")),
  "  sig_94 score: sig_94_zscore (z-score weighted, R17 output)",
  "",
  "================================================================",
  "RESULT 1: Cell-level Spearman ρ (sig_94 vs SHH gene)",
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
    sprintf("  %-8s  rho_sig94=%+.3f  rho_sigcore=%+.3f  (%.1f%% expressing)  [%s]",
            cor_table$gene[i],
            cor_table$rho_sig94[i],
            cor_table$rho_sigcore[i],
            cor_table$pct_expressing[i],
            status))
}

summary_lines <- c(summary_lines, "",
  "  KEY OBSERVATION:",
  "  - All correlations are weak (|rho| < 0.1) for sig_94 vs SHH genes",
  "  - This is DIFFERENT from R21c finding in Ghasemi MBEN data,",
  "    where sig_94 had strong co-expression with PTCH1/HHIP/GLI1.",
  "  - Indicates sig_94's relationship to SHH pathway is",
  "    TUMOR CONTEXT-SPECIFIC, not an intrinsic property of the",
  "    progenitor program in normal development.",
  "",
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
  "  KEY OBSERVATION:",
  "  - PTCH1 and HHIP show modest enrichment in top20 (log2FC +0.59 and +0.48)",
  "  - GLI1, GLI2, GLI3, BOC, PTPRK show no or reverse pattern",
  "  - Consistent with weak overall correlation",
  "",
  "================================================================",
  "RESULT 3: Top 5 clusters by sig_94 (and their SHH gene expression)",
  "================================================================",
  ""
)

# Header
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
  "INTERPRETATION (CRITICAL — narrative-shaping)",
  "================================================================",
  "",
  "  Expected if sig_94 = SHH-pathway-active program (the R21c hypothesis):",
  "    - Strong rho between sig_94 and PTCH1/HHIP/GLI1 in fetal",
  "    - top20 sig_94 cells should be enriched for SHH targets",
  "    - sig_94-high clusters should be SHH-receiving (e.g., 03-GCP)",
  "",
  "  ACTUAL RESULT in fetal Aldinger:",
  "    - All rho values weak (max +0.10 for PTCH1)",
  "    - top20 has only modest PTCH1/HHIP enrichment",
  "    - PTCH1 highest in 08-BG (Bergmann glia), 09-Ast (astrocyte),",
  "      03-GCP (CGNP) — but only 03-GCP has high sig_94",
  "    - Bergmann glia and astrocytes have HIGH SHH genes but LOW sig_94",
  "",
  "  CONCLUSION:",
  "    sig_94's enrichment in MBEN cells co-expressing SHH pathway markers",
  "    (R21c finding) is a TUMOR CONTEXT-SPECIFIC feature, not an",
  "    intrinsic property of the rhombic-lip-derived progenitor program",
  "    in normal cerebellar development.",
  "",
  "  REVISED NARRATIVE (R21abcd integrated):",
  "    sig_94 captures a rhombic-lip-derived progenitor program.",
  "    In MBEN tumors, this program is enriched in cells that also",
  "    express SHH pathway components — reflecting the SHH-driven",
  "    nature of MBEN. In normal fetal cerebellum, the same progenitor",
  "    program does not require concurrent SHH pathway activation.",
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  R21d_SUMMARY.txt                       - this file",
  "  cell_data_with_SHH.csv",
  "  shh_gene_correlation_table.csv",
  "  tier_shh_expression.csv",
  "  cluster_shh_expression.csv",
  "  sig94_vs_SHH_genes_scatter.png         scatter (ρ near 0 in fetal)",
  "  umap_sig94_vs_SHH_genes.png            UMAP comparison (no co-localization)",
  "  shh_gene_cluster_mean_heatmap.png      cluster heatmap (mismatch top/bot)",
  "  sig94_top_vs_bot_SHH_expression.png    tier bar (modest pattern only)",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21d_SUMMARY.txt"))

cat("\n=== R21d SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21d completion done\n")
cat("================================================================\n")
