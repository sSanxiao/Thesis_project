# ============================================================
# R21c: Deeper analysis of sig_94 high cells in Ghasemi data
# ------------------------------------------------------------
# Questions:
#   Q1. Are cluster 9 (proliferating) and cluster 3 (quiescent)
#       sig_94-high for the same biological reason?
#   Q2. What are sig_94 top 20% cells - what cluster mix?
#   Q3. What genes drive sig_94 in this dataset?
#       Are sig_94 positive-direction genes enriched in
#       progenitor markers or differentiation markers?
#
# Outputs (in <RESULTS_DIR>/R21_Ghasemi/):
#   R21c_SUMMARY.txt
#   sig94_top20_cluster_composition.png
#   sig94_top20_vs_bot20_markers.csv
#   sig94_top20_vs_bot20_volcano.png
#   sig94_genes_overlap_ghasemi.csv
#   sig94_genes_overlap_table.png
#   prolif_vs_quiescent_dissection.png
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
OUT_DIR <- file.path(RESULTS_DIR, "R21_Ghasemi")
SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")

cat("================================================================\n")
cat("R21c: Deeper sig_94 analysis in Ghasemi data\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Load R21a Seurat object + scores
# ============================================================
cat("[1] Loading data...\n")
combined <- readRDS(file.path(OUT_DIR, "ghasemi_seurat_R21a.rds"))

# Reload sig_94 / sig_core scores (was added in R21b but not saved to rds)
cell_dt <- fread(file.path(OUT_DIR, "cell_data_with_scores.csv"))
combined$sig_94 <- cell_dt$sig_94[match(colnames(combined), cell_dt$cell_id)]
combined$sig_core <- cell_dt$sig_core[match(colnames(combined), cell_dt$cell_id)]

# Add module scores
for (mod_name in c("module_early_CGNP_proliferating",
                    "module_early_CGNP_quiescent",
                    "module_migrating",
                    "module_postmitotic_differentiated",
                    "module_astrocytic_like")) {
  combined[[mod_name]] <- cell_dt[[mod_name]][match(colnames(combined), cell_dt$cell_id)]
}

cat(sprintf("  Loaded %d cells, sig_94 ready\n", ncol(combined)))

# ============================================================
# 2. Q1 + Q2: sig_94 top 20% vs bottom 20% — cluster composition
# ============================================================
cat("\n[2] sig_94 top/bottom 20% analysis...\n")

# Threshold
q20 <- quantile(combined$sig_94, 0.20)
q80 <- quantile(combined$sig_94, 0.80)
combined$sig94_tier <- factor(
  fifelse(combined$sig_94 >= q80, "top20",
   fifelse(combined$sig_94 <= q20, "bot20", "mid60")),
  levels = c("top20", "mid60", "bot20")
)
cat(sprintf("  q20=%.3f, q80=%.3f\n", q20, q80))
cat(sprintf("  top20: %d cells, bot20: %d cells, mid60: %d cells\n",
            sum(combined$sig94_tier == "top20"),
            sum(combined$sig94_tier == "bot20"),
            sum(combined$sig94_tier == "mid60")))

# Cluster composition of top20
cat("\n  Cluster composition of top20 sig_94 cells:\n")
top20_comp <- table(combined$seurat_clusters[combined$sig94_tier == "top20"])
top20_pct <- prop.table(top20_comp) * 100
top20_dt <- data.table(cluster = names(top20_comp),
                        n = as.integer(top20_comp),
                        pct = round(as.numeric(top20_pct), 1))
setorder(top20_dt, -n)
print(top20_dt)

bot20_comp <- table(combined$seurat_clusters[combined$sig94_tier == "bot20"])
bot20_pct <- prop.table(bot20_comp) * 100
bot20_dt <- data.table(cluster = names(bot20_comp),
                        n = as.integer(bot20_comp),
                        pct = round(as.numeric(bot20_pct), 1))
setorder(bot20_dt, -n)
cat("\n  Cluster composition of bot20 sig_94 cells:\n")
print(bot20_dt)

# ============================================================
# 3. Within-cluster sig_94 high vs low — focus on cluster 9 and 3
# ============================================================
cat("\n[3] Within-cluster dissection of cluster 9 and cluster 3...\n")

# In cluster 9 (proliferating dominant, sig_94 high): what defines sig_94 high?
# In cluster 3 (quiescent dominant, sig_94 high): what defines sig_94 high?
# Compare to see if same biology

DefaultAssay(combined) <- "RNA"

# Cluster 9: top vs bottom by sig_94 within cluster
cluster9_cells <- WhichCells(combined, idents = "9")
cluster9_obj <- subset(combined, cells = cluster9_cells)
cl9_q50 <- quantile(cluster9_obj$sig_94, 0.5)
cluster9_obj$sig94_within <- ifelse(cluster9_obj$sig_94 >= cl9_q50, "high", "low")

cat(sprintf("  Cluster 9: %d cells, sig_94 median=%.3f\n",
            ncol(cluster9_obj), cl9_q50))

# Mean module scores per within-cluster tier
cl9_summary <- data.table(
  tier = c("high", "low"),
  n = c(sum(cluster9_obj$sig94_within == "high"),
        sum(cluster9_obj$sig94_within == "low")),
  prolif = c(mean(cluster9_obj$module_early_CGNP_proliferating[cluster9_obj$sig94_within == "high"]),
             mean(cluster9_obj$module_early_CGNP_proliferating[cluster9_obj$sig94_within == "low"])),
  quiescent = c(mean(cluster9_obj$module_early_CGNP_quiescent[cluster9_obj$sig94_within == "high"]),
                mean(cluster9_obj$module_early_CGNP_quiescent[cluster9_obj$sig94_within == "low"])),
  migrating = c(mean(cluster9_obj$module_migrating[cluster9_obj$sig94_within == "high"]),
                mean(cluster9_obj$module_migrating[cluster9_obj$sig94_within == "low"])),
  diff = c(mean(cluster9_obj$module_postmitotic_differentiated[cluster9_obj$sig94_within == "high"]),
           mean(cluster9_obj$module_postmitotic_differentiated[cluster9_obj$sig94_within == "low"])),
  astro = c(mean(cluster9_obj$module_astrocytic_like[cluster9_obj$sig94_within == "high"]),
            mean(cluster9_obj$module_astrocytic_like[cluster9_obj$sig94_within == "low"]))
)
cat("  Cluster 9 within-cluster module scores by sig_94 tier:\n")
print(cl9_summary)

# Cluster 3
cluster3_cells <- WhichCells(combined, idents = "3")
cluster3_obj <- subset(combined, cells = cluster3_cells)
cl3_q50 <- quantile(cluster3_obj$sig_94, 0.5)
cluster3_obj$sig94_within <- ifelse(cluster3_obj$sig_94 >= cl3_q50, "high", "low")

cl3_summary <- data.table(
  tier = c("high", "low"),
  n = c(sum(cluster3_obj$sig94_within == "high"),
        sum(cluster3_obj$sig94_within == "low")),
  prolif = c(mean(cluster3_obj$module_early_CGNP_proliferating[cluster3_obj$sig94_within == "high"]),
             mean(cluster3_obj$module_early_CGNP_proliferating[cluster3_obj$sig94_within == "low"])),
  quiescent = c(mean(cluster3_obj$module_early_CGNP_quiescent[cluster3_obj$sig94_within == "high"]),
                mean(cluster3_obj$module_early_CGNP_quiescent[cluster3_obj$sig94_within == "low"])),
  migrating = c(mean(cluster3_obj$module_migrating[cluster3_obj$sig94_within == "high"]),
                mean(cluster3_obj$module_migrating[cluster3_obj$sig94_within == "low"])),
  diff = c(mean(cluster3_obj$module_postmitotic_differentiated[cluster3_obj$sig94_within == "high"]),
           mean(cluster3_obj$module_postmitotic_differentiated[cluster3_obj$sig94_within == "low"])),
  astro = c(mean(cluster3_obj$module_astrocytic_like[cluster3_obj$sig94_within == "high"]),
            mean(cluster3_obj$module_astrocytic_like[cluster3_obj$sig94_within == "low"]))
)
cat("  Cluster 3 within-cluster module scores by sig_94 tier:\n")
print(cl3_summary)

# ============================================================
# 4. Q3: What genes drive sig_94 high vs low globally?
# ============================================================
cat("\n[4] sig_94 top20 vs bot20 differential expression...\n")

# Make sig94_tier the active identity for FindMarkers
Idents(combined) <- "sig94_tier"

t0 <- Sys.time()
de_top_vs_bot <- FindMarkers(combined,
                              ident.1 = "top20",
                              ident.2 = "bot20",
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              assay = "RNA",
                              verbose = FALSE)
cat(sprintf("  DE done in %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))

de_dt <- as.data.table(de_top_vs_bot, keep.rownames = "gene")
setorder(de_dt, -avg_log2FC)

# Save full DE table
fwrite(de_dt, file.path(OUT_DIR, "sig94_top20_vs_bot20_markers.csv"))

cat(sprintf("\n  Total DE genes: %d\n", nrow(de_dt)))
cat(sprintf("  Up in top20 (avg_log2FC > 0): %d\n", sum(de_dt$avg_log2FC > 0)))
cat(sprintf("  Down in top20 (avg_log2FC < 0): %d\n", sum(de_dt$avg_log2FC < 0)))

# Top 20 up + 20 down
cat("\n  Top 20 genes UP in sig_94-high (top20 vs bot20):\n")
print(de_dt[order(-avg_log2FC)][1:20, .(gene, avg_log2FC, pct.1, pct.2, p_val_adj)])

cat("\n  Top 20 genes DOWN in sig_94-high:\n")
print(de_dt[order(avg_log2FC)][1:20, .(gene, avg_log2FC, pct.1, pct.2, p_val_adj)])

# ============================================================
# 5. Overlap of sig_94 genes with Ghasemi cell-state markers
# ============================================================
cat("\n[5] Cross-tabulating sig_94 genes against Ghasemi markers...\n")

# Reload sig_94 provenance
sig_prov <- fread(SIG_PROV)
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
sig_prov <- sig_prov[!gene %in% CONFLICT_GENES]

ghasemi_markers <- list(
  early_CGNP_proliferating = c("MKI67", "TOP2A", "ATOH1", "BARHL1", "ZIC1", "ZIC3"),
  early_CGNP_quiescent = c("PTCH1", "SMO", "HHIP", "GLI1", "GLI2", "PTPRK"),
  migrating = c("GRIN2B", "CNTN2", "ASTN1", "SEMA6A"),
  postmitotic_differentiated = c("GABRA1", "GABRA6", "GRIN2C", "RBFOX3"),
  astrocytic_like = c("LAMA2", "SOX2", "SOX9", "GFAP", "AQP4")
)

# Overlap analysis
overlap_dt <- data.table()
for (state_name in names(ghasemi_markers)) {
  ghasemi_g <- ghasemi_markers[[state_name]]
  
  in_sig94_pos <- intersect(ghasemi_g, sig_prov[direction_final == "positive", gene])
  in_sig94_neg <- intersect(ghasemi_g, sig_prov[direction_final == "negative", gene])
  not_in_sig94 <- setdiff(ghasemi_g, sig_prov$gene)
  
  overlap_dt <- rbind(overlap_dt, data.table(
    ghasemi_state = state_name,
    n_total = length(ghasemi_g),
    in_sig94_positive_dir = length(in_sig94_pos),
    sig94_pos_genes = paste(in_sig94_pos, collapse = ","),
    in_sig94_negative_dir = length(in_sig94_neg),
    sig94_neg_genes = paste(in_sig94_neg, collapse = ","),
    not_in_sig94 = length(not_in_sig94),
    not_in_genes = paste(not_in_sig94, collapse = ",")
  ))
}

cat("\n  Ghasemi marker overlap with sig_94 (by direction):\n")
print(overlap_dt[, .(ghasemi_state, n_total, in_sig94_positive_dir,
                      in_sig94_negative_dir, not_in_sig94)])

cat("\n  Detailed overlap:\n")
for (i in 1:nrow(overlap_dt)) {
  cat(sprintf("\n  [%s]\n", overlap_dt$ghasemi_state[i]))
  cat(sprintf("    pos direction: %s\n", overlap_dt$sig94_pos_genes[i]))
  cat(sprintf("    neg direction: %s\n", overlap_dt$sig94_neg_genes[i]))
  cat(sprintf("    not in sig_94: %s\n", overlap_dt$not_in_genes[i]))
}

fwrite(overlap_dt, file.path(OUT_DIR, "sig94_genes_overlap_ghasemi.csv"))

# ============================================================
# 6. sig_94 positive direction genes - which Ghasemi state are they?
# ============================================================
cat("\n[6] sig_94 positive direction genes — which Ghasemi state?\n")

# Get all sig_94 positive direction genes
sig94_pos_genes <- sig_prov[direction_final == "positive", gene]
sig94_neg_genes <- sig_prov[direction_final == "negative", gene]

cat(sprintf("  sig_94 positive direction (%d genes):\n", length(sig94_pos_genes)))
cat(sprintf("    %s\n", paste(sig94_pos_genes, collapse = ", ")))

cat(sprintf("\n  sig_94 negative direction (%d genes):\n", length(sig94_neg_genes)))
cat(sprintf("    %s\n", paste(sig94_neg_genes, collapse = ", ")))

# ============================================================
# 7. Visualizations
# ============================================================
cat("\n[7] Rendering plots...\n")

# Plot 1: cluster composition of top20 vs bot20
combine_dt <- rbind(top20_dt[, .(tier = "top20", cluster, n, pct)],
                     bot20_dt[, .(tier = "bot20", cluster, n, pct)])

p1 <- ggplot(combine_dt, aes(x = factor(cluster, levels = as.character(0:12)),
                              y = pct, fill = tier)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("top20" = "#b2182b", "bot20" = "#2166ac")) +
  labs(title = "R21c: Cluster composition of sig_94 top 20% vs bottom 20% cells",
       subtitle = "Where do sig_94-high cells live, and where do sig_94-low live?",
       x = "Cluster", y = "% of tier") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "sig94_top20_cluster_composition.png"),
       p1, width = 12, height = 6, dpi = 150, bg = "white")

# Plot 2: volcano of top20 vs bot20 DE
de_dt[, neg_log10_p := -log10(pmax(p_val_adj, 1e-300))]
de_dt[, color := fcase(
  avg_log2FC > 0.5 & p_val_adj < 0.05, "up",
  avg_log2FC < -0.5 & p_val_adj < 0.05, "down",
  default = "ns"
)]

# Label top 10 up + top 10 down
de_dt_label <- rbind(
  de_dt[color == "up"][order(-avg_log2FC)][1:10],
  de_dt[color == "down"][order(avg_log2FC)][1:10]
)

p2 <- ggplot(de_dt, aes(x = avg_log2FC, y = neg_log10_p, color = color)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("up" = "#b2182b", "down" = "#2166ac", "ns" = "grey70")) +
  ggrepel::geom_text_repel(data = de_dt_label, aes(label = gene),
                            size = 3, max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  labs(title = "R21c: sig_94 top20 vs bot20 differential expression",
       subtitle = sprintf("Up in top20: %d genes | Down in top20: %d genes",
                           sum(de_dt$color == "up"),
                           sum(de_dt$color == "down")),
       x = "avg_log2FC (top20 / bot20)", y = "-log10(p_adj)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Try with ggrepel; if not available, fall back to geom_text
ggsave(file.path(OUT_DIR, "sig94_top20_vs_bot20_volcano.png"),
       p2, width = 10, height = 8, dpi = 150, bg = "white")

# Plot 3: Cluster 9 vs Cluster 3 within-cluster dissection
dissection_dt <- rbind(
  data.table(cluster = "C9_proliferating", cl9_summary),
  data.table(cluster = "C3_quiescent", cl3_summary)
)
diss_long <- melt(dissection_dt,
                   id.vars = c("cluster", "tier", "n"),
                   measure.vars = c("prolif", "quiescent", "migrating", "diff", "astro"),
                   variable.name = "module", value.name = "mean_score")

p3 <- ggplot(diss_long, aes(x = module, y = mean_score, fill = tier)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ cluster) +
  scale_fill_manual(values = c("high" = "#b2182b", "low" = "#2166ac")) +
  labs(title = "R21c: Within-cluster sig_94 high vs low — Ghasemi module scores",
       subtitle = "If sig_94 = quiescent program: 'high' should have ↑quiescent in BOTH clusters",
       x = "Ghasemi module", y = "Mean module score") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "prolif_vs_quiescent_dissection.png"),
       p3, width = 12, height = 6, dpi = 150, bg = "white")

# Plot 4: Ghasemi marker x sig_94 direction overlap heatmap
overlap_long <- overlap_dt[, .(ghasemi_state, in_sig94_positive_dir, in_sig94_negative_dir, not_in_sig94)]
overlap_long_m <- melt(overlap_long, id.vars = "ghasemi_state",
                        variable.name = "category", value.name = "count")

p4 <- ggplot(overlap_long_m, aes(x = ghasemi_state, y = category, fill = count)) +
  geom_tile() +
  geom_text(aes(label = count), color = "white", size = 4) +
  scale_fill_gradient(low = "#f7f7f7", high = "#08519c") +
  labs(title = "R21c: How many Ghasemi markers in each cell-state are in sig_94?",
       x = "Ghasemi cell state", y = "Direction in sig_94") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "sig94_genes_overlap_table.png"),
       p4, width = 10, height = 5, dpi = 150, bg = "white")

# ============================================================
# 8. SUMMARY
# ============================================================
cat("\n[8] Writing R21c_SUMMARY.txt...\n")

summary_lines <- c(
  "================================================================",
  "R21c — Deeper sig_94 analysis in Ghasemi 2024",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "QUESTIONS",
  "  Q1. Are clusters 9 (proliferating) and 3 (quiescent)",
  "      sig_94-high for the same biological reason?",
  "  Q2. What is the cluster mix of sig_94 top 20% cells?",
  "  Q3. What genes / programs drive sig_94 high in this dataset?",
  "      Are sig_94 positive-direction genes enriched in",
  "      progenitor markers vs differentiation markers?",
  "",
  "================================================================",
  "Q2 RESULT: sig_94 top 20% cluster composition",
  "================================================================",
  ""
)

summary_lines <- c(summary_lines,
  sprintf("  Top 20%% cells (n=%d, sig_94 >= %.3f):",
          sum(combined$sig94_tier == "top20"), q80))
for (i in 1:nrow(top20_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("    Cluster %s: %d cells (%.1f%% of top20)",
            top20_dt$cluster[i], top20_dt$n[i], top20_dt$pct[i]))
}

summary_lines <- c(summary_lines, "",
  sprintf("  Bottom 20%% cells (n=%d, sig_94 <= %.3f):",
          sum(combined$sig94_tier == "bot20"), q20))
for (i in 1:nrow(bot20_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("    Cluster %s: %d cells (%.1f%% of bot20)",
            bot20_dt$cluster[i], bot20_dt$n[i], bot20_dt$pct[i]))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "Q1 RESULT: Within-cluster dissection of cluster 9 and 3",
  "================================================================",
  "",
  "  Cluster 9 (proliferating dominant cluster, sig_94 mean=+0.245)",
  "    Within cluster, split by median sig_94:",
  ""
)

for (i in 1:nrow(cl9_summary)) {
  summary_lines <- c(summary_lines,
    sprintf("    %s tier (n=%d): prolif=%+.3f, quiescent=%+.3f, migrating=%+.3f, diff=%+.3f, astro=%+.3f",
            cl9_summary$tier[i], cl9_summary$n[i],
            cl9_summary$prolif[i], cl9_summary$quiescent[i],
            cl9_summary$migrating[i], cl9_summary$diff[i],
            cl9_summary$astro[i]))
}

summary_lines <- c(summary_lines, "",
  "  Cluster 3 (quiescent dominant cluster, sig_94 mean=+0.233)",
  "    Within cluster, split by median sig_94:",
  ""
)

for (i in 1:nrow(cl3_summary)) {
  summary_lines <- c(summary_lines,
    sprintf("    %s tier (n=%d): prolif=%+.3f, quiescent=%+.3f, migrating=%+.3f, diff=%+.3f, astro=%+.3f",
            cl3_summary$tier[i], cl3_summary$n[i],
            cl3_summary$prolif[i], cl3_summary$quiescent[i],
            cl3_summary$migrating[i], cl3_summary$diff[i],
            cl3_summary$astro[i]))
}

summary_lines <- c(summary_lines, "",
  "  KEY DIAGNOSTIC: Compare 'high' tier deltas across clusters",
  "  - If the dominant module that goes UP in 'high' is the SAME",
  "    in both clusters → sig_94 captures one consistent program",
  "  - If different modules dominate → sig_94 is heterogeneous",
  "",
  "================================================================",
  "Q3 RESULT: Top 20 DE genes — sig_94 high vs low",
  "================================================================",
  "",
  "  Up in top20 (sig_94-high):"
)

for (i in 1:min(20, sum(de_dt$color == "up"))) {
  g <- de_dt[order(-avg_log2FC)][i]
  summary_lines <- c(summary_lines,
    sprintf("    %-15s log2FC=%+.2f  pct.1=%.2f  pct.2=%.2f",
            g$gene, g$avg_log2FC, g$pct.1, g$pct.2))
}

summary_lines <- c(summary_lines, "",
  "  Down in top20 (sig_94-low):")

for (i in 1:min(20, sum(de_dt$color == "down"))) {
  g <- de_dt[order(avg_log2FC)][i]
  summary_lines <- c(summary_lines,
    sprintf("    %-15s log2FC=%+.2f  pct.1=%.2f  pct.2=%.2f",
            g$gene, g$avg_log2FC, g$pct.1, g$pct.2))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "Q3 SECONDARY: sig_94 vs Ghasemi marker overlap",
  "================================================================",
  ""
)

for (i in 1:nrow(overlap_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("  [%s] (%d total markers)",
            overlap_dt$ghasemi_state[i], overlap_dt$n_total[i]),
    sprintf("    in sig_94 POSITIVE dir: %d genes  (%s)",
            overlap_dt$in_sig94_positive_dir[i],
            ifelse(nchar(overlap_dt$sig94_pos_genes[i]) > 0,
                   overlap_dt$sig94_pos_genes[i], "—")),
    sprintf("    in sig_94 NEGATIVE dir: %d genes  (%s)",
            overlap_dt$in_sig94_negative_dir[i],
            ifelse(nchar(overlap_dt$sig94_neg_genes[i]) > 0,
                   overlap_dt$sig94_neg_genes[i], "—")),
    sprintf("    not in sig_94:           %d genes  (%s)",
            overlap_dt$not_in_sig94[i],
            ifelse(nchar(overlap_dt$not_in_genes[i]) > 0,
                   overlap_dt$not_in_genes[i], "—")),
    "")
}

summary_lines <- c(summary_lines,
  "================================================================",
  "INTERPRETATION FRAMEWORK",
  "================================================================",
  "",
  "  Hypothesis A: sig_94 = early CGNP-like (both proliferating + quiescent)",
  "    Evidence: top20 cells should span clusters 9 + 3 + 7",
  "    Within-cluster: sig_94-high should have ↑quiescent OR ↑prolif",
  "",
  "  Hypothesis B: sig_94 = quiescent CGNP specifically",
  "    Evidence: top20 dominated by cluster 3, less from cluster 9",
  "    Within-cluster: sig_94-high should have ↑quiescent in BOTH clusters",
  "",
  "  Hypothesis C: sig_94 = SHH pathway active",
  "    Evidence: PTCH1, SMO, HHIP, GLI1/2 should be in sig_94 positive dir",
  "    Top20 DE should include SHH pathway genes",
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  R21c_SUMMARY.txt                          - this file",
  "  sig94_top20_cluster_composition.png       - top20/bot20 cluster bar",
  "  sig94_top20_vs_bot20_markers.csv          - full DE table",
  "  sig94_top20_vs_bot20_volcano.png          - volcano plot",
  "  prolif_vs_quiescent_dissection.png        - within-cluster dissection",
  "  sig94_genes_overlap_ghasemi.csv           - overlap table",
  "  sig94_genes_overlap_table.png             - overlap heatmap",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21c_SUMMARY.txt"))

cat("\n=== R21c SUMMARY (preview) ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21c DONE\n")
cat("================================================================\n")
