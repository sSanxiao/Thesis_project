# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R17b patch: fix AddModuleScore + complete all downstream
# ------------------------------------------------------------
# Fix: ctrl=50 (not default 100) because only 2000 genes / 25 bins
# Reuses the Seurat object from R17a (aldinger2021_seurat.rds)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cat("================================================================\n")
cat("R17b patch: complete scoring with ctrl=50 for AddModuleScore\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

RES_DIR <- "./results/R17_Aldinger2021"
FIG_DIR <- file.path(RES_DIR, "figures_b")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

INPUT_RDS <- file.path(RES_DIR, "aldinger2021_seurat.rds")
SIG_PROV <- "./results/R12_Gaps/sig_strict_99_provenance.csv"
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
SIG_CORE_GENES <- c("TUBB4A", "APOE", "CCR7", "EOMES", "ST18", "NES", "AQP4", "QKI")

AMS_CTRL <- 50  # key fix: ctrl < bin size (2000/25=80)

cat("[1] Loading Seurat object...\n")
obj <- readRDS(INPUT_RDS)
cat(sprintf("  %d genes x %d cells\n", nrow(obj), ncol(obj)))

sig_prov <- fread(SIG_PROV)
sig_94 <- setdiff(sig_prov$gene, CONFLICT_GENES)
genes_in_obj <- rownames(obj)
match_94 <- intersect(sig_94, genes_in_obj)
match_core <- intersect(SIG_CORE_GENES, genes_in_obj)

dir_table <- sig_prov[gene %in% match_94, .(gene, direction_final)]
setkey(dir_table, gene)
dir_vec <- dir_table[J(match_94), direction_final]
sig_94_pos <- match_94[dir_vec == "positive"]
sig_94_neg <- match_94[dir_vec == "negative"]

core_dir_table <- sig_prov[gene %in% match_core, .(gene, direction_final)]
setkey(core_dir_table, gene)
core_dir_vec <- core_dir_table[J(match_core), direction_final]
sig_core_pos <- match_core[core_dir_vec == "positive"]
sig_core_neg <- match_core[core_dir_vec == "negative"]

cat(sprintf("  sig_94: %d (%d pos + %d neg)\n",
            length(match_94), length(sig_94_pos), length(sig_94_neg)))
cat(sprintf("  sig_core: %d (%d pos + %d neg)\n",
            length(match_core), length(sig_core_pos), length(sig_core_neg)))

cat("\n[2] zscore-weighted scoring...\n")
expr <- GetAssayData(obj, assay = "RNA", layer = "scale.data")

compute_weighted_score <- function(expr, pos_genes, neg_genes) {
  pos_score <- if (length(pos_genes) > 0) {
    colMeans(expr[pos_genes, , drop = FALSE], na.rm = TRUE)
  } else rep(0, ncol(expr))
  neg_score <- if (length(neg_genes) > 0) {
    colMeans(expr[neg_genes, , drop = FALSE], na.rm = TRUE)
  } else rep(0, ncol(expr))
  n_pos <- length(pos_genes); n_neg <- length(neg_genes)
  (pos_score * n_pos - neg_score * n_neg) / (n_pos + n_neg)
}

score_94_zscore <- compute_weighted_score(expr, sig_94_pos, sig_94_neg)
score_core_zscore <- compute_weighted_score(expr, sig_core_pos, sig_core_neg)

cat(sprintf("  sig_94 zscore: [%.3f, %.3f] median=%.3f\n",
            min(score_94_zscore), max(score_94_zscore), median(score_94_zscore)))
cat(sprintf("  sig_core zscore: [%.3f, %.3f] median=%.3f\n",
            min(score_core_zscore), max(score_core_zscore), median(score_core_zscore)))

cat(sprintf("\n[3] AddModuleScore with ctrl=%d...\n", AMS_CTRL))

if (length(sig_94_pos) > 0 && length(sig_94_neg) > 0) {
  obj <- AddModuleScore(obj,
                        features = list(sig94_pos = sig_94_pos,
                                        sig94_neg = sig_94_neg),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sig94_")
  obj$sig_94_AMS <- obj$sig94_1 - obj$sig94_2
} else if (length(sig_94_pos) > 0) {
  obj <- AddModuleScore(obj, features = list(sig94_pos = sig_94_pos),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sig94_pos_")
  obj$sig_94_AMS <- obj$sig94_pos_1
} else {
  obj <- AddModuleScore(obj, features = list(sig94_neg = sig_94_neg),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sig94_neg_")
  obj$sig_94_AMS <- -obj$sig94_neg_1
}

if (length(sig_core_pos) > 0 && length(sig_core_neg) > 0) {
  obj <- AddModuleScore(obj,
                        features = list(sigcore_pos = sig_core_pos,
                                        sigcore_neg = sig_core_neg),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sigcore_")
  obj$sig_core_AMS <- obj$sigcore_1 - obj$sigcore_2
} else if (length(sig_core_neg) > 0) {
  obj <- AddModuleScore(obj, features = list(sigcore_neg = sig_core_neg),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sigcore_neg_")
  obj$sig_core_AMS <- -obj$sigcore_neg_1
} else {
  obj <- AddModuleScore(obj, features = list(sigcore_pos = sig_core_pos),
                        assay = "RNA", ctrl = AMS_CTRL, name = "sigcore_pos_")
  obj$sig_core_AMS <- obj$sigcore_pos_1
}

obj$sig_94_zscore <- score_94_zscore[colnames(obj)]
obj$sig_core_zscore <- score_core_zscore[colnames(obj)]

cat(sprintf("  sig_94 AMS: [%.3f, %.3f]\n",
            min(obj$sig_94_AMS), max(obj$sig_94_AMS)))
cat(sprintf("  sig_core AMS: [%.3f, %.3f]\n",
            min(obj$sig_core_AMS), max(obj$sig_core_AMS)))

cat("\n[4] Method concordance...\n")
cor_94 <- cor(obj$sig_94_zscore, obj$sig_94_AMS, method = "spearman")
cor_core <- cor(obj$sig_core_zscore, obj$sig_core_AMS, method = "spearman")
cat(sprintf("  sig_94 ρ = %.3f\n", cor_94))
cat(sprintf("  sig_core ρ = %.3f\n", cor_core))

cat("\n[5] Per-cluster stats...\n")
md <- obj@meta.data
md_dt <- as.data.table(md)

cluster_stats <- md_dt[, .(
  n = .N,
  sig_94_zscore_mean = mean(sig_94_zscore),
  sig_94_zscore_median = median(sig_94_zscore),
  sig_94_zscore_sd = sd(sig_94_zscore),
  sig_94_AMS_mean = mean(sig_94_AMS),
  sig_core_zscore_mean = mean(sig_core_zscore),
  sig_core_AMS_mean = mean(sig_core_AMS)
), by = cell_type][order(-sig_94_zscore_mean)]

print(cluster_stats)
fwrite(cluster_stats, file.path(RES_DIR, "cluster_score_stats.csv"))

cat("\n[6] ANOVA + Kruskal-Wallis:\n")
stats_dt <- data.table()
for (scorename in c("sig_94_zscore", "sig_94_AMS", "sig_core_zscore", "sig_core_AMS")) {
  kw <- kruskal.test(md_dt[[scorename]] ~ as.factor(md_dt$cell_type))
  anova_fit <- aov(md_dt[[scorename]] ~ as.factor(md_dt$cell_type))
  anova_p <- summary(anova_fit)[[1]][["Pr(>F)"]][1]
  cat(sprintf("  %-18s ANOVA p=%.2e  Kruskal p=%.2e\n",
              scorename, anova_p, kw$p.value))
  stats_dt <- rbind(stats_dt, data.table(score = scorename,
                                          anova_p = anova_p,
                                          kruskal_p = kw$p.value))
}
fwrite(stats_dt, file.path(RES_DIR, "cluster_ANOVA_stats.csv"))

cat("\n[7] Figures...\n")

SHH_RELEVANT <- c("02-RL", "03-GCP", "04-GN", "05-eCN/UBC", "08-BG", "09-Ast")
md_dt[, is_shh_relevant := cell_type %in% SHH_RELEVANT]
cluster_order <- cluster_stats$cell_type
md_dt[, cell_type_ord := factor(cell_type, levels = cluster_order)]

plot_box <- function(score_col, title) {
  ggplot(md_dt, aes(x = cell_type_ord, y = .data[[score_col]],
                     fill = is_shh_relevant)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, alpha = 0.85) +
    scale_fill_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "#999999"),
                      name = "SHH-relevant") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    labs(x = NULL, y = score_col, title = title) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

p1 <- plot_box("sig_94_zscore", "sig_94 (zscore-weighted)")
p2 <- plot_box("sig_94_AMS", "sig_94 (AddModuleScore, ctrl=50)")
p3 <- plot_box("sig_core_zscore", "sig_core (zscore-weighted)")
p4 <- plot_box("sig_core_AMS", "sig_core (AddModuleScore, ctrl=50)")

combined_box <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = sprintf("Aldinger 2021 fetal cerebellum (n=%d, %d cell types)",
                    ncol(obj), length(unique(obj$cell_type))),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "01_cluster_boxplot_4panel.png"),
       combined_box, width = 16, height = 12, dpi = 150)
cat("  Saved: 01_cluster_boxplot_4panel.png\n")

age_cluster_94 <- md_dt[, .(mean_score = mean(sig_94_zscore), n = .N),
                          by = .(cell_type_ord, age_pcw)]
age_cluster_94[, mean_score_show := ifelse(n < 20, NA_real_, mean_score)]

p_heat <- ggplot(age_cluster_94, aes(x = factor(age_pcw), y = cell_type_ord)) +
  geom_tile(aes(fill = mean_score_show)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, name = "sig_94\nscore",
                        na.value = "grey95") +
  geom_text(aes(label = ifelse(n < 20, "", sprintf("%d", n))),
            size = 2.5, color = "grey30") +
  labs(x = "Age (PCW)", y = NULL,
       title = "sig_94 (zscore) mean: cluster × age",
       subtitle = "Numbers = cells; grey = <20 cells") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "02_cluster_age_heatmap_sig94.png"),
       p_heat, width = 10, height = 7, dpi = 150)
cat("  Saved: 02_cluster_age_heatmap_sig94.png\n")

age_cluster_core <- md_dt[, .(mean_score = mean(sig_core_zscore), n = .N),
                            by = .(cell_type_ord, age_pcw)]
age_cluster_core[, mean_score_show := ifelse(n < 20, NA_real_, mean_score)]

p_heat_core <- ggplot(age_cluster_core, aes(x = factor(age_pcw), y = cell_type_ord)) +
  geom_tile(aes(fill = mean_score_show)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, name = "sig_core\nscore",
                        na.value = "grey95") +
  geom_text(aes(label = ifelse(n < 20, "", sprintf("%d", n))),
            size = 2.5, color = "grey30") +
  labs(x = "Age (PCW)", y = NULL,
       title = "sig_core (zscore) mean: cluster × age",
       subtitle = "Numbers = cells; grey = <20 cells") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "03_cluster_age_heatmap_sigcore.png"),
       p_heat_core, width = 10, height = 7, dpi = 150)
cat("  Saved: 03_cluster_age_heatmap_sigcore.png\n")

umap_dt <- data.table(
  UMAP_1 = obj[["umap"]]@cell.embeddings[, 1],
  UMAP_2 = obj[["umap"]]@cell.embeddings[, 2],
  sig_94_zscore = obj$sig_94_zscore,
  sig_core_zscore = obj$sig_core_zscore,
  cell_type = obj$cell_type
)

p_umap_ct <- ggplot(umap_dt, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.6) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "Aldinger 2021 UMAP (cell types)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))

p_umap_94 <- ggplot(umap_dt, aes(x = UMAP_1, y = UMAP_2, color = sig_94_zscore)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b",
                         midpoint = 0, name = "sig_94") +
  labs(title = "sig_94 score on UMAP") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_umap_core <- ggplot(umap_dt, aes(x = UMAP_1, y = UMAP_2, color = sig_core_zscore)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b",
                         midpoint = 0, name = "sig_core") +
  labs(title = "sig_core score on UMAP") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

umap_combined <- p_umap_ct + p_umap_94 + p_umap_core +
  plot_layout(ncol = 3, widths = c(1.5, 1, 1))

ggsave(file.path(FIG_DIR, "04_UMAP_overlay.png"),
       umap_combined, width = 18, height = 6, dpi = 150)
cat("  Saved: 04_UMAP_overlay.png\n")

conc_dt <- data.table(zscore = obj$sig_94_zscore,
                       AMS = obj$sig_94_AMS,
                       cell_type = obj$cell_type)
set.seed(42)
conc_sub <- conc_dt[sample(.N, min(5000, .N))]

p_conc <- ggplot(conc_sub, aes(x = zscore, y = AMS)) +
  geom_point(size = 0.5, alpha = 0.4, color = "#377eb8") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5) +
  labs(x = "sig_94 (zscore-weighted)",
       y = "sig_94 (AddModuleScore)",
       title = "Method concordance (5000-cell subsample)",
       subtitle = sprintf("Spearman ρ = %.3f across all %d cells",
                           cor_94, nrow(conc_dt))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "05_method_concordance.png"),
       p_conc, width = 7, height = 6, dpi = 150)
cat("  Saved: 05_method_concordance.png\n")

cat("\n[8] Saving scored object...\n")
saveRDS(obj, file.path(RES_DIR, "aldinger2021_scored.rds"))
cat(sprintf("  Saved (%.1f MB)\n",
            file.size(file.path(RES_DIR, "aldinger2021_scored.rds")) / 1e6))

summary_lines <- c(
  "================================================================",
  "R17b — Signature Scoring on Aldinger 2021 (fetal cerebellum)",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "COHORT",
  sprintf("  Aldinger 2021 Nat Neurosci, 21 cell types, 10 ages (9-21 PCW)"),
  sprintf("  %d cells x %d genes (UCSC Cell Browser variable genes)",
          ncol(obj), nrow(obj)),
  "",
  "SIGNATURE MATCHING",
  sprintf("  sig_94:   %d/%d genes matched (%d pos + %d neg)",
          length(match_94), length(sig_94),
          length(sig_94_pos), length(sig_94_neg)),
  sprintf("  sig_core: %d/%d genes matched (%d pos + %d neg)",
          length(match_core), length(SIG_CORE_GENES),
          length(sig_core_pos), length(sig_core_neg)),
  sprintf("  NOTE: sig_core matches are all negative-direction genes"),
  sprintf("        (APOE, AQP4, EOMES, QKI, ST18)"),
  "",
  "METHOD CONCORDANCE (zscore-weighted vs AddModuleScore, ctrl=50)",
  sprintf("  sig_94   Spearman ρ = %.3f", cor_94),
  sprintf("  sig_core Spearman ρ = %.3f", cor_core),
  "",
  "CELL-TYPE EFFECT SIZE (ANOVA / Kruskal-Wallis)"
)
for (i in 1:nrow(stats_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-18s ANOVA p=%.2e  Kruskal p=%.2e",
            stats_dt$score[i], stats_dt$anova_p[i], stats_dt$kruskal_p[i]))
}

summary_lines <- c(summary_lines, "",
  "CLUSTER RANKING BY sig_94 zscore (DESCENDING MEAN)")
for (i in 1:nrow(cluster_stats)) {
  is_shh <- cluster_stats$cell_type[i] %in% SHH_RELEVANT
  mark <- ifelse(is_shh, " ★", "")
  summary_lines <- c(summary_lines,
    sprintf("  %2d. %-15s n=%5d  sig_94_mean=%+.3f  sd=%.3f%s",
            i, cluster_stats$cell_type[i], cluster_stats$n[i],
            cluster_stats$sig_94_zscore_mean[i],
            cluster_stats$sig_94_zscore_sd[i], mark))
}
summary_lines <- c(summary_lines, "",
  "  ★ = SHH-relevant (02-RL, 03-GCP, 04-GN, 05-eCN/UBC, 08-BG, 09-Ast)")

summary_lines <- c(summary_lines, "",
  "OUTPUTS",
  sprintf("  Figures: %s/figures_b/", RES_DIR),
  sprintf("  Stats: cluster_score_stats.csv, cluster_ANOVA_stats.csv"),
  sprintf("  Scored rds: aldinger2021_scored.rds"),
  "",
  "READY FOR R17c (neuronal lineage pseudotime)",
  "================================================================"
)

writeLines(summary_lines, file.path(RES_DIR, "R17b_SUMMARY.txt"))
cat(sprintf("\n  SUMMARY: R17b_SUMMARY.txt\n"))

cat("\n================================================================\n")
cat("R17b patch DONE\n")
cat("================================================================\n")
