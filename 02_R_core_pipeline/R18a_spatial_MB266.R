# ============================================================
# R18a: Spatial signature visualization — MB266 probe
# ------------------------------------------------------------
# Single-sample test before scaling to all 4 MB samples.
#
# Goal: render sig_94 and sig_core scores in (x, y) Xenium
# spatial coordinates; verify visualization parameters work.
#
# After this passes review → R18b extends to MB263/295/299.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ---- paths ----
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
SAMPLE_RDS <- file.path(RESULTS_DIR, "R2_Results", "Medulloblastoma_Human", "GSM8840047", "GSM8840047_seurat_R2.rds")
SAMPLE_NAME <- "MB266"

OUT_DIR <- file.path(RESULTS_DIR, "R18_SpatialSignature")
FIG_DIR <- file.path(OUT_DIR, "MB266")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
SIG_CORE_GENES <- c("TUBB4A", "APOE", "CCR7", "EOMES", "ST18", "NES", "AQP4", "QKI")

cat("================================================================\n")
cat(sprintf("R18a: Spatial probe — %s\n", SAMPLE_NAME))
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# 1. Load
# ============================================================
cat("[1] Loading Seurat object...\n")
t0 <- Sys.time()
obj <- readRDS(SAMPLE_RDS)
cat(sprintf("  Loaded in %.1f sec\n", as.numeric(Sys.time() - t0, units = "secs")))
cat(sprintf("  %d genes × %d cells\n", nrow(obj), ncol(obj)))

# Verify required metadata columns
required <- c("x_centroid", "y_centroid", "density_knn_main_piecewise",
              "data_quality_tier", "seurat_clusters")
missing <- setdiff(required, colnames(obj@meta.data))
if (length(missing) > 0) {
  cat(sprintf("  [!] Missing columns: %s\n", paste(missing, collapse=", ")))
  quit(status = 1)
}
cat("  All required metadata columns present\n")

# ============================================================
# 2. Compute signature scores
# ============================================================
cat("\n[2] Computing signature scores (zscore-weighted)...\n")

sig_prov <- fread(SIG_PROV)
sig_94 <- setdiff(sig_prov$gene, CONFLICT_GENES)

# Match against panel
panel_genes <- rownames(obj)
match_94 <- intersect(sig_94, panel_genes)
match_core <- intersect(SIG_CORE_GENES, panel_genes)

cat(sprintf("  sig_94 matched: %d / %d (%.1f%%)\n",
            length(match_94), length(sig_94),
            100 * length(match_94) / length(sig_94)))
cat(sprintf("  sig_core matched: %d / %d\n",
            length(match_core), length(SIG_CORE_GENES)))

# Get scale.data (SCT residuals) for scoring
expr_sct <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
if (is.null(expr_sct) || all(dim(expr_sct) == 0)) {
  cat("  scale.data empty, using normalized data...\n")
  expr_sct <- GetAssayData(obj, assay = "SCT", layer = "data")
}
cat(sprintf("  Expression matrix: %d × %d (range [%.2f, %.2f])\n",
            nrow(expr_sct), ncol(expr_sct),
            min(expr_sct), max(expr_sct)))

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

cat(sprintf("  sig_94: %d pos + %d neg\n",
            length(sig_94_pos), length(sig_94_neg)))
cat(sprintf("  sig_core: %d pos + %d neg\n",
            length(sig_core_pos), length(sig_core_neg)))

compute_score <- function(expr, pos_genes, neg_genes) {
  pos_genes <- intersect(pos_genes, rownames(expr))
  neg_genes <- intersect(neg_genes, rownames(expr))
  pos_score <- if (length(pos_genes) > 0) {
    colMeans(expr[pos_genes, , drop = FALSE], na.rm = TRUE)
  } else rep(0, ncol(expr))
  neg_score <- if (length(neg_genes) > 0) {
    colMeans(expr[neg_genes, , drop = FALSE], na.rm = TRUE)
  } else rep(0, ncol(expr))
  np <- length(pos_genes); nn <- length(neg_genes)
  if (np + nn == 0) return(rep(NA, ncol(expr)))
  (pos_score * np - neg_score * nn) / (np + nn)
}

obj$sig_94 <- compute_score(expr_sct, sig_94_pos, sig_94_neg)
obj$sig_core <- compute_score(expr_sct, sig_core_pos, sig_core_neg)

cat(sprintf("  sig_94 score range: [%.3f, %.3f], median=%.3f\n",
            min(obj$sig_94, na.rm=TRUE),
            max(obj$sig_94, na.rm=TRUE),
            median(obj$sig_94, na.rm=TRUE)))
cat(sprintf("  sig_core score range: [%.3f, %.3f], median=%.3f\n",
            min(obj$sig_core, na.rm=TRUE),
            max(obj$sig_core, na.rm=TRUE),
            median(obj$sig_core, na.rm=TRUE)))

# Filter to high quality cells only (drop tier "low" if any)
md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
cat(sprintf("\n  Quality tier counts:\n"))
print(table(md$data_quality_tier))

# Build plotting dataframe
plot_dt <- md[, .(cell_id, x_centroid, y_centroid,
                   density = density_knn_main_piecewise,
                   sig_94, sig_core,
                   cluster = factor(seurat_clusters),
                   tier = data_quality_tier)]
plot_dt <- plot_dt[!is.na(x_centroid) & !is.na(y_centroid) &
                    !is.na(sig_94) & !is.na(sig_core), ]
cat(sprintf("\n  Plot-ready cells: %d\n", nrow(plot_dt)))

# ============================================================
# 3. Compute spatial extent and aspect ratio
# ============================================================
x_range <- range(plot_dt$x_centroid)
y_range <- range(plot_dt$y_centroid)
x_span <- diff(x_range)
y_span <- diff(y_range)
asp_ratio <- y_span / x_span
cat(sprintf("\n  Spatial extent: X[%.0f, %.0f] (%.0f um), Y[%.0f, %.0f] (%.0f um)\n",
            x_range[1], x_range[2], x_span,
            y_range[1], y_range[2], y_span))
cat(sprintf("  Aspect ratio (Y/X): %.2f\n", asp_ratio))

# Choose figure dimensions to preserve aspect
fig_w <- 10
fig_h <- 10 * asp_ratio
fig_h <- max(min(fig_h, 14), 6)  # clamp to reasonable range
cat(sprintf("  Figure dimensions: %.1f × %.1f inches\n", fig_w, fig_h))

# ============================================================
# 4. Plot helpers
# ============================================================
spatial_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.position = "right",
    aspect.ratio = asp_ratio  # preserve true aspect
  )

# ---- Figure 1: density ----
cat("\n[Figure 1] Density spatial map...\n")
p1 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = density)) +
  geom_point(size = 0.15, alpha = 0.4) +
  scale_color_viridis_c(option = "plasma", name = "density") +
  coord_fixed() +
  labs(title = sprintf("[%s] Cell density (KNN)", SAMPLE_NAME),
       subtitle = sprintf("n=%s cells", format(nrow(plot_dt), big.mark=",")),
       x = "X (μm)", y = "Y (μm)") +
  spatial_theme
ggsave(file.path(FIG_DIR, "01_density_spatial.png"),
       p1, width = fig_w, height = fig_h, dpi = 150)
cat("  Saved\n")

# ---- Figure 2: sig_94 ----
cat("[Figure 2] sig_94 spatial map...\n")
# diverging scale centered on 0
sig94_lim <- max(abs(quantile(plot_dt$sig_94, c(0.01, 0.99), na.rm=TRUE)))
p2 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = sig_94)) +
  geom_point(size = 0.15, alpha = 0.4) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "sig_94",
                         limits = c(-sig94_lim, sig94_lim),
                         oob = scales::squish) +
  coord_fixed() +
  labs(title = sprintf("[%s] sig_94 score (zscore-weighted)", SAMPLE_NAME),
       subtitle = sprintf("range clipped to ±%.2f for color", sig94_lim),
       x = "X (μm)", y = "Y (μm)") +
  spatial_theme
ggsave(file.path(FIG_DIR, "02_sig94_spatial.png"),
       p2, width = fig_w, height = fig_h, dpi = 150)
cat("  Saved\n")

# ---- Figure 3: sig_core ----
cat("[Figure 3] sig_core spatial map...\n")
core_lim <- max(abs(quantile(plot_dt$sig_core, c(0.01, 0.99), na.rm=TRUE)))
p3 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = sig_core)) +
  geom_point(size = 0.15, alpha = 0.4) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "sig_core",
                         limits = c(-core_lim, core_lim),
                         oob = scales::squish) +
  coord_fixed() +
  labs(title = sprintf("[%s] sig_core score (8 genes)", SAMPLE_NAME),
       subtitle = sprintf("range clipped to ±%.2f", core_lim),
       x = "X (μm)", y = "Y (μm)") +
  spatial_theme
ggsave(file.path(FIG_DIR, "03_sigcore_spatial.png"),
       p3, width = fig_w, height = fig_h, dpi = 150)
cat("  Saved\n")

# ---- Figure 4: density vs sig_94 scatter ----
cat("[Figure 4] density vs sig_94 scatter...\n")
# subsample for scatter (60k → 10k for clarity)
n_sub <- min(10000, nrow(plot_dt))
sub <- plot_dt[sample(.N, n_sub)]
sp_cor <- cor(plot_dt$density, plot_dt$sig_94, method = "spearman", use = "complete.obs")

p4 <- ggplot(sub, aes(x = density, y = sig_94)) +
  geom_point(size = 0.3, alpha = 0.3, color = "#377eb8") +
  geom_smooth(method = "loess", color = "red", linewidth = 0.5, se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = sprintf("[%s] density vs sig_94", SAMPLE_NAME),
       subtitle = sprintf("Spearman ρ=%.3f (n=%s, subsample %s shown)",
                           sp_cor,
                           format(nrow(plot_dt), big.mark=","),
                           format(n_sub, big.mark=",")),
       x = "Cell density (KNN)", y = "sig_94 score") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(FIG_DIR, "04_density_vs_sig94_scatter.png"),
       p4, width = 8, height = 6, dpi = 150)
cat("  Saved\n")

# ---- Figure 5: top 25% vs bottom 25% sig_94 cells in space ----
cat("[Figure 5] high vs low sig_94 cells spatial comparison...\n")
q25 <- quantile(plot_dt$sig_94, 0.25, na.rm=TRUE)
q75 <- quantile(plot_dt$sig_94, 0.75, na.rm=TRUE)
plot_dt[, sig94_class := fifelse(sig_94 >= q75, "high",
                                  fifelse(sig_94 <= q25, "low", "mid"))]
plot_dt[, sig94_class := factor(sig94_class, levels = c("low", "mid", "high"))]

# Plot mid in light grey, high in red, low in blue
p5 <- ggplot() +
  geom_point(data = plot_dt[sig94_class == "mid"],
             aes(x = x_centroid, y = y_centroid),
             color = "grey90", size = 0.12, alpha = 0.3) +
  geom_point(data = plot_dt[sig94_class == "low"],
             aes(x = x_centroid, y = y_centroid),
             color = "#2166ac", size = 0.18, alpha = 0.5) +
  geom_point(data = plot_dt[sig94_class == "high"],
             aes(x = x_centroid, y = y_centroid),
             color = "#b2182b", size = 0.18, alpha = 0.5) +
  coord_fixed() +
  labs(title = sprintf("[%s] sig_94 high vs low cell distribution", SAMPLE_NAME),
       subtitle = sprintf("Red = top 25%% (n=%s), Blue = bottom 25%% (n=%s), Grey = middle 50%%",
                           format(sum(plot_dt$sig94_class == "high"), big.mark=","),
                           format(sum(plot_dt$sig94_class == "low"), big.mark=",")),
       x = "X (μm)", y = "Y (μm)") +
  spatial_theme
ggsave(file.path(FIG_DIR, "05_sig94_high_vs_low_spatial.png"),
       p5, width = fig_w, height = fig_h, dpi = 150)
cat("  Saved\n")

# ---- Figure 6: clusters ----
cat("[Figure 6] cluster spatial map...\n")
n_cl <- length(unique(plot_dt$cluster))
p6 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = cluster)) +
  geom_point(size = 0.15, alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                               ncol = if (n_cl > 15) 2 else 1)) +
  coord_fixed() +
  labs(title = sprintf("[%s] Seurat clusters (n=%d)", SAMPLE_NAME, n_cl),
       x = "X (μm)", y = "Y (μm)") +
  spatial_theme +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))
ggsave(file.path(FIG_DIR, "06_clusters_spatial.png"),
       p6, width = fig_w + 1, height = fig_h, dpi = 150)
cat("  Saved\n")

# ============================================================
# 5. Save a summary csv of cell-level data
# ============================================================
cat("\n[Save] cell-level data...\n")
fwrite(plot_dt, file.path(FIG_DIR, "cell_data.csv"))
cat(sprintf("  Saved cell_data.csv (%d cells)\n", nrow(plot_dt)))

# ============================================================
# 6. Summary
# ============================================================
summary_lines <- c(
  "================================================================",
  sprintf("R18a — Spatial Signature Visualization: %s", SAMPLE_NAME),
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "INPUT",
  sprintf("  Seurat: %s", SAMPLE_RDS),
  sprintf("  Cells: %d (after NA removal)", nrow(plot_dt)),
  sprintf("  Genes in panel: %d", nrow(obj)),
  "",
  "SIGNATURE COVERAGE",
  sprintf("  sig_94 panel match: %d/%d (%.1f%%)",
          length(match_94), length(sig_94),
          100 * length(match_94) / length(sig_94)),
  sprintf("    %d positive + %d negative direction genes",
          length(sig_94_pos), length(sig_94_neg)),
  sprintf("  sig_core panel match: %d/%d (%.0f%%)",
          length(match_core), length(SIG_CORE_GENES),
          100 * length(match_core) / length(SIG_CORE_GENES)),
  "",
  "SCORE RANGES",
  sprintf("  sig_94 [%.3f, %.3f] median=%.3f",
          min(plot_dt$sig_94), max(plot_dt$sig_94), median(plot_dt$sig_94)),
  sprintf("  sig_core [%.3f, %.3f] median=%.3f",
          min(plot_dt$sig_core), max(plot_dt$sig_core), median(plot_dt$sig_core)),
  "",
  "DENSITY-SIGNATURE COUPLING (sanity check)",
  sprintf("  Spearman ρ (density vs sig_94): %.3f", sp_cor),
  "  (Should be POSITIVE if signature was correctly derived from this sample)",
  "",
  "SPATIAL EXTENT",
  sprintf("  X span: %.0f μm", x_span),
  sprintf("  Y span: %.0f μm", y_span),
  sprintf("  Aspect ratio (Y/X): %.2f", asp_ratio),
  "",
  "OUTPUTS",
  sprintf("  Figures (6) and cell_data.csv in: %s", FIG_DIR),
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(FIG_DIR, "R18a_SUMMARY.txt"))
cat("\n=== SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R18a DONE — review figures, then approve R18b for 4 samples\n")
cat("================================================================\n")
