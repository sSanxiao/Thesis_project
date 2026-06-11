# ============================================================
# R18b: Spatial signature visualization — extend to all 4 MB
# ------------------------------------------------------------
# - Re-runs MB266 result through unified function (consistency)
# - Adds MB263, MB295, MB299
# - Cross-sample comparison figures
# - MB266 zoom-in to highest sig_94 region
# - Spearman ρ table with quality warnings (ρ < 0.2 flagged)
#
# Builds on R18a; same scoring logic and figure design.
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
R2_MB <- file.path(RESULTS_DIR, "R2_Results", "Medulloblastoma_Human")
MB_SAMPLES <- list(
  MB263 = file.path(R2_MB, "GSM8840046", "GSM8840046_seurat_R2.rds"),
  MB266 = file.path(R2_MB, "GSM8840047", "GSM8840047_seurat_R2.rds"),
  MB295 = file.path(R2_MB, "GSM8840048", "GSM8840048_seurat_R2.rds"),
  MB299 = file.path(R2_MB, "GSM8840049", "GSM8840049_seurat_R2.rds")
)

OUT_ROOT <- file.path(RESULTS_DIR, "R18_SpatialSignature")
COMPARE_DIR <- file.path(OUT_ROOT, "compare_4samples")
dir.create(COMPARE_DIR, showWarnings = FALSE, recursive = TRUE)

SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
SIG_CORE_GENES <- c("TUBB4A", "APOE", "CCR7", "EOMES", "ST18", "NES", "AQP4", "QKI")

# Load signature provenance once
sig_prov <- fread(SIG_PROV)
sig_94_full <- setdiff(sig_prov$gene, CONFLICT_GENES)

cat("================================================================\n")
cat("R18b: Spatial signature — 4 MB samples + comparisons\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# Helper: compute weighted score
# ============================================================
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

# ============================================================
# Helper: render all 6 figures for a sample, return summary
# ============================================================
process_sample <- function(sample_name, rds_path) {
  cat(sprintf("\n========== Processing %s ==========\n", sample_name))
  
  fig_dir <- file.path(OUT_ROOT, sample_name)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load
  cat("  Loading rds...\n")
  t0 <- Sys.time()
  obj <- readRDS(rds_path)
  cat(sprintf("    Loaded in %.1fs (%d cells)\n",
              as.numeric(Sys.time() - t0, units = "secs"), ncol(obj)))
  
  # Match signature genes
  panel_genes <- rownames(obj)
  match_94 <- intersect(sig_94_full, panel_genes)
  match_core <- intersect(SIG_CORE_GENES, panel_genes)
  cat(sprintf("    sig_94 matched: %d / %d (%.1f%%)\n",
              length(match_94), length(sig_94_full),
              100 * length(match_94) / length(sig_94_full)))
  cat(sprintf("    sig_core matched: %d / %d\n",
              length(match_core), length(SIG_CORE_GENES)))
  
  # Get expression
  expr_sct <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
  if (is.null(expr_sct) || all(dim(expr_sct) == 0)) {
    expr_sct <- GetAssayData(obj, assay = "SCT", layer = "data")
  }
  
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
  
  # Score
  cat("    Computing scores...\n")
  obj$sig_94 <- compute_score(expr_sct, sig_94_pos, sig_94_neg)
  obj$sig_core <- compute_score(expr_sct, sig_core_pos, sig_core_neg)
  
  # Plot dataframe
  md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  plot_dt <- md[, .(cell_id, x_centroid, y_centroid,
                     density = density_knn_main_piecewise,
                     sig_94, sig_core,
                     cluster = factor(seurat_clusters),
                     tier = data_quality_tier)]
  plot_dt <- plot_dt[!is.na(x_centroid) & !is.na(y_centroid) &
                      !is.na(sig_94) & !is.na(sig_core)]
  
  # Spearman
  sp_cor <- cor(plot_dt$density, plot_dt$sig_94, method = "spearman", use = "complete.obs")
  sp_cor_core <- cor(plot_dt$density, plot_dt$sig_core, method = "spearman", use = "complete.obs")
  cat(sprintf("    ρ (density vs sig_94): %.3f\n", sp_cor))
  cat(sprintf("    ρ (density vs sig_core): %.3f\n", sp_cor_core))
  
  # Quality warning
  warning_msg <- ""
  if (sp_cor < 0.2) {
    warning_msg <- sprintf("⚠️ WARNING: %s has weak density-sig_94 coupling (ρ=%.3f < 0.2). Signal may be weak in this sample (consistent with R9 AUC if low).",
                            sample_name, sp_cor)
    cat(sprintf("    %s\n", warning_msg))
  }
  
  # Spatial extent
  x_range <- range(plot_dt$x_centroid)
  y_range <- range(plot_dt$y_centroid)
  asp_ratio <- diff(y_range) / diff(x_range)
  fig_w <- 10
  fig_h <- max(min(10 * asp_ratio, 14), 6)
  
  spatial_theme <- theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.position = "right",
          aspect.ratio = asp_ratio)
  
  # ---- Figure 1: density ----
  cat("    [Fig 1/6] density spatial...\n")
  p1 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = density)) +
    geom_point(size = 0.15, alpha = 0.4) +
    scale_color_viridis_c(option = "plasma", name = "density") +
    coord_fixed() +
    labs(title = sprintf("[%s] Cell density (KNN)", sample_name),
         subtitle = sprintf("n=%s cells", format(nrow(plot_dt), big.mark=",")),
         x = "X (μm)", y = "Y (μm)") +
    spatial_theme
  ggsave(file.path(fig_dir, "01_density_spatial.png"),
         p1, width = fig_w, height = fig_h, dpi = 150)
  
  # ---- Figure 2: sig_94 ----
  cat("    [Fig 2/6] sig_94 spatial...\n")
  sig94_lim <- max(abs(quantile(plot_dt$sig_94, c(0.01, 0.99), na.rm=TRUE)))
  p2 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = sig_94)) +
    geom_point(size = 0.15, alpha = 0.4) +
    scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                           midpoint = 0, name = "sig_94",
                           limits = c(-sig94_lim, sig94_lim),
                           oob = scales::squish) +
    coord_fixed() +
    labs(title = sprintf("[%s] sig_94 score (zscore-weighted)", sample_name),
         subtitle = sprintf("color clipped to ±%.2f; ρ(density,sig_94)=%.3f",
                             sig94_lim, sp_cor),
         x = "X (μm)", y = "Y (μm)") +
    spatial_theme
  ggsave(file.path(fig_dir, "02_sig94_spatial.png"),
         p2, width = fig_w, height = fig_h, dpi = 150)
  
  # ---- Figure 3: sig_core ----
  cat("    [Fig 3/6] sig_core spatial...\n")
  core_lim <- max(abs(quantile(plot_dt$sig_core, c(0.01, 0.99), na.rm=TRUE)))
  p3 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = sig_core)) +
    geom_point(size = 0.15, alpha = 0.4) +
    scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                           midpoint = 0, name = "sig_core",
                           limits = c(-core_lim, core_lim),
                           oob = scales::squish) +
    coord_fixed() +
    labs(title = sprintf("[%s] sig_core score (8 genes)", sample_name),
         subtitle = sprintf("color clipped to ±%.2f", core_lim),
         x = "X (μm)", y = "Y (μm)") +
    spatial_theme
  ggsave(file.path(fig_dir, "03_sigcore_spatial.png"),
         p3, width = fig_w, height = fig_h, dpi = 150)
  
  # ---- Figure 4: density vs sig_94 scatter ----
  cat("    [Fig 4/6] density vs sig_94 scatter...\n")
  n_sub <- min(10000, nrow(plot_dt))
  sub <- plot_dt[sample(.N, n_sub)]
  p4 <- ggplot(sub, aes(x = density, y = sig_94)) +
    geom_point(size = 0.3, alpha = 0.3, color = "#377eb8") +
    geom_smooth(method = "loess", color = "red", linewidth = 0.5, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("[%s] density vs sig_94", sample_name),
         subtitle = sprintf("Spearman ρ=%.3f (n=%s, subsample %s shown)",
                             sp_cor,
                             format(nrow(plot_dt), big.mark=","),
                             format(n_sub, big.mark=",")),
         x = "Cell density (KNN)", y = "sig_94 score") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "04_density_vs_sig94_scatter.png"),
         p4, width = 8, height = 6, dpi = 150)
  
  # ---- Figure 5: high vs low ----
  cat("    [Fig 5/6] high vs low sig_94 spatial...\n")
  q25 <- quantile(plot_dt$sig_94, 0.25, na.rm=TRUE)
  q75 <- quantile(plot_dt$sig_94, 0.75, na.rm=TRUE)
  plot_dt[, sig94_class := fifelse(sig_94 >= q75, "high",
                                    fifelse(sig_94 <= q25, "low", "mid"))]
  plot_dt[, sig94_class := factor(sig94_class, levels = c("low", "mid", "high"))]
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
    labs(title = sprintf("[%s] sig_94 high vs low cell distribution", sample_name),
         subtitle = sprintf("Red = top 25%% (n=%s), Blue = bottom 25%% (n=%s), Grey = middle 50%%",
                             format(sum(plot_dt$sig94_class == "high"), big.mark=","),
                             format(sum(plot_dt$sig94_class == "low"), big.mark=",")),
         x = "X (μm)", y = "Y (μm)") +
    spatial_theme
  ggsave(file.path(fig_dir, "05_sig94_high_vs_low_spatial.png"),
         p5, width = fig_w, height = fig_h, dpi = 150)
  
  # ---- Figure 6: clusters ----
  cat("    [Fig 6/6] clusters spatial...\n")
  n_cl <- length(unique(plot_dt$cluster))
  p6 <- ggplot(plot_dt, aes(x = x_centroid, y = y_centroid, color = cluster)) +
    geom_point(size = 0.15, alpha = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                 ncol = if (n_cl > 15) 2 else 1)) +
    coord_fixed() +
    labs(title = sprintf("[%s] Seurat clusters (n=%d)", sample_name, n_cl),
         x = "X (μm)", y = "Y (μm)") +
    spatial_theme +
    theme(legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))
  ggsave(file.path(fig_dir, "06_clusters_spatial.png"),
         p6, width = fig_w + 1, height = fig_h, dpi = 150)
  
  # Save cell data
  fwrite(plot_dt, file.path(fig_dir, "cell_data.csv"))
  
  # Per-sample summary
  summary_lines <- c(
    "================================================================",
    sprintf("R18b — Spatial Signature: %s", sample_name),
    sprintf("Timestamp: %s", Sys.time()),
    "================================================================",
    "",
    sprintf("INPUT: %s", rds_path),
    sprintf("  Cells: %d (after NA removal)", nrow(plot_dt)),
    sprintf("  Genes: %d", nrow(obj)),
    "",
    "SIGNATURE COVERAGE",
    sprintf("  sig_94: %d/%d (%.1f%%) [%d pos + %d neg]",
            length(match_94), length(sig_94_full),
            100 * length(match_94) / length(sig_94_full),
            length(sig_94_pos), length(sig_94_neg)),
    sprintf("  sig_core: %d/%d", length(match_core), length(SIG_CORE_GENES)),
    "",
    "SCORE RANGES",
    sprintf("  sig_94: [%.3f, %.3f] median=%.3f",
            min(plot_dt$sig_94), max(plot_dt$sig_94), median(plot_dt$sig_94)),
    sprintf("  sig_core: [%.3f, %.3f] median=%.3f",
            min(plot_dt$sig_core), max(plot_dt$sig_core), median(plot_dt$sig_core)),
    "",
    "DENSITY-SIGNATURE COUPLING",
    sprintf("  ρ(density, sig_94) = %.3f", sp_cor),
    sprintf("  ρ(density, sig_core) = %.3f", sp_cor_core)
  )
  if (nchar(warning_msg) > 0) {
    summary_lines <- c(summary_lines, "", warning_msg)
  }
  writeLines(summary_lines, file.path(fig_dir, sprintf("R18b_%s_SUMMARY.txt", sample_name)))
  
  # Free memory
  rm(obj, expr_sct, md)
  gc(verbose = FALSE)
  
  # Return summary record
  list(
    sample = sample_name,
    n_cells = nrow(plot_dt),
    sig_94_match = length(match_94),
    sig_94_total = length(sig_94_full),
    sig_94_pos = length(sig_94_pos),
    sig_94_neg = length(sig_94_neg),
    sig_94_min = min(plot_dt$sig_94),
    sig_94_max = max(plot_dt$sig_94),
    sig_94_median = median(plot_dt$sig_94),
    sp_cor_sig94 = sp_cor,
    sp_cor_sigcore = sp_cor_core,
    warning = warning_msg,
    plot_dt = plot_dt  # keep for cross-sample plots
  )
}

# ============================================================
# Run all 4 samples
# ============================================================
results <- list()
for (sn in names(MB_SAMPLES)) {
  results[[sn]] <- process_sample(sn, MB_SAMPLES[[sn]])
}

# ============================================================
# Cross-sample comparison: ρ table
# ============================================================
cat("\n\n========== Cross-sample comparison ==========\n")

cor_table <- data.table(
  sample = sapply(results, function(r) r$sample),
  n_cells = sapply(results, function(r) r$n_cells),
  sig_94_match_pct = sapply(results, function(r) round(100 * r$sig_94_match / r$sig_94_total, 1)),
  sig_94_min = sapply(results, function(r) round(r$sig_94_min, 3)),
  sig_94_max = sapply(results, function(r) round(r$sig_94_max, 3)),
  sig_94_median = sapply(results, function(r) round(r$sig_94_median, 3)),
  rho_sig94 = sapply(results, function(r) round(r$sp_cor_sig94, 3)),
  rho_sigcore = sapply(results, function(r) round(r$sp_cor_sigcore, 3)),
  warning = sapply(results, function(r) if (nchar(r$warning) > 0) "⚠" else "")
)

cat("\nCross-sample ρ summary:\n")
print(cor_table)

fwrite(cor_table, file.path(COMPARE_DIR, "cross_sample_rho_table.csv"))

# ============================================================
# Cross-sample figure: 4-panel sig_94 spatial
# ============================================================
cat("\n[Cross-sample fig 1] 4-panel sig_94 spatial...\n")

panels_sig94 <- list()
for (sn in names(results)) {
  r <- results[[sn]]
  pdt <- r$plot_dt
  asp <- diff(range(pdt$y_centroid)) / diff(range(pdt$x_centroid))
  sig94_lim <- max(abs(quantile(pdt$sig_94, c(0.01, 0.99), na.rm=TRUE)))
  
  p <- ggplot(pdt, aes(x = x_centroid, y = y_centroid, color = sig_94)) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                           midpoint = 0, name = "sig_94",
                           limits = c(-sig94_lim, sig94_lim),
                           oob = scales::squish) +
    coord_fixed() +
    labs(title = sprintf("%s   ρ=%.3f", sn, r$sp_cor_sig94),
         subtitle = sprintf("n=%s cells", format(r$n_cells, big.mark=",")),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 6),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          aspect.ratio = asp)
  panels_sig94[[sn]] <- p
}

p_4panel_sig94 <- wrap_plots(panels_sig94, ncol = 2)
ggsave(file.path(COMPARE_DIR, "01_4samples_sig94.png"),
       p_4panel_sig94, width = 16, height = 14, dpi = 150)

# Same for density
cat("[Cross-sample fig 2] 4-panel density spatial...\n")
panels_density <- list()
for (sn in names(results)) {
  r <- results[[sn]]
  pdt <- r$plot_dt
  asp <- diff(range(pdt$y_centroid)) / diff(range(pdt$x_centroid))
  
  p <- ggplot(pdt, aes(x = x_centroid, y = y_centroid, color = density)) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_viridis_c(option = "plasma", name = "density") +
    coord_fixed() +
    labs(title = sn, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 6),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          aspect.ratio = asp)
  panels_density[[sn]] <- p
}
p_4panel_density <- wrap_plots(panels_density, ncol = 2)
ggsave(file.path(COMPARE_DIR, "02_4samples_density.png"),
       p_4panel_density, width = 16, height = 14, dpi = 150)

# Same for high-vs-low
cat("[Cross-sample fig 3] 4-panel high vs low...\n")
panels_hl <- list()
for (sn in names(results)) {
  r <- results[[sn]]
  pdt <- r$plot_dt
  asp <- diff(range(pdt$y_centroid)) / diff(range(pdt$x_centroid))
  
  p <- ggplot() +
    geom_point(data = pdt[sig94_class == "mid"],
               aes(x = x_centroid, y = y_centroid),
               color = "grey90", size = 0.05, alpha = 0.25) +
    geom_point(data = pdt[sig94_class == "low"],
               aes(x = x_centroid, y = y_centroid),
               color = "#2166ac", size = 0.08, alpha = 0.45) +
    geom_point(data = pdt[sig94_class == "high"],
               aes(x = x_centroid, y = y_centroid),
               color = "#b2182b", size = 0.08, alpha = 0.45) +
    coord_fixed() +
    labs(title = sn, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 6),
          panel.grid = element_blank(),
          aspect.ratio = asp)
  panels_hl[[sn]] <- p
}
p_4panel_hl <- wrap_plots(panels_hl, ncol = 2)
ggsave(file.path(COMPARE_DIR, "03_4samples_high_vs_low.png"),
       p_4panel_hl, width = 16, height = 14, dpi = 150)

# ============================================================
# MB266 zoom-in
# ============================================================
cat("\n[Zoom-in] MB266 highest sig_94 region...\n")

mb266_dt <- results$MB266$plot_dt

# Find center of top 1000 sig_94 cells
top_cells <- mb266_dt[order(-sig_94)][1:1000]
cx <- median(top_cells$x_centroid)
cy <- median(top_cells$y_centroid)
zoom_w <- 1000  # 2000 μm window total (±1000 from center)
cat(sprintf("  Zoom center: (%.0f, %.0f), window 2000×2000 μm\n", cx, cy))

zoom_dt <- mb266_dt[abs(x_centroid - cx) < zoom_w &
                     abs(y_centroid - cy) < zoom_w]
cat(sprintf("  Cells in zoom: %d\n", nrow(zoom_dt)))

# Compute classes specific to zoom
zoom_q25 <- quantile(zoom_dt$sig_94, 0.25)
zoom_q75 <- quantile(zoom_dt$sig_94, 0.75)
zoom_dt[, sig94_class := fifelse(sig_94 >= zoom_q75, "high",
                                  fifelse(sig_94 <= zoom_q25, "low", "mid"))]

# Zoom panels
zoom_lim <- max(abs(quantile(zoom_dt$sig_94, c(0.01, 0.99))))

p_zoom_density <- ggplot(zoom_dt, aes(x = x_centroid, y = y_centroid, color = density)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", name = "density") +
  coord_fixed() +
  labs(title = "[MB266 zoom] Cell density",
       subtitle = sprintf("Center: (%.0f, %.0f), 2000×2000 μm, n=%d cells",
                           cx, cy, nrow(zoom_dt)),
       x = "X (μm)", y = "Y (μm)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

p_zoom_sig94 <- ggplot(zoom_dt, aes(x = x_centroid, y = y_centroid, color = sig_94)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "sig_94",
                         limits = c(-zoom_lim, zoom_lim),
                         oob = scales::squish) +
  coord_fixed() +
  labs(title = "[MB266 zoom] sig_94",
       x = "X (μm)", y = "Y (μm)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

p_zoom_hl <- ggplot() +
  geom_point(data = zoom_dt[sig94_class == "mid"],
             aes(x = x_centroid, y = y_centroid),
             color = "grey90", size = 0.5, alpha = 0.5) +
  geom_point(data = zoom_dt[sig94_class == "low"],
             aes(x = x_centroid, y = y_centroid),
             color = "#2166ac", size = 0.7, alpha = 0.7) +
  geom_point(data = zoom_dt[sig94_class == "high"],
             aes(x = x_centroid, y = y_centroid),
             color = "#b2182b", size = 0.7, alpha = 0.7) +
  coord_fixed() +
  labs(title = "[MB266 zoom] high vs low sig_94",
       x = "X (μm)", y = "Y (μm)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

p_zoom <- (p_zoom_density | p_zoom_sig94 | p_zoom_hl) +
  plot_annotation(title = "MB266 zoom-in to highest sig_94 region",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(COMPARE_DIR, "04_MB266_zoom.png"),
       p_zoom, width = 18, height = 7, dpi = 150)

# ============================================================
# Combined SUMMARY
# ============================================================
cat("\n[Final] writing combined SUMMARY...\n")

summary_lines <- c(
  "================================================================",
  "R18b — Spatial Signature Visualization (4 MB samples)",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "DESIGN",
  "  4 MBEN samples from GSE291688:",
  "    MB263 (GSM8840046)",
  "    MB266 (GSM8840047) - strongest signal sample",
  "    MB295 (GSM8840048)",
  "    MB299 (GSM8840049)",
  "",
  "  Per sample: 6 figures (density, sig_94, sig_core, scatter, high-low, clusters)",
  "  Cross-sample: 4-panel sig_94, 4-panel density, 4-panel high-low, MB266 zoom",
  "",
  "================================================================",
  "PER-SAMPLE RESULTS",
  "================================================================",
  ""
)

for (sn in names(results)) {
  r <- results[[sn]]
  summary_lines <- c(summary_lines,
    sprintf("[%s]", sn),
    sprintf("  n_cells: %d", r$n_cells),
    sprintf("  sig_94 matched: %d/%d (%.1f%%)",
            r$sig_94_match, r$sig_94_total,
            100 * r$sig_94_match / r$sig_94_total),
    sprintf("    %d positive + %d negative direction genes",
            r$sig_94_pos, r$sig_94_neg),
    sprintf("  sig_94 score: [%.3f, %.3f] median=%.3f",
            r$sig_94_min, r$sig_94_max, r$sig_94_median),
    sprintf("  ρ(density, sig_94) = %.3f", r$sp_cor_sig94),
    sprintf("  ρ(density, sig_core) = %.3f", r$sp_cor_sigcore),
    if (nchar(r$warning) > 0) sprintf("  %s", r$warning) else "  (signal strength acceptable)",
    "")
}

summary_lines <- c(summary_lines,
  "================================================================",
  "CROSS-SAMPLE ρ TABLE",
  "================================================================",
  "",
  "  Sample    n_cells   sig_94_match   ρ_sig94   ρ_sigcore   Status",
  "  ----------------------------------------------------------------")
for (i in 1:nrow(cor_table)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-7s   %7s    %6.1f%%       %+.3f    %+.3f       %s",
            cor_table$sample[i],
            format(cor_table$n_cells[i], big.mark=","),
            cor_table$sig_94_match_pct[i],
            cor_table$rho_sig94[i],
            cor_table$rho_sigcore[i],
            cor_table$warning[i]))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "INTERPRETATION",
  "================================================================",
  "",
  "Expected pattern:",
  "  - All 4 samples ρ > 0 (signature derived from density correlation)",
  "  - MB266 strongest (ρ ~0.6, R9 AUC = 0.855)",
  "  - Other samples lower (R9 AUC 0.60-0.65, expect ρ 0.2-0.4)",
  "  - Any sample with ρ < 0.2 has weak signal (flag in Discussion)",
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "",
  "  Per-sample dirs (each with 6 figures + cell_data.csv + per-sample SUMMARY):",
  sprintf("    %s/MB263/", OUT_ROOT),
  sprintf("    %s/MB266/", OUT_ROOT),
  sprintf("    %s/MB295/", OUT_ROOT),
  sprintf("    %s/MB299/", OUT_ROOT),
  "",
  sprintf("  Cross-sample comparison: %s/", COMPARE_DIR),
  "    01_4samples_sig94.png   (★ key figure)",
  "    02_4samples_density.png",
  "    03_4samples_high_vs_low.png",
  "    04_MB266_zoom.png",
  "    cross_sample_rho_table.csv",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(COMPARE_DIR, "R18b_SUMMARY.txt"))
writeLines(summary_lines, file.path(OUT_ROOT, "R18b_SUMMARY.txt"))

cat("\n=== FINAL SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R18b DONE\n")
cat("================================================================\n")
