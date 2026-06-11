# ============================================================
# R17c: Neuronal lineage pseudotime on Aldinger 2021
# ------------------------------------------------------------
# Main analysis: RL lineage 4 clusters (02-RL, 03-GCP, 04-GN, 05-eCN/UBC)
#                = 12243 cells
# Extended:      add 06-iCN → 16384 cells
#
# For each subset:
#   1) Re-embed with subset-specific PCA + UMAP
#   2) Start point = centroid of 9 PCW 02-RL cells
#   3) slingshot to infer trajectory + pseudotime
#   4) For each signature × method × lineage, smooth score vs pseudotime
#
# Outputs: multiple PNG figures + CSV trajectories + SUMMARY
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(slingshot)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

cat("================================================================\n")
cat("R17c: Neuronal lineage pseudotime analysis\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
RES_DIR <- file.path(RESULTS_DIR, "R17_Aldinger2021")
FIG_DIR <- file.path(RES_DIR, "figures_c")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

INPUT_RDS <- file.path(RES_DIR, "aldinger2021_scored.rds")

# ============================================================
# 1. Load + prepare subsets
# ============================================================
cat("[1] Loading scored Seurat object...\n")
obj_full <- readRDS(INPUT_RDS)
cat(sprintf("  Full: %d genes x %d cells\n", nrow(obj_full), ncol(obj_full)))

# MAIN: RL lineage
RL_CLUSTERS <- c("02-RL", "03-GCP", "04-GN", "05-eCN/UBC")
EXT_CLUSTERS <- c(RL_CLUSTERS, "06-iCN")

cat(sprintf("\n  Main RL lineage (4 clusters):\n"))
for (ct in RL_CLUSTERS) {
  cat(sprintf("    %s: %d cells\n", ct, sum(obj_full$cell_type == ct)))
}

cat(sprintf("\n  Extended (5 clusters, +06-iCN):\n"))
for (ct in EXT_CLUSTERS) {
  cat(sprintf("    %s: %d cells\n", ct, sum(obj_full$cell_type == ct)))
}

# ============================================================
# 2. Function: full pipeline for a subset (PCA → UMAP → Slingshot)
# ============================================================
run_trajectory <- function(obj_full, selected_clusters, tag,
                            start_cluster = "02-RL",
                            start_age_pcw = 9) {
  
  cat(sprintf("\n\n========== [%s] subset analysis ==========\n", tag))
  
  # Subset
  keep <- obj_full$cell_type %in% selected_clusters
  obj <- subset(obj_full, cells = colnames(obj_full)[keep])
  cat(sprintf("  Subset: %d cells, %d cell types\n",
              ncol(obj), length(unique(obj$cell_type))))
  
  # Re-identify variable features on subset
  # Data is SCT residuals in scale.data, all 2000 genes are already variable
  # so we use all of them as variable features
  VariableFeatures(obj) <- rownames(obj)
  
  cat("  Running PCA on subset...\n")
  # scale.data is already z-scored; RunPCA directly on it
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
  
  cat("  Running UMAP on subset PCA...\n")
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE,
                  reduction.name = "umap_subset", reduction.key = "UMAPs_")
  
  # -------- Identify start cell --------
  cat(sprintf("  Finding start cell: %s at %d PCW...\n",
              start_cluster, start_age_pcw))
  start_candidates <- which(obj$cell_type == start_cluster &
                             obj$age_pcw == start_age_pcw)
  
  if (length(start_candidates) == 0) {
    # Fallback: just use cluster, any age
    cat(sprintf("    [!] No cells at %d PCW in %s, using youngest age available\n",
                start_age_pcw, start_cluster))
    ages_in_start <- obj$age_pcw[obj$cell_type == start_cluster]
    best_age <- min(ages_in_start, na.rm = TRUE)
    start_candidates <- which(obj$cell_type == start_cluster &
                               obj$age_pcw == best_age)
    cat(sprintf("    Using %d PCW (%d cells)\n", best_age, length(start_candidates)))
  } else {
    cat(sprintf("    Found %d candidate cells\n", length(start_candidates)))
  }
  
  # Use the centroid in PCA space
  pca_emb <- Embeddings(obj, reduction = "pca")
  centroid <- colMeans(pca_emb[start_candidates, , drop = FALSE])
  # Cell closest to centroid
  dists <- sqrt(rowSums(sweep(pca_emb[start_candidates, , drop = FALSE], 2, centroid)^2))
  start_cell_idx <- start_candidates[which.min(dists)]
  start_cell_id <- colnames(obj)[start_cell_idx]
  cat(sprintf("    Start cell ID: %s\n", start_cell_id))
  
  # -------- Slingshot --------
  cat("  Running slingshot...\n")
  t0 <- Sys.time()
  
  # Build SingleCellExperiment
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  reducedDims(sce) <- list(UMAP = Embeddings(obj, "umap_subset"),
                            PCA = Embeddings(obj, "pca"))
  
  # Slingshot uses cluster labels to anchor lineages;
  # start cluster = start_cluster (specified)
  sce <- slingshot(sce,
                    clusterLabels = obj$cell_type,
                    reducedDim = "UMAP",
                    start.clus = start_cluster,
                    approx_points = 150)
  
  t1 <- Sys.time()
  cat(sprintf("    Slingshot done in %.1f sec\n",
              as.numeric(t1 - t0, units = "secs")))
  
  # Extract pseudotime & curves
  sds <- SlingshotDataSet(sce)
  lineages <- slingLineages(sds)
  cat(sprintf("    Identified %d lineage(s):\n", length(lineages)))
  for (i in seq_along(lineages)) {
    cat(sprintf("      Lineage %d: %s\n", i, paste(lineages[[i]], collapse = " → ")))
  }
  
  pseudo_mat <- slingPseudotime(sds)  # cell × lineage
  weights_mat <- slingCurveWeights(sds)
  
  # Join pseudotime back to Seurat obj
  for (i in 1:ncol(pseudo_mat)) {
    obj[[sprintf("slingPseudotime_%d", i)]] <- pseudo_mat[, i]
    obj[[sprintf("slingWeight_%d", i)]] <- weights_mat[, i]
  }
  
  # ============================================================
  # Figures for this subset
  # ============================================================
  cat("  Figures...\n")
  
  umap_df <- data.table(
    UMAP1 = Embeddings(obj, "umap_subset")[, 1],
    UMAP2 = Embeddings(obj, "umap_subset")[, 2],
    cell_type = obj$cell_type,
    age_pcw = obj$age_pcw,
    sig_94_zscore = obj$sig_94_zscore,
    sig_94_AMS = obj$sig_94_AMS,
    sig_core_zscore = obj$sig_core_zscore,
    sig_core_AMS = obj$sig_core_AMS
  )
  for (i in 1:ncol(pseudo_mat)) {
    umap_df[[sprintf("pt_%d", i)]] <- pseudo_mat[, i]
  }
  
  # --- Figure: new UMAP by cell type
  p_ct <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
    geom_point(size = 0.5, alpha = 0.7) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title = sprintf("[%s] subset UMAP (cell types)", tag),
         x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, sprintf("%s_01_UMAP_celltype.png", tag)),
         p_ct, width = 8, height = 6, dpi = 150)
  
  # --- Figure: new UMAP by age
  p_age <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = age_pcw)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", name = "PCW") +
    labs(title = sprintf("[%s] subset UMAP (age)", tag)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, sprintf("%s_02_UMAP_age.png", tag)),
         p_age, width = 8, height = 6, dpi = 150)
  
  # --- Figure: new UMAP with pseudotime (use first lineage)
  p_pt_list <- list()
  for (i in 1:ncol(pseudo_mat)) {
    pt_col <- sprintf("pt_%d", i)
    p <- ggplot(umap_df[!is.na(get(pt_col))], aes(x = UMAP1, y = UMAP2,
                                                    color = .data[[pt_col]])) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_viridis_c(option = "viridis", name = "pseudotime",
                             na.value = "grey80") +
      labs(title = sprintf("[%s] Lineage %d pseudotime\n(%s)",
                            tag, i, paste(lineages[[i]], collapse = " → "))) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    p_pt_list[[i]] <- p
  }
  
  p_pt_combined <- wrap_plots(p_pt_list, ncol = min(length(p_pt_list), 3))
  ggsave(file.path(FIG_DIR, sprintf("%s_03_UMAP_pseudotime.png", tag)),
         p_pt_combined, width = 7 * min(length(p_pt_list), 3),
         height = 6, dpi = 150)
  
  # --- Figure: new UMAP with sig_94 overlay
  p_sig94 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sig_94_zscore)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b",
                          midpoint = 0, name = "sig_94") +
    labs(title = sprintf("[%s] sig_94 (zscore) on subset UMAP", tag)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, sprintf("%s_04_UMAP_sig94.png", tag)),
         p_sig94, width = 8, height = 6, dpi = 150)
  
  # --- Figure: score vs pseudotime, loess smooth, per lineage
  # Prepare long-format: cell × lineage × score
  plot_score_vs_pt <- function(score_col, title_suffix) {
    plot_dt_list <- list()
    for (i in 1:ncol(pseudo_mat)) {
      pt_col <- sprintf("pt_%d", i)
      lineage_name <- paste(lineages[[i]], collapse = "→")
      sub_dt <- umap_df[!is.na(get(pt_col)), .(pt = get(pt_col),
                                                 score = get(score_col),
                                                 cell_type = cell_type,
                                                 lineage = lineage_name)]
      plot_dt_list[[i]] <- sub_dt
    }
    plot_dt <- rbindlist(plot_dt_list)
    
    # Compute Spearman cor per lineage
    cors <- plot_dt[, .(spearman = cor(pt, score, method = "spearman")),
                    by = lineage]
    
    # Subtitle with correlations
    subtitle_txt <- paste(
      sprintf("%s: ρ=%.3f", cors$lineage, cors$spearman),
      collapse = "   |   "
    )
    
    p <- ggplot(plot_dt, aes(x = pt, y = score, color = lineage)) +
      geom_point(size = 0.3, alpha = 0.15) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "pseudotime", y = score_col,
           title = sprintf("[%s] %s vs pseudotime%s",
                            tag, score_col, title_suffix),
           subtitle = subtitle_txt) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")
    p
  }
  
  p_sig94_pt <- plot_score_vs_pt("sig_94_zscore", "")
  p_sig94_ams <- plot_score_vs_pt("sig_94_AMS", " (AMS)")
  p_core_pt <- plot_score_vs_pt("sig_core_zscore", "")
  p_core_ams <- plot_score_vs_pt("sig_core_AMS", " (AMS)")
  
  p_score_combined <- (p_sig94_pt + p_sig94_ams) / (p_core_pt + p_core_ams)
  ggsave(file.path(FIG_DIR, sprintf("%s_05_score_vs_pseudotime_4panel.png", tag)),
         p_score_combined, width = 16, height = 12, dpi = 150)
  
  # --- Figure: cell type distribution along pseudotime
  # (sanity check: RL should be early, GN/UBC should be late)
  p_ct_list <- list()
  for (i in 1:ncol(pseudo_mat)) {
    pt_col <- sprintf("pt_%d", i)
    sub_dt <- umap_df[!is.na(get(pt_col)), .(pt = get(pt_col),
                                              cell_type = cell_type)]
    p <- ggplot(sub_dt, aes(x = pt, y = cell_type, fill = cell_type)) +
      ggridges::geom_density_ridges(alpha = 0.7) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "pseudotime", y = NULL,
           title = sprintf("[%s] Lineage %d cell distribution", tag, i)) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
    p_ct_list[[i]] <- p
  }
  p_ct_ridges <- wrap_plots(p_ct_list, ncol = min(length(p_ct_list), 2))
  ggsave(file.path(FIG_DIR, sprintf("%s_06_celltype_along_pseudotime.png", tag)),
         p_ct_ridges, width = 7 * min(length(p_ct_list), 2),
         height = 5, dpi = 150)
  
  # --- CSV: save pseudotime + scores per cell
  cells_dt <- data.table(
    cell_id = colnames(obj),
    cell_type = obj$cell_type,
    age_pcw = obj$age_pcw,
    sig_94_zscore = obj$sig_94_zscore,
    sig_94_AMS = obj$sig_94_AMS,
    sig_core_zscore = obj$sig_core_zscore,
    sig_core_AMS = obj$sig_core_AMS
  )
  for (i in 1:ncol(pseudo_mat)) {
    cells_dt[[sprintf("pseudotime_lineage_%d", i)]] <- pseudo_mat[, i]
    cells_dt[[sprintf("weight_lineage_%d", i)]] <- weights_mat[, i]
  }
  fwrite(cells_dt, file.path(RES_DIR, sprintf("%s_trajectory_data.csv", tag)))
  
  # --- Compute summary stats: Spearman cor per signature × lineage
  stats_list <- list()
  for (score_col in c("sig_94_zscore", "sig_94_AMS", "sig_core_zscore", "sig_core_AMS")) {
    for (i in 1:ncol(pseudo_mat)) {
      pt_vec <- pseudo_mat[, i]
      score_vec <- obj[[score_col]][, 1]
      valid <- !is.na(pt_vec)
      n_valid <- sum(valid)
      if (n_valid > 10) {
        rho <- cor(pt_vec[valid], score_vec[valid], method = "spearman")
        pval <- cor.test(pt_vec[valid], score_vec[valid], method = "spearman",
                          exact = FALSE)$p.value
      } else {
        rho <- NA; pval <- NA
      }
      stats_list[[length(stats_list) + 1]] <- data.table(
        subset = tag,
        signature = score_col,
        lineage = i,
        lineage_name = paste(lineages[[i]], collapse = "→"),
        n_cells = n_valid,
        spearman_rho = rho,
        p_value = pval
      )
    }
  }
  stats_dt <- rbindlist(stats_list)
  
  list(obj = obj, sds = sds, lineages = lineages,
       pseudotime = pseudo_mat, weights = weights_mat,
       cells_dt = cells_dt, stats_dt = stats_dt)
}

# ============================================================
# 3. Check ggridges
# ============================================================
if (!requireNamespace("ggridges", quietly = TRUE)) {
  cat("  [!] ggridges not installed; cell distribution ridges plot will skip\n")
  # fallback: use boxplot instead
  # let's just install it silently if available from CRAN
  tryCatch(install.packages("ggridges", repos = "https://cloud.r-project.org"),
           error = function(e) cat("  ggridges install failed\n"))
}

# ============================================================
# 4. Run main (RL lineage 4 clusters)
# ============================================================
main_res <- run_trajectory(obj_full, RL_CLUSTERS, tag = "main_RLlineage")

# ============================================================
# 5. Run extended (5 clusters, + iCN)
# ============================================================
ext_res <- run_trajectory(obj_full, EXT_CLUSTERS, tag = "ext_withIcn")

# ============================================================
# 6. Global UMAP with subset highlight (reference)
# ============================================================
cat("\n\n[Global UMAP highlight] ...\n")

global_umap_dt <- data.table(
  UMAP1 = Embeddings(obj_full, "umap")[, 1],
  UMAP2 = Embeddings(obj_full, "umap")[, 2],
  cell_type = obj_full$cell_type,
  is_RLlineage = obj_full$cell_type %in% RL_CLUSTERS,
  is_withIcn = obj_full$cell_type %in% EXT_CLUSTERS
)

p_g1 <- ggplot(global_umap_dt, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = global_umap_dt[!is_RLlineage], color = "grey85",
             size = 0.2, alpha = 0.5) +
  geom_point(data = global_umap_dt[is_RLlineage],
             aes(color = cell_type), size = 0.3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Main analysis subset (RL lineage, 12243 cells) on global UMAP",
       subtitle = "Grey = other cells in Aldinger 2021 dataset") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_g2 <- ggplot(global_umap_dt, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = global_umap_dt[!is_withIcn], color = "grey85",
             size = 0.2, alpha = 0.5) +
  geom_point(data = global_umap_dt[is_withIcn],
             aes(color = cell_type), size = 0.3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Extended subset (+ 06-iCN, 16384 cells) on global UMAP",
       subtitle = "Grey = other cells") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_global <- p_g1 / p_g2
ggsave(file.path(FIG_DIR, "00_global_UMAP_subset_highlight.png"),
       p_global, width = 10, height = 14, dpi = 150)
cat("  Saved: 00_global_UMAP_subset_highlight.png\n")

# ============================================================
# 7. Save combined stats + SUMMARY
# ============================================================
all_stats <- rbindlist(list(main_res$stats_dt, ext_res$stats_dt))
fwrite(all_stats, file.path(RES_DIR, "R17c_trajectory_stats.csv"))
cat("\n[Stats] R17c_trajectory_stats.csv saved\n")
print(all_stats)

# Save subset Seurat objects
saveRDS(main_res$obj, file.path(RES_DIR, "main_RLlineage_scored.rds"))
saveRDS(ext_res$obj, file.path(RES_DIR, "ext_withIcn_scored.rds"))

# SUMMARY
summary_lines <- c(
  "================================================================",
  "R17c — Neuronal lineage pseudotime analysis",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "MAIN ANALYSIS: RL lineage (4 clusters)",
  "  Clusters: 02-RL, 03-GCP, 04-GN, 05-eCN/UBC",
  sprintf("  Total cells: %d", ncol(main_res$obj)),
  sprintf("  Lineages detected: %d", length(main_res$lineages))
)
for (i in seq_along(main_res$lineages)) {
  summary_lines <- c(summary_lines,
    sprintf("    Lineage %d: %s", i,
            paste(main_res$lineages[[i]], collapse = " → ")))
}

summary_lines <- c(summary_lines, "",
  "EXTENDED ANALYSIS: +06-iCN (5 clusters)",
  "  Clusters: 02-RL, 03-GCP, 04-GN, 05-eCN/UBC, 06-iCN",
  sprintf("  Total cells: %d", ncol(ext_res$obj)),
  sprintf("  Lineages detected: %d", length(ext_res$lineages)))
for (i in seq_along(ext_res$lineages)) {
  summary_lines <- c(summary_lines,
    sprintf("    Lineage %d: %s", i,
            paste(ext_res$lineages[[i]], collapse = " → ")))
}

summary_lines <- c(summary_lines, "",
  "SPEARMAN ρ: signature score vs pseudotime (main)")
for (i in 1:nrow(main_res$stats_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-18s  Lineage %d (%s)  n=%d  ρ=%+.3f  p=%.2e",
            main_res$stats_dt$signature[i],
            main_res$stats_dt$lineage[i],
            main_res$stats_dt$lineage_name[i],
            main_res$stats_dt$n_cells[i],
            main_res$stats_dt$spearman_rho[i],
            main_res$stats_dt$p_value[i]))
}

summary_lines <- c(summary_lines, "",
  "SPEARMAN ρ: signature score vs pseudotime (extended)")
for (i in 1:nrow(ext_res$stats_dt)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-18s  Lineage %d (%s)  n=%d  ρ=%+.3f  p=%.2e",
            ext_res$stats_dt$signature[i],
            ext_res$stats_dt$lineage[i],
            ext_res$stats_dt$lineage_name[i],
            ext_res$stats_dt$n_cells[i],
            ext_res$stats_dt$spearman_rho[i],
            ext_res$stats_dt$p_value[i]))
}

summary_lines <- c(summary_lines, "",
  "OUTPUTS",
  sprintf("  Figures: %s/figures_c/", RES_DIR),
  sprintf("  Trajectory data (main): main_RLlineage_trajectory_data.csv"),
  sprintf("  Trajectory data (ext):  ext_withIcn_trajectory_data.csv"),
  sprintf("  Seurat objects: main_RLlineage_scored.rds, ext_withIcn_scored.rds"),
  sprintf("  Combined stats: R17c_trajectory_stats.csv"),
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(RES_DIR, "R17c_SUMMARY.txt"))
cat(sprintf("\nSUMMARY: %s/R17c_SUMMARY.txt\n", RES_DIR))

cat("\n================================================================\n")
cat("R17c DONE\n")
cat("================================================================\n")
