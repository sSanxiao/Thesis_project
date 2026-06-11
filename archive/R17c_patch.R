# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R17c patch v2: complete stats + SUMMARY from CSV outputs
# ------------------------------------------------------------
# Note: original R17c crashed before saving subset rds files,
# but trajectory_data.csv files DO exist with all needed info.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cat("================================================================\n")
cat("R17c patch v2: complete from CSVs\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

RES_DIR <- "./results/R17_Aldinger2021"
FIG_DIR <- file.path(RES_DIR, "figures_c")

MAIN_CSV <- file.path(RES_DIR, "main_RLlineage_trajectory_data.csv")
EXT_CSV <- file.path(RES_DIR, "ext_withIcn_trajectory_data.csv")
FULL_RDS <- file.path(RES_DIR, "aldinger2021_scored.rds")

RL_CLUSTERS <- c("02-RL", "03-GCP", "04-GN", "05-eCN/UBC")
EXT_CLUSTERS <- c(RL_CLUSTERS, "06-iCN")

cat("[1] Loading trajectory CSVs...\n")
main_dt <- fread(MAIN_CSV)
ext_dt <- fread(EXT_CSV)

cat(sprintf("  main: %d cells × %d cols\n", nrow(main_dt), ncol(main_dt)))
cat(sprintf("    Cols: %s\n", paste(names(main_dt), collapse = ", ")))
cat(sprintf("  ext:  %d cells × %d cols\n", nrow(ext_dt), ncol(ext_dt)))

cat("\n[2] Spearman correlations...\n")

compute_traj_stats <- function(dt, tag, lineage_labels) {
  stats_list <- list()
  score_cols <- c("sig_94_zscore", "sig_94_AMS", "sig_core_zscore", "sig_core_AMS")
  pt_cols <- grep("^pseudotime_lineage_", names(dt), value = TRUE)
  for (scol in score_cols) {
    for (i in seq_along(pt_cols)) {
      pt_vec <- dt[[pt_cols[i]]]
      score_vec <- dt[[scol]]
      valid <- !is.na(pt_vec)
      n_valid <- sum(valid)
      if (n_valid > 10) {
        rho <- cor(pt_vec[valid], score_vec[valid], method = "spearman")
        pval <- cor.test(pt_vec[valid], score_vec[valid], method = "spearman",
                          exact = FALSE)$p.value
      } else { rho <- NA; pval <- NA }
      stats_list[[length(stats_list) + 1]] <- data.table(
        subset = tag, signature = scol, lineage = i,
        lineage_name = lineage_labels[i],
        n_cells = n_valid, spearman_rho = rho, p_value = pval)
    }
  }
  rbindlist(stats_list)
}

main_lineages <- c("02-RL→03-GCP→04-GN", "02-RL→05-eCN/UBC")
ext_lineages <- c("02-RL→03-GCP→04-GN→06-iCN", "02-RL→05-eCN/UBC")

main_stats <- compute_traj_stats(main_dt, "main_RLlineage", main_lineages)
ext_stats <- compute_traj_stats(ext_dt, "ext_withIcn", ext_lineages)

all_stats <- rbindlist(list(main_stats, ext_stats))
fwrite(all_stats, file.path(RES_DIR, "R17c_trajectory_stats.csv"))
print(all_stats)

cat("\n[3] Per-cluster mean pt + score (main)...\n")

summary_by_ct_main <- list()
for (i in 1:2) {
  pt_col <- sprintf("pseudotime_lineage_%d", i)
  sub_dt <- main_dt[!is.na(get(pt_col)),
                    .(cell_type, sig_94_zscore, sig_94_AMS,
                      sig_core_zscore, sig_core_AMS, pt = get(pt_col))]
  agg <- sub_dt[, .(n = .N,
                     mean_pt = mean(pt, na.rm = TRUE),
                     mean_sig94 = mean(sig_94_zscore, na.rm = TRUE),
                     mean_sigcore = mean(sig_core_zscore, na.rm = TRUE)),
                 by = cell_type][order(mean_pt)]
  agg$lineage <- i
  summary_by_ct_main[[i]] <- agg
}
ct_main_dt <- rbindlist(summary_by_ct_main)
fwrite(ct_main_dt, file.path(RES_DIR, "R17c_cluster_pseudotime_main.csv"))
print(ct_main_dt)

cat("\n[4] Global UMAP highlight plot...\n")

if (file.exists(FULL_RDS)) {
  obj_full <- readRDS(FULL_RDS)
  
  global_umap_dt <- data.table(
    UMAP1 = Embeddings(obj_full, "umap")[, 1],
    UMAP2 = Embeddings(obj_full, "umap")[, 2],
    cell_type = obj_full$cell_type,
    is_RLlineage = obj_full$cell_type %in% RL_CLUSTERS,
    is_withIcn = obj_full$cell_type %in% EXT_CLUSTERS
  )
  
  p_g1 <- ggplot(global_umap_dt, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = global_umap_dt[is_RLlineage == FALSE],
               color = "grey85", size = 0.2, alpha = 0.5) +
    geom_point(data = global_umap_dt[is_RLlineage == TRUE],
               aes(color = cell_type), size = 0.3, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title = "Main subset (RL lineage, 12243 cells) on global UMAP",
         subtitle = "Grey = other cells in Aldinger 2021") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  p_g2 <- ggplot(global_umap_dt, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = global_umap_dt[is_withIcn == FALSE],
               color = "grey85", size = 0.2, alpha = 0.5) +
    geom_point(data = global_umap_dt[is_withIcn == TRUE],
               aes(color = cell_type), size = 0.3, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title = "Extended (+ 06-iCN, 16384 cells) on global UMAP",
         subtitle = "Grey = other cells") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  p_global <- p_g1 / p_g2
  ggsave(file.path(FIG_DIR, "00_global_UMAP_subset_highlight.png"),
         p_global, width = 10, height = 14, dpi = 150)
  cat("  Saved: 00_global_UMAP_subset_highlight.png\n")
  rm(obj_full); gc(verbose = FALSE)
} else {
  cat("  [!] aldinger2021_scored.rds not found, skipping global plot\n")
}

cat("\n[5] Writing SUMMARY...\n")

summary_lines <- c(
  "================================================================",
  "R17c — Neuronal lineage pseudotime on Aldinger 2021",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "MAIN ANALYSIS: RL lineage (4 clusters, 12243 cells)",
  "  Start: 02-RL at 9 PCW (centroid-closest)",
  "  Lineages detected: 2",
  sprintf("    Lineage 1: %s  (granule)", main_lineages[1]),
  sprintf("    Lineage 2: %s  (UBC)", main_lineages[2]),
  "",
  "EXTENDED: +06-iCN (5 clusters, 16384 cells)",
  sprintf("    Lineage 1: %s", ext_lineages[1]),
  sprintf("    Lineage 2: %s", ext_lineages[2]),
  "  NOTE: iCN at end of L1 is UMAP-geometric artifact",
  "        (iCN is VZ-origin, not rhombic lip). Use main analysis.",
  "",
  "================================================================",
  "SPEARMAN ρ (score vs pseudotime)",
  "================================================================",
  "",
  "[Main: RL lineage]"
)
for (i in 1:nrow(main_stats)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-18s L%d  ρ=%+.3f  p=%.2e  n=%d",
            main_stats$signature[i], main_stats$lineage[i],
            main_stats$spearman_rho[i], main_stats$p_value[i],
            main_stats$n_cells[i]))
}
summary_lines <- c(summary_lines, "",
  sprintf("  L1 = %s", main_lineages[1]),
  sprintf("  L2 = %s", main_lineages[2]))

summary_lines <- c(summary_lines, "", "[Extended: +iCN]")
for (i in 1:nrow(ext_stats)) {
  summary_lines <- c(summary_lines,
    sprintf("  %-18s L%d  ρ=%+.3f  p=%.2e  n=%d",
            ext_stats$signature[i], ext_stats$lineage[i],
            ext_stats$spearman_rho[i], ext_stats$p_value[i],
            ext_stats$n_cells[i]))
}
summary_lines <- c(summary_lines, "",
  sprintf("  L1 = %s", ext_lineages[1]),
  sprintf("  L2 = %s", ext_lineages[2]))

summary_lines <- c(summary_lines, "",
  "================================================================",
  "CLUSTER MEAN PSEUDOTIME + SCORE (main analysis)",
  "================================================================")
for (ln in 1:2) {
  summary_lines <- c(summary_lines, "",
    sprintf("Lineage %d (%s):", ln, main_lineages[ln]),
    "  cell_type       n     mean_pt   sig_94   sig_core")
  sub <- ct_main_dt[lineage == ln]
  for (i in 1:nrow(sub)) {
    summary_lines <- c(summary_lines,
      sprintf("  %-15s %5d  %8.2f  %+.3f   %+.3f",
              sub$cell_type[i], sub$n[i], sub$mean_pt[i],
              sub$mean_sig94[i], sub$mean_sigcore[i]))
  }
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  sprintf("  Figures: %s/figures_c/ (12 panels + global highlight)", FIG_DIR),
  sprintf("  Data: main_RLlineage_trajectory_data.csv (%d cells)", nrow(main_dt)),
  sprintf("        ext_withIcn_trajectory_data.csv (%d cells)", nrow(ext_dt)),
  sprintf("        R17c_trajectory_stats.csv"),
  sprintf("        R17c_cluster_pseudotime_main.csv"),
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(RES_DIR, "R17c_SUMMARY.txt"))
cat("  Saved: R17c_SUMMARY.txt\n")

cat("\n================================================================\n")
cat("R17c patch v2 DONE\n")
cat("================================================================\n")
