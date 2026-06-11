# ============================================================
# R21a v4: Ghasemi 2024 snRNA-seq pipeline
# Fix: Harmony API (dims.use → dims, use harmony 1.2+ syntax)
# Also: save intermediate Seurat object before Harmony for safety
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(harmony)
  library(patchwork)
})

set.seed(42)

# ---- paths ----
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
GHASEMI_DIR <- file.path(EXTDATA_DIR, "Ghasemi_Rademacher_2024")
SAMPLE_PATHS <- list(
  MB263 = file.path(GHASEMI_DIR, "GSM6604617_ICGC_MB263_10X_raw_counts.tsv.gz"),
  MB266 = file.path(GHASEMI_DIR, "GSM6604618_ICGC_MB266_10X_raw_counts.tsv.gz"),
  MB295 = file.path(GHASEMI_DIR, "GSM6604620_ICGC_MB295_10X_raw_counts.tsv.gz"),
  MB299 = file.path(GHASEMI_DIR, "GSM6604621_ICGC_MB299_10X_raw_counts.tsv.gz")
)

OUT_DIR <- file.path(RESULTS_DIR, "R21_Ghasemi")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Print harmony version for diagnosis
cat(sprintf("Harmony version: %s\n", packageVersion("harmony")))
cat(sprintf("Seurat version: %s\n\n", packageVersion("Seurat")))

cat("================================================================\n")
cat("R21a v4: Ghasemi 2024 snRNA-seq processing\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# Helper: read Ghasemi TSV
# ============================================================
read_ghasemi_tsv <- function(path) {
  cat(sprintf("    Reading %s...\n", basename(path)))
  
  con <- gzfile(path, "rt")
  header_line <- readLines(con, n = 1)
  close(con)
  cell_ids <- strsplit(header_line, "\t")[[1]]
  
  cat(sprintf("    Header: %d cell IDs\n", length(cell_ids)))
  
  dt <- fread(path, header = FALSE, sep = "\t", skip = 1,
              col.names = c("gene", cell_ids))
  
  gene_names <- dt$gene
  if (sum(duplicated(gene_names)) > 0) {
    gene_names <- make.unique(gene_names)
  }
  
  count_mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(count_mat) <- gene_names
  rm(dt); gc(verbose = FALSE)
  
  cat(sprintf("    Matrix: %d genes × %d cells\n",
              nrow(count_mat), ncol(count_mat)))
  count_mat
}

# ============================================================
# 1. Load + Seurat objects + QC
# ============================================================
cat("[1] Loading 4 samples...\n")

seurat_list <- list()
qc_summary <- data.table()

for (sn in names(SAMPLE_PATHS)) {
  cat(sprintf("\n  ==== %s ====\n", sn))
  count_mat <- read_ghasemi_tsv(SAMPLE_PATHS[[sn]])
  
  obj <- CreateSeuratObject(counts = count_mat,
                             project = sn,
                             min.cells = 3,
                             min.features = 200)
  obj$sample <- sn
  obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  cat(sprintf("    After CreateSeuratObject: %d cells, %d genes\n",
              ncol(obj), nrow(obj)))
  
  cells_before <- ncol(obj)
  obj <- subset(obj,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 8000 &
                         percent_mito < 20)
  cells_after <- ncol(obj)
  cat(sprintf("    After QC: %d → %d cells (%.1f%% kept)\n",
              cells_before, cells_after,
              100 * cells_after / cells_before))
  
  qc_summary <- rbind(qc_summary, data.table(
    sample = sn,
    cells_before_qc = cells_before,
    cells_after_qc = cells_after,
    pct_kept = round(100 * cells_after / cells_before, 1),
    median_nFeature = median(obj$nFeature_RNA),
    median_nCount = median(obj$nCount_RNA),
    median_percent_mito = round(median(obj$percent_mito), 2)
  ))
  
  seurat_list[[sn]] <- obj
  rm(count_mat); gc(verbose = FALSE)
}

cat("\n=== QC summary ===\n")
print(qc_summary)

# ============================================================
# 2. SCTransform per sample
# ============================================================
cat("\n[2] SCTransform v2 per sample...\n")

for (sn in names(seurat_list)) {
  cat(sprintf("  SCTransform %s...\n", sn))
  seurat_list[[sn]] <- SCTransform(seurat_list[[sn]],
                                    vst.flavor = "v2",
                                    verbose = FALSE,
                                    return.only.var.genes = FALSE)
}

# ============================================================
# 3. Integration features + Merge
# ============================================================
cat("\n[3] SelectIntegrationFeatures + Merge...\n")

integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)
cat(sprintf("  %d integration features selected\n", length(integration_features)))

combined <- merge(seurat_list[[1]], y = seurat_list[2:4],
                   add.cell.ids = names(seurat_list),
                   project = "Ghasemi_MBEN",
                   merge.data = TRUE)

DefaultAssay(combined) <- "SCT"
VariableFeatures(combined[["SCT"]]) <- integration_features

cat(sprintf("  Combined: %d cells, %d genes, %d var features\n",
            ncol(combined), nrow(combined), length(VariableFeatures(combined))))

rm(seurat_list); gc(verbose = FALSE)

# ============================================================
# 4. PCA
# ============================================================
cat("\n[4] PCA...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE,
                    features = integration_features)
cat(sprintf("  PCA done. Variance in first 5 PCs: %.1f%%\n",
            100 * sum(combined@reductions$pca@stdev[1:5]^2) /
              sum(combined@reductions$pca@stdev^2)))

# Save intermediate
cat("\n  [Checkpoint] Saving intermediate object before Harmony...\n")
saveRDS(combined, file.path(OUT_DIR, "ghasemi_seurat_pre_harmony.rds"))

# ============================================================
# 5. Harmony integration (with new API)
# ============================================================
cat("\n[5] Running Harmony (batch = sample)...\n")

# Detect Harmony version
harmony_ver <- as.character(packageVersion("harmony"))
cat(sprintf("  Harmony version detected: %s\n", harmony_ver))

# harmony 1.0+: param is 'dims', no longer 'dims.use'
# Try the modern API call
combined <- tryCatch({
  RunHarmony(combined,
             group.by.vars = "sample",
             reduction.use = "pca",
             dims.use = 1:30,
             max_iter = 20,
             verbose = TRUE)
}, error = function(e1) {
  cat(sprintf("  First attempt failed: %s\n", conditionMessage(e1)))
  cat("  Trying alternative API...\n")
  
  tryCatch({
    RunHarmony(combined,
               group.by.vars = "sample",
               reduction = "pca",
               dims = 1:30,
               max.iter.harmony = 20,
               verbose = TRUE)
  }, error = function(e2) {
    cat(sprintf("  Second attempt failed: %s\n", conditionMessage(e2)))
    cat("  Trying minimal API...\n")
    
    RunHarmony(combined,
               group.by.vars = "sample")
  })
})

# Verify Harmony worked
if (!"harmony" %in% Reductions(combined)) {
  stop("Harmony reduction not created. Check API call.")
}

cat(sprintf("  Harmony done. Reductions available: %s\n",
            paste(Reductions(combined), collapse = ", ")))

# Save another checkpoint
cat("  [Checkpoint] Saving post-Harmony object...\n")
saveRDS(combined, file.path(OUT_DIR, "ghasemi_seurat_post_harmony.rds"))

# ============================================================
# 6. UMAP + Clustering
# ============================================================
cat("\n[6] UMAP + Louvain clustering (resolution=0.5)...\n")

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

n_clusters <- length(unique(combined$seurat_clusters))
cat(sprintf("  Found %d clusters\n", n_clusters))

cluster_sizes <- table(combined$seurat_clusters)
cat("  Cluster sizes:\n")
print(cluster_sizes)

cat("\n  Cluster × Sample composition:\n")
comp <- table(combined$seurat_clusters, combined$sample)
print(comp)
cat("\n  Cluster composition (proportions per cluster):\n")
print(round(prop.table(comp, margin = 1), 3))

# ============================================================
# 7. Visualizations
# ============================================================
cat("\n[7] Rendering UMAPs...\n")

p_cluster <- DimPlot(combined, reduction = "umap",
                      group.by = "seurat_clusters",
                      label = TRUE, label.size = 4,
                      pt.size = 0.3) +
  labs(title = sprintf("R21a: %d clusters from 4 MBEN snRNA-seq samples",
                        n_clusters),
       subtitle = sprintf("n=%d cells, Harmony-integrated, resolution=0.5",
                           ncol(combined))) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "umap_clusters.png"),
       p_cluster, width = 10, height = 8, dpi = 150, bg = "white")

p_sample <- DimPlot(combined, reduction = "umap",
                     group.by = "sample",
                     pt.size = 0.3) +
  labs(title = "R21a: UMAP colored by sample (batch effect check)",
       subtitle = "Good integration = clusters NOT dominated by single sample") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "umap_by_sample.png"),
       p_sample, width = 10, height = 8, dpi = 150, bg = "white")

p_combined <- p_cluster + p_sample
ggsave(file.path(OUT_DIR, "umap_combined.png"),
       p_combined, width = 18, height = 8, dpi = 150, bg = "white")

# ============================================================
# 8. Cluster markers
# ============================================================
cat("\n[8] Finding cluster markers...\n")

DefaultAssay(combined) <- "RNA"
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
combined <- NormalizeData(combined, verbose = FALSE)

t0 <- Sys.time()
markers <- FindAllMarkers(combined,
                           only.pos = TRUE,
                           min.pct = 0.25,
                           logfc.threshold = 0.5,
                           verbose = FALSE)
cat(sprintf("  FindAllMarkers done in %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))

markers_dt <- as.data.table(markers)
top_markers <- markers_dt[, .SD[order(-avg_log2FC)][1:10],
                            by = cluster]
fwrite(top_markers, file.path(OUT_DIR, "cluster_top_markers.csv"))

cat("\n  Top 5 markers per cluster (preview):\n")
preview <- markers_dt[, .SD[order(-avg_log2FC)][1:5], by = cluster][, .(cluster, gene, avg_log2FC, pct.1, pct.2)]
print(preview)

DefaultAssay(combined) <- "SCT"

# ============================================================
# 9. Ghasemi marker DotPlot
# ============================================================
cat("\n[9] Ghasemi 2024 marker dotplot...\n")

ghasemi_markers <- list(
  early_CGNP_proliferating = c("MKI67", "TOP2A", "ATOH1", "BARHL1", "ZIC1", "ZIC3"),
  early_CGNP_quiescent = c("PTCH1", "SMO", "HHIP", "GLI1", "GLI2", "PTPRK"),
  migrating = c("GRIN2B", "CNTN2", "ASTN1", "SEMA6A"),
  postmitotic_differentiated = c("GABRA1", "GABRA6", "GRIN2C", "RBFOX3"),
  astrocytic_like = c("LAMA2", "SOX2", "SOX9", "GFAP", "AQP4")
)

all_markers_flat <- unlist(ghasemi_markers)
present_markers <- intersect(all_markers_flat, rownames(combined))
absent_markers <- setdiff(all_markers_flat, rownames(combined))

cat(sprintf("  Coverage: %d / %d (%.1f%%)\n",
            length(present_markers), length(all_markers_flat),
            100 * length(present_markers) / length(all_markers_flat)))

if (length(present_markers) > 0) {
  p_dot <- DotPlot(combined, features = present_markers,
                    group.by = "seurat_clusters",
                    cols = c("lightgrey", "red")) +
    RotatedAxis() +
    labs(title = "R21a: Ghasemi 2024 cell-state markers across clusters",
         x = "Marker gene", y = "Cluster") +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))
  
  ggsave(file.path(OUT_DIR, "ghasemi_marker_expression.png"),
         p_dot, width = 14, height = 7, dpi = 150, bg = "white")
}

# ============================================================
# 10. Save final object + SUMMARY
# ============================================================
cat("\n[10] Saving final Seurat object...\n")
saveRDS(combined, file.path(OUT_DIR, "ghasemi_seurat_R21a.rds"))

# Cleanup intermediate checkpoints once final object is saved
file.remove(file.path(OUT_DIR, "ghasemi_seurat_pre_harmony.rds"))
file.remove(file.path(OUT_DIR, "ghasemi_seurat_post_harmony.rds"))

cat("[11] Writing R21a_SUMMARY.txt...\n")

summary_lines <- c(
  "================================================================",
  "R21a — Ghasemi 2024 snRNA-seq pipeline + clustering",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "DATA",
  "  4 samples from GSE239854 matching our Xenium cohort:",
  "    MB263 (GSM6604617)",
  "    MB266 (GSM6604618)",
  "    MB295 (GSM6604620)",
  "    MB299 (GSM6604621)",
  "",
  "================================================================",
  "QC SUMMARY",
  "================================================================",
  ""
)

for (i in 1:nrow(qc_summary)) {
  summary_lines <- c(summary_lines,
    sprintf("[%s]  %d → %d cells (%.1f%% kept)",
            qc_summary$sample[i],
            qc_summary$cells_before_qc[i],
            qc_summary$cells_after_qc[i],
            qc_summary$pct_kept[i]),
    sprintf("    median nFeature: %d, nCount: %d, %%mito: %.2f",
            qc_summary$median_nFeature[i],
            qc_summary$median_nCount[i],
            qc_summary$median_percent_mito[i]),
    "")
}

summary_lines <- c(summary_lines,
  "================================================================",
  "INTEGRATION + CLUSTERING",
  "================================================================",
  sprintf("  Total cells after QC: %d", ncol(combined)),
  sprintf("  Integration features: %d", length(integration_features)),
  "  Method: SCT v2 per sample → SelectIntegrationFeatures → Merge → PCA → Harmony → UMAP",
  sprintf("  Clusters: %d (Louvain, resolution=0.5)", n_clusters),
  "",
  "  Cluster sizes:")

for (cl in names(cluster_sizes)) {
  summary_lines <- c(summary_lines,
    sprintf("    Cluster %s: %d cells (%.1f%%)",
            cl, cluster_sizes[[cl]],
            100 * cluster_sizes[[cl]] / ncol(combined)))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "CLUSTER × SAMPLE COMPOSITION (proportions)",
  "================================================================",
  "")

prop_mat <- prop.table(comp, margin = 1)
header <- sprintf("  %-12s %s", "Cluster",
                   paste(sprintf("%-8s", colnames(prop_mat)), collapse = ""))
summary_lines <- c(summary_lines, header)
for (i in 1:nrow(prop_mat)) {
  row_str <- sprintf("  %-12s %s", paste0("Cluster_", rownames(prop_mat)[i]),
                      paste(sprintf("%-8.3f", prop_mat[i, ]), collapse = ""))
  summary_lines <- c(summary_lines, row_str)
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "GHASEMI 2024 MARKER COVERAGE",
  "================================================================",
  sprintf("  Coverage: %d / %d (%.1f%%)",
          length(present_markers), length(all_markers_flat),
          100 * length(present_markers) / length(all_markers_flat))
)
if (length(absent_markers) > 0) {
  summary_lines <- c(summary_lines,
    sprintf("  Missing: %s", paste(absent_markers, collapse = ", ")))
}

summary_lines <- c(summary_lines, "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  ghasemi_seurat_R21a.rds          - integrated Seurat object",
  "  umap_clusters.png                - UMAP by cluster",
  "  umap_by_sample.png               - UMAP by sample",
  "  umap_combined.png                - both side by side",
  "  cluster_top_markers.csv          - top 10 markers per cluster",
  "  ghasemi_marker_expression.png    - DotPlot Ghasemi markers x clusters",
  "  R21a_SUMMARY.txt                 - this file",
  "",
  "================================================================"
)

writeLines(summary_lines, file.path(OUT_DIR, "R21a_SUMMARY.txt"))

cat("\n=== R21a SUMMARY ===\n")
cat(paste(summary_lines, collapse="\n"), "\n")

cat("\n================================================================\n")
cat("R21a DONE\n")
cat("================================================================\n")
