# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R21a: Ghasemi 2024 snRNA-seq pipeline (4 samples)
# v2: fix TSV reading - manually parse header (cell IDs) and data
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
SAMPLE_PATHS <- list(
  MB263 = "./external_data/Ghasemi_Rademacher_2024/GSM6604617_ICGC_MB263_10X_raw_counts.tsv.gz",
  MB266 = "./external_data/Ghasemi_Rademacher_2024/GSM6604618_ICGC_MB266_10X_raw_counts.tsv.gz",
  MB295 = "./external_data/Ghasemi_Rademacher_2024/GSM6604620_ICGC_MB295_10X_raw_counts.tsv.gz",
  MB299 = "./external_data/Ghasemi_Rademacher_2024/GSM6604621_ICGC_MB299_10X_raw_counts.tsv.gz"
)

OUT_DIR <- "./results/R21_Ghasemi"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("R21a v2: Ghasemi 2024 snRNA-seq processing\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# Helper: read Ghasemi TSV format
# Format: line 1 = cell IDs (N entries, no gene column header)
#         lines 2+ = gene_name TAB val1 TAB val2 ... TAB valN  (N+1 entries)
# ============================================================
read_ghasemi_tsv <- function(path) {
  cat(sprintf("    Reading %s...\n", basename(path)))
  
  # Step 1: read header line manually to get cell IDs
  con <- gzfile(path, "rt")
  header_line <- readLines(con, n = 1)
  close(con)
  cell_ids <- strsplit(header_line, "\t")[[1]]
  n_cells <- length(cell_ids)
  cat(sprintf("    Header: %d cell IDs (first: %s, last: %s)\n",
              n_cells, cell_ids[1], cell_ids[n_cells]))
  
  # Step 2: read data rows with fread, skipping header
  # All data rows have N+1 columns (gene + N values)
  dt <- fread(path, header = FALSE, sep = "\t", skip = 1,
              col.names = c("gene", cell_ids))
  
  cat(sprintf("    Data rows: %d, columns: %d\n", nrow(dt), ncol(dt)))
  
  # Step 3: extract gene names, dedup, build matrix
  gene_names <- dt$gene
  dup_count <- sum(duplicated(gene_names))
  if (dup_count > 0) {
    cat(sprintf("    %d duplicate gene names found, applying make.unique\n", dup_count))
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
# 1. Load each sample, build Seurat object, basic QC
# ============================================================
cat("[1] Loading 4 samples...\n")

seurat_list <- list()
qc_summary <- data.table()

for (sn in names(SAMPLE_PATHS)) {
  cat(sprintf("\n  ==== %s ====\n", sn))
  
  count_mat <- read_ghasemi_tsv(SAMPLE_PATHS[[sn]])
  
  # Build Seurat object
  obj <- CreateSeuratObject(counts = count_mat,
                             project = sn,
                             min.cells = 3,
                             min.features = 200)
  obj$sample <- sn
  obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  cat(sprintf("    After CreateSeuratObject: %d cells, %d genes\n",
              ncol(obj), nrow(obj)))
  cat(sprintf("    QC stats: median nFeature=%.0f, nCount=%.0f, %%mito=%.2f\n",
              median(obj$nFeature_RNA),
              median(obj$nCount_RNA),
              median(obj$percent_mito)))
  
  # QC filter
  cells_before <- ncol(obj)
  obj <- subset(obj,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 8000 &
                         percent_mito < 20)
  cells_after <- ncol(obj)
  cat(sprintf("    QC filter: %d → %d cells (%.1f%% kept)\n",
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
# 2. SCTransform per sample, then merge
# ============================================================
cat("\n[2] SCTransform v2 per sample...\n")

for (sn in names(seurat_list)) {
  cat(sprintf("  SCTransform %s...\n", sn))
  seurat_list[[sn]] <- SCTransform(seurat_list[[sn]],
                                    vst.flavor = "v2",
                                    verbose = FALSE,
                                    return.only.var.genes = FALSE)
}

cat("\n  Merging 4 SCT-normalized objects...\n")
combined <- merge(seurat_list[[1]], y = seurat_list[2:4],
                   add.cell.ids = names(seurat_list),
                   project = "Ghasemi_MBEN")

DefaultAssay(combined) <- "SCT"
combined <- FindVariableFeatures(combined, nfeatures = 3000, verbose = FALSE)

var_features_union <- unique(unlist(lapply(seurat_list, VariableFeatures)))
VariableFeatures(combined) <- intersect(var_features_union, rownames(combined))

cat(sprintf("  Combined object: %d cells, %d genes, %d variable features\n",
            ncol(combined), nrow(combined), length(VariableFeatures(combined))))

rm(seurat_list); gc(verbose = FALSE)

# ============================================================
# 3. PCA + Harmony integration
# ============================================================
cat("\n[3] PCA + Harmony integration...\n")

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
cat(sprintf("  PCA done. Variance in first 5 PCs: %.1f%%\n",
            100 * sum(combined@reductions$pca@stdev[1:5]^2) /
              sum(combined@reductions$pca@stdev^2)))

cat("  Running Harmony (batch = sample)...\n")
combined <- RunHarmony(combined,
                        group.by.vars = "sample",
                        reduction = "pca",
                        dims.use = 1:30,
                        max.iter.harmony = 20,
                        verbose = FALSE)

# ============================================================
# 4. UMAP + Clustering
# ============================================================
cat("\n[4] UMAP + Louvain clustering (resolution=0.5)...\n")

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

n_clusters <- length(unique(combined$seurat_clusters))
cat(sprintf("  Found %d clusters\n", n_clusters))

cluster_sizes <- table(combined$seurat_clusters)
cat("  Cluster sizes:\n")
print(cluster_sizes)

cat("\n  Cluster composition by sample:\n")
comp <- table(combined$seurat_clusters, combined$sample)
print(comp)
cat("\n  Cluster composition (proportions per cluster):\n")
print(round(prop.table(comp, margin = 1), 3))

# ============================================================
# 5. Visualizations
# ============================================================
cat("\n[5] Rendering UMAPs...\n")

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
# 6. Top markers per cluster
# ============================================================
cat("\n[6] Finding cluster markers (top 10 per cluster)...\n")

DefaultAssay(combined) <- "RNA"

# Need to join layers in Seurat v5 before normalization
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
# 7. Sanity check: Ghasemi 2024 known markers
# ============================================================
cat("\n[7] Ghasemi 2024 marker dotplot...\n")

ghasemi_markers <- list(
  early_CGNP_proliferating = c("MKI67", "TOP2A", "ATOH1", "BARHL1", "ZIC1", "ZIC3"),
  early_CGNP_quiescent = c("PTCH1", "SMO", "HHIP", "GLI1", "GLI2", "PTPRK"),
  migrating = c("GRIN2B", "CNTN2", "ASTN1", "SEMA6A"),
  postmitotic_differentiated = c("GABRA1", "GABRA6", "GRIN2C", "RBFOX3"),
  astrocytic_like = c("LAMA2", "SOX2", "SOX9", "GFAP", "AQP4")
)

all_markers_flat <- unlist(ghasemi_markers)
present_markers <- intersect(all_markers_flat, rownames(combined))
cat(sprintf("  Ghasemi marker coverage: %d / %d (%.1f%%)\n",
            length(present_markers), length(all_markers_flat),
            100 * length(present_markers) / length(all_markers_flat)))

absent_markers <- setdiff(all_markers_flat, rownames(combined))
if (length(absent_markers) > 0) {
  cat(sprintf("  Missing: %s\n", paste(absent_markers, collapse = ", ")))
}

p_dot <- DotPlot(combined, features = present_markers,
                  group.by = "seurat_clusters",
                  cols = c("lightgrey", "red")) +
  RotatedAxis() +
  labs(title = "R21a: Ghasemi 2024 cell-state markers across clusters",
       subtitle = "Each cluster's expression of canonical MBEN cell-state markers",
       x = "Marker gene", y = "Cluster") +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))

ggsave(file.path(OUT_DIR, "ghasemi_marker_expression.png"),
       p_dot, width = 14, height = 7, dpi = 150, bg = "white")

# ============================================================
# 8. Save + SUMMARY
# ============================================================
cat("\n[8] Saving Seurat object...\n")
saveRDS(combined, file.path(OUT_DIR, "ghasemi_seurat_R21a.rds"))

cat("[9] Writing R21a_SUMMARY.txt...\n")

summary_lines <- c(
  "================================================================",
  "R21a — Ghasemi 2024 snRNA-seq pipeline + clustering",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "DATA",
  "  4 samples from GSE239854 (Ghasemi 2024):",
  "    MB263  -  GSM6604617",
  "    MB266  -  GSM6604618",
  "    MB295  -  GSM6604620",
  "    MB299  -  GSM6604621",
  "",
  "  Same patients as our 4 Xenium spatial samples (GSE291688).",
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
  sprintf("  Variable features: %d", length(VariableFeatures(combined))),
  "  Method: SCTransform v2 → PCA(30) → Harmony(batch=sample) → UMAP",
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
  sprintf("  Coverage: %d / %d markers present in data",
          length(present_markers), length(all_markers_flat)),
  ""
)

if (length(absent_markers) > 0) {
  summary_lines <- c(summary_lines,
    sprintf("  Missing: %s", paste(absent_markers, collapse = ", ")))
}

summary_lines <- c(summary_lines, "",
  "  Marker sets:",
  sprintf("    early_CGNP_proliferating: %s", paste(ghasemi_markers$early_CGNP_proliferating, collapse = ", ")),
  sprintf("    early_CGNP_quiescent: %s", paste(ghasemi_markers$early_CGNP_quiescent, collapse = ", ")),
  sprintf("    migrating: %s", paste(ghasemi_markers$migrating, collapse = ", ")),
  sprintf("    postmitotic_differentiated: %s", paste(ghasemi_markers$postmitotic_differentiated, collapse = ", ")),
  sprintf("    astrocytic_like: %s", paste(ghasemi_markers$astrocytic_like, collapse = ", ")),
  "",
  "================================================================",
  "OUTPUTS",
  "================================================================",
  "  ghasemi_seurat_R21a.rds          - integrated Seurat object",
  "  umap_clusters.png                - UMAP by cluster",
  "  umap_by_sample.png               - UMAP by sample (batch check)",
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
cat("R21a DONE — review UMAP and cluster markers, then approve R21b\n")
cat("================================================================\n")
