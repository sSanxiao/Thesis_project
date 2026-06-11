# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
suppressPackageStartupMessages({library(Seurat); library(data.table); library(jsonlite)})
registry <- fromJSON("./data/sample_registry.json")
R2_DIR <- "./results/R2_Results"
qc_list <- list()
for (sname in names(registry)) {
  parts <- strsplit(sname, "/")[[1]]
  rds_path <- file.path(R2_DIR, parts[1], parts[2], paste0(parts[2], "_seurat_R2.rds"))
  if (!file.exists(rds_path)) { cat("Missing:", sname, "\n"); next }
  cat("Reading:", sname, "...")
  obj <- readRDS(rds_path)
  sct_data <- GetAssayData(obj, assay="SCT", layer="data")
  pc_val <- tryCatch(obj$n_pcs_selected[1], error=function(e) NA)
  qc_list[[sname]] <- data.frame(
    sample_name=sname, dataset=parts[1], n_genes=nrow(obj), n_cells=ncol(obj),
    n_var_features=length(VariableFeatures(obj)), n_pcs_selected=pc_val,
    n_clusters=length(levels(Idents(obj))),
    residual_mean=round(mean(sct_data@x),4), residual_sd=round(sd(sct_data@x),4),
    rds_size_mb=round(file.size(rds_path)/1024/1024,1), stringsAsFactors=FALSE)
  rm(obj, sct_data); gc(verbose=FALSE)
  cat(" done\n")
}
qc_df <- do.call(rbind, qc_list)
fwrite(qc_df, file.path(R2_DIR, "ALL_SAMPLES_R2_QC_FULL.csv"))
cat("\nSaved full QC:", nrow(qc_df), "samples\n")
print(qc_df[, c("sample_name","n_genes","n_cells","n_pcs_selected","n_clusters","rds_size_mb")])
