# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# 在R2脚本开头设置future限制，然后source它
options(future.globals.maxSize = 8000 * 1024^2)

# 读取原始脚本内容
lines <- readLines("./R2_sctransform.R")

# 在循环中的"创建输出目录"之前插入跳过逻辑
skip_code <- '
  # 跳过已完成的样本
  out_rds_check <- file.path(OUTPUT_DIR, dataset_name, sample_subname,
                             paste0(sample_subname, "_seurat_R2.rds"))
  if (file.exists(out_rds_check)) {
    cat("  已完成, 跳过\n\n")
    tryCatch({
      prev_obj <- readRDS(out_rds_check)
      pc_val <- tryCatch(prev_obj$n_pcs_selected[1], error=function(e) NA)
      qc_list[[sample_name]] <<- data.frame(
        sample_name = sample_name, dataset = dataset_name,
        n_genes = nrow(prev_obj), n_cells = ncol(prev_obj),
        n_var_features = length(VariableFeatures(prev_obj)),
        n_pcs_selected = pc_val,
        n_pcs_elbow = NA, n_pcs_threshold = NA,
        cum_var_pct = NA,
        n_clusters = length(levels(Idents(prev_obj))),
        residual_mean = NA, residual_sd = NA,
        time_seconds = 0,
        rds_size_mb = round(file.size(out_rds_check)/1024/1024, 1),
        stringsAsFactors = FALSE)
      rm(prev_obj); gc(verbose=FALSE)
    }, error = function(e) cat("  无法读取已有rds\n"))
    next
  }
'

# 找到插入位置
idx <- grep("# --- 创建输出目录 ---", lines)
if (length(idx) == 1) {
  new_lines <- c(lines[1:idx-1], skip_code, lines[idx:length(lines)])
  tmp <- tempfile(fileext = ".R")
  writeLines(new_lines, tmp)
  source(tmp)
} else {
  cat("找不到插入位置，直接运行原始脚本\n")
  source("./R2_sctransform.R")
}
