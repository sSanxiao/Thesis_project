############################################################
# R2_normalize_cluster.R
# 功能: 批量 SCTransform 归一化 + PCA + UMAP + 聚类
# 输入: seurat_raw.rds (R1 输出)
# 输出: seurat_umap.rds (每样本一个)
#
# 注意: SCTransform 是计算密集型操作, 21 个样本串行可能
#       需要数小时。建议 nohup 后台运行。
#       ATRT/33 (45万细胞) 和 ATRT/32 (34万细胞) 会特别慢。
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R2_normalize_cluster.R > R2_run.log 2>&1 &
#   tail -f R2_run.log
############################################################

library(Seurat)
library(future)

# 并行 + 大内存
plan("multicore", workers = 8)
options(future.globals.maxSize = 200 * 1024^3)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"

# 默认参数 (后续 R5 会验证这些值是否合理)
DIMS       <- 1:30
RESOLUTION <- 0.3

# ===========================================================
# 样本注册表 (与 R1 完全一致)
# ===========================================================

sample_list <- list(
  list(project="Alzheimer_Mouse", sample="Wild_13_4",     species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_5_7",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_2_5",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_17_9",  species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_5_7",   species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_2_5",   species="mouse", condition="TgCRND8_AD"),
  list(project="Brain_Human",     sample="Alz",           species="human", condition="alzheimer"),
  list(project="Brain_Human",     sample="Gilo",          species="human", condition="glioblastoma"),
  list(project="Brain_Human",     sample="Healthy",       species="human", condition="healthy_brain"),
  list(project="Brain_Mouse",     sample="single",        species="mouse", condition="normal_brain"),
  list(project="ATRT_Human",      sample="28",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="29",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="30",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="31",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="32",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="33",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="34",            species="human", condition="ATRT"),
  list(project="Medulloblastoma_Human", sample="GSM8840046_MB263", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840047_MB266", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840048_MB295", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840049_MB299", species="human", condition="medulloblastoma")
)

# ===========================================================
# 主函数
# ===========================================================

normalize_one_sample <- function(s) {

  project <- s$project
  sname   <- s$sample

  sample_dir <- file.path(RESULTS_ROOT, project, sname)
  rds_in     <- file.path(sample_dir, "seurat_raw.rds")
  rds_out    <- file.path(sample_dir, "seurat_umap.rds")

  # 跳过已完成的样本 (断点恢复)
  if (file.exists(rds_out)) {
    fsize <- round(file.size(rds_out) / 1024^2, 1)
    cat("  [SKIP] Already exists:", rds_out, "(", fsize, "MB)\n")
    return(data.frame(
      project = project, sample = sname,
      status = "skipped", time_min = 0,
      n_clusters = NA, rds_size_MB = fsize,
      stringsAsFactors = FALSE
    ))
  }

  if (!file.exists(rds_in)) {
    cat("  [ERROR] seurat_raw.rds not found\n")
    return(data.frame(
      project = project, sample = sname,
      status = "missing_input", time_min = 0,
      n_clusters = NA, rds_size_MB = NA,
      stringsAsFactors = FALSE
    ))
  }

  t0 <- Sys.time()

  # --- 读取 ---
  cat("  Loading seurat_raw.rds...\n")
  so <- readRDS(rds_in)
  cat("  Cells:", ncol(so), " Genes:", nrow(so), "\n")

  # --- SCTransform ---
  cat("  SCTransform...\n")
  so <- SCTransform(so, assay = "spatial", verbose = FALSE)

  # --- PCA ---
  cat("  PCA...\n")
  so <- RunPCA(so, assay = "SCT", verbose = FALSE)

  # --- UMAP ---
  cat("  UMAP (dims=1:", max(DIMS), ")...\n", sep = "")
  so <- RunUMAP(so, reduction = "pca", dims = DIMS, verbose = FALSE)

  # --- Neighbors + Clusters ---
  cat("  FindNeighbors + FindClusters (res=", RESOLUTION, ")...\n", sep = "")
  so <- FindNeighbors(so, reduction = "pca", dims = DIMS, verbose = FALSE)
  so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)

  n_cl <- length(unique(Idents(so)))
  cat("  Clusters found:", n_cl, "\n")

  # --- 保存 ---
  saveRDS(so, rds_out)
  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  fsize   <- round(file.size(rds_out) / 1024^2, 1)

  cat("  Saved:", rds_out, "(", fsize, "MB,", elapsed, "min)\n")

  return(data.frame(
    project     = project,
    sample      = sname,
    status      = "done",
    time_min    = as.numeric(elapsed),
    n_clusters  = n_cl,
    rds_size_MB = fsize,
    stringsAsFactors = FALSE
  ))
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R2: NORMALIZE + CLUSTER\n")
cat("  Samples:", length(sample_list), "\n")
cat("  Params:  dims=1:", max(DIMS), ", resolution=", RESOLUTION, "\n", sep = "")
cat("  Output:  ", RESULTS_ROOT, "\n")
cat("  断点恢复: 已有 seurat_umap.rds 的样本会自动跳过\n")
cat("============================================================\n\n")

t_total <- Sys.time()
all_qc  <- list()

for (i in seq_along(sample_list)) {
  s     <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("\n[%d/%d] %s\n", i, length(sample_list), label))

  qc <- normalize_one_sample(s)
  all_qc[[i]] <- qc
}

# ===========================================================
# 汇总 QC
# ===========================================================

qc_dir <- file.path(RESULTS_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

qc_all <- do.call(rbind, all_qc)
write.csv(qc_all,
          file.path(qc_dir, "R2_normalize_cluster_qc.csv"),
          row.names = FALSE)

total_time <- round(difftime(Sys.time(), t_total, units = "mins"), 1)
done_count <- sum(qc_all$status == "done")
skip_count <- sum(qc_all$status == "skipped")

cat("\n============================================================\n")
cat("  R2 COMPLETE\n")
cat("  New:", done_count, " Skipped:", skip_count, "\n")
cat("  Total time:", total_time, "min\n")
cat("  QC saved:", file.path(qc_dir, "R2_normalize_cluster_qc.csv"), "\n")
cat("============================================================\n")

print(qc_all)
