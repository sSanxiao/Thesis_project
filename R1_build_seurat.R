############################################################
# R1_build_seurat.R
# 功能: 批量构建 21 个样本的 Seurat 对象
# 输入: xenium_count_matrix.csv + xenium_cell_metadata.csv
# 输出: seurat_raw.rds (每样本一个)
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R1_build_seurat.R > R1_run.log 2>&1 &
#   tail -f R1_run.log
############################################################

library(Seurat)
library(data.table)

# ===========================================================
# 配置
# ===========================================================

INPUT_ROOT  <- "/home/disk/wangqilu/Stage_2_Results_v2"
OUTPUT_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"

# ===========================================================
# 样本注册表
# ===========================================================

sample_list <- list(
  # Alzheimer_Mouse (6)
  list(project="Alzheimer_Mouse", sample="Wild_13_4",     species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_5_7",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="Wild_2_5",      species="mouse", condition="wild_type"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_17_9",  species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_5_7",   species="mouse", condition="TgCRND8_AD"),
  list(project="Alzheimer_Mouse", sample="TgCRND8_2_5",   species="mouse", condition="TgCRND8_AD"),
  # Brain_Human (3)
  list(project="Brain_Human",     sample="Alz",           species="human", condition="alzheimer"),
  list(project="Brain_Human",     sample="Gilo",          species="human", condition="glioblastoma"),
  list(project="Brain_Human",     sample="Healthy",       species="human", condition="healthy_brain"),
  # Brain_Mouse (1) — 特殊路径: single/
  list(project="Brain_Mouse",     sample="single",        species="mouse", condition="normal_brain"),
  # ATRT_Human (7)
  list(project="ATRT_Human",      sample="28",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="29",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="30",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="31",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="32",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="33",            species="human", condition="ATRT"),
  list(project="ATRT_Human",      sample="34",            species="human", condition="ATRT"),
  # Medulloblastoma_Human (4)
  list(project="Medulloblastoma_Human", sample="GSM8840046_MB263", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840047_MB266", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840048_MB295", species="human", condition="medulloblastoma"),
  list(project="Medulloblastoma_Human", sample="GSM8840049_MB299", species="human", condition="medulloblastoma")
)

# ===========================================================
# 主函数: 处理单个样本
# ===========================================================

build_one_seurat <- function(s) {

  project <- s$project
  sname   <- s$sample
  species <- s$species
  cond    <- s$condition

  # --- 路径 ---
  input_dir  <- file.path(INPUT_ROOT, project, sname)
  output_dir <- file.path(OUTPUT_ROOT, project, sname)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cm_file   <- file.path(input_dir, "xenium_count_matrix.csv")
  meta_file <- file.path(input_dir, "xenium_cell_metadata.csv")

  # --- 检查文件 ---
  if (!file.exists(cm_file)) {
    cat("  [SKIP] count matrix not found:", cm_file, "\n")
    return(NULL)
  }
  if (!file.exists(meta_file)) {
    cat("  [SKIP] metadata not found:", meta_file, "\n")
    return(NULL)
  }

  # --- 读取 count matrix ---
  cat("  Reading count matrix...\n")
  dt       <- fread(cm_file, header = TRUE)
  cell_ids <- as.character(dt$cell_id)
  gnames   <- setdiff(colnames(dt), "cell_id")
  mat      <- as.matrix(dt[, ..gnames])
  rownames(mat) <- cell_ids
  mat      <- t(mat)  # genes × cells

  # --- 读取 metadata ---
  cat("  Reading metadata...\n")
  metadata <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)

  # --- cell_id 对齐 ---
  common <- intersect(colnames(mat), rownames(metadata))
  if (length(common) == 0) {
    cat("  [ERROR] No shared cell IDs between matrix and metadata\n")
    return(NULL)
  }
  mat      <- mat[, common, drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  metadata$x_centroid <- as.numeric(metadata$x_centroid)
  metadata$y_centroid <- as.numeric(metadata$y_centroid)

  cat("  Matrix:", nrow(mat), "genes x", ncol(mat), "cells\n")

  # --- 创建 Seurat 对象 ---
  cat("  Building Seurat object...\n")
  so <- CreateSeuratObject(counts = mat, meta.data = metadata)

  # 注册空间坐标到 metadata (供后续直接取用)
  so$x_centroid <- metadata$x_centroid
  so$y_centroid <- metadata$y_centroid

  # 注册样本元信息
  so$project_id <- project
  so$sample_id  <- sname
  so$species    <- species
  so$condition  <- cond

  # 创建 spatial assay
  so[["spatial"]] <- CreateAssayObject(counts = mat)
  DefaultAssay(so) <- "spatial"

  # 注册空间坐标到 images 槽 (SlideSeq 格式, 供 SpatialFeaturePlot 等使用)
  coords_df <- metadata[, c("x_centroid", "y_centroid")]
  so@images <- list()
  so@images[["image"]] <- new(
    Class       = "SlideSeq",
    assay       = "spatial",
    coordinates = coords_df,
    key         = "image_"
  )

  # --- 保存 ---
  rds_path <- file.path(output_dir, "seurat_raw.rds")
  saveRDS(so, rds_path)

  # --- QC 记录 ---
  qc <- data.frame(
    project   = project,
    sample    = sname,
    species   = species,
    condition = cond,
    n_genes   = nrow(mat),
    n_cells   = ncol(mat),
    rds_size_MB = round(file.size(rds_path) / 1024^2, 1),
    stringsAsFactors = FALSE
  )

  cat("  Saved:", rds_path, "\n")
  return(qc)
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R1: BUILD SEURAT OBJECTS\n")
cat("  Samples:", length(sample_list), "\n")
cat("  Input root: ", INPUT_ROOT, "\n")
cat("  Output root:", OUTPUT_ROOT, "\n")
cat("============================================================\n\n")

t_start <- Sys.time()
all_qc  <- list()

for (i in seq_along(sample_list)) {
  s <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("[%d/%d] %s (%s, %s)\n", i, length(sample_list), label, s$species, s$condition))

  qc <- build_one_seurat(s)
  if (!is.null(qc)) {
    all_qc[[i]] <- qc
    cat("  OK.\n\n")
  } else {
    cat("  FAILED.\n\n")
  }
}

# ===========================================================
# 汇总 QC
# ===========================================================

qc_dir <- file.path(OUTPUT_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

qc_all <- do.call(rbind, all_qc)
write.csv(qc_all,
          file.path(qc_dir, "R1_build_seurat_qc.csv"),
          row.names = FALSE)

elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)

cat("============================================================\n")
cat("  R1 COMPLETE\n")
cat("  Samples processed:", nrow(qc_all), "/", length(sample_list), "\n")
cat("  Total time:", elapsed, "min\n")
cat("  QC saved:", file.path(qc_dir, "R1_build_seurat_qc.csv"), "\n")
cat("============================================================\n")

# 打印汇总表
print(qc_all)
