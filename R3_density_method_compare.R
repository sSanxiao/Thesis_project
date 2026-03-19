############################################################
# R3_density_method_compare.R
# 功能: 读取三种 density, 各跑一次全基因组 Spearman 相关,
#       比较三种方法的信号强度和基因重叠度, 选定主方法
# 输入: seurat_umap.rds (R2) + cell_density_three_methods.csv (P2)
# 输出: 每样本 density_results_{method}.csv ×3
#       每样本 method_comparison.csv
#       全局  METHOD_SELECTION_REPORT.csv
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R3_density_method_compare.R > R3_run.log 2>&1 &
#   tail -f R3_run.log
############################################################

library(Seurat)
library(data.table)
library(dplyr)

# ===========================================================
# 配置
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"
DENSITY_ROOT <- "/home/disk/wangqilu/Density_Caculation/Results"

FDR_THRESHOLD <- 0.05
# 使用一个宽松的初始阈值, R4 会精细验证
COR_THRESHOLD <- 0.05

METHODS <- c("density_knn", "density_voronoi", "density_delaunay")
METHOD_LABELS <- c("KNN", "Voronoi", "Delaunay")

# ===========================================================
# P2 → R 路径映射
# P2 的输出目录名和 R 管线的样本名不完全一致
# ===========================================================

# key = R管线中的 "project/sample"
# value = P2 输出中的子目录名
P2_DIR_MAP <- list(
  "Alzheimer_Mouse/Wild_13_4"    = "Alzheimer_Mouse/Wild_13_4",
  "Alzheimer_Mouse/Wild_5_7"     = "Alzheimer_Mouse/Wild_5_7",
  "Alzheimer_Mouse/Wild_2_5"     = "Alzheimer_Mouse/Wild_2_5",
  "Alzheimer_Mouse/TgCRND8_17_9" = "Alzheimer_Mouse/TgCRND8_17_9",
  "Alzheimer_Mouse/TgCRND8_5_7"  = "Alzheimer_Mouse/TgCRND8_5_7",
  "Alzheimer_Mouse/TgCRND8_2_5"  = "Alzheimer_Mouse/TgCRND8_2_5",
  "Brain_Human/Alz"              = "Brain_Human/Alz",
  "Brain_Human/Gilo"             = "Brain_Human/Glio",
  "Brain_Human/Healthy"          = "Brain_Human/Healthy",
  "Brain_Mouse/single"           = "Brain_Mouse/Normal",
  "ATRT_Human/28"                = "ATRT_Human/28",
  "ATRT_Human/29"                = "ATRT_Human/29",
  "ATRT_Human/30"                = "ATRT_Human/30",
  "ATRT_Human/31"                = "ATRT_Human/31",
  "ATRT_Human/32"                = "ATRT_Human/32",
  "ATRT_Human/33"                = "ATRT_Human/33",
  "ATRT_Human/34"                = "ATRT_Human/34",
  "Medulloblastoma_Human/GSM8840046_MB263" = "Medulloblastoma_Human/MB263",
  "Medulloblastoma_Human/GSM8840047_MB266" = "Medulloblastoma_Human/MB266",
  "Medulloblastoma_Human/GSM8840048_MB295" = "Medulloblastoma_Human/MB295",
  "Medulloblastoma_Human/GSM8840049_MB299" = "Medulloblastoma_Human/MB299"
)

# ===========================================================
# 样本注册表 (与 R1/R2 一致)
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
# 向量化 Spearman (来自 stage2_step2_analysis.R)
# ===========================================================

compute_spearman_vectorized <- function(expr_mat, density_vec) {
  expr_dense <- as.matrix(expr_mat)
  n          <- ncol(expr_dense)
  rank_d     <- rank(density_vec, ties.method = "average")
  rd_c       <- rank_d - mean(rank_d)
  ss_d       <- sum(rd_c^2)
  rank_e     <- t(apply(expr_dense, 1, rank, ties.method = "average"))
  re_c       <- rank_e - rowMeans(rank_e)
  ss_e       <- rowSums(re_c^2)
  cross      <- as.numeric(re_c %*% rd_c)
  denom      <- sqrt(ss_e * ss_d)
  denom[denom == 0] <- 1e-30
  rho        <- pmax(pmin(cross / denom, 1), -1)
  tstat      <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-10))
  pval       <- 2 * pt(-abs(tstat), df = n - 2)
  data.frame(
    gene         = rownames(expr_mat),
    spearman_cor = rho,
    p_value      = pval,
    stringsAsFactors = FALSE
  )
}

# ===========================================================
# Jaccard 相似系数
# ===========================================================

jaccard <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  return(inter / uni)
}

# ===========================================================
# 主函数: 单样本三方法比较
# ===========================================================

process_one_sample <- function(s) {

  project <- s$project
  sname   <- s$sample
  r_label <- paste0(project, "/", sname)

  sample_dir <- file.path(RESULTS_ROOT, project, sname)
  rds_path   <- file.path(sample_dir, "seurat_umap.rds")

  # P2 density 文件路径
  p2_label    <- P2_DIR_MAP[[r_label]]
  density_csv <- file.path(DENSITY_ROOT, p2_label, "cell_density_three_methods.csv")

  # --- 检查输入 ---
  if (!file.exists(rds_path)) {
    cat("  [SKIP] seurat_umap.rds not found\n")
    return(NULL)
  }
  if (!file.exists(density_csv)) {
    cat("  [SKIP] density CSV not found:", density_csv, "\n")
    return(NULL)
  }

  t0 <- Sys.time()

  # --- 读取 Seurat ---
  cat("  Loading Seurat...\n")
  so <- readRDS(rds_path)

  # --- 读取 density ---
  cat("  Loading density CSV...\n")
  dens_df <- fread(density_csv)
  dens_df$cell_id <- as.character(dens_df$cell_id)

  # --- cell_id 对齐 ---
  common <- intersect(colnames(so), dens_df$cell_id)
  if (length(common) < 100) {
    cat("  [ERROR] Only", length(common), "shared cell IDs. Trying string cleanup...\n")
    # 尝试去掉 bytes 前缀残留
    dens_df$cell_id <- gsub("^b'|'$", "", dens_df$cell_id)
    common <- intersect(colnames(so), dens_df$cell_id)
    if (length(common) < 100) {
      cat("  [FAIL] Still only", length(common), "shared cells. Skipping.\n")
      return(NULL)
    }
  }
  cat("  Matched cells:", length(common), "/", ncol(so), "\n")

  # 对齐
  so_sub   <- so[, common]
  dens_sub <- dens_df[match(common, dens_df$cell_id), ]

  # --- 获取 SCT 表达矩阵 ---
  expr_mat <- GetAssayData(so_sub, assay = "SCT", layer = "data")

  # 去除零方差基因
  gvar  <- apply(expr_mat, 1, var)
  valid <- gvar > 0
  expr_v <- expr_mat[valid, , drop = FALSE]
  cat("  Valid genes:", sum(valid), "/", nrow(expr_mat), "\n")

  # --- 三种方法各跑一次 Spearman ---
  all_results <- list()
  method_summary <- list()

  for (m_idx in seq_along(METHODS)) {
    method_col   <- METHODS[m_idx]
    method_label <- METHOD_LABELS[m_idx]

    d_vec <- dens_sub[[method_col]]

    # 处理 NaN (Voronoi/Delaunay 的边界细胞)
    valid_cells <- !is.na(d_vec)
    n_valid     <- sum(valid_cells)

    if (n_valid < 100) {
      cat("  [WARN]", method_label, ": only", n_valid, "valid cells, skipping\n")
      next
    }

    cat("  Spearman [", method_label, "] (", n_valid, " cells)...", sep = "")

    expr_sub <- expr_v[, valid_cells, drop = FALSE]
    d_sub    <- d_vec[valid_cells]

    res <- compute_spearman_vectorized(expr_sub, d_sub)
    res$FDR <- p.adjust(res$p_value, method = "BH")
    res <- res[order(-res$spearman_cor), ]

    # 保存
    write.csv(res,
              file.path(sample_dir, paste0("density_results_", method_label, ".csv")),
              row.names = FALSE)

    # 统计
    sig_genes  <- res$gene[res$FDR < FDR_THRESHOLD & abs(res$spearman_cor) > COR_THRESHOLD]
    pos_genes  <- res$gene[res$FDR < FDR_THRESHOLD & res$spearman_cor > COR_THRESHOLD]
    neg_genes  <- res$gene[res$FDR < FDR_THRESHOLD & res$spearman_cor < -COR_THRESHOLD]

    method_summary[[method_label]] <- list(
      method      = method_label,
      n_valid     = n_valid,
      n_sig       = length(sig_genes),
      n_pos       = length(pos_genes),
      n_neg       = length(neg_genes),
      median_abs_rho = median(abs(res$spearman_cor)),
      max_rho     = max(res$spearman_cor),
      min_rho     = min(res$spearman_cor),
      sig_genes   = sig_genes,
      top50       = head(res$gene, 50)
    )

    all_results[[method_label]] <- res
    cat(" ", length(sig_genes), "sig genes\n")
  }

  if (length(method_summary) < 2) {
    cat("  [WARN] Fewer than 2 methods completed, cannot compare\n")
    return(NULL)
  }

  # --- 方法间比较 ---
  pairs <- list(
    c("KNN", "Voronoi"), c("KNN", "Delaunay"), c("Voronoi", "Delaunay")
  )

  compare_rows <- list()
  for (p in pairs) {
    if (is.null(method_summary[[p[1]]]) || is.null(method_summary[[p[2]]])) next

    jac_sig  <- jaccard(method_summary[[p[1]]]$sig_genes,
                        method_summary[[p[2]]]$sig_genes)
    jac_top50 <- jaccard(method_summary[[p[1]]]$top50,
                         method_summary[[p[2]]]$top50)

    # 排名相关: 共有基因的 spearman_cor 排名一致性
    res_a <- all_results[[p[1]]]
    res_b <- all_results[[p[2]]]
    shared_genes <- intersect(res_a$gene, res_b$gene)
    if (length(shared_genes) > 10) {
      rho_a <- res_a$spearman_cor[match(shared_genes, res_a$gene)]
      rho_b <- res_b$spearman_cor[match(shared_genes, res_b$gene)]
      rank_cor <- cor(rho_a, rho_b, method = "spearman", use = "complete.obs")
    } else {
      rank_cor <- NA
    }

    compare_rows[[length(compare_rows) + 1]] <- data.frame(
      pair           = paste0(p[1], "_vs_", p[2]),
      jaccard_sig    = round(jac_sig, 4),
      jaccard_top50  = round(jac_top50, 4),
      rank_spearman  = round(rank_cor, 4),
      stringsAsFactors = FALSE
    )
  }

  compare_df <- do.call(rbind, compare_rows)
  write.csv(compare_df,
            file.path(sample_dir, "method_comparison.csv"),
            row.names = FALSE)

  # --- 方法汇总 (每样本) ---
  summ_rows <- lapply(method_summary, function(ms) {
    data.frame(
      method         = ms$method,
      n_valid_cells  = ms$n_valid,
      n_sig_genes    = ms$n_sig,
      n_pos          = ms$n_pos,
      n_neg          = ms$n_neg,
      median_abs_rho = round(ms$median_abs_rho, 5),
      max_rho        = round(ms$max_rho, 4),
      min_rho        = round(ms$min_rho, 4),
      stringsAsFactors = FALSE
    )
  })
  summ_df <- do.call(rbind, summ_rows)
  write.csv(summ_df,
            file.path(sample_dir, "method_summary.csv"),
            row.names = FALSE)

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  cat("  Done (", elapsed, " min)\n", sep = "")

  # 返回全局汇总数据
  return(list(
    sample     = r_label,
    species    = s$species,
    condition  = s$condition,
    summ       = summ_df,
    compare    = compare_df,
    elapsed    = elapsed
  ))
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R3: THREE-METHOD DENSITY COMPARISON\n")
cat("  Samples:", length(sample_list), "\n")
cat("  Methods:", paste(METHOD_LABELS, collapse = ", "), "\n")
cat("  Seurat from:", RESULTS_ROOT, "\n")
cat("  Density from:", DENSITY_ROOT, "\n")
cat("============================================================\n\n")

t_total  <- Sys.time()
all_info <- list()

for (i in seq_along(sample_list)) {
  s     <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("\n[%d/%d] %s (%s)\n", i, length(sample_list), label, s$species))

  info <- process_one_sample(s)
  if (!is.null(info)) {
    all_info[[length(all_info) + 1]] <- info
  }
}

# ===========================================================
# 全局方法选择报告
# ===========================================================

cat("\n============================================================\n")
cat("  METHOD SELECTION ANALYSIS\n")
cat("============================================================\n\n")

qc_dir <- file.path(RESULTS_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# 1. 汇总每个样本每种方法的 sig 基因数
global_summ <- list()
for (info in all_info) {
  for (r in seq_len(nrow(info$summ))) {
    global_summ[[length(global_summ) + 1]] <- data.frame(
      sample    = info$sample,
      species   = info$species,
      condition = info$condition,
      method    = info$summ$method[r],
      n_sig     = info$summ$n_sig_genes[r],
      n_pos     = info$summ$n_pos[r],
      n_neg     = info$summ$n_neg[r],
      median_abs_rho = info$summ$median_abs_rho[r],
      stringsAsFactors = FALSE
    )
  }
}
gs_df <- do.call(rbind, global_summ)

write.csv(gs_df,
          file.path(qc_dir, "R3_all_samples_method_summary.csv"),
          row.names = FALSE)

# 2. 汇总方法间比较
global_compare <- list()
for (info in all_info) {
  comp <- info$compare
  comp$sample <- info$sample
  global_compare[[length(global_compare) + 1]] <- comp
}
gc_df <- do.call(rbind, global_compare)

write.csv(gc_df,
          file.path(qc_dir, "R3_all_samples_method_comparison.csv"),
          row.names = FALSE)

# 3. 方法选择决策
cat("  Per-method average significant genes across samples:\n")
method_avg <- gs_df %>%
  group_by(method) %>%
  summarise(
    mean_n_sig     = round(mean(n_sig), 1),
    median_n_sig   = median(n_sig),
    mean_abs_rho   = round(mean(median_abs_rho), 5),
    n_samples      = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_n_sig))

print(as.data.frame(method_avg))

cat("\n  Pairwise consistency (mean Jaccard of sig genes):\n")
pair_avg <- gc_df %>%
  group_by(pair) %>%
  summarise(
    mean_jaccard_sig   = round(mean(jaccard_sig, na.rm = TRUE), 3),
    mean_jaccard_top50 = round(mean(jaccard_top50, na.rm = TRUE), 3),
    mean_rank_cor      = round(mean(rank_spearman, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(as.data.frame(pair_avg))

# 决策逻辑:
# 如果三方法高度一致 (mean Jaccard > 0.7), 选 KNN (文献标准)
# 否则选 mean_n_sig 最高的方法

knn_vor_jac <- pair_avg$mean_jaccard_sig[pair_avg$pair == "KNN_vs_Voronoi"]
knn_del_jac <- pair_avg$mean_jaccard_sig[pair_avg$pair == "KNN_vs_Delaunay"]

if (length(knn_vor_jac) > 0 && length(knn_del_jac) > 0 &&
    !is.na(knn_vor_jac) && !is.na(knn_del_jac) &&
    knn_vor_jac > 0.7 && knn_del_jac > 0.7) {
  selected_method <- "KNN"
  selection_reason <- paste0(
    "Three methods highly consistent (Jaccard KNN-Vor=", knn_vor_jac,
    ", KNN-Del=", knn_del_jac, "). KNN selected as literature standard."
  )
} else {
  selected_method <- as.character(method_avg$method[1])
  selection_reason <- paste0(
    "Methods show moderate divergence. ", selected_method,
    " selected for highest mean significant gene count (",
    method_avg$mean_n_sig[1], ")."
  )
}

cat("\n  >>> SELECTED METHOD:", selected_method, "\n")
cat("  >>> Reason:", selection_reason, "\n")

# 保存决策
decision_df <- data.frame(
  selected_method  = selected_method,
  reason           = selection_reason,
  knn_vor_jaccard  = ifelse(length(knn_vor_jac) > 0, knn_vor_jac, NA),
  knn_del_jaccard  = ifelse(length(knn_del_jac) > 0, knn_del_jac, NA),
  stringsAsFactors = FALSE
)

write.csv(decision_df,
          file.path(qc_dir, "R3_METHOD_SELECTION.csv"),
          row.names = FALSE)

# 保存完整的汇总
write.csv(as.data.frame(method_avg),
          file.path(qc_dir, "R3_method_avg_summary.csv"),
          row.names = FALSE)
write.csv(as.data.frame(pair_avg),
          file.path(qc_dir, "R3_pairwise_consistency.csv"),
          row.names = FALSE)

total_time <- round(difftime(Sys.time(), t_total, units = "mins"), 1)

cat("\n============================================================\n")
cat("  R3 COMPLETE\n")
cat("  Selected method:", selected_method, "\n")
cat("  Total time:", total_time, "min\n")
cat("  QC dir:", qc_dir, "\n")
cat("============================================================\n")
