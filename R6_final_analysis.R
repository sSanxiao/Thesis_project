############################################################
# R6_final_analysis.R
# 功能: 基于 R3-R5 确定的参数, 执行最终分析:
#   1. 挂载 KNN density + tier 分类到 Seurat 对象
#   2. 计算 Density Signature (方向感知, |ρ| 加权)
#   3. Cell State Coupling (Task 2): 细胞周期, 炎症,
#      PC1/PC2, cluster identity 与 density 的相关性
#   4. 导出完整 cell metadata
#   5. 保存 seurat_final.rds
#
# 输入: seurat_umap.rds (R2)
#       cell_density_three_methods.csv (P2)
#       density_results_KNN.csv (R3)
#       gene_tier_classification.csv (R5)
# 输出: seurat_final.rds, cell_metadata_full.csv,
#       density_signature_vector.csv, density_state_coupling.csv
#
# 运行:
#   cd /home/disk/wangqilu/Stage2_new/Scripts/
#   nohup Rscript R6_final_analysis.R > R6_run.log 2>&1 &
#   tail -f R6_run.log
############################################################

library(Seurat)
library(data.table)
library(dplyr)

# ===========================================================
# 配置 (所有参数均经过 R3-R5 验证)
# ===========================================================

RESULTS_ROOT <- "/home/disk/wangqilu/Stage2_new/Results"
DENSITY_ROOT <- "/home/disk/wangqilu/Density_Caculation/Results"

# 确定的参数
SELECTED_METHOD <- "KNN"
FDR_THRESHOLD   <- 0.05
COR_THRESHOLD   <- 0.05
MIN_SIG_GENES   <- 5

# P2 路径映射
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
# 样本注册表
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
# Cell cycle gene lists (Tirosh et al. 2016, uppercase)
# ===========================================================

s_genes <- c(
  "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2",
  "MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2",
  "RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7",
  "POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
  "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B",
  "BRIP1","E2F8"
)

g2m_genes <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
  "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A",
  "SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E",
  "TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK",
  "CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8",
  "ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5",
  "CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
)

# 炎症 marker (保守的, 跨物种)
inflammation_markers <- c("PTPRC","AIF1","CD68","GFAP","TYROBP","CX3CR1")

# ===========================================================
# 主函数
# ===========================================================

process_one_sample <- function(s) {

  project <- s$project
  sname   <- s$sample
  species <- s$species
  cond    <- s$condition
  r_label <- paste0(project, "/", sname)
  sample_dir <- file.path(RESULTS_ROOT, project, sname)

  rds_path    <- file.path(sample_dir, "seurat_umap.rds")
  p2_label    <- P2_DIR_MAP[[r_label]]
  density_csv <- file.path(DENSITY_ROOT, p2_label, "cell_density_three_methods.csv")
  global_csv  <- file.path(sample_dir, "density_results_KNN.csv")
  tier_csv    <- file.path(sample_dir, "gene_tier_classification.csv")

  if (!file.exists(rds_path)) { cat("  [SKIP] no seurat_umap.rds\n"); return(NULL) }
  if (!file.exists(density_csv)) { cat("  [SKIP] no density CSV\n"); return(NULL) }

  t0 <- Sys.time()

  # --- 加载 ---
  cat("  Loading Seurat + density...\n")
  so <- readRDS(rds_path)

  dens_df <- fread(density_csv)
  dens_df$cell_id <- as.character(dens_df$cell_id)
  if (nrow(dens_df) > 0 && grepl("^b'", dens_df$cell_id[1])) {
    dens_df$cell_id <- gsub("^b'|'$", "", dens_df$cell_id)
  }

  common <- intersect(colnames(so), dens_df$cell_id)
  if (length(common) < 100) { cat("  [FAIL] <100 shared cells\n"); return(NULL) }

  so <- so[, common]
  dens_sub <- dens_df[match(common, dens_df$cell_id), ]

  # =========================================================
  # 1. 挂载 density
  # =========================================================
  cat("  Mounting density values...\n")
  so$density_knn      <- dens_sub$density_knn
  so$density_voronoi  <- dens_sub$density_voronoi
  so$density_delaunay <- dens_sub$density_delaunay

  # =========================================================
  # 2. Density Signature (方向感知, |ρ| 加权)
  # =========================================================
  cat("  Computing density signature...\n")
  expr_mat <- GetAssayData(so, assay = "SCT", layer = "data")

  if (file.exists(global_csv)) {
    res <- read.csv(global_csv, stringsAsFactors = FALSE)
  } else {
    cat("    [WARN] No R3 results, skipping signature\n")
    so$Density_Signature <- NA_real_
    res <- NULL
  }

  if (!is.null(res)) {
    pos_g <- res$gene[res$FDR < FDR_THRESHOLD & res$spearman_cor > COR_THRESHOLD]
    neg_g <- res$gene[res$FDR < FDR_THRESHOLD & res$spearman_cor < -COR_THRESHOLD]

    pos_g <- intersect(pos_g, rownames(expr_mat))
    neg_g <- intersect(neg_g, rownames(expr_mat))

    if (length(pos_g) >= MIN_SIG_GENES && length(neg_g) >= MIN_SIG_GENES) {
      # |ρ| 加权 signature: 正向基因加权均值 - 负向基因加权均值
      pos_weights <- abs(res$spearman_cor[match(pos_g, res$gene)])
      neg_weights <- abs(res$spearman_cor[match(neg_g, res$gene)])

      pos_weights <- pos_weights / sum(pos_weights)
      neg_weights <- neg_weights / sum(neg_weights)

      pos_score <- as.numeric(pos_weights %*% as.matrix(expr_mat[pos_g, ]))
      neg_score <- as.numeric(neg_weights %*% as.matrix(expr_mat[neg_g, ]))

      dsig <- pos_score - neg_score

    } else if (length(c(pos_g, neg_g)) >= MIN_SIG_GENES) {
      all_sig <- c(pos_g, neg_g)
      signs   <- ifelse(res$spearman_cor[match(all_sig, res$gene)] > 0, 1, -1)
      weights <- abs(res$spearman_cor[match(all_sig, res$gene)])
      weights <- weights / sum(weights)

      dsig <- as.numeric((weights * signs) %*% as.matrix(expr_mat[all_sig, ]))
    } else {
      cat("    [WARN] <", MIN_SIG_GENES, "sig genes, signature = NA\n")
      dsig <- rep(NA_real_, ncol(so))
    }

    so$Density_Signature <- dsig

    write.csv(
      data.frame(
        cell              = colnames(so),
        density_signature = dsig,
        local_density     = so$density_knn
      ),
      file.path(sample_dir, "density_signature_vector.csv"),
      row.names = FALSE
    )

    cat("    Pos genes:", length(pos_g), " Neg genes:", length(neg_g), "\n")
  }

  # =========================================================
  # 3. Cell State Coupling (Task 2)
  # =========================================================
  cat("  Cell state coupling...\n")
  coupling <- list()
  ld <- so$density_knn

  # 3a. Cell cycle scoring
  # 转为大写匹配 (跨物种兼容)
  gene_panel <- toupper(rownames(so))
  names(gene_panel) <- rownames(so)

  s_avail   <- intersect(s_genes, gene_panel)
  g2m_avail <- intersect(g2m_genes, gene_panel)

  # 映射回原始基因名
  s_orig   <- names(gene_panel)[gene_panel %in% s_avail]
  g2m_orig <- names(gene_panel)[gene_panel %in% g2m_avail]

  if (length(s_orig) >= 3 && length(g2m_orig) >= 3) {
    cc_ok <- tryCatch({
      so <- CellCycleScoring(so, s.features = s_orig, g2m.features = g2m_orig,
                             set.ident = FALSE)
      TRUE
    }, error = function(e) {
      cat("    [WARN] CellCycleScoring failed:", conditionMessage(e), "\n")
      FALSE
    })
    if (cc_ok && "S.Score" %in% colnames(so@meta.data)) {
      coupling$rho_S_Score   <- cor(ld, so$S.Score, method = "spearman", use = "complete.obs")
      coupling$rho_G2M_Score <- cor(ld, so$G2M.Score, method = "spearman", use = "complete.obs")
      cat("    Cell cycle: S.Score rho=", round(coupling$rho_S_Score, 3),
          " G2M.Score rho=", round(coupling$rho_G2M_Score, 3), "\n")
    } else {
      coupling$rho_S_Score   <- NA
      coupling$rho_G2M_Score <- NA
    }
  } else {
    coupling$rho_S_Score   <- NA
    coupling$rho_G2M_Score <- NA
    cat("    [WARN] Not enough cell cycle genes in panel (",
        length(s_orig), "S,", length(g2m_orig), "G2M)\n")
  }

  # 3b. Inflammation score
  infl_panel <- toupper(rownames(so))
  infl_avail <- inflammation_markers[inflammation_markers %in% infl_panel]
  infl_orig  <- names(gene_panel)[gene_panel %in% infl_avail]

  if (length(infl_orig) >= 2) {
    # ctrl must be < total genes in panel; use min(ctrl, ngenes/2)
    n_panel <- nrow(so)
    ctrl_size <- min(50, max(5, floor(n_panel / 5)))

    infl_ok <- tryCatch({
      so <- AddModuleScore(so, features = list(infl_orig), name = "Inflammation",
                           ctrl = ctrl_size, verbose = FALSE)
      TRUE
    }, error = function(e) {
      cat("    [WARN] AddModuleScore failed:", conditionMessage(e), "\n")
      FALSE
    })

    if (infl_ok && "Inflammation1" %in% colnames(so@meta.data)) {
      coupling$rho_inflammation <- cor(ld, so$Inflammation1, method = "spearman", use = "complete.obs")
      cat("    Inflammation rho=", round(coupling$rho_inflammation, 3),
          " (", length(infl_orig), "markers, ctrl=", ctrl_size, ")\n")
    } else {
      coupling$rho_inflammation <- NA
    }
  } else {
    coupling$rho_inflammation <- NA
    cat("    [WARN] <2 inflammation markers in panel (",
        length(infl_orig), "found)\n")
  }

  # 3c. PC1, PC2
  pca_emb <- Embeddings(so, reduction = "pca")
  coupling$rho_PC1 <- cor(ld, pca_emb[, 1], method = "spearman", use = "complete.obs")
  coupling$rho_PC2 <- cor(ld, pca_emb[, 2], method = "spearman", use = "complete.obs")
  cat("    PC1 rho=", round(coupling$rho_PC1, 3),
      " PC2 rho=", round(coupling$rho_PC2, 3), "\n")

  # 3d. Cluster mean density
  clusters <- as.character(Idents(so))
  cl_means <- tapply(ld, clusters, mean)
  cell_cl_mean <- cl_means[clusters]
  coupling$rho_cluster_mean <- cor(ld, cell_cl_mean, method = "spearman", use = "complete.obs")
  cat("    Cluster mean density rho=", round(coupling$rho_cluster_mean, 3), "\n")

  # 3e. Density signature (if computed)
  if (!all(is.na(so$Density_Signature))) {
    coupling$rho_signature <- cor(ld, so$Density_Signature, method = "spearman", use = "complete.obs")
  } else {
    coupling$rho_signature <- NA
  }

  # 保存 coupling
  coupling_df <- data.frame(
    sample    = r_label,
    species   = species,
    condition = cond,
    as.data.frame(coupling),
    stringsAsFactors = FALSE
  )
  write.csv(coupling_df,
            file.path(sample_dir, "density_state_coupling.csv"),
            row.names = FALSE)

  # =========================================================
  # 4. 导出完整 cell metadata
  # =========================================================
  cat("  Exporting cell metadata...\n")
  cmf <- data.frame(
    cell_id           = colnames(so),
    project           = project,
    sample            = sname,
    species           = species,
    condition         = cond,
    x_centroid        = so$x_centroid,
    y_centroid        = so$y_centroid,
    cluster           = as.character(Idents(so)),
    density_knn       = so$density_knn,
    density_voronoi   = so$density_voronoi,
    density_delaunay  = so$density_delaunay,
    density_signature = so$Density_Signature,
    PC1               = pca_emb[, 1],
    PC2               = pca_emb[, 2],
    nCount_RNA        = so$nCount_RNA,
    nFeature_RNA      = so$nFeature_RNA,
    stringsAsFactors  = FALSE
  )

  # 加入 cell cycle 和 inflammation (如果有)
  if ("S.Score" %in% colnames(so@meta.data)) {
    cmf$S_Score   <- so$S.Score
    cmf$G2M_Score <- so$G2M.Score
    cmf$Phase     <- so$Phase
  }
  if ("Inflammation1" %in% colnames(so@meta.data)) {
    cmf$Inflammation <- so$Inflammation1
  }

  write.csv(cmf,
            file.path(sample_dir, "cell_metadata_full.csv"),
            row.names = FALSE)

  # =========================================================
  # 5. 保存最终 Seurat
  # =========================================================
  cat("  Saving seurat_final.rds...\n")
  saveRDS(so, file.path(sample_dir, "seurat_final.rds"))

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  cat("  Done (", elapsed, " min)\n", sep = "")

  return(coupling_df)
}

# ===========================================================
# 批量运行
# ===========================================================

cat("============================================================\n")
cat("  R6: FINAL ANALYSIS\n")
cat("  Method:", SELECTED_METHOD, "\n")
cat("  Threshold: FDR<", FDR_THRESHOLD, " & |ρ|>", COR_THRESHOLD, "\n")
cat("  Signature: direction-aware, |ρ|-weighted\n")
cat("  Task 2: cell cycle + inflammation + PC + cluster coupling\n")
cat("  Samples:", length(sample_list), "\n")
cat("============================================================\n\n")

t_total <- Sys.time()
all_coupling <- list()

for (i in seq_along(sample_list)) {
  s     <- sample_list[[i]]
  label <- paste0(s$project, "/", s$sample)
  cat(sprintf("\n[%d/%d] %s (%s, %s)\n", i, length(sample_list), label, s$species, s$condition))

  cp <- process_one_sample(s)
  if (!is.null(cp)) {
    all_coupling[[length(all_coupling) + 1]] <- cp
  }
}

# ===========================================================
# 汇总 Coupling
# ===========================================================

qc_dir <- file.path(RESULTS_ROOT, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

coupling_all <- do.call(rbind, all_coupling)
write.csv(coupling_all,
          file.path(qc_dir, "R6_all_coupling.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("  DENSITY-STATE COUPLING SUMMARY\n")
cat("============================================================\n\n")

# 按 condition 汇总
coupling_cols <- c("rho_S_Score", "rho_G2M_Score", "rho_inflammation",
                   "rho_PC1", "rho_PC2", "rho_cluster_mean", "rho_signature")

cat("  Mean coupling by condition:\n")
cond_coupling <- coupling_all %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    across(all_of(coupling_cols), ~ round(mean(.x, na.rm = TRUE), 3)),
    .groups = "drop"
  )
print(as.data.frame(cond_coupling))

write.csv(as.data.frame(cond_coupling),
          file.path(qc_dir, "R6_coupling_by_condition.csv"),
          row.names = FALSE)

# 全样本均值
cat("\n  Overall mean coupling:\n")
for (col in coupling_cols) {
  vals <- coupling_all[[col]]
  cat("    ", col, ":", round(mean(vals, na.rm = TRUE), 3),
      " (range:", round(min(vals, na.rm = TRUE), 3), "to",
      round(max(vals, na.rm = TRUE), 3), ")\n")
}

total_time <- round(difftime(Sys.time(), t_total, units = "mins"), 1)

cat("\n============================================================\n")
cat("  R6 COMPLETE\n")
cat("  Total time:", total_time, "min\n")
cat("  Output: seurat_final.rds + cell_metadata_full.csv per sample\n")
cat("  Coupling summary:", file.path(qc_dir, "R6_all_coupling.csv"), "\n")
cat("============================================================\n")
