#!/usr/bin/env Rscript
# ============================================================================
# R05_visualization.R
# 功能：生成 density-gene 关联的论文级图
# 输入：R2的rds + R4的filtered_density_genes.csv
# 输出：散点图（hex+box）、跨样本热力图、棒棒糖图、代表性空间图
# Run (EN): Rscript R05_visualization.R
#   Purpose: publication-grade figures of density-gene associations.
#   Paths configured via env vars (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================================

options(future.globals.maxSize = Inf)  # 无内存限制

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(jsonlite); library(ggplot2); library(viridis)
  library(pheatmap); library(ggrepel)
})

# Configurable roots (see config/paths.R); override via environment variables.
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")

REGISTRY_PATH <- file.path(DATA_DIR, "sample_registry.json")
R2_DIR        <- file.path(RESULTS_DIR, "R2_Results")
R4_DIR        <- file.path(RESULTS_DIR, "R4_Results")
OUTPUT_DIR    <- file.path(RESULTS_DIR, "R5_Results")
DENSITY_COL   <- "density_knn_main_piecewise"

# 代表性空间图：(sample_name, gene_name)
SPATIAL_TARGETS <- list(
  c("Alzheimer_Mouse/TgCRND8_17_9", "Cst3"),
  c("ATRT_Human/GSM8672834",        "OTX2"),
  c("Glioblastoma_Human/Single_Sample", "CCND1"),
  c("Brain_Human_Preview/GBM",      "ENC1"),
  c("Medulloblastoma_Human/GSM8840047", "TUBB4A"),
  c("Brain_Mouse/Single_Sample",    "Ryr2")
)

# 创建输出子目录
for (sub in c("scatter_plots", "heatmaps", "lollipop", "spatial")) {
  dir.create(file.path(OUTPUT_DIR, sub), recursive=TRUE, showWarnings=FALSE)
}

cat("============================================\n")
cat("R5: 可视化\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)

# ============================================================================
# 辅助函数
# ============================================================================

# 从rds提取单个基因的残差+密度+坐标
extract_gene_data <- function(seurat_obj, gene_name) {
  sct_data <- GetAssayData(seurat_obj, assay="SCT", layer="data")
  if (!(gene_name %in% rownames(sct_data))) return(NULL)
  data.frame(
    residual = as.numeric(sct_data[gene_name, ]),
    density  = seurat_obj@meta.data[[DENSITY_COL]],
    x        = seurat_obj$x_centroid,
    y        = seurat_obj$y_centroid
  )
}

# 画 hex + box 并排散点图
plot_scatter_hex_box <- function(df, gene, sample_name, out_path) {
  df <- df[!is.na(df$density) & !is.na(df$residual), ]
  if (nrow(df) < 100) return(FALSE)

  # hex bin
  p1 <- ggplot(df, aes(x=density, y=residual)) +
    geom_hex(bins=60) + scale_fill_viridis(option="plasma") +
    labs(title=paste0(sample_name, " — ", gene),
         x="KNN density (main)", y="SCT residual") +
    theme_minimal() + theme(plot.title=element_text(size=9))

  # 分箱箱线图
  df$density_bin <- cut(df$density, breaks=quantile(df$density, probs=seq(0,1,0.1),
                        na.rm=TRUE), include.lowest=TRUE, labels=1:10)
  df2 <- df[!is.na(df$density_bin), ]
  p2 <- ggplot(df2, aes(x=density_bin, y=residual)) +
    geom_boxplot(fill="steelblue", alpha=0.6, outlier.size=0.3) +
    labs(title="Residual by density decile",
         x="Density decile (1=low, 10=high)", y="SCT residual") +
    theme_minimal() + theme(plot.title=element_text(size=9))

  png(out_path, width=1400, height=500, res=140)
  tryCatch({
    if (requireNamespace("gridExtra", quietly=TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol=2)
    } else { print(p1) }
  }, finally = dev.off())
  TRUE
}

# 棒棒糖图
plot_lollipop <- function(r4_df, sample_name, out_path) {
  tier1 <- r4_df[r4_df$tier == "tier1_strict", ]
  if (nrow(tier1) == 0) return(FALSE)
  tier1 <- tier1[order(tier1$rho_knn_main), ]
  tier1$gene <- factor(tier1$gene, levels=tier1$gene)
  tier1$direction <- ifelse(tier1$rho_knn_main > 0, "positive", "negative")

  p <- ggplot(tier1, aes(x=gene, y=rho_knn_main, color=direction)) +
    geom_segment(aes(xend=gene, yend=0)) + geom_point(size=2) +
    scale_color_manual(values=c("positive"="#D62728","negative"="#1F77B4")) +
    coord_flip() + geom_hline(yintercept=0, color="grey50") +
    labs(title=paste0(sample_name, " — Tier1 Density Genes (n=", nrow(tier1), ")"),
         x=NULL, y="ρ (Spearman)") +
    theme_minimal() + theme(plot.title=element_text(size=9),
                            axis.text.y=element_text(size=max(4, 10-nrow(tier1)/20)))

  h <- max(5, min(20, nrow(tier1) * 0.15))
  ggsave(out_path, p, width=8, height=h, dpi=130, limitsize=FALSE)
  TRUE
}

# 空间图：并排展示密度和基因表达
plot_spatial <- function(df, gene, sample_name, out_path) {
  df <- df[!is.na(df$density) & !is.na(df$residual), ]
  if (nrow(df) < 100) return(FALSE)
  pt_size <- ifelse(nrow(df) > 200000, 0.05, ifelse(nrow(df) > 50000, 0.1, 0.3))

  p1 <- ggplot(df, aes(x=x, y=y, color=density)) +
    geom_point(size=pt_size) + scale_color_viridis(option="viridis") +
    coord_fixed() + labs(title=paste0(sample_name, " — Density"),
                         x="x (µm)", y="y (µm)") +
    theme_minimal() + theme(plot.title=element_text(size=9))

  # 残差用对称色标
  lim <- quantile(abs(df$residual), 0.98, na.rm=TRUE)
  p2 <- ggplot(df, aes(x=x, y=y, color=residual)) +
    geom_point(size=pt_size) +
    scale_color_gradient2(low="#1F77B4", mid="white", high="#D62728",
                          midpoint=0, limits=c(-lim, lim), oob=scales::squish) +
    coord_fixed() + labs(title=paste0(gene, " — SCT residual"),
                         x="x (µm)", y="y (µm)") +
    theme_minimal() + theme(plot.title=element_text(size=9))

  png(out_path, width=1600, height=700, res=140)
  tryCatch({
    if (requireNamespace("gridExtra", quietly=TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol=2)
    } else { print(p1) }
  }, finally = dev.off())
  TRUE
}

# ============================================================================
# 第一阶段：逐样本处理（散点图 + 棒棒糖图）
# ============================================================================

cat("----- 阶段1: 散点图 + 棒棒糖图 -----\n\n")

for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  parts <- strsplit(sample_name, "/")[[1]]
  dataset_name <- parts[1]; sample_subname <- parts[2]
  cat(sprintf("[%d/22] %s\n", i, sample_name))

  r4_path <- file.path(R4_DIR, dataset_name, sample_subname, "filtered_density_genes.csv")
  rds_path <- file.path(R2_DIR, dataset_name, sample_subname, paste0(sample_subname, "_seurat_R2.rds"))
  if (!file.exists(r4_path) || !file.exists(rds_path)) { cat("  ✗ 缺文件\n\n"); next }

  r4_df <- fread(r4_path)
  tier1_top3 <- head(r4_df[r4_df$tier == "tier1_strict", ], 3)

  # 棒棒糖图
  lolli_path <- file.path(OUTPUT_DIR, "lollipop",
                          paste0(dataset_name, "__", sample_subname, "_lollipop.png"))
  tryCatch({
    if (plot_lollipop(r4_df, sample_name, lolli_path)) cat("  ✓ 棒棒糖图\n")
    else cat("  ⚠ 无tier1基因, 跳过棒棒糖\n")
  }, error=function(e) cat(sprintf("  ⚠ 棒棒糖失败: %s\n", e$message)))

  # 散点图（需读rds）
  if (nrow(tier1_top3) > 0) {
    cat(sprintf("  读取 rds (top3: %s) ...\n", paste(tier1_top3$gene, collapse=",")))
    t0 <- Sys.time()
    seurat_obj <- readRDS(rds_path)
    cat(sprintf("    rds加载耗时 %.1f秒\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))

    for (g in tier1_top3$gene) {
      df <- extract_gene_data(seurat_obj, g)
      if (is.null(df)) { cat(sprintf("    ⚠ %s 不在对象中\n", g)); next }
      out_path <- file.path(OUTPUT_DIR, "scatter_plots",
                            paste0(dataset_name, "__", sample_subname, "__", g, "_scatter.png"))
      tryCatch({
        if (plot_scatter_hex_box(df, g, sample_name, out_path))
          cat(sprintf("    ✓ 散点图 %s\n", g))
      }, error=function(e) cat(sprintf("    ⚠ %s 散点图失败: %s\n", g, e$message)))
    }
    rm(seurat_obj); gc(verbose=FALSE)
  }
  cat("\n")
}

# ============================================================================
# 第二阶段：跨样本热力图（多样本数据集）
# ============================================================================

cat("----- 阶段2: 跨样本热力图 -----\n\n")

multi_datasets <- list(
  Alzheimer_Mouse     = grep("^Alzheimer_Mouse/", sample_names, value=TRUE),
  ATRT_Human          = grep("^ATRT_Human/", sample_names, value=TRUE),
  Brain_Human_Preview = grep("^Brain_Human_Preview/", sample_names, value=TRUE),
  Medulloblastoma_Human = grep("^Medulloblastoma_Human/", sample_names, value=TRUE)
)

for (ds_name in names(multi_datasets)) {
  samples <- multi_datasets[[ds_name]]
  cat(sprintf("[%s] %d 样本\n", ds_name, length(samples)))

  # 收集所有样本的tier1基因并集
  tier1_union <- c()
  sample_dfs <- list()
  for (s in samples) {
    parts <- strsplit(s, "/")[[1]]
    r4_path <- file.path(R4_DIR, parts[1], parts[2], "filtered_density_genes.csv")
    if (!file.exists(r4_path)) next
    df <- fread(r4_path)
    sample_dfs[[s]] <- df
    tier1_union <- union(tier1_union, df[df$tier == "tier1_strict", ]$gene)
  }

  if (length(tier1_union) == 0) { cat("  ⚠ 无tier1基因\n\n"); next }

  # 构建基因×样本矩阵
  mat <- matrix(NA_real_, nrow=length(tier1_union), ncol=length(samples),
                dimnames=list(tier1_union, basename(samples)))
  for (s in samples) {
    if (is.null(sample_dfs[[s]])) next
    df <- sample_dfs[[s]]
    matched <- match(tier1_union, df$gene)
    mat[, basename(s)] <- df$rho_knn_main[matched]
  }

  # 按行均值排序（|ρ|均值高的在上）
  row_order <- order(-rowMeans(abs(mat), na.rm=TRUE))
  mat <- mat[row_order, , drop=FALSE]

  # 限制行数（最多50行）
  if (nrow(mat) > 50) mat <- mat[1:50, , drop=FALSE]

  out_path <- file.path(OUTPUT_DIR, "heatmaps", paste0(ds_name, "_heatmap.png"))
  tryCatch({
    png(out_path, width=1200, height=max(800, nrow(mat)*25), res=130)
    pheatmap(mat, cluster_rows=FALSE, cluster_cols=FALSE,
             color=colorRampPalette(c("#1F77B4","white","#D62728"))(100),
             breaks=seq(-0.5, 0.5, length.out=101),
             main=paste0(ds_name, " — Tier1 genes (top ", nrow(mat), ")"),
             fontsize_row=7, fontsize_col=8, na_col="grey90")
    dev.off()
    cat(sprintf("  ✓ 热力图 (%d genes × %d samples)\n\n", nrow(mat), ncol(mat)))
  }, error=function(e) { try(dev.off(), silent=TRUE); cat(sprintf("  ⚠ 失败: %s\n\n", e$message)) })
}

# ============================================================================
# 第三阶段：代表性空间图
# ============================================================================

cat("----- 阶段3: 代表性空间图 -----\n\n")

for (k in seq_along(SPATIAL_TARGETS)) {
  target <- SPATIAL_TARGETS[[k]]
  sample_name <- target[1]; gene <- target[2]
  cat(sprintf("[%d/%d] %s — %s\n", k, length(SPATIAL_TARGETS), sample_name, gene))

  parts <- strsplit(sample_name, "/")[[1]]
  rds_path <- file.path(R2_DIR, parts[1], parts[2], paste0(parts[2], "_seurat_R2.rds"))
  if (!file.exists(rds_path)) { cat("  ✗ rds缺失\n\n"); next }

  t0 <- Sys.time()
  seurat_obj <- readRDS(rds_path)
  cat(sprintf("  rds加载 %.1f秒\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))

  df <- extract_gene_data(seurat_obj, gene)
  if (is.null(df)) { cat(sprintf("  ⚠ %s 不在对象中\n\n", gene)); rm(seurat_obj); gc(verbose=FALSE); next }

  out_path <- file.path(OUTPUT_DIR, "spatial",
                        paste0(parts[1], "__", parts[2], "__", gene, "_spatial.png"))
  tryCatch({
    if (plot_spatial(df, gene, sample_name, out_path)) cat("  ✓ 空间图\n")
  }, error=function(e) cat(sprintf("  ⚠ 失败: %s\n", e$message)))

  rm(seurat_obj, df); gc(verbose=FALSE)
  cat("\n")
}

cat("============================================\n")
cat("R5 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
