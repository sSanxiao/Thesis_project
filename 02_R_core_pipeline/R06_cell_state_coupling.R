#!/usr/bin/env Rscript
# ============================================================================
# R06_cell_state_coupling.R
# 功能：拆解 density-gene 关联的来源（composition vs regulation effect）
# 输入：R2 的 .rds + R4 的 filtered_density_genes.csv
# 输出：分层分析结果 + 细胞周期关联 + 诊断图
# Run (EN): Rscript R06_cell_state_coupling.R
#   Purpose: decompose density-gene association into composition vs regulation effects.
#   Paths configured via env vars (DATA_DIR/RESULTS_DIR); see config/paths.R.
# ============================================================================

options(future.globals.maxSize = Inf)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(jsonlite); library(ggplot2); library(viridis)
})

# Configurable roots (see config/paths.R); override via environment variables.
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")

REGISTRY_PATH <- file.path(DATA_DIR, "sample_registry.json")
R2_DIR        <- file.path(RESULTS_DIR, "R2_Results")
R4_DIR        <- file.path(RESULTS_DIR, "R4_Results")
OUTPUT_DIR    <- file.path(RESULTS_DIR, "R6_Results")
DENSITY_COL   <- "density_knn_main_piecewise"
MIN_CELLS_PER_CLUSTER <- 100
COMPOSITION_THRESHOLD <- 0.5  # cluster内|ρ| < 全局|ρ|×0.5 → composition_driven
MIN_CC_GENES <- 5             # 细胞周期基因匹配数<5则跳过

cat("============================================\n")
cat("R6: Cell State Coupling 分析\n")
cat("============================================\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("参数: MIN_CELLS=%d, COMP_THRESHOLD=%.1f, MIN_CC_GENES=%d\n\n",
            MIN_CELLS_PER_CLUSTER, COMPOSITION_THRESHOLD, MIN_CC_GENES))

registry <- fromJSON(REGISTRY_PATH)
sample_names <- names(registry)
n_samples <- length(sample_names)
cat("样本数:", n_samples, "\n\n")

# Tirosh细胞周期基因（Seurat内置）
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

summary_list <- list()

for (i in seq_along(sample_names)) {

  sample_name <- sample_names[i]
  parts <- strsplit(sample_name, "/")[[1]]
  dataset_name <- parts[1]; sample_subname <- parts[2]

  cat("========================================\n")
  cat(sprintf("[%d/%d] %s\n", i, n_samples, sample_name))
  cat("========================================\n")

  t_start <- Sys.time()

  rds_path <- file.path(R2_DIR, dataset_name, sample_subname,
                        paste0(sample_subname, "_seurat_R2.rds"))
  r4_path <- file.path(R4_DIR, dataset_name, sample_subname, "filtered_density_genes.csv")

  if (!file.exists(rds_path) || !file.exists(r4_path)) {
    cat("  ✗ 缺文件, 跳过\n\n"); next
  }

  out_dir <- file.path(OUTPUT_DIR, dataset_name, sample_subname)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------
  # 第一步：读取数据
  # ---------------------------------------------------------------

  cat("  [1/6] 读取数据 ...\n")
  seurat_obj <- readRDS(rds_path)
  r4_df <- fread(r4_path)
  tier1_genes <- r4_df[r4_df$tier == "tier1_strict", ]$gene
  n_tier1 <- length(tier1_genes)
  n_cells <- ncol(seurat_obj)
  clusters <- levels(Idents(seurat_obj))
  n_clusters <- length(clusters)
  cat(sprintf("        %d 细胞, %d clusters, %d tier1 基因\n", n_cells, n_clusters, n_tier1))

  # 提取密度和残差
  density_vec <- seurat_obj@meta.data[[DENSITY_COL]]
  cluster_vec <- as.character(Idents(seurat_obj))
  sct_data <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")

  # ---------------------------------------------------------------
  # 第二步：每个cluster的密度分布
  # ---------------------------------------------------------------

  cat("  [2/6] Cluster 密度分布 ...\n")

  cluster_profile <- data.frame(
    cluster = character(), n_cells = integer(),
    density_median = numeric(), density_mean = numeric(),
    density_q25 = numeric(), density_q75 = numeric(),
    density_category = character(), stringsAsFactors = FALSE
  )

  global_q25 <- quantile(density_vec, 0.25, na.rm = TRUE)
  global_q75 <- quantile(density_vec, 0.75, na.rm = TRUE)

  for (cl in clusters) {
    mask <- cluster_vec == cl
    d <- density_vec[mask]
    d <- d[!is.na(d)]
    if (length(d) == 0) next
    med <- median(d)
    cat_label <- ifelse(med > global_q75, "high_density",
                   ifelse(med < global_q25, "low_density", "medium_density"))
    cluster_profile <- rbind(cluster_profile, data.frame(
      cluster = cl, n_cells = sum(mask),
      density_median = round(med, 6), density_mean = round(mean(d), 6),
      density_q25 = round(quantile(d, 0.25), 6), density_q75 = round(quantile(d, 0.75), 6),
      density_category = cat_label, stringsAsFactors = FALSE
    ))
  }

  fwrite(cluster_profile, file.path(out_dir, "cluster_density_profile.csv"))

  n_high <- sum(cluster_profile$density_category == "high_density")
  n_low <- sum(cluster_profile$density_category == "low_density")
  cat(sprintf("        高密度cluster=%d, 低密度=%d, 中等=%d\n", n_high, n_low,
              n_clusters - n_high - n_low))

  # 箱线图
  tryCatch({
    box_df <- data.frame(cluster = cluster_vec, density = density_vec)
    box_df <- box_df[!is.na(box_df$density), ]
    # 按中位密度排序
    cl_order <- cluster_profile$cluster[order(cluster_profile$density_median)]
    box_df$cluster <- factor(box_df$cluster, levels = cl_order)

    p_box <- ggplot(box_df, aes(x = cluster, y = density, fill = cluster)) +
      geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
      labs(title = paste0(sample_name, " — Density by Cluster"),
           x = "Cluster (sorted by median density)", y = "KNN density (main)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10), legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    ggsave(file.path(out_dir, "cluster_density_boxplot.png"), p_box,
           width = max(8, n_clusters * 0.3), height = 5, dpi = 130, limitsize = FALSE)
    cat("        ✓ 密度箱线图\n")
  }, error = function(e) cat(sprintf("        ⚠ 箱线图失败: %s\n", e$message)))

  # ---------------------------------------------------------------
  # 第三步：cluster内Spearman分层分析
  # ---------------------------------------------------------------

  cat("  [3/6] 分层 Spearman 分析 ...\n")

  # 哪些cluster有足够细胞
  valid_clusters <- cluster_profile$cluster[cluster_profile$n_cells >= MIN_CELLS_PER_CLUSTER]
  cat(sprintf("        %d/%d clusters 有 ≥%d 细胞可用于分层分析\n",
              length(valid_clusters), n_clusters, MIN_CELLS_PER_CLUSTER))

  coupling_list <- list()

  if (n_tier1 > 0 && length(valid_clusters) > 0) {

    for (g in tier1_genes) {
      if (!(g %in% rownames(sct_data))) next

      gene_resid <- as.numeric(sct_data[g, ])

      # 全局ρ（从R4结果中取）
      r4_row <- r4_df[r4_df$gene == g, ]
      global_rho <- r4_row$rho_knn_main[1]

      # 各cluster内ρ
      within_rhos <- numeric(length(valid_clusters))
      within_ps <- numeric(length(valid_clusters))
      names(within_rhos) <- valid_clusters
      names(within_ps) <- valid_clusters

      for (cl in valid_clusters) {
        mask <- cluster_vec == cl
        d_sub <- density_vec[mask]
        g_sub <- gene_resid[mask]
        valid_sub <- !is.na(d_sub) & !is.na(g_sub)
        if (sum(valid_sub) < 20) {
          within_rhos[cl] <- NA; within_ps[cl] <- NA; next
        }
        sp <- suppressWarnings(cor.test(g_sub[valid_sub], d_sub[valid_sub], method = "spearman"))
        within_rhos[cl] <- sp$estimate
        within_ps[cl] <- sp$p.value
      }

      # 判定：composition vs regulation
      valid_within <- within_rhos[!is.na(within_rhos)]
      if (length(valid_within) == 0) {
        effect_label <- "insufficient_data"
      } else {
        all_below <- all(abs(valid_within) < abs(global_rho) * COMPOSITION_THRESHOLD)
        any_opposite <- any(sign(valid_within) != sign(global_rho) & abs(valid_within) > 0.03)
        any_strong <- any(abs(valid_within) >= abs(global_rho) * COMPOSITION_THRESHOLD)

        if (all_below) {
          effect_label <- "composition_driven"
        } else if (any_opposite) {
          effect_label <- "cluster_heterogeneous"
        } else if (any_strong) {
          effect_label <- "regulation_present"
        } else {
          effect_label <- "mixed"
        }
      }

      median_within_rho <- round(median(valid_within, na.rm = TRUE), 4)
      max_within_abs_rho <- round(max(abs(valid_within), na.rm = TRUE), 4)
      regulation_ratio <- ifelse(abs(global_rho) > 0,
                                  round(abs(median_within_rho) / abs(global_rho), 3), NA)

      coupling_list[[g]] <- data.frame(
        gene = g,
        global_rho = round(global_rho, 4),
        median_within_rho = median_within_rho,
        max_within_abs_rho = max_within_abs_rho,
        regulation_ratio = regulation_ratio,
        n_clusters_tested = length(valid_within),
        n_clusters_sig = sum(within_ps[!is.na(within_ps)] < 0.05),
        effect_label = effect_label,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(coupling_list) > 0) {
    coupling_df <- do.call(rbind, coupling_list)
    rownames(coupling_df) <- NULL
    coupling_df <- coupling_df[order(-abs(coupling_df$global_rho)), ]
    fwrite(coupling_df, file.path(out_dir, "cell_state_coupling.csv"))

    # 统计
    n_comp <- sum(coupling_df$effect_label == "composition_driven")
    n_reg <- sum(coupling_df$effect_label == "regulation_present")
    n_het <- sum(coupling_df$effect_label == "cluster_heterogeneous")
    n_mix <- sum(coupling_df$effect_label == "mixed")

    cat(sprintf("        composition_driven=%d, regulation_present=%d, heterogeneous=%d, mixed=%d\n",
                n_comp, n_reg, n_het, n_mix))
  } else {
    coupling_df <- data.frame()
    n_comp <- 0; n_reg <- 0; n_het <- 0; n_mix <- 0
    cat("        ⚠ 无tier1基因或无有效cluster, 跳过分层分析\n")
  }

  # ---------------------------------------------------------------
  # 第四步：ρ分解对比图
  # ---------------------------------------------------------------

  if (nrow(coupling_df) > 0) {
    tryCatch({
      plot_df <- coupling_df[1:min(30, nrow(coupling_df)), ]
      plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene))

      p_decomp <- ggplot(plot_df) +
        geom_segment(aes(x = gene, xend = gene, y = 0, yend = global_rho),
                     color = "grey60", linewidth = 0.5) +
        geom_point(aes(x = gene, y = global_rho, color = "Global ρ"), size = 2.5) +
        geom_point(aes(x = gene, y = median_within_rho, color = "Median within-cluster ρ"),
                   size = 2, shape = 17) +
        scale_color_manual(values = c("Global ρ" = "#D62728",
                                       "Median within-cluster ρ" = "#1F77B4"),
                           name = NULL) +
        geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
        coord_flip() +
        labs(title = paste0(sample_name, " — Effect Decomposition (top ", nrow(plot_df), ")"),
             x = NULL, y = "ρ (Spearman)") +
        theme_minimal() +
        theme(plot.title = element_text(size = 10), legend.position = "bottom")

      h <- max(5, nrow(plot_df) * 0.25)
      ggsave(file.path(out_dir, "rho_decomposition.png"), p_decomp,
             width = 9, height = h, dpi = 130, limitsize = FALSE)
      cat("        ✓ ρ分解图\n")
    }, error = function(e) cat(sprintf("        ⚠ 分解图失败: %s\n", e$message)))
  }

  # ---------------------------------------------------------------
  # 第五步：细胞周期分析
  # ---------------------------------------------------------------

  cat("  [5/6] 细胞周期分析 ...\n")

  # 检查Tirosh基因匹配率
  all_genes <- rownames(sct_data)
  # 处理Brain_Mouse下划线→短横线问题
  s_match <- intersect(s_genes, all_genes)
  g2m_match <- intersect(g2m_genes, all_genes)

  cc_done <- FALSE
  if (length(s_match) >= MIN_CC_GENES && length(g2m_match) >= MIN_CC_GENES) {
    cat(sprintf("        S期匹配=%d, G2M匹配=%d → 执行细胞周期评分\n",
                length(s_match), length(g2m_match)))

    tryCatch({
      seurat_obj <- CellCycleScoring(seurat_obj,
                                      s.features = s_match,
                                      g2m.features = g2m_match,
                                      set.ident = FALSE)

      s_score <- seurat_obj$S.Score
      g2m_score <- seurat_obj$G2M.Score
      phase <- seurat_obj$Phase

      # 密度 vs S.Score
      sp_s <- suppressWarnings(cor.test(density_vec, s_score, method = "spearman",
                                         use = "complete.obs"))
      # 密度 vs G2M.Score
      sp_g2m <- suppressWarnings(cor.test(density_vec, g2m_score, method = "spearman",
                                           use = "complete.obs"))

      # 各phase的密度分布
      phase_density <- data.frame(phase = phase, density = density_vec)
      phase_density <- phase_density[!is.na(phase_density$density), ]
      phase_stats <- aggregate(density ~ phase, phase_density, function(x) {
        c(median = median(x), mean = mean(x), n = length(x))
      })

      cc_result <- data.frame(
        sample_name = sample_name,
        n_s_genes_matched = length(s_match),
        n_g2m_genes_matched = length(g2m_match),
        rho_density_s_score = round(sp_s$estimate, 4),
        p_density_s_score = sp_s$p.value,
        rho_density_g2m_score = round(sp_g2m$estimate, 4),
        p_density_g2m_score = sp_g2m$p.value,
        pct_G1 = round(sum(phase == "G1") / length(phase) * 100, 1),
        pct_S = round(sum(phase == "S") / length(phase) * 100, 1),
        pct_G2M = round(sum(phase == "G2M") / length(phase) * 100, 1),
        stringsAsFactors = FALSE
      )

      fwrite(cc_result, file.path(out_dir, "cell_cycle_density.csv"))

      cat(sprintf("        density↔S.Score: ρ=%.4f (p=%.2e)\n",
                  sp_s$estimate, sp_s$p.value))
      cat(sprintf("        density↔G2M.Score: ρ=%.4f (p=%.2e)\n",
                  sp_g2m$estimate, sp_g2m$p.value))
      cat(sprintf("        Phase分布: G1=%.1f%%, S=%.1f%%, G2M=%.1f%%\n",
                  cc_result$pct_G1, cc_result$pct_S, cc_result$pct_G2M))
      cc_done <- TRUE

    }, error = function(e) {
      cat(sprintf("        ⚠ 细胞周期评分失败: %s\n", e$message))
    })
  } else {
    cat(sprintf("        ⚠ 匹配不足 (S=%d, G2M=%d, 需≥%d), 跳过\n",
                length(s_match), length(g2m_match), MIN_CC_GENES))
  }

  # ---------------------------------------------------------------
  # 第六步：汇总
  # ---------------------------------------------------------------

  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  cat(sprintf("  完成! 耗时 %.1f 秒\n\n", t_elapsed))

  summary_list[[sample_name]] <- data.frame(
    sample_name = sample_name,
    dataset = dataset_name,
    n_cells = n_cells,
    n_clusters = n_clusters,
    n_clusters_valid = length(valid_clusters),
    n_tier1 = n_tier1,
    n_composition_driven = n_comp,
    n_regulation_present = n_reg,
    n_cluster_heterogeneous = n_het,
    n_mixed = n_mix,
    cc_analyzed = cc_done,
    cc_s_genes_matched = ifelse(cc_done, length(s_match), NA),
    cc_g2m_genes_matched = ifelse(cc_done, length(g2m_match), NA),
    cc_rho_s = ifelse(cc_done, round(sp_s$estimate, 4), NA),
    cc_rho_g2m = ifelse(cc_done, round(sp_g2m$estimate, 4), NA),
    time_seconds = t_elapsed,
    stringsAsFactors = FALSE
  )

  rm(seurat_obj, sct_data); gc(verbose = FALSE)
}

# ============================================================================
# 全局汇总
# ============================================================================

cat("============================================\n")
cat("全局汇总\n")
cat("============================================\n")

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL
fwrite(summary_df, file.path(OUTPUT_DIR, "ALL_SAMPLES_R6_SUMMARY.csv"))

cat(sprintf("%-45s %5s %5s %5s %5s %5s %5s %5s %6s %6s %7s\n",
            "sample", "tier1", "comp", "reg", "het", "mix", "cc?", "S_g", "ρ_S", "ρ_G2M", "time"))
cat(paste(rep("-", 130), collapse = ""), "\n")

for (j in 1:nrow(summary_df)) {
  r <- summary_df[j, ]
  cc_str <- ifelse(r$cc_analyzed, "Y", "N")
  s_g_str <- ifelse(is.na(r$cc_s_genes_matched), "-", as.character(r$cc_s_genes_matched))
  rho_s_str <- ifelse(is.na(r$cc_rho_s), "  -", sprintf("%.3f", r$cc_rho_s))
  rho_g_str <- ifelse(is.na(r$cc_rho_g2m), "  -", sprintf("%.3f", r$cc_rho_g2m))
  cat(sprintf("%-45s %5d %5d %5d %5d %5d %5s %5s %6s %6s %7.1f\n",
              r$sample_name, r$n_tier1, r$n_composition_driven,
              r$n_regulation_present, r$n_cluster_heterogeneous, r$n_mixed,
              cc_str, s_g_str, rho_s_str, rho_g_str, r$time_seconds))
}

cat("\n--- 按数据集分组 ---\n")
for (ds in unique(summary_df$dataset)) {
  sub <- summary_df[summary_df$dataset == ds, ]
  cat(sprintf("  %-30s %d 样本, composition=%d-%d, regulation=%d-%d, cc=%d/%d 样本可分析\n",
              ds, nrow(sub),
              min(sub$n_composition_driven), max(sub$n_composition_driven),
              min(sub$n_regulation_present), max(sub$n_regulation_present),
              sum(sub$cc_analyzed), nrow(sub)))
}

cat("\n============================================\n")
cat("R6 全部完成!\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
