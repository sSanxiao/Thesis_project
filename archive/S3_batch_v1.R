# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# S3_batch_v1.R
# Stage 3 准备工作的一体化脚本：第一波 + 第二波任务
# ============================================================================
# 任务清单：
#   Part A: 修复 R7 effect_label bug 并重新生成 ALL_SAMPLES_R7_PROFILE.csv
#   Part B: 为 EOMES/HHIP/GLI2/PTCH1 等关键基因生成 MB 4 样本的空间图
#   Part C: 提取 R8 关键数字到单一汇总 csv
#   Part D: 查 EOMES 在所有数据集 panel 的状态
#   Part E: EOMES+ cluster 识别（4 MB 样本）
#   Part F: MB 47 的增强版 ρ 分解图（确保 EOMES/HHIP/GLI2 在图上）
#   Part G: 生成汇总 Markdown 报告
#
# 设计原则：
#   - 不修改任何原 pipeline 代码
#   - 所有输出到 ./results/20260420/
#   - 每个 Part 独立 tryCatch 包裹，失败不影响其他 Part
#   - 生成人类可读的 run.log 和 summary_report.md
#
# 运行：
#   cd ./
#   nohup Rscript S3_batch_v1.R > S3_batch_v1.log 2>&1 &
#
# 预期耗时：2-3 小时（主要在 Part B 和 Part E 的 Seurat rds 加载）
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(viridis)
  library(scales)  # for scales::squish in spatial plots
})

options(future.globals.maxSize = Inf)

# ============================================================================
# 全局配置
# ============================================================================

REGISTRY_PATH <- "./data/sample_registry.json"
R2_DIR        <- "./results/R2_Results"
R3_DIR        <- "./results/R3_Results"
R4_DIR        <- "./results/R4_Results"
R6_DIR        <- "./results/R6_Results"
R7_DIR        <- "./results/R7_Results"
R8_DIR        <- "./results/R8_Results"

# 新输出目录
OUTPUT_DIR <- "./results/20260420"

DENSITY_COL   <- "density_knn_main_piecewise"

# MB 样本列表（Stage 3 主要关注）
MB_SAMPLES <- c(
  "Medulloblastoma_Human/GSM8840046",
  "Medulloblastoma_Human/GSM8840047",
  "Medulloblastoma_Human/GSM8840048",
  "Medulloblastoma_Human/GSM8840049"
)

# Stage 3 关键基因
STAGE3_TARGET_GENES <- c("EOMES", "HHIP", "GLI2", "PTCH1", "OTX2", "MEIS2",
                          "PAX6", "LMX1A", "TBR1", "MYCN", "MKI67", "TOP2A")

# 空间图目标：MB 4 样本 × (EOMES, HHIP, GLI2, PTCH1)
# MB 47 是主样本，跑全部；其他 3 个只跑 EOMES
SPATIAL_TARGETS <- list(
  c("Medulloblastoma_Human/GSM8840047", "EOMES"),
  c("Medulloblastoma_Human/GSM8840047", "HHIP"),
  c("Medulloblastoma_Human/GSM8840047", "GLI2"),
  c("Medulloblastoma_Human/GSM8840047", "PTCH1"),
  c("Medulloblastoma_Human/GSM8840046", "EOMES"),
  c("Medulloblastoma_Human/GSM8840048", "EOMES"),
  c("Medulloblastoma_Human/GSM8840049", "EOMES")
)

# 创建输出子目录
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (sub in c("partA_R7_fixed", "partB_spatial_plots", "partC_R8_summary",
              "partD_gene_presence", "partE_eomes_clusters", "partF_rho_decomp")) {
  dir.create(file.path(OUTPUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# Logging utilities
# ============================================================================

LOG_FILE <- file.path(OUTPUT_DIR, "run.log")
REPORT_FILE <- file.path(OUTPUT_DIR, "summary_report.md")

# 初始化 log
cat("", file = LOG_FILE)  # 清空
cat("", file = REPORT_FILE)  # 清空

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", ts, msg)
  cat(line)
  cat(line, file = LOG_FILE, append = TRUE)
}

log_section <- function(title) {
  sep <- paste(rep("=", 70), collapse = "")
  log_msg(sep)
  log_msg(title)
  log_msg(sep)
}

report_append <- function(content) {
  cat(content, file = REPORT_FILE, append = TRUE)
  cat("\n", file = REPORT_FILE, append = TRUE)
}

log_section("S3 批处理脚本启动")
log_msg(sprintf("输出目录: %s", OUTPUT_DIR))
log_msg(sprintf("开始时间: %s", format(Sys.time())))

# 初始化报告头部
report_append("# Stage 3 批处理结果报告\n")
report_append(sprintf("**运行时间**: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
report_append(sprintf("**输出目录**: `%s`\n", OUTPUT_DIR))
report_append("---\n")

# 全局状态追踪
task_status <- list()

# ============================================================================
# Part A: 修复 R7 bug 并重新生成 ALL_SAMPLES_R7_PROFILE.csv
# ============================================================================
# 原 bug: R7 找 "effect_class" 但 R6 产出 "effect_label"
# 修复方案: 重新读 R4 和 R6 的 csv，用正确列名做 merge，regenerate profile
# ============================================================================

log_section("Part A: 修复 R7 effect_label bug")

tryCatch({
  registry <- fromJSON(REGISTRY_PATH)
  sample_names <- names(registry)
  log_msg(sprintf("样本数: %d", length(sample_names)))
  
  profile_list <- list()
  
  for (i in seq_along(sample_names)) {
    sname <- sample_names[i]
    parts <- strsplit(sname, "/")[[1]]
    dataset_name <- parts[1]
    sample_subname <- parts[2]
    
    r4_path <- file.path(R4_DIR, dataset_name, sample_subname, "filtered_density_genes.csv")
    r6_path <- file.path(R6_DIR, dataset_name, sample_subname, "cell_state_coupling.csv")
    
    if (!file.exists(r4_path)) {
      log_msg(sprintf("[%d/%d] %s  ⚠ R4 缺失", i, length(sample_names), sname))
      next
    }
    
    r4 <- fread(r4_path)
    r6_exists <- file.exists(r6_path)
    
    # 合并 R4 和 R6（用正确列名）
    merged <- copy(r4)
    if (r6_exists) {
      r6 <- fread(r6_path)
      if (nrow(r6) > 0) {
        # 正确的 R6 列名
        r6_desired_cols <- c("gene", "effect_label", "median_within_rho",
                              "max_within_abs_rho", "regulation_ratio",
                              "n_clusters_tested", "n_clusters_sig", "global_rho")
        r6_cols <- intersect(names(r6), r6_desired_cols)
        if (length(r6_cols) > 1) {
          r6_sub <- r6[, ..r6_cols]
          # 重命名 global_rho 避免和 R4 的 rho_knn_main 混淆（如果有）
          if ("global_rho" %in% names(r6_sub)) {
            setnames(r6_sub, "global_rho", "r6_global_rho")
          }
          merged <- merge(merged, r6_sub, by = "gene", all.x = TRUE)
        }
      }
    }
    
    # 计算样本级统计
    n_total <- nrow(merged)
    n_tier1 <- sum(merged$tier == "tier1_strict", na.rm = TRUE)
    n_tier2 <- sum(merged$tier == "tier2_moderate", na.rm = TRUE)
    n_tier3 <- sum(merged$tier == "tier3_lenient", na.rm = TRUE)
    n_ns <- sum(merged$tier == "not_significant", na.rm = TRUE)
    
    tier1_rows <- merged[tier == "tier1_strict"]
    
    if (nrow(tier1_rows) > 0) {
      n_pos <- sum(tier1_rows$rho_knn_main > 0, na.rm = TRUE)
      n_neg <- sum(tier1_rows$rho_knn_main < 0, na.rm = TRUE)
      med_rho <- round(median(abs(tier1_rows$rho_knn_main), na.rm = TRUE), 4)
      max_rho <- round(max(abs(tier1_rows$rho_knn_main), na.rm = TRUE), 4)
      
      # top5
      top5 <- tier1_rows[order(-abs(rho_knn_main))][1:min(5, nrow(tier1_rows))]
      top5_str <- paste(sprintf("%s(%.3f)", top5$gene, top5$rho_knn_main),
                        collapse = ", ")
      
      # effect_label 统计（这次用正确列名！）
      if ("effect_label" %in% names(tier1_rows)) {
        n_comp <- sum(tier1_rows$effect_label == "composition_driven", na.rm = TRUE)
        n_reg <- sum(tier1_rows$effect_label == "regulation_present", na.rm = TRUE)
        n_het <- sum(tier1_rows$effect_label == "cluster_heterogeneous", na.rm = TRUE)
        n_mix <- sum(tier1_rows$effect_label == "mixed", na.rm = TRUE)
        reg_pct <- round(n_reg / n_tier1 * 100, 1)
      } else {
        n_comp <- NA; n_reg <- NA; n_het <- NA; n_mix <- NA; reg_pct <- NA
      }
    } else {
      n_pos <- 0; n_neg <- 0; med_rho <- NA; max_rho <- NA; top5_str <- ""
      n_comp <- NA; n_reg <- NA; n_het <- NA; n_mix <- NA; reg_pct <- NA
    }
    
    profile_list[[sname]] <- data.frame(
      sample_name = sname,
      dataset = dataset_name,
      n_genes_total = n_total,
      n_tier1 = n_tier1,
      n_tier2 = n_tier2,
      n_tier3 = n_tier3,
      n_not_sig = n_ns,
      tier1_pos = n_pos,
      tier1_neg = n_neg,
      tier1_neg_pct = ifelse(n_tier1 > 0, round(n_neg / n_tier1 * 100, 1), NA),
      tier1_median_abs_rho = med_rho,
      tier1_max_abs_rho = max_rho,
      tier1_composition = ifelse(is.na(n_comp), NA, n_comp),
      tier1_regulation = ifelse(is.na(n_reg), NA, n_reg),
      tier1_heterogeneous = ifelse(is.na(n_het), NA, n_het),
      tier1_mixed = ifelse(is.na(n_mix), NA, n_mix),
      tier1_regulation_pct = reg_pct,
      top5_genes = top5_str,
      stringsAsFactors = FALSE
    )
    
    log_msg(sprintf("[%2d/%d] %-45s tier1=%3d (reg=%s comp=%s het=%s mix=%s)",
                    i, length(sample_names), sname, n_tier1,
                    ifelse(is.na(n_reg), "NA", n_reg),
                    ifelse(is.na(n_comp), "NA", n_comp),
                    ifelse(is.na(n_het), "NA", n_het),
                    ifelse(is.na(n_mix), "NA", n_mix)))
  }
  
  profile_df <- do.call(rbind, profile_list)
  rownames(profile_df) <- NULL
  
  out_path <- file.path(OUTPUT_DIR, "partA_R7_fixed", "ALL_SAMPLES_R7_PROFILE_fixed.csv")
  fwrite(profile_df, out_path)
  log_msg(sprintf("✓ 保存: %s", out_path))
  
  # 对比原始 bug 版本（如果存在）
  orig_path <- file.path(R7_DIR, "ALL_SAMPLES_R7_PROFILE.csv")
  n_orig_na <- NA_integer_
  n_fixed_na <- sum(is.na(profile_df$tier1_regulation_pct))
  
  if (file.exists(orig_path)) {
    orig <- fread(orig_path)
    if ("tier1_regulation_pct" %in% names(orig)) {
      n_orig_na <- sum(is.na(orig$tier1_regulation_pct))
    }
    log_msg(sprintf("对比：原始 regulation_pct NA 数 = %s，修复后 NA 数 = %d",
                    ifelse(is.na(n_orig_na), "N/A", as.character(n_orig_na)),
                    n_fixed_na))
  }
  
  # 报告
  report_append("## Part A: R7 bug 修复 ✓\n")
  report_append(sprintf("修复前 `tier1_regulation_pct` NA 数: **%s/22**\n",
                        ifelse(is.na(n_orig_na), "N/A", as.character(n_orig_na))))
  report_append(sprintf("修复后 NA 数: **%d/22**\n", n_fixed_na))
  report_append("\n**修复后 profile（关键列）：**\n\n")
  report_append("| sample | tier1 | reg | comp | het | mix | reg_pct |\n")
  report_append("|--------|-------|-----|------|-----|-----|--------|\n")
  for (i in 1:nrow(profile_df)) {
    r <- profile_df[i, ]
    report_append(sprintf("| %s | %d | %s | %s | %s | %s | %s |\n",
                          r$sample_name, r$n_tier1,
                          ifelse(is.na(r$tier1_regulation), "—", r$tier1_regulation),
                          ifelse(is.na(r$tier1_composition), "—", r$tier1_composition),
                          ifelse(is.na(r$tier1_heterogeneous), "—", r$tier1_heterogeneous),
                          ifelse(is.na(r$tier1_mixed), "—", r$tier1_mixed),
                          ifelse(is.na(r$tier1_regulation_pct), "—",
                                 sprintf("%.1f%%", r$tier1_regulation_pct))))
  }
  report_append("\n")
  
  task_status$partA <<- "success"
}, error = function(e) {
  log_msg(sprintf("✗ Part A 失败: %s", e$message))
  report_append("## Part A: R7 bug 修复 ✗\n")
  report_append(sprintf("**失败**: %s\n", e$message))
  task_status$partA <<- "failed"
})

# ============================================================================
# Part D（先做，不依赖 rds 加载）: 查 EOMES 等关键基因在各数据集 panel 的状态
# ============================================================================

log_section("Part D: EOMES 等关键基因在各 panel 的状态")

tryCatch({
  registry <- fromJSON(REGISTRY_PATH)
  sample_names <- names(registry)
  
  # 对每个样本，查关键基因在 R3 csv 中的情况
  presence_list <- list()
  
  for (gene in STAGE3_TARGET_GENES) {
    log_msg(sprintf("查询基因: %s", gene))
    
    # 大小写都查（mouse vs human）
    gene_variants <- unique(c(gene, toupper(gene),
                               paste0(substr(gene, 1, 1),
                                      tolower(substr(gene, 2, nchar(gene))))))
    
    for (sname in sample_names) {
      parts <- strsplit(sname, "/")[[1]]
      r3_path <- file.path(R3_DIR, parts[1], parts[2], "density_gene_correlations.csv")
      r4_path <- file.path(R4_DIR, parts[1], parts[2], "filtered_density_genes.csv")
      
      if (!file.exists(r3_path)) next
      
      r3 <- fread(r3_path)
      r4 <- if (file.exists(r4_path)) fread(r4_path) else NULL
      
      # 查这个基因的任何变体
      found_row <- NULL
      found_gene_name <- NA
      for (gv in gene_variants) {
        if (gv %in% r3$gene) {
          found_row <- r3[gene == gv][1]
          found_gene_name <- gv
          break
        }
      }
      
      if (is.null(found_row)) {
        # 不在 panel
        presence_list[[paste0(gene, "__", sname)]] <- data.frame(
          query_gene = gene,
          actual_gene = NA,
          sample_name = sname,
          dataset = parts[1],
          in_panel = FALSE,
          rho_knn_main = NA_real_,
          q_knn_main = NA_real_,
          convergence = NA_character_,
          tier = NA_character_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # 找到了 - 查 R4 的 tier
      tier_value <- NA_character_
      if (!is.null(r4)) {
        r4_row <- r4[gene == found_gene_name]
        if (nrow(r4_row) > 0) tier_value <- as.character(r4_row$tier[1])
      }
      
      presence_list[[paste0(gene, "__", sname)]] <- data.frame(
        query_gene = gene,
        actual_gene = found_gene_name,
        sample_name = sname,
        dataset = parts[1],
        in_panel = TRUE,
        rho_knn_main = as.numeric(found_row$rho_knn_main),
        q_knn_main = as.numeric(found_row$q_knn_main),
        convergence = as.character(found_row$convergence),
        tier = tier_value,
        stringsAsFactors = FALSE
      )
    }
  }
  
  presence_df <- do.call(rbind, presence_list)
  rownames(presence_df) <- NULL
  
  out_path <- file.path(OUTPUT_DIR, "partD_gene_presence", "key_genes_panel_status.csv")
  fwrite(presence_df, out_path)
  log_msg(sprintf("✓ 保存: %s", out_path))
  
  # EOMES 专门报告
  eomes_rows <- presence_df[presence_df$query_gene == "EOMES", ]
  in_panel_count <- sum(eomes_rows$in_panel)
  tier1_count <- sum(!is.na(eomes_rows$tier) & eomes_rows$tier == "tier1_strict")
  
  log_msg(sprintf("EOMES 在 %d/%d 样本的 panel 中", in_panel_count, nrow(eomes_rows)))
  log_msg(sprintf("EOMES 在 %d 样本达到 tier1", tier1_count))
  
  # 报告
  report_append("## Part D: 关键基因在 panel 的状态\n")
  report_append("\n**EOMES 在所有 22 样本的状态**：\n\n")
  report_append("| sample | in_panel | ρ | q | convergence | tier |\n")
  report_append("|--------|----------|---|---|-------------|------|\n")
  for (i in 1:nrow(eomes_rows)) {
    r <- eomes_rows[i, ]
    report_append(sprintf("| %s | %s | %s | %s | %s | %s |\n",
                          r$sample_name,
                          ifelse(r$in_panel, "✓", "✗"),
                          ifelse(is.na(r$rho_knn_main), "—",
                                 sprintf("%+.4f", r$rho_knn_main)),
                          ifelse(is.na(r$q_knn_main), "—",
                                 sprintf("%.2e", r$q_knn_main)),
                          ifelse(is.na(r$convergence), "—", r$convergence),
                          ifelse(is.na(r$tier), "—", r$tier)))
  }
  report_append("\n")
  
  # 其他关键基因汇总
  report_append("**其他关键基因的 panel 覆盖汇总**：\n\n")
  report_append("| gene | in_panel | tier1_samples | datasets_with |\n")
  report_append("|------|----------|---------------|---------------|\n")
  for (gene in STAGE3_TARGET_GENES) {
    sub <- presence_df[presence_df$query_gene == gene, ]
    n_in_panel <- sum(sub$in_panel)
    n_tier1 <- sum(!is.na(sub$tier) & sub$tier == "tier1_strict")
    ds_list <- unique(sub$dataset[sub$in_panel])
    report_append(sprintf("| %s | %d/%d | %d | %s |\n",
                          gene, n_in_panel, nrow(sub), n_tier1,
                          paste(ds_list, collapse = ", ")))
  }
  report_append("\n")
  
  task_status$partD <<- "success"
}, error = function(e) {
  log_msg(sprintf("✗ Part D 失败: %s", e$message))
  report_append("## Part D: 关键基因在 panel 的状态 ✗\n")
  report_append(sprintf("**失败**: %s\n", e$message))
  task_status$partD <<- "failed"
})

# ============================================================================
# Part C: 提取 R8 关键数字到单一汇总
# ============================================================================

log_section("Part C: R8 关键数字汇总")

tryCatch({
  # 读 R8 的主要输出
  comparisons_path <- file.path(R8_DIR, "ALL_COMPARISONS_R8_SUMMARY.csv")
  gene_summary_path <- file.path(R8_DIR, "global_gene_summary.csv")
  landscape_path <- file.path(R8_DIR, "global_density_gene_landscape.csv")
  
  report_append("## Part C: R8 关键数字汇总\n\n")
  
  if (file.exists(comparisons_path)) {
    comp <- fread(comparisons_path)
    log_msg(sprintf("R8 comparisons 行数: %d", nrow(comp)))
    
    # 重点行：MB 相关的比较
    mb_comparisons <- comp[dataset_A == "Medulloblastoma_Human" | 
                           dataset_B == "Medulloblastoma_Human"]
    
    if (nrow(mb_comparisons) > 0) {
      # 保存 MB 相关比较
      fwrite(mb_comparisons, 
             file.path(OUTPUT_DIR, "partC_R8_summary", "MB_related_comparisons.csv"))
      log_msg(sprintf("MB 相关比较: %d 对", nrow(mb_comparisons)))
      
      report_append("### MB 相关的跨数据集比较\n\n")
      report_append("| vs | type | shared | rank_ρ | Fisher_p | t1_A | t1_B | t1_ov | dir% |\n")
      report_append("|----|------|--------|--------|----------|------|------|-------|------|\n")
      for (i in 1:nrow(mb_comparisons)) {
        r <- mb_comparisons[i, ]
        other_ds <- ifelse(r$dataset_A == "Medulloblastoma_Human",
                           r$dataset_B, r$dataset_A)
        report_append(sprintf("| %s | %s | %d | %.3f | %.2e | %d | %d | %d | %s |\n",
                              other_ds,
                              substr(r$comparison_type, 1, 5),
                              r$n_shared_genes,
                              r$rank_correlation_rho,
                              ifelse(is.na(r$tier1_fisher_p), 1, r$tier1_fisher_p),
                              r$tier1_A, r$tier1_B, r$tier1_overlap,
                              ifelse(is.na(r$direction_consistency_pct), "—",
                                     sprintf("%.1f%%", r$direction_consistency_pct))))
      }
      report_append("\n")
    }
    
    # 全部比较汇总表
    fwrite(comp, file.path(OUTPUT_DIR, "partC_R8_summary",
                           "ALL_COMPARISONS_copy.csv"))
    
    report_append("### 全部跨数据集比较汇总\n\n")
    report_append("| A | B | type | shared | rank_ρ | Fisher_p | t1_A | t1_B | t1_ov | dir% |\n")
    report_append("|---|---|------|--------|--------|----------|------|------|-------|------|\n")
    for (i in 1:nrow(comp)) {
      r <- comp[i, ]
      report_append(sprintf("| %s | %s | %s | %d | %.3f | %.2e | %d | %d | %d | %s |\n",
                            r$dataset_A, r$dataset_B,
                            substr(r$comparison_type, 1, 5),
                            r$n_shared_genes, r$rank_correlation_rho,
                            ifelse(is.na(r$tier1_fisher_p), 1, r$tier1_fisher_p),
                            r$tier1_A, r$tier1_B, r$tier1_overlap,
                            ifelse(is.na(r$direction_consistency_pct), "—",
                                   sprintf("%.1f%%", r$direction_consistency_pct))))
    }
    report_append("\n")
  } else {
    log_msg("⚠ R8 comparisons 文件不存在")
    report_append("⚠ `ALL_COMPARISONS_R8_SUMMARY.csv` 不存在\n\n")
  }
  
  # EOMES 在 global_gene_summary 的位置
  if (file.exists(gene_summary_path)) {
    gs <- fread(gene_summary_path)
    log_msg(sprintf("global_gene_summary 行数: %d", nrow(gs)))
    
    eomes_row <- gs[gene_upper == "EOMES"]
    
    report_append("### EOMES 在 global_gene_summary 的状态\n\n")
    if (nrow(eomes_row) > 0) {
      r <- eomes_row[1, ]
      report_append(sprintf("- **gene_upper**: %s\n", r$gene_upper))
      report_append(sprintf("- **n_datasets**: %d\n", r$n_datasets))
      report_append(sprintf("- **datasets**: %s\n", r$datasets))
      report_append(sprintf("- **species_list**: %s\n", r$species_list))
      report_append(sprintf("- **directions**: %s\n", r$directions))
      report_append(sprintf("- **median_rhos**: %s\n", r$median_rhos))
    } else {
      report_append("⚠ EOMES 不在 global_gene_summary！\n")
    }
    report_append("\n")
    
    # n_datasets 分布
    report_append("### global_gene_summary 的 n_datasets 分布\n\n")
    report_append("| n_datasets | gene count |\n")
    report_append("|-----------|------------|\n")
    for (n in sort(unique(gs$n_datasets), decreasing = TRUE)) {
      report_append(sprintf("| %d | %d |\n", n, sum(gs$n_datasets == n)))
    }
    report_append(sprintf("\n**总基因数**: %d\n\n", nrow(gs)))
    
    # n_datasets >= 2 的基因列表
    multi_ds <- gs[n_datasets >= 2]
    fwrite(multi_ds, file.path(OUTPUT_DIR, "partC_R8_summary",
                               "genes_in_multiple_datasets.csv"))
    
    report_append(sprintf("### 在 ≥2 数据集中 consistent 的基因（n=%d）\n\n",
                          nrow(multi_ds)))
    report_append("（完整列表见 `partC_R8_summary/genes_in_multiple_datasets.csv`）\n\n")
    if (nrow(multi_ds) > 0) {
      report_append("**Top 20**：\n\n")
      report_append("| gene | n_datasets | datasets | directions |\n")
      report_append("|------|-----------|----------|------------|\n")
      for (i in 1:min(20, nrow(multi_ds))) {
        r <- multi_ds[i, ]
        report_append(sprintf("| %s | %d | %s | %s |\n",
                              r$gene_upper, r$n_datasets, r$datasets,
                              r$directions))
      }
      report_append("\n")
    }
  }
  
  task_status$partC <<- "success"
}, error = function(e) {
  log_msg(sprintf("✗ Part C 失败: %s", e$message))
  report_append(sprintf("**Part C 失败**: %s\n", e$message))
  task_status$partC <<- "failed"
})

# ============================================================================
# Part F: 增强版 ρ 分解图（MB 47，确保 EOMES/HHIP/GLI2 出现）
# ============================================================================

log_section("Part F: MB 47 增强版 ρ 分解图")

tryCatch({
  coupling_path <- file.path(R6_DIR, "Medulloblastoma_Human", "GSM8840047",
                              "cell_state_coupling.csv")
  
  if (!file.exists(coupling_path)) {
    log_msg("⚠ MB 47 cell_state_coupling.csv 不存在")
    report_append("## Part F: ρ 分解图 ⚠\n**跳过**: MB 47 coupling csv 不存在\n\n")
  } else {
    coupling_df <- fread(coupling_path)
    log_msg(sprintf("MB 47 coupling 行数: %d", nrow(coupling_df)))
    
    # 按 |global_rho| 降序
    coupling_df <- coupling_df[order(-abs(global_rho))]
    
    # 检查 EOMES 等基因的位置
    for (g in STAGE3_TARGET_GENES) {
      rank_pos <- which(coupling_df$gene == g)
      if (length(rank_pos) > 0) {
        log_msg(sprintf("  %s: rank=%d, global_rho=%.4f, effect_label=%s",
                        g, rank_pos[1], coupling_df$global_rho[rank_pos[1]],
                        coupling_df$effect_label[rank_pos[1]]))
      } else {
        log_msg(sprintf("  %s: 不在 coupling_df（不是 tier1）", g))
      }
    }
    
    # 增强版 plot_df：top30 + target genes 并集
    top30 <- coupling_df[1:min(30, nrow(coupling_df))]
    target_rows <- coupling_df[gene %in% STAGE3_TARGET_GENES]
    plot_df <- unique(rbind(top30, target_rows))
    plot_df <- plot_df[order(-abs(global_rho))]
    
    # 标记 Stage 3 目标基因（在 factor 化之前）
    plot_df[, is_target := gene %in% STAGE3_TARGET_GENES]
    
    # 因子化 gene 保持顺序（coord_flip 后图顶部 = 列表前面）
    plot_df[, gene := factor(gene, levels = rev(as.character(gene)))]
    
    # 预计算 x 位置的常量
    text_x_pos <- max(abs(plot_df$global_rho)) * 1.15
    
    p_decomp <- ggplot(plot_df) +
      geom_segment(aes(x = gene, xend = gene, y = 0, yend = global_rho),
                   color = "grey60", linewidth = 0.5) +
      geom_point(aes(x = gene, y = global_rho, color = "Global ρ"), size = 2.5) +
      geom_point(aes(x = gene, y = median_within_rho,
                     color = "Median within-cluster ρ"),
                 size = 2, shape = 17) +
      # 高亮 target genes 的基因名
      geom_text(data = plot_df[is_target == TRUE],
                aes(x = gene, y = text_x_pos, label = as.character(gene)),
                hjust = 0, size = 3, color = "darkred", fontface = "bold",
                inherit.aes = FALSE) +
      scale_color_manual(values = c("Global ρ" = "#D62728",
                                     "Median within-cluster ρ" = "#1F77B4"),
                         name = NULL) +
      geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
      coord_flip() +
      labs(title = sprintf("MB GSM8840047 — Effect Decomposition (n=%d, Stage 3 targets highlighted)",
                            nrow(plot_df)),
           x = NULL, y = "ρ (Spearman)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            legend.position = "bottom",
            axis.text.y = element_text(size = 8))
    
    h <- max(6, nrow(plot_df) * 0.25)
    out_path <- file.path(OUTPUT_DIR, "partF_rho_decomp",
                           "MB47_rho_decomposition_enhanced.png")
    ggsave(out_path, p_decomp, width = 10, height = h, dpi = 130, limitsize = FALSE)
    log_msg(sprintf("✓ 保存: %s", out_path))
    
    # 导出 plot_df 供你查
    fwrite(plot_df, file.path(OUTPUT_DIR, "partF_rho_decomp",
                               "MB47_rho_decomp_data.csv"))
    
    report_append("## Part F: MB 47 增强版 ρ 分解图 ✓\n\n")
    report_append(sprintf("图含 **%d** 基因（top 30 + %d Stage 3 targets 并集去重）\n\n",
                          nrow(plot_df), nrow(target_rows)))
    report_append("**Stage 3 target 基因的位置**：\n\n")
    report_append("| gene | rank (by |global_ρ|) | global_ρ | median_within | ratio | effect_label |\n")
    report_append("|------|----------------------|----------|---------------|-------|--------------|\n")
    
    for (g in STAGE3_TARGET_GENES) {
      rank_pos <- which(coupling_df$gene == g)
      if (length(rank_pos) > 0) {
        r <- coupling_df[rank_pos[1]]
        report_append(sprintf("| %s | %d | %+.4f | %+.4f | %.3f | %s |\n",
                              g, rank_pos[1], r$global_rho, r$median_within_rho,
                              ifelse(is.na(r$regulation_ratio), "—", r$regulation_ratio),
                              r$effect_label))
      } else {
        report_append(sprintf("| %s | — | — | — | — | not in tier1 |\n", g))
      }
    }
    report_append("\n")
    
    task_status$partF <<- "success"
  }
}, error = function(e) {
  log_msg(sprintf("✗ Part F 失败: %s", e$message))
  report_append(sprintf("**Part F 失败**: %s\n", e$message))
  task_status$partF <<- "failed"
})

# ============================================================================
# Part E + Part B（合并）: MB 样本处理 - EOMES+ cluster 识别 + 空间图
# ----------------------------------------------------------------------------
# 对每个 MB 样本加载 rds 一次，同时完成 Part E 和 Part B 的工作
# 节省 ~15-20 分钟（避免重复加载同一个 rds）
# ============================================================================

log_section("Part E + B (合并): MB 样本处理")

# 先按样本分组 SPATIAL_TARGETS
targets_by_sample <- list()
for (t in SPATIAL_TARGETS) {
  sname <- t[1]
  gene <- t[2]
  if (is.null(targets_by_sample[[sname]])) targets_by_sample[[sname]] <- c()
  targets_by_sample[[sname]] <- c(targets_by_sample[[sname]], gene)
}

# 辅助函数：画空间图
plot_spatial_fn <- function(df, gene, sample_name, out_path) {
  df <- df[!is.na(df$density) & !is.na(df$residual), ]
  if (nrow(df) < 100) return(FALSE)
  pt_size <- ifelse(nrow(df) > 200000, 0.05,
                    ifelse(nrow(df) > 50000, 0.1, 0.3))
  
  p1 <- ggplot(df, aes(x = x, y = y, color = density)) +
    geom_point(size = pt_size) + scale_color_viridis(option = "viridis") +
    coord_fixed() +
    labs(title = paste0(sample_name, " — Density"),
         x = "x (µm)", y = "y (µm)") +
    theme_minimal() + theme(plot.title = element_text(size = 9))
  
  lim <- quantile(abs(df$residual), 0.98, na.rm = TRUE)
  p2 <- ggplot(df, aes(x = x, y = y, color = residual)) +
    geom_point(size = pt_size) +
    scale_color_gradient2(low = "#1F77B4", mid = "white", high = "#D62728",
                          midpoint = 0, limits = c(-lim, lim),
                          oob = scales::squish) +
    coord_fixed() +
    labs(title = paste0(gene, " — SCT residual"),
         x = "x (µm)", y = "y (µm)") +
    theme_minimal() + theme(plot.title = element_text(size = 9))
  
  png(out_path, width = 1600, height = 700, res = 140)
  tryCatch({
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    } else {
      print(p1)
    }
  }, finally = dev.off())
  TRUE
}

# 辅助函数：画散点图（hex + boxplot）
plot_scatter_fn <- function(df, gene, sample_name, out_path) {
  df <- df[!is.na(df$density) & !is.na(df$residual), ]
  if (nrow(df) < 100) return(FALSE)
  
  p1 <- ggplot(df, aes(x = density, y = residual)) +
    geom_hex(bins = 60) + scale_fill_viridis(option = "plasma") +
    labs(title = paste0(sample_name, " — ", gene),
         x = "KNN density (main)", y = "SCT residual") +
    theme_minimal() + theme(plot.title = element_text(size = 9))
  
  df$density_bin <- cut(df$density,
                         breaks = quantile(df$density,
                                            probs = seq(0, 1, 0.1),
                                            na.rm = TRUE),
                         include.lowest = TRUE, labels = 1:10)
  df2 <- df[!is.na(df$density_bin), ]
  p2 <- ggplot(df2, aes(x = density_bin, y = residual)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.3) +
    labs(title = "Residual by density decile",
         x = "Density decile (1=low, 10=high)", y = "SCT residual") +
    theme_minimal() + theme(plot.title = element_text(size = 9))
  
  png(out_path, width = 1400, height = 500, res = 140)
  tryCatch({
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    } else {
      print(p1)
    }
  }, finally = dev.off())
  TRUE
}

eomes_cluster_results <- list()
partB_success_count <- 0
partB_total <- length(SPATIAL_TARGETS)

for (sname in MB_SAMPLES) {
  tryCatch({
    parts <- strsplit(sname, "/")[[1]]
    rds_path <- file.path(R2_DIR, parts[1], parts[2],
                          paste0(parts[2], "_seurat_R2.rds"))
    profile_path <- file.path(R6_DIR, parts[1], parts[2],
                              "cluster_density_profile.csv")
    
    if (!file.exists(rds_path)) {
      log_msg(sprintf("⚠ %s rds 不存在", sname))
      next
    }
    
    log_msg(sprintf("========== 开始处理 %s ==========", sname))
    log_msg(sprintf("读取 rds..."))
    t0 <- Sys.time()
    obj <- readRDS(rds_path)
    log_msg(sprintf("  加载耗时 %.1fs, %d 细胞",
                    as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    ncol(obj)))
    
    # 一次性提取所有要的数据
    sct_data <- GetAssayData(obj, assay = "SCT", layer = "data")
    density_vec <- obj@meta.data[[DENSITY_COL]]
    cluster_vec <- as.character(Idents(obj))
    x_coord <- obj$x_centroid
    y_coord <- obj$y_centroid
    
    # ---------- Part E: EOMES+ cluster 识别 ----------
    if ("EOMES" %in% rownames(sct_data)) {
      log_msg("  [Part E] 计算 EOMES+ cluster...")
      
      eomes_expr <- as.numeric(sct_data["EOMES", ])
      
      eomes_mean_by_cluster <- tapply(eomes_expr, cluster_vec, mean)
      eomes_median_by_cluster <- tapply(eomes_expr, cluster_vec, median)
      eomes_pct_pos_by_cluster <- tapply(eomes_expr, cluster_vec,
                                          function(x) mean(x > 0) * 100)
      n_cells_by_cluster <- table(cluster_vec)
      
      cluster_table <- data.table(
        cluster = names(eomes_mean_by_cluster),
        n_cells = as.integer(n_cells_by_cluster[names(eomes_mean_by_cluster)]),
        eomes_mean = round(as.numeric(eomes_mean_by_cluster), 4),
        eomes_median = round(as.numeric(eomes_median_by_cluster), 4),
        eomes_pct_pos = round(as.numeric(eomes_pct_pos_by_cluster), 1)
      )
      cluster_table <- cluster_table[order(-eomes_mean)]
      
      # 合并 density_category
      if (file.exists(profile_path)) {
        profile <- fread(profile_path)
        cluster_table <- merge(cluster_table,
                                profile[, .(cluster, density_median,
                                            density_category)],
                                by = "cluster", all.x = TRUE)
        cluster_table <- cluster_table[order(-eomes_mean)]
      }
      
      fwrite(cluster_table, file.path(OUTPUT_DIR, "partE_eomes_clusters",
                                       sprintf("%s_eomes_clusters.csv", parts[2])))
      
      eomes_cluster_results[[sname]] <- cluster_table
      
      log_msg(sprintf("    Top 3 EOMES+ clusters:"))
      for (i in 1:min(3, nrow(cluster_table))) {
        r <- cluster_table[i]
        dens_cat <- if ("density_category" %in% names(r)) r$density_category else NA
        log_msg(sprintf("      cluster %s: eomes_mean=%.3f, n=%d, density=%s",
                        r$cluster, r$eomes_mean, r$n_cells,
                        ifelse(is.na(dens_cat), "?", dens_cat)))
      }
    } else {
      log_msg(sprintf("  ⚠ %s panel 中无 EOMES - 跳过 Part E", sname))
    }
    
    # ---------- Part B: 空间图 ----------
    genes_for_this_sample <- targets_by_sample[[sname]]
    if (!is.null(genes_for_this_sample)) {
      log_msg(sprintf("  [Part B] 生成 %d 个基因的图...", length(genes_for_this_sample)))
      
      for (gene in genes_for_this_sample) {
        if (!(gene %in% rownames(sct_data))) {
          log_msg(sprintf("    ⚠ %s 不在 panel", gene))
          next
        }
        
        df <- data.frame(
          residual = as.numeric(sct_data[gene, ]),
          density = density_vec,
          x = x_coord,
          y = y_coord
        )
        
        # 空间图
        sp_path <- file.path(OUTPUT_DIR, "partB_spatial_plots",
                              sprintf("%s__%s__%s_spatial.png",
                                      parts[1], parts[2], gene))
        tryCatch({
          if (plot_spatial_fn(df, gene, sname, sp_path)) {
            log_msg(sprintf("    ✓ 空间图: %s", gene))
            partB_success_count <- partB_success_count + 1
          }
        }, error = function(e) {
          log_msg(sprintf("    ✗ 空间图 %s 失败: %s", gene, e$message))
        })
        
        # 散点图（额外，便宜）
        sc_path <- file.path(OUTPUT_DIR, "partB_spatial_plots",
                              sprintf("%s__%s__%s_scatter.png",
                                      parts[1], parts[2], gene))
        tryCatch({
          if (plot_scatter_fn(df, gene, sname, sc_path)) {
            log_msg(sprintf("    ✓ 散点图: %s", gene))
          }
        }, error = function(e) {
          log_msg(sprintf("    ✗ 散点图 %s 失败: %s", gene, e$message))
        })
        
        rm(df); gc(verbose = FALSE)
      }
    }
    
    rm(obj, sct_data, density_vec, cluster_vec, x_coord, y_coord)
    gc(verbose = FALSE)
    log_msg(sprintf("========== 完成 %s ==========\n", sname))
    
  }, error = function(e) {
    log_msg(sprintf("✗ %s 处理失败: %s", sname, e$message))
  })
}

# ---------- 汇总 Part E 到报告 ----------
if (length(eomes_cluster_results) > 0) {
  report_append("## Part E: EOMES+ cluster 识别 ✓\n\n")
  
  for (sname in names(eomes_cluster_results)) {
    tbl <- eomes_cluster_results[[sname]]
    parts <- strsplit(sname, "/")[[1]]
    
    report_append(sprintf("### %s\n\n", parts[2]))
    report_append("**按 EOMES 平均表达排序的 cluster（top 10）**：\n\n")
    
    has_density <- "density_category" %in% names(tbl)
    if (has_density) {
      report_append("| cluster | n_cells | EOMES_mean | EOMES_median | pct_pos | density_median | density_cat |\n")
      report_append("|---------|---------|------------|--------------|---------|----------------|-------------|\n")
    } else {
      report_append("| cluster | n_cells | EOMES_mean | EOMES_median | pct_pos |\n")
      report_append("|---------|---------|------------|--------------|---------|\n")
    }
    
    for (i in 1:min(10, nrow(tbl))) {
      r <- tbl[i]
      if (has_density) {
        report_append(sprintf("| %s | %d | %.3f | %.3f | %.1f%% | %s | %s |\n",
                              r$cluster, r$n_cells, r$eomes_mean,
                              r$eomes_median, r$eomes_pct_pos,
                              ifelse(is.na(r$density_median), "—",
                                     sprintf("%.2f", r$density_median)),
                              ifelse(is.na(r$density_category), "—",
                                     r$density_category)))
      } else {
        report_append(sprintf("| %s | %d | %.3f | %.3f | %.1f%% |\n",
                              r$cluster, r$n_cells, r$eomes_mean,
                              r$eomes_median, r$eomes_pct_pos))
      }
    }
    report_append("\n")
    
    # 关键观察：top EOMES+ cluster 的 density category
    if (has_density) {
      top3_cats <- tbl$density_category[1:min(3, nrow(tbl))]
      n_low_in_top3 <- sum(top3_cats == "low_density", na.rm = TRUE)
      report_append(sprintf("**观察**：top 3 EOMES+ cluster 中有 **%d/3** 个是 `low_density`\n",
                            n_low_in_top3))
      if (n_low_in_top3 >= 2) {
        report_append("→ 与 EOMES-密度负相关的全局模式一致 ✓\n\n")
      } else {
        report_append("→ 与 EOMES-密度负相关预期不完全一致，需要进一步解读\n\n")
      }
    }
  }
  
  task_status$partE <- "success"
} else {
  report_append("## Part E: EOMES+ cluster 识别 ✗\n\n**失败**: 没有成功的样本\n\n")
  task_status$partE <- "failed"
}

# ---------- 汇总 Part B 到报告 ----------
report_append(sprintf("## Part B: 空间图生成 (%d/%d 成功)\n\n",
                      partB_success_count, partB_total))
report_append("**生成的空间图清单**：\n\n")
for (t in SPATIAL_TARGETS) {
  parts <- strsplit(t[1], "/")[[1]]
  expected <- sprintf("%s__%s__%s_spatial.png", parts[1], parts[2], t[2])
  exists_now <- file.exists(file.path(OUTPUT_DIR, "partB_spatial_plots", expected))
  report_append(sprintf("- %s %s: `%s`\n",
                        ifelse(exists_now, "✓", "✗"),
                        t[1], expected))
}
report_append("\n")

task_status$partB <- ifelse(partB_success_count > 0, "success", "failed")

# ============================================================================
# Part G: 最终报告生成
# ============================================================================

log_section("Part G: 最终报告")

# 汇总
report_prefix <- paste0(
  "# Stage 3 批处理结果报告\n\n",
  sprintf("**运行完成时间**: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**输出目录**: `%s`\n\n", OUTPUT_DIR),
  "## 任务完成状态\n\n",
  "| Part | 任务 | 状态 |\n",
  "|------|------|------|\n",
  sprintf("| A | R7 effect_label bug 修复 | %s |\n",
          ifelse(is.null(task_status$partA), "未运行",
                 ifelse(task_status$partA == "success", "✓ 成功", "✗ 失败"))),
  sprintf("| B | MB 空间图生成 | %s |\n",
          ifelse(is.null(task_status$partB), "未运行",
                 ifelse(task_status$partB == "success", "✓ 成功", "✗ 失败"))),
  sprintf("| C | R8 关键数字汇总 | %s |\n",
          ifelse(is.null(task_status$partC), "未运行",
                 ifelse(task_status$partC == "success", "✓ 成功", "✗ 失败"))),
  sprintf("| D | 关键基因 panel 状态 | %s |\n",
          ifelse(is.null(task_status$partD), "未运行",
                 ifelse(task_status$partD == "success", "✓ 成功", "✗ 失败"))),
  sprintf("| E | EOMES+ cluster 识别 | %s |\n",
          ifelse(is.null(task_status$partE), "未运行",
                 ifelse(task_status$partE == "success", "✓ 成功", "✗ 失败"))),
  sprintf("| F | MB 47 增强 ρ 分解图 | %s |\n",
          ifelse(is.null(task_status$partF), "未运行",
                 ifelse(task_status$partF == "success", "✓ 成功", "✗ 失败"))),
  "\n---\n\n"
)

# 读取原报告内容，删除原 header 部分
existing_content <- readLines(REPORT_FILE)
existing_text <- paste(existing_content, collapse = "\n")

# 找 "## Part A" 的位置，以后才是真正的内容
part_a_pos <- regexpr("## Part A", existing_text, fixed = TRUE)
if (part_a_pos > 0) {
  body_text <- substring(existing_text, part_a_pos)
} else {
  body_text <- existing_text
}

final_content <- paste0(report_prefix, body_text)
cat(final_content, file = REPORT_FILE)

log_msg(sprintf("✓ 最终报告: %s", REPORT_FILE))
log_msg("========== 全部完成 ==========")
log_msg(sprintf("结束时间: %s", format(Sys.time())))

# 同时打印到 console
cat("\n\n========================================\n")
cat("Stage 3 批处理完成\n")
cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
cat(sprintf("报告文件: %s\n", REPORT_FILE))
cat(sprintf("日志文件: %s\n", LOG_FILE))
cat("========================================\n")
