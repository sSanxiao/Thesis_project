# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# R10_precheck.R
# 功能: 在写 R10 主脚本前, 侦查所有字段和数据结构
# 不做任何分析, 不下载表达矩阵, 只拉元数据
#
# 输出: ./results/R10_Cavalli/PRECHECK_REPORT.txt
#        把这个文件贴给 Claude, 再写 R10 主脚本
#
# 运行: cd .
#       Rscript R10_precheck.R 2>&1 | tee R10_precheck.log
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

R4_MB_DIR    <- "./results/R4_Results/Medulloblastoma_Human"
R7_DIR       <- "./results/R7_Results"
OUT_DIR      <- "./results/R10_Cavalli"
GEO_DATA_DIR <- "./external_data/Cavalli_GSE85217"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GEO_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

REPORT <- file.path(OUT_DIR, "PRECHECK_REPORT.txt")
sink(REPORT, split = TRUE)   # 同时写文件和屏幕

cat("================================================================\n")
cat("  R10 PRECHECK REPORT\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")

# ============================================================================
# Part 1: R4 的 4 个 MB 样本 filtered_density_genes.csv 侦查
# ============================================================================

cat("================================================================\n")
cat("Part 1: R4 MB 样本 filtered_density_genes.csv 结构\n")
cat("================================================================\n\n")

mb_sample_dirs <- list.dirs(R4_MB_DIR, recursive = FALSE)
cat("MB 样本目录数:", length(mb_sample_dirs), "\n")
for (d in mb_sample_dirs) cat("  -", basename(d), "\n")
cat("\n")

for (d in mb_sample_dirs) {
  sname <- basename(d)
  f <- file.path(d, "filtered_density_genes.csv")
  
  cat("----------------------------------------------------------------\n")
  cat(sprintf("[%s]\n", sname))
  cat("----------------------------------------------------------------\n")
  
  if (!file.exists(f)) {
    cat("  ✗ 文件不存在\n\n"); next
  }
  
  df <- fread(f)
  cat(sprintf("  文件: %s\n", f))
  cat(sprintf("  行数: %d, 列数: %d\n", nrow(df), ncol(df)))
  cat(sprintf("  所有列名:\n"))
  cat(sprintf("    %s\n", paste(names(df), collapse = ", ")))
  
  # tier 列的所有值和频次
  if ("tier" %in% names(df)) {
    cat("\n  tier 列的值频次:\n")
    print(table(df$tier, useNA = "ifany"))
  }
  
  # direction 列的所有值
  if ("direction" %in% names(df)) {
    cat("\n  direction 列的值频次:\n")
    print(table(df$direction, useNA = "ifany"))
  }
  
  # convergence 列 (如果有)
  if ("convergence" %in% names(df)) {
    cat("\n  convergence 列的值频次:\n")
    print(table(df$convergence, useNA = "ifany"))
  }
  
  # rho_knn_main 列的基本统计
  if ("rho_knn_main" %in% names(df)) {
    cat("\n  rho_knn_main 统计:\n")
    cat(sprintf("    range: [%.4f, %.4f]\n",
                min(df$rho_knn_main, na.rm=TRUE),
                max(df$rho_knn_main, na.rm=TRUE)))
    cat(sprintf("    median |rho|: %.4f\n",
                median(abs(df$rho_knn_main), na.rm=TRUE)))
  }
  
  # tier1_strict 基因的完整列表和 rho
  if ("tier" %in% names(df) && "tier1_strict" %in% df$tier) {
    t1 <- df[tier == "tier1_strict"]
    cat(sprintf("\n  tier1_strict 基因 (%d 个):\n", nrow(t1)))
    if (nrow(t1) > 0) {
      show_cols <- intersect(c("gene", "rho_knn_main", "direction", "convergence"),
                             names(t1))
      print(head(t1[, ..show_cols], 20))
      if (nrow(t1) > 20) cat(sprintf("    ... (还有 %d 个)\n", nrow(t1) - 20))
    }
  }
  
  cat("\n")
}

# ============================================================================
# Part 2: R7 的数据集一致性表侦查 (用它看交集)
# ============================================================================

cat("================================================================\n")
cat("Part 2: R7 MB 数据集一致性表\n")
cat("================================================================\n\n")

r7_consist <- file.path(R7_DIR, "ALL_DATASETS_R7_CONSISTENCY.csv")
if (file.exists(r7_consist)) {
  r7 <- fread(r7_consist)
  cat("列名:", paste(names(r7), collapse = ", "), "\n\n")
  cat("MB 行:\n")
  mb_row <- r7[dataset == "Medulloblastoma_Human"]
  print(mb_row)
  if (nrow(mb_row) > 0 && "consistent_genes" %in% names(mb_row)) {
    cat("\nconsistent_genes 完整内容:\n")
    cat(mb_row$consistent_genes[1], "\n")
  }
} else {
  cat("⚠ 文件不存在:", r7_consist, "\n")
}

# ============================================================================
# Part 3: R7 ALL_SAMPLES_R7_PROFILE.csv (另一个可能的数据源)
# ============================================================================

cat("\n================================================================\n")
cat("Part 3: R7 ALL_SAMPLES_R7_PROFILE.csv (如果存在)\n")
cat("================================================================\n\n")

r7_profile <- file.path(R7_DIR, "ALL_SAMPLES_R7_PROFILE.csv")
if (file.exists(r7_profile)) {
  p <- fread(r7_profile)
  cat("列名:", paste(names(p), collapse = ", "), "\n")
  cat("MB 行:\n")
  print(p[grepl("Medulloblastoma", sample)])
} else {
  cat("文件不存在\n")
}

# ============================================================================
# Part 4: GSE85217 元数据侦查 (不下表达矩阵, 只看结构)
# ============================================================================

cat("\n================================================================\n")
cat("Part 4: GSE85217 (Cavalli 2017) 元数据侦查\n")
cat("================================================================\n\n")

cat("开始拉 GSE85217 metadata (不下表达矩阵, 预计 1-3 分钟)...\n")
t0 <- Sys.time()

# 关键: GSEMatrix = TRUE 但指定 destdir 本地缓存后续会重用
# 如果只要元数据, 可以用 GSEMatrix=FALSE (更快, 但不给表达矩阵结构)
# 这里用 GSEMatrix=TRUE 但只看元数据部分, 表达矩阵拉下来但不用
gse <- tryCatch({
  getGEO("GSE85217", destdir = GEO_DATA_DIR, getGPL = FALSE)
}, error = function(e) {
  cat("\n✗ 拉取失败:", conditionMessage(e), "\n"); NULL
})

if (is.null(gse)) {
  cat("\n无法继续 Part 4, 请检查网络\n")
  sink(); quit(save = "no", status = 1)
}

cat(sprintf("完成, 耗时 %.1f 分钟\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# 可能返回 list (多 platform) 或单个 eset
if (is.list(gse) && !is(gse, "ExpressionSet")) {
  cat(sprintf("Platform 数: %d\n", length(gse)))
  for (i in seq_along(gse)) {
    cat(sprintf("  [%d] %s: %d genes × %d samples\n",
                i, names(gse)[i], nrow(gse[[i]]), ncol(gse[[i]])))
  }
  eset <- gse[[1]]  # 取第一个 platform
} else {
  eset <- gse
  cat(sprintf("单 platform: %d genes × %d samples\n", nrow(eset), ncol(eset)))
}

# ---- 4.1 Phenotype data ----
cat("\n----------------------------------------------------------------\n")
cat("4.1 Phenotype data (pData)\n")
cat("----------------------------------------------------------------\n")
pd <- pData(eset)
cat(sprintf("  维度: %d 样本 × %d 字段\n\n", nrow(pd), ncol(pd)))
cat("  所有列名:\n")
for (i in seq_along(names(pd))) {
  cat(sprintf("    [%3d] %s\n", i, names(pd)[i]))
}

cat("\n  关键列的唯一值 (前 20):\n")
interesting <- grep("subgroup|subtype|subclass|group|cluster|class|tumor|histology|survival|os|overall|dead|death|status|vital|event|time|month|year|age|sex|gender",
                    names(pd), value = TRUE, ignore.case = TRUE)
for (col in interesting) {
  vals <- pd[[col]]
  uniq <- unique(vals)
  cat(sprintf("\n    [%s]\n", col))
  cat(sprintf("      n_unique = %d, n_missing = %d\n",
              length(uniq), sum(is.na(vals) | vals == "")))
  if (length(uniq) <= 20) {
    cat("      values:\n")
    print(table(vals, useNA = "ifany"))
  } else {
    cat("      first 20 unique:\n")
    print(head(uniq, 20))
    cat(sprintf("      ... (%d more)\n", length(uniq) - 20))
  }
}

# 前 3 个样本的完整 phenotype
cat("\n  前 3 个样本完整 phenotype (转置显示):\n")
print(t(head(pd, 3)))

# ---- 4.2 Feature data ----
cat("\n----------------------------------------------------------------\n")
cat("4.2 Feature data (fData) — probe to gene symbol 映射\n")
cat("----------------------------------------------------------------\n")
fd <- fData(eset)
cat(sprintf("  维度: %d probes × %d 字段\n", nrow(fd), ncol(fd)))

if (ncol(fd) > 0) {
  cat("\n  所有列名:\n")
  for (i in seq_along(names(fd))) {
    cat(sprintf("    [%3d] %s\n", i, names(fd)[i]))
  }
  
  cat("\n  前 5 行:\n")
  print(head(fd, 5))
  
  # 尝试探测 symbol 列
  sym_candidates <- grep("symbol|gene.name|gene_symbol|hgnc",
                         names(fd), value = TRUE, ignore.case = TRUE)
  if (length(sym_candidates) > 0) {
    cat(sprintf("\n  可能的 symbol 列: %s\n",
                paste(sym_candidates, collapse = ", ")))
    for (sc in sym_candidates) {
      n_nonempty <- sum(!is.na(fd[[sc]]) & fd[[sc]] != "")
      cat(sprintf("    [%s] 非空: %d/%d, 前 10 个: %s\n",
                  sc, n_nonempty, nrow(fd),
                  paste(head(fd[[sc]][!is.na(fd[[sc]]) & fd[[sc]] != ""], 10),
                        collapse = ", ")))
    }
  } else {
    cat("\n  ⚠ 未探到 symbol 列, 可能需要单独下 GPL 注释\n")
  }
} else {
  cat("\n  ⚠ fData 为空, 需要通过 GPL 平台文件获取 probe->symbol 映射\n")
}

# ---- 4.3 Expression matrix 维度和样例 ----
cat("\n----------------------------------------------------------------\n")
cat("4.3 表达矩阵基本信息\n")
cat("----------------------------------------------------------------\n")
expr <- exprs(eset)
cat(sprintf("  维度: %d × %d\n", nrow(expr), ncol(expr)))
cat(sprintf("  值范围: [%.2f, %.2f]\n",
            min(expr, na.rm=TRUE), max(expr, na.rm=TRUE)))
cat(sprintf("  中位数: %.2f\n", median(expr, na.rm=TRUE)))
cat("  (如果是 log2 转换的, 典型范围是 2-14)\n")
cat(sprintf("  前 5 个 rownames (probe ID 或 gene symbol?):\n"))
cat("    ", paste(head(rownames(expr), 5), collapse = ", "), "\n")
cat(sprintf("  前 3 个 colnames:\n"))
cat("    ", paste(head(colnames(expr), 3), collapse = ", "), "\n")
cat(sprintf("  左上角 3x3:\n"))
print(expr[1:3, 1:3])

# ---- 4.4 关键验证: 我打算用的 signature 基因能不能匹配 ----
cat("\n----------------------------------------------------------------\n")
cat("4.4 先验 MB 基因在 bulk 里的匹配情况\n")
cat("----------------------------------------------------------------\n")
# 从 R7 日志里我已经知道 MB consistent 基因列表前 10:
# RBFOX3, TUBB4A, CCR7, EOMES, ST18, SLC17A7, SV2B, NRXN3, DCN, APOE
prior_genes <- c("TUBB4A", "CCR7", "EOMES", "RBFOX3", "ST18",
                 "SLC17A7", "SV2B", "NRXN3", "DCN", "APOE",
                 "TENM1", "CXCR4", "GLI2", "BOC", "CRYM",
                 "HHIP", "OTX2", "PTCH1", "MYC")

# 在 rownames (probe ID) 里找
direct_match <- intersect(prior_genes, rownames(expr))
cat(sprintf("  直接按 rownames 匹配 (%d/%d): %s\n",
            length(direct_match), length(prior_genes),
            paste(direct_match, collapse = ", ")))

# 如果 fData 有 symbol 列, 也试一下
if (ncol(fd) > 0) {
  sym_cands <- grep("symbol|gene.name|gene_symbol",
                    names(fd), value = TRUE, ignore.case = TRUE)
  for (sc in sym_cands) {
    m <- intersect(prior_genes, fd[[sc]])
    cat(sprintf("  按 fData[[%s]] 匹配 (%d/%d): %s\n",
                sc, length(m), length(prior_genes),
                paste(m, collapse = ", ")))
  }
}

# ---- 4.5 Platform 信息 ----
cat("\n----------------------------------------------------------------\n")
cat("4.5 Platform / Annotation 信息\n")
cat("----------------------------------------------------------------\n")
if ("annotation" %in% slotNames(eset)) {
  cat("  annotation slot:", annotation(eset), "\n")
}
exp_info <- experimentData(eset)
cat("  Title:", exp_info@title, "\n")
cat("  URL:", exp_info@url, "\n")

# ============================================================================
# Part 5: 磁盘空间最后确认
# ============================================================================

cat("\n================================================================\n")
cat("Part 5: 磁盘空间\n")
cat("================================================================\n\n")
system("df -h .")
system(sprintf("du -sh %s", GEO_DATA_DIR))

cat("\n================================================================\n")
cat("PRECHECK 完成\n")
cat("================================================================\n")

sink()

cat("\n报告已保存: ", REPORT, "\n")
cat("把这个文件的完整内容贴给 Claude, 然后写最终 R10 脚本\n")
