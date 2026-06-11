# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# R10_precheck_part4.R
# 功能: 跑 Part 4 (GSE85217 元数据侦查)
# Part 1/2 已经从上一次跑的日志里拿到足够信息, 不再重复
#
# 运行: cd .
#       Rscript R10_precheck_part4.R 2>&1 | tee R10_precheck_part4.log
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

OUT_DIR      <- "./results/R10_Cavalli"
GEO_DATA_DIR <- "./external_data/Cavalli_GSE85217"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GEO_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

REPORT <- file.path(OUT_DIR, "PRECHECK_PART4_REPORT.txt")
sink(REPORT, split = TRUE)

cat("================================================================\n")
cat("  R10 PRECHECK Part 4: GSE85217 元数据侦查\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("拉 GSE85217 metadata (预计 2-5 分钟)...\n")
t0 <- Sys.time()

gse <- tryCatch({
  getGEO("GSE85217", destdir = GEO_DATA_DIR, getGPL = FALSE)
}, error = function(e) {
  cat("\n✗ 拉取失败:", conditionMessage(e), "\n"); NULL
})

if (is.null(gse)) { sink(); quit(save = "no", status = 1) }
cat(sprintf("完成, 耗时 %.1f 分钟\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Platform
if (is.list(gse) && !is(gse, "ExpressionSet")) {
  cat(sprintf("Platform 数: %d\n", length(gse)))
  for (i in seq_along(gse)) {
    cat(sprintf("  [%d] %s: %d genes × %d samples\n",
                i, names(gse)[i], nrow(gse[[i]]), ncol(gse[[i]])))
  }
  eset <- gse[[1]]
} else {
  eset <- gse
}
cat(sprintf("\n使用 eset: %d genes × %d samples\n", nrow(eset), ncol(eset)))

# ---- 4.1 pData ----
cat("\n----------------------------------------------------------------\n")
cat("4.1 Phenotype data (pData)\n")
cat("----------------------------------------------------------------\n")
pd <- pData(eset)
cat(sprintf("维度: %d 样本 × %d 字段\n\n", nrow(pd), ncol(pd)))

cat("所有列名:\n")
for (i in seq_along(names(pd))) {
  cat(sprintf("  [%3d] %s\n", i, names(pd)[i]))
}

cat("\n关键列枚举:\n")
interesting <- grep("subgroup|subtype|subclass|group|cluster|class|tumor|histology|survival|os|overall|dead|death|status|vital|event|time|month|year|age|sex|gender|pathology",
                    names(pd), value = TRUE, ignore.case = TRUE)
for (col in interesting) {
  vals <- pd[[col]]
  uniq <- unique(vals)
  cat(sprintf("\n[%s]\n", col))
  cat(sprintf("  n_unique=%d, n_missing=%d\n",
              length(uniq), sum(is.na(vals) | vals == "")))
  if (length(uniq) <= 25) {
    print(table(vals, useNA = "ifany"))
  } else {
    cat("  first 20 unique values:\n")
    print(head(uniq, 20))
    cat(sprintf("  ... (%d more)\n", length(uniq) - 20))
  }
}

cat("\n\n前 3 个样本完整 phenotype (转置):\n")
print(t(head(pd, 3)))

# ---- 4.2 fData ----
cat("\n\n----------------------------------------------------------------\n")
cat("4.2 Feature data (fData) — probe to gene symbol 映射\n")
cat("----------------------------------------------------------------\n")
fd <- fData(eset)
cat(sprintf("维度: %d probes × %d 字段\n", nrow(fd), ncol(fd)))

if (ncol(fd) > 0) {
  cat("\n所有列名:\n")
  for (i in seq_along(names(fd))) {
    cat(sprintf("  [%3d] %s\n", i, names(fd)[i]))
  }
  
  cat("\n前 5 行:\n")
  print(head(fd, 5))
  
  sym_candidates <- grep("symbol|gene.name|gene_symbol|hgnc|geneassignment|gene_assignment",
                         names(fd), value = TRUE, ignore.case = TRUE)
  if (length(sym_candidates) > 0) {
    cat(sprintf("\n可能的 symbol 列: %s\n\n",
                paste(sym_candidates, collapse = ", ")))
    for (sc in sym_candidates) {
      n_nonempty <- sum(!is.na(fd[[sc]]) & fd[[sc]] != "" & fd[[sc]] != "---")
      cat(sprintf("[%s] 非空: %d/%d (%.1f%%)\n",
                  sc, n_nonempty, nrow(fd), 100*n_nonempty/nrow(fd)))
      sample_vals <- head(fd[[sc]][!is.na(fd[[sc]]) & fd[[sc]] != ""], 10)
      cat("  前 10 个值 (看是否是 gene symbol):\n")
      for (v in sample_vals) cat(sprintf("    %s\n", v))
    }
  } else {
    cat("\n⚠ 未探到 symbol 列\n")
  }
} else {
  cat("\n⚠ fData 为空\n")
}

# ---- 4.3 Expression matrix ----
cat("\n\n----------------------------------------------------------------\n")
cat("4.3 表达矩阵\n")
cat("----------------------------------------------------------------\n")
expr <- exprs(eset)
cat(sprintf("维度: %d × %d\n", nrow(expr), ncol(expr)))
cat(sprintf("值范围: [%.2f, %.2f]\n",
            min(expr, na.rm=TRUE), max(expr, na.rm=TRUE)))
cat(sprintf("中位数: %.2f\n", median(expr, na.rm=TRUE)))
cat(sprintf("前 10 个 rownames:\n"))
for (r in head(rownames(expr), 10)) cat(sprintf("  %s\n", r))
cat(sprintf("左上角 3x3:\n"))
print(expr[1:3, 1:3])

# ---- 4.4 关键基因匹配 ----
cat("\n\n----------------------------------------------------------------\n")
cat("4.4 MB signature 先验基因在 bulk 里的匹配\n")
cat("----------------------------------------------------------------\n")

# 从 R7 一致性 + MB266 强信号合起来的扩展集
prior_genes <- c(
  # R7 跨样本一致的 15 个
  "RBFOX3", "TUBB4A", "CCR7", "EOMES", "ST18",
  "SLC17A7", "SV2B", "NRXN3", "DCN", "APOE",
  "TENM1", "CXCR4", "GLI2", "BOC", "CRYM",
  # MB266 独有强信号 (|rho|>0.3)
  "HHIP", "THSD7B", "PLCH1", "MYO16", "NES", "NNAT",
  "COL25A1", "CDH12", "C1QL3", "SOX9",
  # GSM8840046 tier1
  "POU5F1B", "IGFBP5",
  # GSM8840049 tier1
  "DNER", "SPHKAP", "IGFBP4", "TLL1", "RELN", "CA10", "SOX11"
)
prior_genes <- unique(prior_genes)
cat(sprintf("先验基因总数: %d\n", length(prior_genes)))

# 按 rownames 直接匹配
direct_match <- intersect(prior_genes, rownames(expr))
cat(sprintf("\n按 rownames 直接匹配: %d / %d\n",
            length(direct_match), length(prior_genes)))
cat("  匹配的:", paste(direct_match, collapse=", "), "\n")

# 按 fData 里每个 symbol 候选列匹配
if (ncol(fd) > 0) {
  sym_cands <- grep("symbol|gene.name|gene_symbol|hgnc|gene_assignment|geneassignment",
                    names(fd), value = TRUE, ignore.case = TRUE)
  for (sc in sym_cands) {
    m <- intersect(prior_genes, fd[[sc]])
    cat(sprintf("\n按 fData[[%s]] 匹配: %d / %d\n",
                sc, length(m), length(prior_genes)))
    if (length(m) > 0) {
      cat("  匹配的:", paste(m, collapse=", "), "\n")
    }
  }
}

# ---- 4.5 Platform ----
cat("\n\n----------------------------------------------------------------\n")
cat("4.5 Platform / Annotation\n")
cat("----------------------------------------------------------------\n")
if (length(annotation(eset)) > 0) cat("annotation:", annotation(eset), "\n")
exp_info <- experimentData(eset)
if (length(exp_info@title) > 0) cat("Title:", exp_info@title, "\n")
if (length(exp_info@url) > 0) cat("URL:", exp_info@url, "\n")

cat("\n================================================================\n")
cat("Part 4 完成\n")
cat("================================================================\n")

sink()
cat("\n报告已保存:", REPORT, "\n")
