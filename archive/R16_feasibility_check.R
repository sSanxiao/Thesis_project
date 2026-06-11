# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R16 Feasibility Check: GSE124814 as independent validation
# ------------------------------------------------------------
# 只读 metadata, 不下载完整 expression data (先确认能不能用)
#
# 核查:
# 1) GSE124814 元数据能否获取
# 2) 样本总数 + 平台信息
# 3) 与 Cavalli GSE85217 的样本重叠
# 4) 有无生存数据字段
# 5) 有无分子亚型字段 (SHH/WNT/Group3/Group4)
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

OUT_DIR <- "./results/R16_Feasibility"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CAVALLI_CLIN <- "./external_data/Cavalli_GSE85217/cavalli2017_mmc2_TableS1_clinical.csv"
CAVALLI_SERIES <- "./external_data/Cavalli_GSE85217/GSE85217_series_matrix.txt.gz"

cat("================================================================\n")
cat("R16 Feasibility Check: GSE124814 vs Cavalli\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ------------------------------------------------------------
# Step 1: Load Cavalli sample IDs (for overlap check)
# ------------------------------------------------------------
cat("[Step 1] Loading Cavalli sample IDs...\n")

# 从 Cavalli series matrix 提取 GSM + Study_ID
cav_gse <- getGEO(filename = CAVALLI_SERIES, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(cav_gse)) cav_gse <- cav_gse[[1]]
cav_pdata <- pData(cav_gse)
cav_samples <- data.table(
  GSM = rownames(cav_pdata),
  Study_ID = cav_pdata$title,
  source = "Cavalli_GSE85217"
)
cat(sprintf("  Cavalli GSE85217 samples: %d\n", nrow(cav_samples)))

# 从 clinical csv 也拉一份 (含 Subgroup, OS)
cav_clin <- fread(CAVALLI_CLIN)
cav_clin[, .N, by = Subgroup]
cat(sprintf("  Cavalli clinical rows: %d with survival: %d\n",
            nrow(cav_clin), sum(!is.na(as.numeric(cav_clin$`OS (years)`)))))

# ------------------------------------------------------------
# Step 2: Try to fetch GSE124814 metadata (light-weight)
# ------------------------------------------------------------
cat("\n[Step 2] Fetching GSE124814 metadata (GEOquery, no GPL)...\n")

# 用 destdir 避免污染 HOME 的 cache
tmp_dir <- file.path(OUT_DIR, "geoquery_cache")
dir.create(tmp_dir, showWarnings = FALSE)

gse_124 <- tryCatch({
  getGEO("GSE124814", GSEMatrix = TRUE, getGPL = FALSE, destdir = tmp_dir)
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  NULL
})

if (is.null(gse_124)) {
  cat("\n  [!] GEOquery FAILED. Trying raw FTP...\n")
  # 备用: 直接 curl GEO FTP
  ftp_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/"
  cat(sprintf("  Checking FTP availability: %s\n", ftp_url))
  system(sprintf("curl -sI --max-time 30 %s | head -3", ftp_url))
  
  cat("\n  [?] If FTP works but GEOquery doesn't, we can download raw series matrix.\n")
  cat("      But first, try re-running this script — GEOquery sometimes\n")
  cat("      transient-fails. If still broken, paste this log.\n")
  quit(status = 1)
}

if (is.list(gse_124)) {
  cat(sprintf("  Returned list of %d elements (sub-platforms):\n", length(gse_124)))
  for (i in seq_along(gse_124)) {
    cat(sprintf("    [%d] %s: %d samples, platform %s\n",
                i, names(gse_124)[i], ncol(exprs(gse_124[[i]])),
                annotation(gse_124[[i]])))
  }
  # 主矩阵通常是第一个
  gse_main <- gse_124[[1]]
} else {
  gse_main <- gse_124
}

pdata_124 <- pData(gse_main)
cat(sprintf("\n  Main matrix dims: %d probes x %d samples\n",
            nrow(exprs(gse_main)), ncol(exprs(gse_main))))
cat(sprintf("  pData columns (first 15): %s\n",
            paste(head(names(pdata_124), 15), collapse = ", ")))

# ------------------------------------------------------------
# Step 3: Scan pData for key fields
# ------------------------------------------------------------
cat("\n[Step 3] Scanning pData for relevant fields...\n")

# 查看所有 characteristics_ch1.* 字段 (GEO 的亚型/临床通常在这里)
char_cols <- grep("characteristics|subgroup|subtype|survival|os|time|dead|event|age",
                  names(pdata_124), ignore.case = TRUE, value = TRUE)
cat(sprintf("  Potentially relevant fields (%d):\n", length(char_cols)))
for (cc in char_cols) {
  vals <- unique(pdata_124[[cc]])
  if (length(vals) <= 15) {
    cat(sprintf("    %s: %s\n", cc, paste(vals[1:min(10, length(vals))], collapse = " | ")))
  } else {
    cat(sprintf("    %s: [%d unique values, first 5:] %s\n", cc, length(vals),
                paste(head(vals, 5), collapse = " | ")))
  }
}

# ------------------------------------------------------------
# Step 4: Overlap check — GSE124814 vs GSE85217 GSMs
# ------------------------------------------------------------
cat("\n[Step 4] Sample overlap analysis...\n")

gsm_124 <- rownames(pdata_124)
gsm_cav <- cav_samples$GSM

overlap_gsm <- intersect(gsm_124, gsm_cav)
cat(sprintf("  GSE124814 total samples:  %d\n", length(gsm_124)))
cat(sprintf("  Cavalli samples:          %d\n", length(gsm_cav)))
cat(sprintf("  Overlap (same GSM IDs):   %d (%.1f%% of GSE124814)\n",
            length(overlap_gsm), 100 * length(overlap_gsm) / length(gsm_124)))
cat(sprintf("  Non-overlap (independent): %d\n",
            length(gsm_124) - length(overlap_gsm)))

# GSE124814 的关键在于它是 meta-analysis — 样本来自多个 series
# 看 pData 有没有 original_series 字段
series_cols <- grep("series|source|batch|dataset|cohort",
                    names(pdata_124), ignore.case = TRUE, value = TRUE)
cat("\n  Checking for source-series fields:\n")
for (sc in series_cols) {
  vals <- unique(pdata_124[[sc]])
  if (length(vals) <= 30) {
    cat(sprintf("    %s [%d unique]: %s\n", sc, length(vals),
                paste(vals, collapse = " | ")))
  }
}

# ------------------------------------------------------------
# Step 5: Summary + Go/No-Go
# ------------------------------------------------------------
cat("\n================================================================\n")
cat("FEASIBILITY SUMMARY\n")
cat("================================================================\n")

# 保存 pData 全量
fwrite(as.data.table(pdata_124, keep.rownames = "GSM"),
       file.path(OUT_DIR, "GSE124814_pdata.csv"))
cat(sprintf("  Saved: GSE124814_pdata.csv (%d samples)\n", nrow(pdata_124)))

# 保存重叠样本列表
fwrite(data.table(GSM = overlap_gsm),
       file.path(OUT_DIR, "overlap_with_cavalli.csv"))
cat(sprintf("  Saved: overlap_with_cavalli.csv\n"))

cat("\nKEY QUESTIONS FOR GO/NO-GO:\n")
cat(sprintf("  Q1. Can we fetch GSE124814? ......... %s\n",
            ifelse(!is.null(gse_124), "YES", "NO")))
cat(sprintf("  Q2. Independent samples (non-Cavalli)? %d\n",
            length(gsm_124) - length(overlap_gsm)))

# 初步判断有无生存数据
has_surv <- any(grepl("survival|os|time|dead|event",
                      unlist(lapply(pdata_124, as.character)),
                      ignore.case = TRUE))
cat(sprintf("  Q3. Survival data present? .......... %s (need manual inspection below)\n",
            ifelse(has_surv, "PROBABLY YES", "LIKELY NO")))

has_subgroup <- any(grepl("SHH|WNT|group.?3|group.?4",
                          unlist(lapply(pdata_124, as.character)),
                          ignore.case = TRUE))
cat(sprintf("  Q4. Subgroup labels present? ........ %s\n",
            ifelse(has_subgroup, "YES", "NO")))

cat("\n  >>> Next step: examine pData fields above to decide.\n")
cat("  >>> If Q2>=200 AND Q3=YES AND Q4=YES: R16 is GO.\n")
cat("  >>> If Q3 or Q4 NO: fall back to Northcott 2017 or Schwalbe 2017.\n")

cat("\n================================================================\n")
cat("R16 FEASIBILITY DONE\n")
cat("================================================================\n")
