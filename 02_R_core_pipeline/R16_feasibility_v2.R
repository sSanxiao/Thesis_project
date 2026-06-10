# ============================================================
# R16 Feasibility v2 — network-robust version
# ------------------------------------------------------------
# 避免 GEOquery 对本地文件的隐式网络调用 (v1 在 Step 1 就挂了)
# 做法:
# 1) Cavalli sample IDs 直接从 series matrix 文件 grep 出来
# 2) GSE124814 series matrix 用 curl 下载到本地再解析
# 3) 如果下载失败, 给出备用诊断信息
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(RESULTS_DIR, "R16_Feasibility")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CAVALLI_SERIES <- file.path(EXTDATA_DIR, "Cavalli_GSE85217", "GSE85217_series_matrix.txt.gz")
CAVALLI_CLIN <- file.path(EXTDATA_DIR, "Cavalli_GSE85217", "cavalli2017_mmc2_TableS1_clinical.csv")

GSE124_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/matrix/GSE124814_series_matrix.txt.gz"
GSE124_LOCAL <- file.path(OUT_DIR, "GSE124814_series_matrix.txt.gz")

cat("================================================================\n")
cat("R16 Feasibility v2 (network-robust)\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ============================================================
# Helper: extract !Key lines from series matrix .gz
# ============================================================
extract_series_keys <- function(gz_file, keys) {
  # keys: vector of keys like c("Sample_geo_accession", "Sample_title")
  # returns: named list of vectors (each vector = values in that line)
  
  result <- list()
  for (k in keys) {
    pattern <- sprintf("^!%s", k)
    # zgrep the line, split by tab
    line <- system(sprintf("zgrep -m1 '%s' %s", pattern, shQuote(gz_file)),
                   intern = TRUE, ignore.stderr = TRUE)
    if (length(line) == 0) {
      result[[k]] <- character(0)
      next
    }
    # parse: remove !Key_name, strip quotes, split by tab
    line <- sub(sprintf("^!%s\\s*", k), "", line[1])
    values <- strsplit(line, "\t")[[1]]
    values <- gsub('^"|"$', "", values)
    values <- trimws(values)
    values <- values[values != ""]
    result[[k]] <- values
  }
  result
}

extract_all_char_lines <- function(gz_file) {
  # 返回所有 !Sample_characteristics_ch1 行 (可能多行)
  lines <- system(sprintf("zgrep '^!Sample_characteristics' %s", shQuote(gz_file)),
                  intern = TRUE, ignore.stderr = TRUE)
  if (length(lines) == 0) return(list())
  
  parsed <- list()
  for (i in seq_along(lines)) {
    raw <- lines[i]
    # 去掉前缀
    raw <- sub("^!Sample_characteristics[^\t]*\t", "", raw)
    vals <- strsplit(raw, "\t")[[1]]
    vals <- gsub('^"|"$', "", vals)
    vals <- trimws(vals)
    parsed[[i]] <- vals
  }
  parsed
}

# ============================================================
# Step 1: Cavalli sample IDs (from local file, no network)
# ============================================================
cat("[Step 1] Extracting Cavalli sample IDs (local file only)...\n")

if (!file.exists(CAVALLI_SERIES)) {
  stop("Cavalli series matrix not found: ", CAVALLI_SERIES)
}

cav_keys <- extract_series_keys(CAVALLI_SERIES,
                                 c("Sample_geo_accession", "Sample_title"))

cav_samples <- data.table(
  GSM = cav_keys[["Sample_geo_accession"]],
  title = cav_keys[["Sample_title"]]
)
cat(sprintf("  Cavalli samples: %d\n", nrow(cav_samples)))
cat(sprintf("  First 3 GSMs: %s\n",
            paste(head(cav_samples$GSM, 3), collapse = ", ")))

# ============================================================
# Step 2: Download GSE124814 series matrix (via curl)
# ============================================================
cat("\n[Step 2] Downloading GSE124814 series matrix via curl...\n")
cat(sprintf("  URL: %s\n", GSE124_URL))
cat(sprintf("  Local: %s\n", GSE124_LOCAL))

if (file.exists(GSE124_LOCAL) && file.size(GSE124_LOCAL) > 1e6) {
  cat(sprintf("  Already exists (%.1f MB), skipping download.\n",
              file.size(GSE124_LOCAL) / 1e6))
} else {
  cat("  Downloading (may take 1-5 minutes)...\n")
  cmd <- sprintf("curl -L --max-time 600 --retry 3 -o %s %s",
                 shQuote(GSE124_LOCAL), shQuote(GSE124_URL))
  rc <- system(cmd)
  
  if (rc != 0 || !file.exists(GSE124_LOCAL) || file.size(GSE124_LOCAL) < 1e5) {
    cat("\n  [!] Download FAILED. Diagnostics:\n")
    system(sprintf("curl -sI --max-time 30 %s | head -5", shQuote(GSE124_URL)))
    cat("\n  Possible causes:\n")
    cat("    - Network interruption (retry later)\n")
    cat("    - GSE124814 filename might differ (check FTP directory)\n")
    cat(sprintf("    - FTP dir: %s\n",
                dirname(dirname(GSE124_URL))))
    quit(status = 1)
  }
  cat(sprintf("  Downloaded: %.1f MB\n", file.size(GSE124_LOCAL) / 1e6))
}

# ============================================================
# Step 3: Extract GSE124814 sample info (no GEOquery)
# ============================================================
cat("\n[Step 3] Parsing GSE124814 sample info...\n")

gse124_keys <- extract_series_keys(GSE124_LOCAL,
                                    c("Sample_geo_accession", "Sample_title",
                                      "Sample_source_name_ch1",
                                      "Sample_platform_id"))

n_124 <- length(gse124_keys[["Sample_geo_accession"]])
cat(sprintf("  GSE124814 samples: %d\n", n_124))
cat(sprintf("  Platform(s): %s\n",
            paste(unique(gse124_keys[["Sample_platform_id"]]), collapse = ", ")))

# characteristics_ch1 (通常有多行, 每行不同字段如 subgroup/os/age)
cat("\n  Parsing !Sample_characteristics_ch1 lines...\n")
char_lines_raw <- system(sprintf("zgrep '^!Sample_characteristics' %s",
                                  shQuote(GSE124_LOCAL)),
                          intern = TRUE, ignore.stderr = TRUE)
cat(sprintf("  Found %d characteristics lines\n", length(char_lines_raw)))

# 每行第一个非空值做 preview (看字段名)
if (length(char_lines_raw) > 0) {
  for (i in seq_along(char_lines_raw)) {
    raw <- char_lines_raw[i]
    raw <- sub("^!Sample_characteristics[^\t]*\t", "", raw)
    vals <- strsplit(raw, "\t")[[1]]
    vals <- gsub('^"|"$', "", vals)
    vals <- trimws(vals)
    vals <- vals[vals != ""]
    
    # 提取字段名 (key:value pattern)
    sample_val <- vals[1]
    if (grepl(":", sample_val)) {
      field_name <- trimws(sub(":.*$", "", sample_val))
      # 取前 5 个 unique values 看看
      all_vals <- unique(sub("^[^:]+:\\s*", "", vals))
      preview <- head(all_vals, 8)
      cat(sprintf("    Line %d [%s] (%d unique): %s\n",
                  i, field_name, length(all_vals),
                  paste(preview, collapse = " | ")))
    } else {
      cat(sprintf("    Line %d [no key]: %s\n",
                  i, paste(head(unique(vals), 5), collapse = " | ")))
    }
  }
}

# ============================================================
# Step 4: Overlap with Cavalli
# ============================================================
cat("\n[Step 4] Overlap analysis...\n")

gsm_124 <- gse124_keys[["Sample_geo_accession"]]
gsm_cav <- cav_samples$GSM

overlap_gsm <- intersect(gsm_124, gsm_cav)

cat(sprintf("  GSE124814 total: %d\n", length(gsm_124)))
cat(sprintf("  Cavalli total:   %d\n", length(gsm_cav)))
cat(sprintf("  Overlapping GSMs: %d (%.1f%% of GSE124814, %.1f%% of Cavalli)\n",
            length(overlap_gsm),
            100 * length(overlap_gsm) / length(gsm_124),
            100 * length(overlap_gsm) / length(gsm_cav)))
cat(sprintf("  Independent (non-Cavalli) GSMs in GSE124814: %d\n",
            length(gsm_124) - length(overlap_gsm)))

# ============================================================
# Step 5: Save outputs
# ============================================================
cat("\n[Step 5] Saving outputs...\n")

# GSE124814 sample summary
gse124_dt <- data.table(
  GSM = gse124_keys[["Sample_geo_accession"]],
  title = gse124_keys[["Sample_title"]],
  source = if (length(gse124_keys[["Sample_source_name_ch1"]]) == n_124)
             gse124_keys[["Sample_source_name_ch1"]] else NA,
  platform = if (length(gse124_keys[["Sample_platform_id"]]) == n_124)
               gse124_keys[["Sample_platform_id"]] else NA
)
gse124_dt[, in_cavalli := GSM %in% gsm_cav]
fwrite(gse124_dt, file.path(OUT_DIR, "GSE124814_samples.csv"))
cat(sprintf("  Saved: GSE124814_samples.csv\n"))

fwrite(data.table(GSM = overlap_gsm),
       file.path(OUT_DIR, "overlap_with_cavalli.csv"))
cat(sprintf("  Saved: overlap_with_cavalli.csv\n"))

# ============================================================
# Step 6: Go/No-Go judgment
# ============================================================
cat("\n================================================================\n")
cat("FEASIBILITY SUMMARY\n")
cat("================================================================\n\n")

indep_n <- length(gsm_124) - length(overlap_gsm)

cat(sprintf("  Q1. Can we fetch GSE124814? ............ YES (curl OK)\n"))
cat(sprintf("  Q2. Independent samples (non-Cavalli): %d\n", indep_n))
cat(sprintf("  Q3. Survival / subgroup fields:        see Step 3 above\n"))

go_no <- if (indep_n >= 200) "GO (proceed to R16 main)" else
         if (indep_n >= 100) "MARGINAL (proceed but underpowered)" else
         "NO-GO (use Northcott 2017 instead)"

cat(sprintf("\n  Verdict: %s\n", go_no))
cat("\n  Next: inspect Step 3 output for survival/subgroup fields.\n")
cat("        Paste full log to decide R16 design.\n")

cat("\n================================================================\n")
cat("R16 FEASIBILITY v2 DONE\n")
cat("================================================================\n")
