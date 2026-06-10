# ============================================================
# R16a: GSE124814 Expression Matrix — Download & Inspect
# ------------------------------------------------------------
# 目标:
# 1) 下载 109MB HW_expr_matrix
# 2) 探查矩阵结构 (行是 probe 还是 symbol? log2 还是 raw?)
# 3) 验证 sig_94 / sig_core 基因能否匹配到矩阵
# 4) 输出诊断报告, 不做分析
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(EXTDATA_DIR, "GSE124814")
RES_DIR <- file.path(RESULTS_DIR, "R16_GSE124814")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

EXPR_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/suppl/GSE124814_HW_expr_matrix.tsv.gz"
EXPR_LOCAL <- file.path(OUT_DIR, "GSE124814_HW_expr_matrix.tsv.gz")
SAMPLE_CSV <- file.path(OUT_DIR, "sample_descriptions.csv")
SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")

cat("================================================================\n")
cat("R16a: GSE124814 expression matrix inspection\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

# ------------------------------------------------------------
# Step 1: Download 109MB expression matrix
# ------------------------------------------------------------
cat("[Step 1] Downloading GSE124814_HW_expr_matrix.tsv.gz (~109 MB)...\n")

if (file.exists(EXPR_LOCAL) && file.size(EXPR_LOCAL) > 1e8) {
  cat(sprintf("  Already exists (%.1f MB), skipping.\n",
              file.size(EXPR_LOCAL) / 1e6))
} else {
  cat("  Downloading... may take 3-20 minutes depending on network.\n")
  t0 <- Sys.time()
  cmd <- sprintf("curl -L --max-time 1800 --retry 3 -o %s %s",
                 shQuote(EXPR_LOCAL), shQuote(EXPR_URL))
  rc <- system(cmd)
  t1 <- Sys.time()
  
  if (rc != 0 || !file.exists(EXPR_LOCAL) || file.size(EXPR_LOCAL) < 1e7) {
    cat("\n  [!] DOWNLOAD FAILED.\n")
    cat(sprintf("  Elapsed: %.1f min\n", as.numeric(t1 - t0, units = "mins")))
    cat("  Size: ", file.size(EXPR_LOCAL), " bytes\n")
    quit(status = 1)
  }
  cat(sprintf("  Downloaded: %.1f MB in %.1f min\n",
              file.size(EXPR_LOCAL) / 1e6,
              as.numeric(t1 - t0, units = "mins")))
}

# ------------------------------------------------------------
# Step 2: Inspect structure (top rows, dim, data ranges)
# ------------------------------------------------------------
cat("\n[Step 2] Inspecting matrix structure...\n")

# 看前 3 行, 先判断 header + 行名结构
cat("\n  First 3 lines (peek):\n")
peek <- system(sprintf("zcat %s | head -3 | cut -c1-300", EXPR_LOCAL),
               intern = TRUE)
for (l in peek) {
  # 截断显示避免太长
  cat(sprintf("    %s%s\n", substr(l, 1, 250),
              if (nchar(l) > 250) "..." else ""))
}

# 看 header (第 1 行) 的列数和前几个列名
hdr_line <- system(sprintf("zcat %s | head -1", EXPR_LOCAL), intern = TRUE)[1]
header <- strsplit(hdr_line, "\t")[[1]]
cat(sprintf("\n  Header: %d columns\n", length(header)))
cat(sprintf("  First 5 col names: %s\n", paste(head(header, 5), collapse = " | ")))
cat(sprintf("  Last 3 col names: %s\n", paste(tail(header, 3), collapse = " | ")))

# 行数 (不含 header)
cat("\n  Counting total lines...\n")
nline <- as.integer(system(sprintf("zcat %s | wc -l", EXPR_LOCAL), intern = TRUE))
cat(sprintf("  Total lines: %d (= %d rows + 1 header)\n", nline, nline - 1))

# ------------------------------------------------------------
# Step 3: Read full matrix and inspect
# ------------------------------------------------------------
cat("\n[Step 3] Loading full matrix (this may take 1-2 min)...\n")
t0 <- Sys.time()
expr_dt <- fread(cmd = sprintf("zcat %s", EXPR_LOCAL),
                 sep = "\t", header = TRUE, data.table = TRUE)
t1 <- Sys.time()
cat(sprintf("  Loaded in %.1f sec\n", as.numeric(t1 - t0, units = "secs")))
cat(sprintf("  Dimensions: %d rows x %d cols\n", nrow(expr_dt), ncol(expr_dt)))

# 第一列通常是 gene/probe ID
id_col <- names(expr_dt)[1]
cat(sprintf("  ID column: '%s'\n", id_col))
cat(sprintf("  First 5 IDs: %s\n",
            paste(head(expr_dt[[id_col]], 5), collapse = " | ")))
cat(sprintf("  Last 3 IDs: %s\n",
            paste(tail(expr_dt[[id_col]], 3), collapse = " | ")))

# 检查 ID 类型: 是 ENSG / SYMBOL / probe ID?
sample_ids <- head(expr_dt[[id_col]], 20)
id_type <- "UNKNOWN"
if (any(grepl("^ENSG", sample_ids))) id_type <- "ENSEMBL"
if (all(sample_ids %in% c("TP53", "MYC", "APP", "APOE", "NES", "GFAP",
                           "MBP", "OLIG2", "TUBB3"))) id_type <- "SYMBOL"
if (any(grepl("^[0-9]+_at$", sample_ids))) id_type <- "AFFY_PROBE"
cat(sprintf("  Inferred ID type: %s\n", id_type))

# 数据列 (排除第一列)
data_cols <- setdiff(names(expr_dt), id_col)
n_samples <- length(data_cols)
cat(sprintf("  Sample columns: %d\n", n_samples))
cat(sprintf("  First 3 sample IDs: %s\n",
            paste(head(data_cols, 3), collapse = " | ")))

# 数据范围 (从前 1000 行取 sample)
cat("\n  Checking data value ranges...\n")
sample_data <- as.matrix(expr_dt[1:min(1000, nrow(expr_dt)),
                                  data_cols[1:min(100, length(data_cols))],
                                  with = FALSE])
cat(sprintf("  Sample (1000 rows x 100 cols) summary:\n"))
cat(sprintf("    min:    %.4f\n", min(sample_data, na.rm = TRUE)))
cat(sprintf("    max:    %.4f\n", max(sample_data, na.rm = TRUE)))
cat(sprintf("    mean:   %.4f\n", mean(sample_data, na.rm = TRUE)))
cat(sprintf("    median: %.4f\n", median(sample_data, na.rm = TRUE)))
cat(sprintf("    NA fraction: %.2f%%\n",
            100 * sum(is.na(sample_data)) / length(sample_data)))

# 判断 log scale
is_log_scale <- max(sample_data, na.rm = TRUE) < 50
cat(sprintf("  Likely log-scale: %s (max=%.2f)\n",
            ifelse(is_log_scale, "YES", "NO"),
            max(sample_data, na.rm = TRUE)))

# ------------------------------------------------------------
# Step 4: Load signature + sample descriptions
# ------------------------------------------------------------
cat("\n[Step 4] Loading signature and sample descriptions...\n")

sig_prov <- fread(SIG_PROV)
sig_99 <- sig_prov$gene
sig_94 <- setdiff(sig_99, CONFLICT_GENES)
sig_core <- setdiff(sig_prov[n_samples >= 2, gene], CONFLICT_GENES)

cat(sprintf("  sig_99: %d genes\n", length(sig_99)))
cat(sprintf("  sig_94 (cleaned): %d genes\n", length(sig_94)))
cat(sprintf("  sig_core: %d genes (%s)\n",
            length(sig_core), paste(sig_core, collapse = ", ")))

samples_dt <- fread(SAMPLE_CSV)
cat(sprintf("  Samples from xlsx: %d rows\n", nrow(samples_dt)))

# ------------------------------------------------------------
# Step 5: Try matching signature genes to matrix IDs
# ------------------------------------------------------------
cat("\n[Step 5] Matching signature genes to matrix IDs...\n")

# 先直接尝试 symbol 匹配
matrix_ids <- expr_dt[[id_col]]
direct_match_94 <- intersect(sig_94, matrix_ids)
direct_match_core <- intersect(sig_core, matrix_ids)

cat(sprintf("  Direct symbol match sig_94: %d/%d (%.1f%%)\n",
            length(direct_match_94), length(sig_94),
            100 * length(direct_match_94) / length(sig_94)))
cat(sprintf("  Direct symbol match sig_core: %d/%d\n",
            length(direct_match_core), length(sig_core)))

# 如果 ID 是 ENSG, 尝试 ENSG -> SYMBOL
ensembl_match_n <- 0
ensembl_tried <- FALSE
if (id_type == "ENSEMBL" || any(grepl("^ENSG", matrix_ids[1:100]))) {
  ensembl_tried <- TRUE
  cat("\n  IDs look like ENSG — trying ENSG -> SYMBOL mapping...\n")
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  
  clean_ids <- gsub("\\..*$", "", matrix_ids)  # strip version
  sym_map <- mapIds(org.Hs.eg.db,
                    keys = clean_ids,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")
  
  mapped_syms <- as.character(sym_map)
  names(mapped_syms) <- matrix_ids
  
  ensembl_match_94 <- intersect(sig_94, mapped_syms)
  ensembl_match_core <- intersect(sig_core, mapped_syms)
  ensembl_match_n <- length(ensembl_match_94)
  
  cat(sprintf("  ENSG mapped to %d unique symbols\n",
              length(unique(na.omit(mapped_syms)))))
  cat(sprintf("  sig_94 match via ENSG: %d/%d\n",
              length(ensembl_match_94), length(sig_94)))
  cat(sprintf("  sig_core match via ENSG: %d/%d\n",
              length(ensembl_match_core), length(sig_core)))
}

# 显示 sig_core 里每个基因的匹配状态 (因为只有 8 个, 好看)
cat("\n  sig_core gene-by-gene match check:\n")
for (g in sig_core) {
  direct <- g %in% matrix_ids
  marker <- if (direct) "OK-direct" else "MISS"
  cat(sprintf("    %-10s %s\n", g, marker))
}

# ------------------------------------------------------------
# Step 6: Cross-check with xlsx sample IDs
# ------------------------------------------------------------
cat("\n[Step 6] Cross-checking matrix sample IDs with xlsx metadata...\n")

# xlsx 里 "Sample name" 列 (Sample_1, Sample_2, ...)
xlsx_sample_ids <- samples_dt$`Sample name`
cat(sprintf("  xlsx sample IDs (first 3): %s\n",
            paste(head(xlsx_sample_ids, 3), collapse = " | ")))
cat(sprintf("  matrix sample IDs (first 3): %s\n",
            paste(head(data_cols, 3), collapse = " | ")))

overlap_samples <- intersect(data_cols, xlsx_sample_ids)
cat(sprintf("  Direct overlap: %d / %d (matrix) / %d (xlsx)\n",
            length(overlap_samples), n_samples, length(xlsx_sample_ids)))

# 如果 overlap 低, 可能 matrix 用别的 ID (比如 GSM, title)
if (length(overlap_samples) < 100) {
  cat("\n  Low overlap — trying alt IDs from xlsx...\n")
  for (col in names(samples_dt)) {
    vals <- as.character(samples_dt[[col]])
    m <- length(intersect(data_cols, vals))
    if (m > 100) {
      cat(sprintf("    Col '%s' matches %d matrix columns!\n", col, m))
    }
  }
}

# ------------------------------------------------------------
# Step 7: Save inspection report
# ------------------------------------------------------------
cat("\n[Step 7] Saving inspection report...\n")

report_lines <- c(
  "================================================================",
  "R16a: GSE124814 Expression Matrix Inspection Report",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "FILE",
  sprintf("  Path: %s", EXPR_LOCAL),
  sprintf("  Size: %.1f MB", file.size(EXPR_LOCAL) / 1e6),
  "",
  "MATRIX STRUCTURE",
  sprintf("  Rows: %d (gene/probe level)", nrow(expr_dt)),
  sprintf("  Columns: %d (1 ID + %d samples)", ncol(expr_dt), n_samples),
  sprintf("  ID column: '%s'", id_col),
  sprintf("  ID type (inferred): %s", id_type),
  sprintf("  Data range: [%.2f, %.2f], median=%.2f",
          min(sample_data, na.rm = TRUE),
          max(sample_data, na.rm = TRUE),
          median(sample_data, na.rm = TRUE)),
  sprintf("  Log-scale: %s", ifelse(is_log_scale, "YES", "NO")),
  sprintf("  NA fraction (sample): %.2f%%",
          100 * sum(is.na(sample_data)) / length(sample_data)),
  "",
  "SIGNATURE MATCHING",
  sprintf("  sig_94 direct match: %d/%d (%.1f%%)",
          length(direct_match_94), length(sig_94),
          100 * length(direct_match_94) / length(sig_94)),
  sprintf("  sig_core direct match: %d/%d",
          length(direct_match_core), length(sig_core)),
  if (ensembl_tried) sprintf("  sig_94 via ENSG mapping: %d/%d",
                              ensembl_match_n, length(sig_94)) else "",
  "",
  "SAMPLE ID ALIGNMENT",
  sprintf("  Matrix columns: %d", n_samples),
  sprintf("  xlsx rows: %d", nrow(samples_dt)),
  sprintf("  Overlap via 'Sample name': %d", length(overlap_samples)),
  "",
  "VERDICT"
)

# 自动判断
if (length(direct_match_94) > 70 || ensembl_match_n > 70) {
  report_lines <- c(report_lines,
    "  [PASS] sig_94 > 70/94 genes found on matrix",
    "  >>> R16b can proceed")
} else if (length(direct_match_94) > 40 || ensembl_match_n > 40) {
  report_lines <- c(report_lines,
    sprintf("  [MARGINAL] only ~%d/94 signature genes found",
            max(length(direct_match_94), ensembl_match_n)),
    "  >>> R16b possible but with reduced signature")
} else {
  report_lines <- c(report_lines,
    sprintf("  [FAIL] only %d/94 genes found",
            max(length(direct_match_94), ensembl_match_n)),
    "  >>> R16b cannot proceed — consider alternative data")
}

if (length(overlap_samples) < 100) {
  report_lines <- c(report_lines,
    "  [!] Sample ID alignment NOT via 'Sample name' — check 'title' or other cols")
} else {
  report_lines <- c(report_lines,
    sprintf("  [OK] Sample ID alignment via 'Sample name' (%d overlap)",
            length(overlap_samples)))
}

report_lines <- c(report_lines, "",
                  "================================================================")

writeLines(report_lines, file.path(RES_DIR, "R16a_INSPECTION.txt"))
cat(sprintf("  Saved: %s\n", file.path(RES_DIR, "R16a_INSPECTION.txt")))

cat("\n================================================================\n")
cat("R16a DONE — review INSPECTION.txt before running R16b\n")
cat("================================================================\n")
