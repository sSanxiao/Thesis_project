# ============================================================
# R15: Signature Cleanup — remove direction-conflicting genes
# ------------------------------------------------------------
# 目标:
# 1) 从 R12 的 99 基因 signature 中剔除 5 个方向冲突基因
#    (TENM1, SLC17A7, NRXN3, DCN, SV2B) → 94 基因 cleaned
# 2) 额外构建 8 基因 "multi-sample core" (跨 ≥2 样本且方向完全一致)
# 3) 在 Cavalli 763 例 bulk 上重跑:
#    - Stage A: Subgroup-stratified Cox (ALL/SHH/WNT/Group3/Group4)
#    - Stage D: 12 subtype boxplot
#    - 对 94-gene 和 8-gene core 分别输出
# 4) 核心对比: 94 vs 99 vs 8 — 验证清理不伤害 signal, 且 core 单独也可用
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(GEOquery)
  library(Biobase)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(42)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(RESULTS_DIR, "R15_Cleanup")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)

PROV_FILE <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")
CAVALLI_DIR <- file.path(EXTDATA_DIR, "Cavalli_GSE85217")
CLIN_CSV <- file.path(CAVALLI_DIR, "cavalli2017_mmc2_TableS1_clinical.csv")

# 5 个方向冲突基因 (来自审计)
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")

# ------------------------------------------------------------
# Stage 0: Load signature provenance & construct 3 signatures
# ------------------------------------------------------------
cat("================================================================\n")
cat("R15: Signature Cleanup & Validation\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

cat("[Stage 0] Loading R12 signature provenance...\n")
prov <- fread(PROV_FILE)
cat(sprintf("  Total genes in provenance: %d\n", nrow(prov)))
cat(sprintf("  Columns: %s\n", paste(names(prov), collapse=", ")))

# 找到 gene symbol 列和 direction 列
gene_col <- grep("^gene$|symbol", names(prov), ignore.case=TRUE, value=TRUE)[1]
dir_col <- grep("direction|final_dir", names(prov), ignore.case=TRUE, value=TRUE)[1]
nsamp_col <- grep("n_samples|n_sample", names(prov), ignore.case=TRUE, value=TRUE)[1]

cat(sprintf("  Gene col: %s | Direction col: %s | n_samples col: %s\n",
            gene_col, dir_col, nsamp_col))

# 构建 3 个 signature
sig_99 <- prov[[gene_col]]
sig_94 <- setdiff(sig_99, CONFLICT_GENES)

# core = n_samples >= 2 AND direction 完全一致 (根据审计发现是 8 个基因)
# 从 provenance 里我们只能看 n_samples 和 final direction, 不能直接看 per-sample direction
# 所以先挑 n_samples >= 2 的基因, 再从 CONFLICT_GENES 里排除
core_candidate <- prov[get(nsamp_col) >= 2, ][[gene_col]]
sig_core <- setdiff(core_candidate, CONFLICT_GENES)

cat(sprintf("\n  sig_99 (original):          %d genes\n", length(sig_99)))
cat(sprintf("  sig_94 (cleaned):           %d genes\n", length(sig_94)))
cat(sprintf("  sig_core (multi-sample):    %d genes\n", length(sig_core)))
cat(sprintf("  Conflict genes removed:     %s\n", paste(CONFLICT_GENES, collapse=", ")))
cat(sprintf("  Core genes:                 %s\n", paste(sig_core, collapse=", ")))

# 保存基因列表 + 每个 signature 的方向
for (sig_name in c("sig_99", "sig_94", "sig_core")) {
  genes <- get(sig_name)
  sub <- prov[get(gene_col) %in% genes, ]
  fwrite(sub, file.path(OUT_DIR, sprintf("%s_genes.csv", sig_name)))
}

# ------------------------------------------------------------
# Stage 0.5: Load Cavalli data
# ------------------------------------------------------------
cat("\n[Stage 0.5] Loading Cavalli GSE85217 (cached)...\n")

# series matrix (本地缓存)
gse_file <- list.files(CAVALLI_DIR, pattern = "series_matrix", full.names = TRUE)[1]
cat(sprintf("  Using: %s\n", gse_file))

gse <- getGEO(filename = gse_file, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]

exprs_mat <- exprs(gse)
pdata <- pData(gse)
cat(sprintf("  Expression matrix: %d x %d\n", nrow(exprs_mat), ncol(exprs_mat)))

# ENSG -> SYMBOL mapping
rn <- rownames(exprs_mat)
ensg_ids <- gsub("_at$", "", rn)

cat("  Mapping ENSG -> SYMBOL...\n")
sym_map <- mapIds(org.Hs.eg.db,
                  keys = ensg_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# 聚合到 symbol 级别 (取最大方差 probe)
df_map <- data.table(ensg = ensg_ids, symbol = as.character(sym_map), probe = rn)
df_map <- df_map[!is.na(symbol) & symbol != "", ]
cat(sprintf("  After symbol mapping: %d probes with symbols\n", nrow(df_map)))

# 对于一个 symbol 多个 probe, 保留方差最大的那个
probe_var <- apply(exprs_mat, 1, var)
df_map[, var := probe_var[probe]]
df_map_best <- df_map[order(-var), .SD[1], by = symbol]

expr_sym <- exprs_mat[df_map_best$probe, ]
rownames(expr_sym) <- df_map_best$symbol
cat(sprintf("  After aggregation: %d unique symbols\n", nrow(expr_sym)))

# ------------------------------------------------------------
# Load clinical
# ------------------------------------------------------------
clin <- fread(CLIN_CSV)
cat(sprintf("\n  Clinical rows: %d\n", nrow(clin)))

# pData 的 title -> Study_ID 匹配
pdata_dt <- as.data.table(pdata, keep.rownames = "GSM")
pdata_dt[, Study_ID := title]

bridge <- merge(pdata_dt[, .(GSM, Study_ID)],
                clin, by = "Study_ID", all.x = TRUE)
cat(sprintf("  Bridge matched: %d / %d samples\n",
            sum(!is.na(bridge$Subgroup)), nrow(bridge)))

# 清理生存字段
bridge[, OS_years := as.numeric(`OS (years)`)]
bridge[, event := as.numeric(Dead)]
bridge[, age := as.numeric(Age)]
bridge[, subgroup := as.character(Subgroup)]
bridge[, subtype := as.character(Subtype)]

# 过滤有生存数据的样本
bridge_surv <- bridge[!is.na(OS_years) & !is.na(event) & OS_years > 0, ]
cat(sprintf("  Samples with survival: %d (events: %d)\n",
            nrow(bridge_surv), sum(bridge_surv$event, na.rm = TRUE)))

# ------------------------------------------------------------
# Helper: compute signature score
# ------------------------------------------------------------
compute_score <- function(expr, genes, directions) {
  # genes: vector of gene symbols
  # directions: vector of "positive" / "negative", same length as genes
  
  genes_avail <- intersect(genes, rownames(expr))
  if (length(genes_avail) == 0) return(NULL)
  
  # 对齐方向
  dir_map <- setNames(directions, genes)[genes_avail]
  
  expr_sub <- expr[genes_avail, , drop = FALSE]
  
  # z-score each gene across samples
  z <- t(scale(t(expr_sub)))
  
  # apply direction sign: positive -> +z, negative -> -z
  sign_vec <- ifelse(dir_map == "positive", 1, -1)
  z_signed <- z * sign_vec
  
  # mean across genes
  score <- colMeans(z_signed, na.rm = TRUE)
  
  list(score = score, n_genes_used = length(genes_avail),
       genes_used = genes_avail)
}

# ------------------------------------------------------------
# Helper: run Cox + log-rank
# ------------------------------------------------------------
run_cox <- function(df, score_col) {
  # df: data.table with OS_years, event, age
  # score_col: column name in df
  
  formula_cont <- as.formula(sprintf("Surv(OS_years, event) ~ %s", score_col))
  formula_age <- as.formula(sprintf("Surv(OS_years, event) ~ %s + age", score_col))
  
  cox_cont <- tryCatch(coxph(formula_cont, data = df), error = function(e) NULL)
  cox_age <- tryCatch(coxph(formula_age, data = df), error = function(e) NULL)
  
  # median split log-rank
  med <- median(df[[score_col]], na.rm = TRUE)
  df[, score_bin := ifelse(get(score_col) >= med, "High", "Low")]
  
  sdf <- tryCatch(survdiff(Surv(OS_years, event) ~ score_bin, data = df),
                  error = function(e) NULL)
  
  lr_p <- if (!is.null(sdf)) 1 - pchisq(sdf$chisq, length(sdf$n) - 1) else NA
  
  list(
    n = nrow(df),
    n_events = sum(df$event, na.rm = TRUE),
    cox_HR = if (!is.null(cox_cont)) exp(coef(cox_cont))[1] else NA,
    cox_p = if (!is.null(cox_cont)) summary(cox_cont)$coefficients[1, "Pr(>|z|)"] else NA,
    cox_age_HR = if (!is.null(cox_age)) exp(coef(cox_age))[1] else NA,
    cox_age_p = if (!is.null(cox_age)) summary(cox_age)$coefficients[1, "Pr(>|z|)"] else NA,
    logrank_p = lr_p
  )
}

# ------------------------------------------------------------
# Stage A: 3 signatures × 5 strata Cox
# ------------------------------------------------------------
cat("\n================================================================\n")
cat("STAGE A: Subgroup-stratified Cox for 3 signatures\n")
cat("================================================================\n")

results_all <- list()

for (sig_name in c("sig_99", "sig_94", "sig_core")) {
  cat(sprintf("\n--- %s ---\n", sig_name))
  genes <- get(sig_name)
  sub <- prov[get(gene_col) %in% genes, ]
  
  score_info <- compute_score(expr_sym, sub[[gene_col]], sub[[dir_col]])
  if (is.null(score_info)) {
    cat(sprintf("  SKIP: no genes matched on bulk\n"))
    next
  }
  
  cat(sprintf("  Genes matched on bulk: %d/%d\n",
              score_info$n_genes_used, length(genes)))
  
  # 加 score 到 bridge_surv
  score_vec <- score_info$score
  bridge_scored <- copy(bridge_surv)
  bridge_scored[, score := score_vec[GSM]]
  bridge_scored <- bridge_scored[!is.na(score), ]
  
  for (sg in c("ALL", "SHH", "WNT", "Group3", "Group4")) {
    df <- if (sg == "ALL") bridge_scored else bridge_scored[subgroup == sg, ]
    if (nrow(df) < 20 || sum(df$event) < 5) {
      cat(sprintf("  %s: n=%d events=%d -- SKIP (too few)\n",
                  sg, nrow(df), sum(df$event)))
      next
    }
    
    res <- run_cox(df, "score")
    cat(sprintf("  %s: n=%d ev=%d  HR=%.3f p=%.4f  age-adj HR=%.3f p=%.4f  LR p=%.4f\n",
                sg, res$n, res$n_events, res$cox_HR, res$cox_p,
                res$cox_age_HR, res$cox_age_p, res$logrank_p))
    
    results_all[[length(results_all) + 1]] <- data.table(
      signature = sig_name,
      n_genes = score_info$n_genes_used,
      stratum = sg,
      n = res$n, n_events = res$n_events,
      cox_HR = res$cox_HR, cox_p = res$cox_p,
      cox_age_HR = res$cox_age_HR, cox_age_p = res$cox_age_p,
      logrank_p = res$logrank_p
    )
  }
}

results_dt <- rbindlist(results_all, fill = TRUE)
fwrite(results_dt, file.path(OUT_DIR, "Stage_A_cox_results.csv"))
cat(sprintf("\n  Saved: %s\n", file.path(OUT_DIR, "Stage_A_cox_results.csv")))

# ------------------------------------------------------------
# Stage A Forest plot (subgroup HR for 94-gene)
# ------------------------------------------------------------
cat("\n[Stage A Forest plot for sig_94]\n")
forest_dt <- results_dt[signature == "sig_94" & stratum != "ALL", ]
if (nrow(forest_dt) > 0) {
  forest_dt[, stratum := factor(stratum, levels = c("SHH", "WNT", "Group3", "Group4"))]
  # 近似 95% CI
  forest_dt[, cox_logHR := log(cox_HR)]
  forest_dt[, se := abs(cox_logHR / qnorm(1 - cox_p / 2))]
  forest_dt[, lower := exp(cox_logHR - 1.96 * se)]
  forest_dt[, upper := exp(cox_logHR + 1.96 * se)]
  
  p1 <- ggplot(forest_dt, aes(x = cox_HR, y = stratum)) +
    geom_point(size = 4, color = "#2c7fb8") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "#2c7fb8") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_x_log10(limits = c(0.05, 5), breaks = c(0.1, 0.3, 0.5, 1, 2, 5)) +
    labs(x = "Cox HR (log scale)", y = "Subgroup",
         title = "sig_94 (cleaned): subgroup-stratified HR",
         subtitle = "Point = HR, whiskers = ~95% CI (from Wald p)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUT_DIR, "figures", "stage_A_forest_sig94.png"),
         p1, width = 7, height = 4, dpi = 150)
  cat(sprintf("  Saved: stage_A_forest_sig94.png\n"))
}

# ------------------------------------------------------------
# Stage D: 12 subtype boxplot for 94-gene signature
# ------------------------------------------------------------
cat("\n================================================================\n")
cat("STAGE D: 12-subtype signature score distribution (sig_94)\n")
cat("================================================================\n")

# 重新算 sig_94 score 给全部 bridge 样本 (不只是有生存的)
sub94 <- prov[get(gene_col) %in% sig_94, ]
score_info_94 <- compute_score(expr_sym, sub94[[gene_col]], sub94[[dir_col]])
bridge_all <- copy(bridge)
bridge_all[, score_94 := score_info_94$score[GSM]]
bridge_all <- bridge_all[!is.na(score_94) & !is.na(subtype), ]

# 12 subtype 统计
subt_stats <- bridge_all[, .(
  n = .N,
  mean_score = mean(score_94),
  sd_score = sd(score_94)
), by = subtype][order(-mean_score)]

fwrite(subt_stats, file.path(OUT_DIR, "Stage_D_subtype_stats.csv"))
print(subt_stats)

anova_res <- aov(score_94 ~ subtype, data = bridge_all)
anova_p <- summary(anova_res)[[1]][["Pr(>F)"]][1]
cat(sprintf("\n  12-subtype ANOVA p = %.3e\n", anova_p))

# Boxplot
bridge_all[, subtype := factor(subtype, levels = subt_stats$subtype)]
bridge_all[, subgroup_color := substring(subtype, 1, regexpr("_", subtype) - 1)]
bridge_all[subgroup_color == "", subgroup_color := as.character(subtype)]

p2 <- ggplot(bridge_all, aes(x = subtype, y = score_94, fill = subgroup_color)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  scale_fill_manual(values = c(SHH = "#e41a1c", WNT = "#377eb8",
                               Group3 = "#4daf4a", Group4 = "#984ea3")) +
  labs(x = NULL, y = "sig_94 score",
       title = "sig_94 (94 genes) score across 12 MB subtypes",
       subtitle = sprintf("ANOVA p = %.2e", anova_p),
       fill = "Subgroup") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "figures", "stage_D_12subtype_sig94.png"),
       p2, width = 10, height = 5, dpi = 150)
cat(sprintf("  Saved: stage_D_12subtype_sig94.png\n"))

# ------------------------------------------------------------
# Comparison Table: 99 vs 94 vs core, ALL & SHH
# ------------------------------------------------------------
cat("\n================================================================\n")
cat("KEY COMPARISON TABLE: sig_99 vs sig_94 vs sig_core\n")
cat("================================================================\n")

key_dt <- results_dt[stratum %in% c("ALL", "SHH"), .(
  signature, stratum, n_genes, n, n_events,
  cox_HR = round(cox_HR, 3),
  cox_p = sprintf("%.4f", cox_p),
  age_adj_p = sprintf("%.4f", cox_age_p),
  logrank_p = sprintf("%.4f", logrank_p)
)]
print(key_dt)
fwrite(key_dt, file.path(OUT_DIR, "key_comparison.csv"))

# ------------------------------------------------------------
# Write SUMMARY
# ------------------------------------------------------------
cat("\n[Writing SUMMARY...]\n")

summary_lines <- c(
  "================================================================",
  "R15 SIGNATURE CLEANUP — SUMMARY",
  sprintf("Timestamp: %s", Sys.time()),
  "================================================================",
  "",
  "BACKGROUND",
  "----------",
  "Audit identified 5 direction-conflicting genes in the 99-gene sig:",
  sprintf("  %s", paste(CONFLICT_GENES, collapse = ", ")),
  "These genes' mean_rho < 0.06 due to cancellation between samples.",
  "R15 removes them, producing sig_94 (cleaned).",
  "Additionally constructs sig_core (multi-sample consensus).",
  "",
  "SIGNATURE COMPOSITION",
  "---------------------",
  sprintf("  sig_99 (original):  %d genes", length(sig_99)),
  sprintf("  sig_94 (cleaned):   %d genes", length(sig_94)),
  sprintf("  sig_core:           %d genes (%s)",
          length(sig_core), paste(sig_core, collapse = ", ")),
  "",
  "KEY RESULTS (see key_comparison.csv for full table)",
  "---------------------------------------------------"
)

# append table
summary_lines <- c(summary_lines,
                   capture.output(print(key_dt)))

summary_lines <- c(summary_lines, "",
                   "STAGE D — 12-subtype ANOVA",
                   "--------------------------",
                   sprintf("  sig_94 ANOVA p = %.3e", anova_p),
                   "",
                   "FILES",
                   "-----",
                   sprintf("  Gene lists:  sig_99_genes.csv, sig_94_genes.csv, sig_core_genes.csv"),
                   sprintf("  Cox table:   Stage_A_cox_results.csv"),
                   sprintf("  Subtype:     Stage_D_subtype_stats.csv"),
                   sprintf("  Figures:     figures/stage_A_forest_sig94.png"),
                   sprintf("               figures/stage_D_12subtype_sig94.png"),
                   "",
                   "INTERPRETATION",
                   "--------------")

# 自动判断
sig94_shh <- results_dt[signature == "sig_94" & stratum == "SHH", ]
sig94_all <- results_dt[signature == "sig_94" & stratum == "ALL", ]
sigc_shh <- results_dt[signature == "sig_core" & stratum == "SHH", ]
sigc_all <- results_dt[signature == "sig_core" & stratum == "ALL", ]

if (nrow(sig94_shh) > 0 && !is.na(sig94_shh$cox_age_p) && sig94_shh$cox_age_p < 0.05) {
  summary_lines <- c(summary_lines,
    sprintf("  [OK] sig_94 retains SHH signal: HR=%.3f, age-adj p=%.4f",
            sig94_shh$cox_HR, sig94_shh$cox_age_p))
} else if (nrow(sig94_shh) > 0) {
  summary_lines <- c(summary_lines,
    sprintf("  [!] sig_94 SHH age-adj p=%.4f (not <0.05, review)",
            sig94_shh$cox_age_p))
}

if (nrow(sigc_all) > 0 && !is.na(sigc_all$cox_p) && sigc_all$cox_p < 0.05) {
  summary_lines <- c(summary_lines,
    sprintf("  [OK] sig_core (%d genes) ALL Cox p=%.4f HR=%.3f — robustness anchor confirmed",
            length(sig_core), sigc_all$cox_p, sigc_all$cox_HR))
} else if (nrow(sigc_all) > 0) {
  summary_lines <- c(summary_lines,
    sprintf("  [?] sig_core (%d genes) ALL Cox p=%.4f — see notes",
            length(sig_core), sigc_all$cox_p))
}

summary_lines <- c(summary_lines, "",
                   "================================================================")

writeLines(summary_lines, file.path(OUT_DIR, "R15_SUMMARY.txt"))
cat(sprintf("  SUMMARY: %s\n", file.path(OUT_DIR, "R15_SUMMARY.txt")))

cat("\n================================================================\n")
cat("R15 DONE\n")
cat("================================================================\n")
