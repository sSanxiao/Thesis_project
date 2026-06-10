#!/usr/bin/env Rscript
# ============================================================================
# R14 PATCH: re-lock signature to R12's 99 genes and re-run R13 stages A/C/D
#
# Why: R13 reconstructed the strict signature from scratch and produced 131
#      genes, which diluted the SHH-specific survival signal that R11/R12
#      reported (SHH HR=0.22, p=0.007 -> R13 got HR=0.51, p=0.16).
#      Paper narrative needs one signature used consistently from R10 onward.
#
# Plan:
#   - Load the canonical 99-gene signature from R12 (with per-gene direction)
#   - Re-score Cavalli bulk
#   - Re-run Stage A (subgroup-stratified Cox) and Stage C (leave-one-out)
#     and Stage D (12-subtype boxplot) under the 99-gene definition
#   - Stage B (Xenium fingerprint match) is signature-independent -> skip
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(data.table)
  library(survival); library(ggplot2)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(matrixStats)
})

options(stringsAsFactors = FALSE)
set.seed(42)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ---- paths ----
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
R12_DIR  <- file.path(RESULTS_DIR, "R12_Gaps")
R4_DIR   <- file.path(RESULTS_DIR, "R4_Results", "Medulloblastoma_Human")
OUT_DIR  <- file.path(RESULTS_DIR, "R14_Patch")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

GEO_CACHE    <- file.path(EXTDATA_DIR, "Cavalli_GSE85217")
CLINICAL_CSV <- file.path(GEO_CACHE, "cavalli2017_mmc2_TableS1_clinical.csv")

MB_SAMPLES <- c("GSM8840046", "GSM8840047", "GSM8840048", "GSM8840049")
MB_NAMES   <- c("MB263",      "MB266",      "MB295",      "MB299")
names(MB_NAMES) <- MB_SAMPLES

cat("================================================================\n")
cat("R14 PATCH: re-lock to R12 99-gene signature\n")
cat(format(Sys.time()), "\n")
cat("================================================================\n\n")

# ============================================================================
# STAGE 0: Load the canonical 99-gene signature from R12
# ============================================================================
cat("[Stage 0] Loading canonical 99-gene signature from R12...\n")

sig99_file <- file.path(R12_DIR, "sig_strict_99_provenance.csv")
if (!file.exists(sig99_file)) {
  sig99_file_slim <- file.path(R12_DIR, "sig_strict_99_provenance_slim.csv")
  if (file.exists(sig99_file_slim)) {
    sig99_file <- sig99_file_slim
    cat("  Using slim version: ", sig99_file_slim, "\n", sep = "")
  } else {
    stop("Cannot find 99-gene signature from R12. Expected one of:\n",
         "  ", file.path(R12_DIR, "sig_strict_99_provenance.csv"), "\n",
         "  ", file.path(R12_DIR, "sig_strict_99_provenance_slim.csv"))
  }
}

sig99 <- fread(sig99_file)
cat("  Read: ", nrow(sig99), " rows from ", basename(sig99_file), "\n", sep = "")
cat("  Columns: ", paste(names(sig99), collapse = ", "), "\n", sep = "")
cat("  First few rows:\n")
print(head(sig99, 5))

# Detect gene column and direction column
gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol"), names(sig99))[1]
dir_col  <- intersect(c("direction", "Direction"), names(sig99))[1]

if (is.na(gene_col)) stop("Cannot find gene column in ", sig99_file)

# If direction is not there, we need to reconstruct from R4 (strict = tier1_strict, any sample)
if (is.na(dir_col)) {
  cat("  No direction column in provenance file; reconstructing from R4...\n")
  # Read R4 tier1_strict genes
  all_mb <- rbindlist(lapply(seq_along(MB_SAMPLES), function(i) {
    f <- file.path(R4_DIR, MB_SAMPLES[i], "filtered_density_genes.csv")
    dt <- fread(f)
    dt[, sample_id := MB_SAMPLES[i]]
    dt
  }), fill = TRUE)
  t1 <- all_mb[tier == "tier1_strict"]
  # For each gene, pick direction from sample with max |rho|
  dir_df <- t1[, .(direction = direction[which.max(abs(rho_knn_main))],
                   mean_rho  = mean(rho_knn_main, na.rm = TRUE)),
               by = gene]
  sig99_final <- merge(sig99[, .SD, .SDcols = gene_col], dir_df,
                       by.x = gene_col, by.y = "gene", all.x = TRUE)
  setnames(sig99_final, gene_col, "gene")
} else {
  sig99_final <- data.table(gene = sig99[[gene_col]],
                            direction = sig99[[dir_col]])
  if ("mean_rho" %in% names(sig99)) sig99_final$mean_rho <- sig99$mean_rho
  setnames(sig99_final, gene_col, "gene", skip_absent = TRUE)
  # ensure gene column exists
  if (!"gene" %in% names(sig99_final)) sig99_final$gene <- sig99[[gene_col]]
}

# If the file had 99 rows but some are duplicates, deduplicate
sig99_final <- unique(sig99_final, by = "gene")
sig99_final <- sig99_final[!is.na(direction)]
cat("  Final locked signature: ", nrow(sig99_final), " unique genes\n", sep = "")
cat("  Direction breakdown:\n")
print(table(sig99_final$direction))

# Write out the locked signature for reference
fwrite(sig99_final, file.path(OUT_DIR, "LOCKED_signature_genes.csv"))
cat("  Saved: LOCKED_signature_genes.csv\n\n")

# ============================================================================
# STAGE 0.1: Load Cavalli data
# ============================================================================
cat("[Stage 0.1] Loading Cavalli GSE85217 from local cache...\n")
gse <- getGEO("GSE85217", GSEMatrix = TRUE, destdir = GEO_CACHE, AnnotGPL = FALSE)
eset <- gse[[1]]
exprs_mat <- exprs(eset)
pdata <- pData(eset)
cat("  Expression matrix: ", nrow(exprs_mat), " x ", ncol(exprs_mat), "\n", sep = "")

# ENSG -> SYMBOL
rn <- rownames(exprs_mat)
ensg <- sub("_at$", "", rn)
sym_map <- suppressMessages(mapIds(
  org.Hs.eg.db, keys = ensg, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"
))
valid <- !is.na(sym_map) & sym_map != ""
exprs_mat <- exprs_mat[valid, , drop = FALSE]
sym_map <- sym_map[valid]
rv <- matrixStats::rowVars(exprs_mat)
o <- order(sym_map, -rv)
exprs_mat <- exprs_mat[o, ]
sym_map <- sym_map[o]
keep <- !duplicated(sym_map)
exprs_mat <- exprs_mat[keep, ]
rownames(exprs_mat) <- sym_map[keep]
cat("  After symbol mapping: ", nrow(exprs_mat), " unique gene symbols\n", sep = "")

# Clinical bridge
clin <- fread(CLINICAL_CSV)
title_map <- setNames(pdata$title, rownames(pdata))
gsm_to_title <- data.table(GSM = names(title_map), Study_ID = unname(title_map))
merged_clin <- merge(gsm_to_title, clin, by = "Study_ID", all.x = TRUE)
merged_clin <- merged_clin[!is.na(Subgroup)]
setnames(merged_clin, "OS (years)", "OS_years", skip_absent = TRUE)
setnames(merged_clin, "Dead", "Dead", skip_absent = TRUE)
merged_clin[, OS_years := as.numeric(OS_years)]
merged_clin[, Dead := as.integer(Dead)]
merged_clin[, Age := as.numeric(Age)]

cat("  Clinical bridge: ", nrow(merged_clin), " samples with subgroup\n", sep = "")

# ============================================================================
# HELPERS
# ============================================================================
score_signature <- function(sig_df, exprs_mat) {
  common <- intersect(sig_df$gene, rownames(exprs_mat))
  if (length(common) < 3) return(NULL)
  sub_sig <- sig_df[match(common, sig_df$gene), ]
  sub_expr <- exprs_mat[common, , drop = FALSE]
  z <- t(scale(t(sub_expr)))
  sign_vec <- ifelse(sub_sig$direction == "negative", -1, 1)
  z_signed <- z * sign_vec
  score <- colMeans(z_signed, na.rm = TRUE)
  list(score = score, n_matched = length(common), genes_used = common)
}

run_5_tests <- function(score, clin_sub, label = "") {
  common_s <- intersect(names(score), clin_sub$GSM)
  if (length(common_s) < 20) return(NULL)
  s_aligned <- score[common_s]
  c_aligned <- clin_sub[match(common_s, GSM)]
  c_aligned[, score := s_aligned]

  out <- list(label = label, n = nrow(c_aligned))

  if (length(unique(c_aligned$Subgroup)) > 1) {
    out$A1_p <- tryCatch(summary(aov(score ~ Subgroup, data = c_aligned))[[1]][1, "Pr(>F)"],
                         error = function(e) NA)
  } else { out$A1_p <- NA }

  if (length(unique(c_aligned$Subtype)) > 1) {
    out$A2_p <- tryCatch(summary(aov(score ~ Subtype, data = c_aligned))[[1]][1, "Pr(>F)"],
                         error = function(e) NA)
  } else { out$A2_p <- NA }

  surv_ok <- !is.na(c_aligned$OS_years) & !is.na(c_aligned$Dead)
  c_surv <- c_aligned[surv_ok]
  if (nrow(c_surv) >= 20 && sum(c_surv$Dead) >= 5) {
    c_surv[, score_grp := ifelse(score >= median(score), "High", "Low")]
    out$B_logrank_p <- tryCatch({
      sdiff <- survdiff(Surv(OS_years, Dead) ~ score_grp, data = c_surv)
      1 - pchisq(sdiff$chisq, df = 1)
    }, error = function(e) NA)

    cox_fit <- tryCatch(coxph(Surv(OS_years, Dead) ~ score, data = c_surv),
                        error = function(e) NULL)
    if (!is.null(cox_fit)) {
      s <- summary(cox_fit)
      out$Cox_HR <- s$coefficients[1, "exp(coef)"]
      out$Cox_p  <- s$coefficients[1, "Pr(>|z|)"]
    }

    if ("Age" %in% names(c_surv) && sum(!is.na(c_surv$Age)) > 20) {
      cox2 <- tryCatch(coxph(Surv(OS_years, Dead) ~ score + Age, data = c_surv),
                       error = function(e) NULL)
      if (!is.null(cox2)) {
        s2 <- summary(cox2)
        out$CoxAge_HR <- s2$coefficients["score", "exp(coef)"]
        out$CoxAge_p  <- s2$coefficients["score", "Pr(>|z|)"]
      }
    }
    out$n_surv <- nrow(c_surv)
    out$n_events <- sum(c_surv$Dead)
  }
  out
}

# ============================================================================
# STAGE A (locked 99-gene): subgroup-stratified validation
# ============================================================================
cat("\n================================================================\n")
cat("STAGE A (99-gene): subgroup-stratified validation\n")
cat("================================================================\n")

score_99 <- score_signature(sig99_final, exprs_mat)
cat("  99-gene signature matched ", score_99$n_matched, "/",
    nrow(sig99_final), " genes on bulk\n", sep = "")

stage_a_results <- list()
for (sg in c("SHH", "WNT", "Group3", "Group4", "ALL")) {
  c_sub <- if (sg == "ALL") merged_clin else merged_clin[Subgroup == sg]
  if (nrow(c_sub) < 20) next
  res <- run_5_tests(score_99$score, c_sub, label = sg)
  if (!is.null(res)) stage_a_results[[sg]] <- res
  cat(sprintf("  %-8s n=%d  A1=%.3g  B_logrank=%.3g  Cox_HR=%.3g (p=%.3g)  CoxAge_p=%.3g\n",
              sg, res$n %||% 0,
              res$A1_p %||% NA, res$B_logrank_p %||% NA,
              res$Cox_HR %||% NA, res$Cox_p %||% NA, res$CoxAge_p %||% NA))
}

stage_a_tbl <- rbindlist(lapply(stage_a_results, function(x) {
  data.table(
    Subgroup = x$label,
    n = x$n,
    n_surv = x$n_surv %||% NA,
    n_events = x$n_events %||% NA,
    A1_Subgroup_p = x$A1_p %||% NA,
    A2_Subtype_p = x$A2_p %||% NA,
    B_logrank_p = x$B_logrank_p %||% NA,
    Cox_HR = x$Cox_HR %||% NA,
    Cox_p = x$Cox_p %||% NA,
    CoxAge_HR = x$CoxAge_HR %||% NA,
    CoxAge_p = x$CoxAge_p %||% NA
  )
}), fill = TRUE)
fwrite(stage_a_tbl, file.path(OUT_DIR, "STAGE_A_subgroup_stratified_99gene.csv"))

# Forest plot with CIs
fa <- stage_a_tbl[Subgroup != "ALL" & !is.na(Cox_HR)]
if (nrow(fa) >= 2) {
  fa[, CI_low := NA_real_]; fa[, CI_high := NA_real_]
  for (i in seq_len(nrow(fa))) {
    sg <- fa$Subgroup[i]
    c_sub <- merged_clin[Subgroup == sg & !is.na(OS_years) & !is.na(Dead)]
    s <- score_99$score[c_sub$GSM]
    c_sub[, score := s]
    c_sub <- c_sub[!is.na(score)]
    if (nrow(c_sub) < 20 || sum(c_sub$Dead) < 5) next
    fit <- tryCatch(coxph(Surv(OS_years, Dead) ~ score, data = c_sub),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      ci <- confint(fit)
      fa$CI_low[i]  <- exp(ci[1, 1])
      fa$CI_high[i] <- exp(ci[1, 2])
    }
  }

  p_forest <- ggplot(fa, aes(x = Cox_HR, y = Subgroup)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    geom_text(aes(label = sprintf("HR=%.2f\np=%.3f", Cox_HR, Cox_p)),
              hjust = -0.15, size = 3) +
    scale_x_log10(limits = c(0.05, 10)) +
    labs(x = "Cox HR (log scale)", y = NULL,
         title = "Stage A (locked 99-gene): subgroup-stratified Cox HR",
         subtitle = "MBEN-derived signature, median-split") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(PLOT_DIR, "stage_A_forest_subgroup_99gene.png"),
         p_forest, width = 7, height = 4, dpi = 200)
  cat("  Saved plot: stage_A_forest_subgroup_99gene.png\n")
}

# ============================================================================
# STAGE C (locked 99-gene): leave-one-out stability
# ============================================================================
cat("\n================================================================\n")
cat("STAGE C (99-gene): leave-one-out stability\n")
cat("================================================================\n")

# Build LOO signatures: for each excluded sample, use tier1_strict genes
# from remaining samples, but restricted to the 99-gene locked set
# (i.e., a gene is in the LOO signature only if it was in the locked 99
#  AND appears as tier1_strict in at least one non-excluded sample)

# Re-read R4 tier1_strict for provenance
all_mb_t1 <- rbindlist(lapply(seq_along(MB_SAMPLES), function(i) {
  f <- file.path(R4_DIR, MB_SAMPLES[i], "filtered_density_genes.csv")
  dt <- fread(f)
  dt[, sample_id := MB_SAMPLES[i]]
  dt[tier == "tier1_strict"]
}), fill = TRUE)

# Which locked-99 genes appear in which samples as tier1_strict
locked_genes <- sig99_final$gene
provenance_t1 <- all_mb_t1[gene %in% locked_genes,
                           .(samples = list(unique(sample_id))), by = gene]

stage_c_results <- list()
# Full signature baseline
stage_c_results[["full"]] <- list(
  excluded = "none",
  n_genes = nrow(sig99_final),
  n_matched = score_99$n_matched,
  all = run_5_tests(score_99$score, merged_clin, label = "full"),
  shh = run_5_tests(score_99$score, merged_clin[Subgroup == "SHH"], label = "full_SHH")
)

for (excl in MB_SAMPLES) {
  # Keep only genes that are tier1_strict in at least one non-excluded sample
  loo_genes <- provenance_t1[sapply(samples, function(s) any(s != excl)), gene]
  loo_sig <- sig99_final[gene %in% loo_genes]

  if (nrow(loo_sig) < 3) {
    cat(sprintf("  excl %s: only %d genes remain; skipping\n", MB_NAMES[excl], nrow(loo_sig)))
    next
  }

  sc_loo <- score_signature(loo_sig, exprs_mat)
  res_all <- run_5_tests(sc_loo$score, merged_clin, label = paste0("excl_", MB_NAMES[excl]))
  res_shh <- run_5_tests(sc_loo$score, merged_clin[Subgroup == "SHH"],
                         label = paste0("excl_", MB_NAMES[excl], "_SHH"))
  stage_c_results[[excl]] <- list(
    excluded = MB_NAMES[excl],
    n_genes = nrow(loo_sig),
    n_matched = sc_loo$n_matched,
    all = res_all, shh = res_shh
  )
}

stage_c_tbl <- rbindlist(lapply(stage_c_results, function(x) {
  rbind(
    data.table(excluded = x$excluded, scope = "ALL",
               n_genes = x$n_genes, n_matched = x$n_matched,
               n = x$all$n %||% NA,
               A1_p = x$all$A1_p %||% NA,
               B_logrank_p = x$all$B_logrank_p %||% NA,
               Cox_HR = x$all$Cox_HR %||% NA,
               Cox_p  = x$all$Cox_p %||% NA,
               CoxAge_HR = x$all$CoxAge_HR %||% NA,
               CoxAge_p = x$all$CoxAge_p %||% NA),
    data.table(excluded = x$excluded, scope = "SHH",
               n_genes = x$n_genes, n_matched = x$n_matched,
               n = x$shh$n %||% NA,
               A1_p = NA,
               B_logrank_p = x$shh$B_logrank_p %||% NA,
               Cox_HR = x$shh$Cox_HR %||% NA,
               Cox_p  = x$shh$Cox_p %||% NA,
               CoxAge_HR = x$shh$CoxAge_HR %||% NA,
               CoxAge_p = x$shh$CoxAge_p %||% NA)
  )
}))
fwrite(stage_c_tbl, file.path(OUT_DIR, "STAGE_C_leave_one_out_99gene.csv"))

cat("\n  ALL cohort (99-gene LOO):\n")
print(stage_c_tbl[scope == "ALL", .(excluded, n_genes, Cox_HR, Cox_p, CoxAge_p)])
cat("\n  SHH subgroup (99-gene LOO):\n")
print(stage_c_tbl[scope == "SHH", .(excluded, n_genes, Cox_HR, Cox_p, B_logrank_p)])

# ============================================================================
# STAGE D (locked 99-gene): 12-subtype profile
# ============================================================================
cat("\n================================================================\n")
cat("STAGE D (99-gene): 12-subtype profile\n")
cat("================================================================\n")

df_score <- data.table(GSM = names(score_99$score),
                       score = as.numeric(score_99$score))
df_score <- merge(df_score, merged_clin[, .(GSM, Subgroup, Subtype, Age, OS_years, Dead)],
                  by = "GSM")

df_score[, Subgroup := factor(Subgroup, levels = c("WNT", "SHH", "Group3", "Group4"))]
subtype_order <- df_score[, .(mean_score = mean(score, na.rm = TRUE)), by = Subtype][order(-mean_score), Subtype]
df_score[, Subtype := factor(Subtype, levels = subtype_order)]

sub_stats <- df_score[, .(n = .N, mean_score = mean(score, na.rm = TRUE),
                          median_score = median(score, na.rm = TRUE),
                          sd_score = sd(score, na.rm = TRUE)),
                      by = .(Subgroup, Subtype)]
fwrite(sub_stats[order(-mean_score)], file.path(OUT_DIR, "STAGE_D_subtype_stats_99gene.csv"))
cat("  Subtype score rankings (99-gene):\n")
print(sub_stats[order(-mean_score)])

subgroup_colors <- c(WNT = "#66c2a5", SHH = "#fc8d62", Group3 = "#8da0cb", Group4 = "#e78ac3")
p_12 <- ggplot(df_score, aes(x = Subtype, y = score, fill = Subgroup)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  scale_fill_manual(values = subgroup_colors) +
  labs(x = "MB molecular subtype (sorted by mean score)",
       y = "Density signature score (z-score)",
       title = "Stage D (locked 99-gene): signature score across 12 MB subtypes",
       subtitle = sprintf("n=%d; MBEN-derived 99-gene signature", nrow(df_score))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
ggsave(file.path(PLOT_DIR, "stage_D_12subtype_boxplot_99gene.png"),
       p_12, width = 10, height = 5, dpi = 200)
cat("  Saved plot: stage_D_12subtype_boxplot_99gene.png\n")

a_subtype <- summary(aov(score ~ Subtype, data = df_score))[[1]][1, "Pr(>F)"]
cat(sprintf("  12-subtype ANOVA p = %.3g\n", a_subtype))

# Compare SHH_gamma variance between R13 (131) and R14 (99)
shh_gamma_sd_99 <- sub_stats[Subtype == "SHH_gamma", sd_score]
cat(sprintf("  SHH_gamma sd (99-gene): %.3f  (R13 with 131 genes was 0.365)\n",
            shh_gamma_sd_99 %||% NA))

# ============================================================================
# FINAL SUMMARY
# ============================================================================
sink(file.path(OUT_DIR, "R14_SUMMARY.txt"))
cat("================================================================\n")
cat("  R14 PATCH: locked 99-gene signature re-run of Stages A/C/D\n")
cat("  ", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("LOCKED SIGNATURE\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  %d genes (from R12 sig_strict_99_provenance.csv)\n", nrow(sig99_final)))
cat(sprintf("  Direction: %d negative, %d positive\n",
            sum(sig99_final$direction == "negative"),
            sum(sig99_final$direction == "positive")))
cat(sprintf("  Matched in Cavalli bulk: %d / %d\n",
            score_99$n_matched, nrow(sig99_final)))

cat("\nSTAGE A (99-gene): subgroup-stratified validation\n")
cat("----------------------------------------------------------------\n")
print(stage_a_tbl)
sh <- stage_a_tbl[Subgroup == "SHH"]
if (nrow(sh) == 1) {
  cat(sprintf("\n  Key SHH numbers: n=%d, n_surv=%d, events=%d\n",
              sh$n, sh$n_surv, sh$n_events))
  cat(sprintf("                   Cox HR=%.3f p=%.4f, CoxAge p=%.4f\n",
              sh$Cox_HR, sh$Cox_p, sh$CoxAge_p))
  cat(sprintf("                   log-rank p=%.4f\n", sh$B_logrank_p))
}

cat("\n  Comparison vs R13 (131-gene) and R11/R12 expectations:\n")
cat("  Test              R14 (99)          R13 (131)       R11/R12 reported\n")
cat(sprintf("  SHH Cox HR        %.3f             0.508           0.22\n",
            sh$Cox_HR %||% NA))
cat(sprintf("  SHH Cox p         %.4f            0.157           0.007\n",
            sh$Cox_p %||% NA))
cat(sprintf("  SHH log-rank p    %.4f            0.032           0.0073\n",
            sh$B_logrank_p %||% NA))

cat("\nSTAGE C (99-gene): leave-one-out stability\n")
cat("----------------------------------------------------------------\n")
cat("  ALL cohort:\n")
print(stage_c_tbl[scope == "ALL", .(excluded, n_genes, Cox_HR, Cox_p, CoxAge_p)])
cat("\n  Within SHH subgroup:\n")
print(stage_c_tbl[scope == "SHH", .(excluded, n_genes, Cox_HR, Cox_p, B_logrank_p)])

cat("\nSTAGE D (99-gene): 12-subtype profile\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  Overall 12-subtype ANOVA: p = %.3g\n", a_subtype))
cat(sprintf("  SHH_gamma sd: %.3f (R13 131-gene version was 0.365)\n",
            shh_gamma_sd_99 %||% NA))
cat("\n  Full rankings:\n")
print(sub_stats[order(-mean_score)])

cat("\n================================================================\n")
cat("  DECISION LOG for paper\n")
cat("================================================================\n")
cat("  - Main signature for paper:  LOCKED 99-gene (this R14 run)\n")
cat("  - Sensitivity reference:     131-gene R13 (retain as SI)\n")
cat("  - Stage B (fingerprint):     signature-independent, use R13 output\n")
cat("  - Key numbers to report:\n")
cat(sprintf("    * SHH Cox HR = %.3f (p=%.4f) in n=%d\n",
            sh$Cox_HR %||% NA, sh$Cox_p %||% NA, sh$n_surv %||% NA))
cat(sprintf("    * 12-subtype ANOVA p = %.3g\n", a_subtype))
sink()

cat("\n\nR14_SUMMARY.txt saved to: ", file.path(OUT_DIR, "R14_SUMMARY.txt"), "\n", sep = "")
cat("\nAll outputs:\n")
print(list.files(OUT_DIR))
print(list.files(PLOT_DIR))
cat("\nR14 patch complete.\n")
