#!/usr/bin/env Rscript
# ============================================================================
# R13: MBEN-SHH subtype validation — full ABCD package
#
# Background: All 4 Xenium MB samples (MB263/266/295/299) are MBEN histology,
#             which per WHO classification is SHH-activated molecular subgroup.
#
# Plan:
#   A. Subgroup-stratified validation  (SHH/WNT/G3/G4 独立跑 5 tests)
#   B. Xenium -> Cavalli reverse matching  (4 MB samples 最像哪个亚型)
#   C. Leave-one-out signature stability  (排 MB266 后还剩多少信号)
#   D. 12-subtype signature score profile  (细分亚型梯度 + rhombic lip 对接)
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(data.table)
  library(survival); library(ggplot2)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

options(stringsAsFactors = FALSE)
set.seed(42)

# Helper: null-coalesce operator (must be defined before first use)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ---- paths ----
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
OUT_DIR <- file.path(RESULTS_DIR, "R13_MBEN")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

R4_DIR <- file.path(RESULTS_DIR, "R4_Results", "Medulloblastoma_Human")
GEO_CACHE <- file.path(EXTDATA_DIR, "Cavalli_GSE85217")
CLINICAL_CSV <- file.path(GEO_CACHE, "cavalli2017_mmc2_TableS1_clinical.csv")

MB_SAMPLES <- c("GSM8840046", "GSM8840047", "GSM8840048", "GSM8840049")
MB_NAMES   <- c("MB263",      "MB266",      "MB295",      "MB299")
names(MB_NAMES) <- MB_SAMPLES

# ============================================================================
# STAGE 0: Load data (reuses R10/R11/R12 cached outputs)
# ============================================================================
cat("================================================================\n")
cat("R13: MBEN-SHH full validation (A+B+C+D)\n")
cat(format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("[Stage 0] Loading cached data...\n")

# --- 0.1 R4 gene lists per MB sample ---
mb_gene_data <- list()
for (i in seq_along(MB_SAMPLES)) {
  f <- file.path(R4_DIR, MB_SAMPLES[i], "filtered_density_genes.csv")
  if (!file.exists(f)) stop("Missing: ", f)
  dt <- fread(f)
  dt[, sample_id := MB_SAMPLES[i]]
  dt[, sample_name := MB_NAMES[i]]
  mb_gene_data[[MB_SAMPLES[i]]] <- dt
}
all_mb <- rbindlist(mb_gene_data, fill = TRUE)
cat("  R4 gene data: ", nrow(all_mb), " records across ", length(MB_SAMPLES), " samples\n", sep = "")

# --- 0.2 Signature strategies (from R11/R12) ---
build_sig_strict <- function(dt) {
  # tier1_strict, any sample (>=1)
  t1 <- dt[tier == "tier1_strict"]
  # keep unique gene x direction; if conflicting direction, use max |rho| sample's direction
  t1 <- t1[, .(direction = direction[which.max(abs(rho_knn_main))],
               mean_rho  = mean(rho_knn_main, na.rm = TRUE),
               n_samples = .N),
           by = gene]
  t1
}

sig_strict_full <- build_sig_strict(all_mb)
cat("  strict signature (all 4 samples): ", nrow(sig_strict_full), " genes\n", sep = "")

# Leave-one-out: build signature excluding each sample
sig_loo <- list()
for (excl in MB_SAMPLES) {
  dt_sub <- all_mb[sample_id != excl]
  sig_loo[[excl]] <- build_sig_strict(dt_sub)
  cat("    LOO excluding ", MB_NAMES[excl], " (", excl, "): ",
      nrow(sig_loo[[excl]]), " genes\n", sep = "")
}

# --- 0.3 Cavalli data ---
cat("\n[Stage 0.3] Loading Cavalli GSE85217...\n")
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
# collapse duplicated symbols: keep max-variance probe
rv <- matrixStats::rowVars(exprs_mat)
o <- order(sym_map, -rv)
exprs_mat <- exprs_mat[o, ]
sym_map <- sym_map[o]
keep <- !duplicated(sym_map)
exprs_mat <- exprs_mat[keep, ]
rownames(exprs_mat) <- sym_map[keep]
cat("  After symbol mapping: ", nrow(exprs_mat), " unique gene symbols\n", sep = "")

# --- 0.4 Clinical data bridge ---
clin <- fread(CLINICAL_CSV)
title_map <- setNames(pdata$title, rownames(pdata))
gsm_to_title <- data.table(GSM = names(title_map), Study_ID = unname(title_map))
merged_clin <- merge(gsm_to_title, clin, by = "Study_ID", all.x = TRUE)
merged_clin <- merged_clin[!is.na(Subgroup)]
rownames(merged_clin) <- merged_clin$GSM

setnames(merged_clin, "OS (years)", "OS_years", skip_absent = TRUE)
setnames(merged_clin, "Dead", "Dead", skip_absent = TRUE)
merged_clin[, OS_years := as.numeric(OS_years)]
merged_clin[, Dead := as.integer(Dead)]
merged_clin[, Age := as.numeric(Age)]

cat("  Clinical bridge: ", nrow(merged_clin), " samples with subgroup\n", sep = "")
cat("  Subgroup distribution:\n")
print(table(merged_clin$Subgroup, useNA = "ifany"))
cat("  Subtype distribution:\n")
print(table(merged_clin$Subtype, useNA = "ifany"))

# ============================================================================
# Helper: compute signature score on bulk
# ============================================================================
score_signature <- function(sig_df, exprs_mat) {
  # sig_df: gene, direction
  # exprs_mat: rownames = gene symbols, cols = samples
  common <- intersect(sig_df$gene, rownames(exprs_mat))
  if (length(common) < 3) return(NULL)
  sub_sig <- sig_df[match(common, sig_df$gene), ]
  sub_expr <- exprs_mat[common, , drop = FALSE]
  # z-score per gene
  z <- t(scale(t(sub_expr)))
  # flip negative-direction genes
  sign_vec <- ifelse(sub_sig$direction == "negative", -1, 1)
  z_signed <- z * sign_vec
  score <- colMeans(z_signed, na.rm = TRUE)
  list(score = score, n_matched = length(common))
}

# Helper: run 5 tests on a given score + clinical subset
run_5_tests <- function(score, clin_sub, label = "") {
  # align
  common_s <- intersect(names(score), clin_sub$GSM)
  if (length(common_s) < 20) return(NULL)
  s_aligned <- score[common_s]
  c_aligned <- clin_sub[match(common_s, GSM)]
  c_aligned[, score := s_aligned]

  out <- list(label = label, n = nrow(c_aligned))

  # A1: subgroup ANOVA (if multiple subgroups)
  if (length(unique(c_aligned$Subgroup)) > 1) {
    a1 <- tryCatch(summary(aov(score ~ Subgroup, data = c_aligned))[[1]][1, "Pr(>F)"],
                   error = function(e) NA)
    out$A1_p <- a1
  } else {
    out$A1_p <- NA
  }

  # A2: subtype ANOVA
  if (length(unique(c_aligned$Subtype)) > 1) {
    a2 <- tryCatch(summary(aov(score ~ Subtype, data = c_aligned))[[1]][1, "Pr(>F)"],
                   error = function(e) NA)
    out$A2_p <- a2
  } else {
    out$A2_p <- NA
  }

  # Survival subset
  surv_ok <- !is.na(c_aligned$OS_years) & !is.na(c_aligned$Dead)
  c_surv <- c_aligned[surv_ok]
  if (nrow(c_surv) >= 20 && sum(c_surv$Dead) >= 5) {
    # log-rank (median split)
    c_surv[, score_grp := ifelse(score >= median(score), "High", "Low")]
    lr <- tryCatch({
      sdiff <- survdiff(Surv(OS_years, Dead) ~ score_grp, data = c_surv)
      1 - pchisq(sdiff$chisq, df = 1)
    }, error = function(e) NA)
    out$B_logrank_p <- lr

    # Cox
    cox_fit <- tryCatch(coxph(Surv(OS_years, Dead) ~ score, data = c_surv),
                        error = function(e) NULL)
    if (!is.null(cox_fit)) {
      s <- summary(cox_fit)
      out$Cox_HR <- s$coefficients[1, "exp(coef)"]
      out$Cox_p  <- s$coefficients[1, "Pr(>|z|)"]
    }

    # Cox age-adjusted (if age available)
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
# STAGE A: Subgroup-stratified validation
# ============================================================================
cat("\n================================================================\n")
cat("STAGE A: Subgroup-stratified validation (strict 99-gene signature)\n")
cat("================================================================\n")

score_strict <- score_signature(sig_strict_full, exprs_mat)
cat("  strict signature matched ", score_strict$n_matched, "/",
    nrow(sig_strict_full), " genes on bulk\n", sep = "")

stage_a_results <- list()
for (sg in c("SHH", "WNT", "Group3", "Group4", "ALL")) {
  if (sg == "ALL") {
    c_sub <- merged_clin
  } else {
    c_sub <- merged_clin[Subgroup == sg]
  }
  if (nrow(c_sub) < 20) next
  res <- run_5_tests(score_strict$score, c_sub, label = sg)
  if (!is.null(res)) stage_a_results[[sg]] <- res
  cat(sprintf("  %-8s n=%d  A1=%.3g  B_logrank=%.3g  Cox_HR=%.3g (p=%.3g)  CoxAge_p=%.3g\n",
              sg, res$n %||% 0,
              res$A1_p %||% NA, res$B_logrank_p %||% NA,
              res$Cox_HR %||% NA, res$Cox_p %||% NA, res$CoxAge_p %||% NA))
}

# Save Stage A table
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
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
fwrite(stage_a_tbl, file.path(OUT_DIR, "STAGE_A_subgroup_stratified.csv"))
cat("  Saved: STAGE_A_subgroup_stratified.csv\n")

# Stage A Forest plot: Cox HR by subgroup
fa <- stage_a_tbl[Subgroup != "ALL" & !is.na(Cox_HR)]
if (nrow(fa) >= 2) {
  # compute 95% CI from HR and p (approximation via Wald: se = log(HR)/qnorm)
  # use coxph re-fit for proper CI
  fa[, CI_low := NA_real_]
  fa[, CI_high := NA_real_]
  for (i in seq_len(nrow(fa))) {
    sg <- fa$Subgroup[i]
    c_sub <- merged_clin[Subgroup == sg & !is.na(OS_years) & !is.na(Dead)]
    s <- score_strict$score[c_sub$GSM]
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
         title = "Stage A: Subgroup-stratified Cox HR",
         subtitle = "strict signature (99 genes), median-split") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(PLOT_DIR, "stage_A_forest_subgroup.png"),
         p_forest, width = 7, height = 4, dpi = 200)
  cat("  Saved plot: stage_A_forest_subgroup.png\n")
}

# ============================================================================
# STAGE B: Xenium -> Cavalli reverse matching
# ============================================================================
cat("\n================================================================\n")
cat("STAGE B: Xenium sample fingerprints match to Cavalli subgroups\n")
cat("================================================================\n")

# For each of 4 MB samples, build a density-correlation "fingerprint"
# using rho_knn_main across all genes; then compute per-subgroup mean
# expression profile in Cavalli restricted to those genes;
# correlate the two.

panel_genes <- unique(all_mb$gene)
panel_genes_on_bulk <- intersect(panel_genes, rownames(exprs_mat))
cat("  Xenium panel genes found in Cavalli bulk: ",
    length(panel_genes_on_bulk), " / ", length(panel_genes), "\n", sep = "")

# Per-subgroup mean expression profiles on panel genes
subgroup_profiles <- list()
for (sg in c("SHH", "WNT", "Group3", "Group4")) {
  gsms_in_sg <- merged_clin[Subgroup == sg, GSM]
  gsms_avail <- intersect(gsms_in_sg, colnames(exprs_mat))
  if (length(gsms_avail) < 5) next
  prof <- rowMeans(exprs_mat[panel_genes_on_bulk, gsms_avail, drop = FALSE], na.rm = TRUE)
  subgroup_profiles[[sg]] <- prof
}

# For each MB sample, get rho fingerprint on same gene set
fingerprint_match <- data.table(MB_sample = character(), Subgroup = character(),
                                spearman_rho = numeric(), spearman_p = numeric())
for (s in MB_SAMPLES) {
  dt_s <- mb_gene_data[[s]]
  dt_s <- dt_s[gene %in% panel_genes_on_bulk]
  setkey(dt_s, gene)
  fp_vals <- dt_s[panel_genes_on_bulk, rho_knn_main]
  for (sg in names(subgroup_profiles)) {
    ct <- tryCatch(cor.test(fp_vals, subgroup_profiles[[sg]], method = "spearman"),
                   error = function(e) NULL)
    if (!is.null(ct)) {
      fingerprint_match <- rbind(fingerprint_match,
        data.table(MB_sample = MB_NAMES[s], Subgroup = sg,
                   spearman_rho = ct$estimate, spearman_p = ct$p.value))
    }
  }
}
fwrite(fingerprint_match, file.path(OUT_DIR, "STAGE_B_fingerprint_match.csv"))

cat("  Fingerprint match results:\n")
print(dcast(fingerprint_match, MB_sample ~ Subgroup, value.var = "spearman_rho"))

# Heatmap
fp_wide <- dcast(fingerprint_match, MB_sample ~ Subgroup, value.var = "spearman_rho")
fp_mat <- as.matrix(fp_wide[, -1])
rownames(fp_mat) <- fp_wide$MB_sample

library(reshape2)
fp_long <- melt(fp_mat)
colnames(fp_long) <- c("MB_sample", "Subgroup", "rho")

p_heat <- ggplot(fp_long, aes(x = Subgroup, y = MB_sample, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 4) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(title = "Stage B: Xenium MB -> Cavalli subgroup match",
       subtitle = "Spearman rho of rho_knn_main fingerprint vs subgroup mean expression") +
  theme_minimal(base_size = 12)
ggsave(file.path(PLOT_DIR, "stage_B_fingerprint_heatmap.png"),
       p_heat, width = 6, height = 3.5, dpi = 200)
cat("  Saved plot: stage_B_fingerprint_heatmap.png\n")

# ============================================================================
# STAGE C: Leave-one-out signature stability
# ============================================================================
cat("\n================================================================\n")
cat("STAGE C: Leave-one-out signature stability\n")
cat("================================================================\n")

stage_c_results <- list()
for (excl in MB_SAMPLES) {
  sig_df <- sig_loo[[excl]]
  sc <- score_signature(sig_df, exprs_mat)
  if (is.null(sc)) next
  res_all <- run_5_tests(sc$score, merged_clin, label = paste0("excl_", MB_NAMES[excl]))
  res_shh <- run_5_tests(sc$score, merged_clin[Subgroup == "SHH"],
                         label = paste0("excl_", MB_NAMES[excl], "_SHH"))
  stage_c_results[[excl]] <- list(
    excluded = MB_NAMES[excl],
    n_genes = nrow(sig_df),
    n_matched = sc$n_matched,
    all = res_all,
    shh = res_shh
  )
}

# Also include full signature as baseline
sc_full <- score_strict
res_full_all <- run_5_tests(sc_full$score, merged_clin, label = "full")
res_full_shh <- run_5_tests(sc_full$score, merged_clin[Subgroup == "SHH"], label = "full_SHH")

stage_c_tbl <- rbindlist(lapply(c(list(full = list(
    excluded = "none",
    n_genes = nrow(sig_strict_full),
    n_matched = sc_full$n_matched,
    all = res_full_all, shh = res_full_shh
  )), stage_c_results), function(x) {
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
               A1_p = NA,  # single subgroup, ANOVA不适用
               B_logrank_p = x$shh$B_logrank_p %||% NA,
               Cox_HR = x$shh$Cox_HR %||% NA,
               Cox_p  = x$shh$Cox_p %||% NA,
               CoxAge_HR = x$shh$CoxAge_HR %||% NA,
               CoxAge_p = x$shh$CoxAge_p %||% NA)
  )
}))
fwrite(stage_c_tbl, file.path(OUT_DIR, "STAGE_C_leave_one_out.csv"))

cat("\n  Leave-one-out results:\n")
print(stage_c_tbl[scope == "ALL", .(excluded, n_genes, Cox_HR, Cox_p, CoxAge_p)])
cat("\n  Within SHH subgroup:\n")
print(stage_c_tbl[scope == "SHH", .(excluded, n_genes, Cox_HR, Cox_p, B_logrank_p)])

# ============================================================================
# STAGE D: 12-subtype signature profile
# ============================================================================
cat("\n================================================================\n")
cat("STAGE D: 12-subtype score profile\n")
cat("================================================================\n")

df_score <- data.table(GSM = names(score_strict$score),
                       score = as.numeric(score_strict$score))
df_score <- merge(df_score, merged_clin[, .(GSM, Subgroup, Subtype, Age, OS_years, Dead)],
                  by = "GSM")

# Box plot by subgroup (4 groups) - already have in R10, but include here with consistent style
df_score[, Subgroup := factor(Subgroup, levels = c("WNT", "SHH", "Group3", "Group4"))]
subtype_order <- df_score[, .(mean_score = mean(score, na.rm = TRUE)), by = Subtype][order(-mean_score), Subtype]
df_score[, Subtype := factor(Subtype, levels = subtype_order)]

# stats per subtype
sub_stats <- df_score[, .(n = .N, mean_score = mean(score, na.rm = TRUE),
                          median_score = median(score, na.rm = TRUE),
                          sd_score = sd(score, na.rm = TRUE)),
                      by = .(Subgroup, Subtype)]
fwrite(sub_stats[order(-mean_score)], file.path(OUT_DIR, "STAGE_D_subtype_stats.csv"))
cat("  Subtype score rankings:\n")
print(sub_stats[order(-mean_score)])

# plot: 12 subtype box plot, colored by subgroup
subgroup_colors <- c(WNT = "#66c2a5", SHH = "#fc8d62", Group3 = "#8da0cb", Group4 = "#e78ac3")
p_12 <- ggplot(df_score, aes(x = Subtype, y = score, fill = Subgroup)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  scale_fill_manual(values = subgroup_colors) +
  labs(x = "MB molecular subtype (sorted by mean score)",
       y = "Density signature score (z-score)",
       title = "Stage D: Signature score across 12 MB subtypes",
       subtitle = sprintf("n=%d; strict signature, 99 genes", nrow(df_score))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
ggsave(file.path(PLOT_DIR, "stage_D_12subtype_boxplot.png"),
       p_12, width = 10, height = 5, dpi = 200)
cat("  Saved plot: stage_D_12subtype_boxplot.png\n")

# ANOVA on 12 subtypes (already done in R10 but re-report)
a_subtype <- summary(aov(score ~ Subtype, data = df_score))[[1]][1, "Pr(>F)"]
cat(sprintf("  12-subtype ANOVA p = %.3g\n", a_subtype))

# ============================================================================
# FINAL SUMMARY
# ============================================================================
sink(file.path(OUT_DIR, "R13_SUMMARY.txt"))
cat("================================================================\n")
cat("  R13 MBEN-SHH Full Validation: A+B+C+D\n")
cat("  ", format(Sys.time()), "\n")
cat("================================================================\n\n")

cat("BACKGROUND\n")
cat("----------------------------------------------------------------\n")
cat("All 4 Xenium MB samples are MBEN histology (WHO: SHH-activated).\n")
cat("This R13 tests whether the signature's validation pattern is\n")
cat("consistent with this a priori SHH-enriched origin.\n\n")

cat("STAGE A: Subgroup-stratified validation (strict 99-gene signature)\n")
cat("----------------------------------------------------------------\n")
print(stage_a_tbl)
cat("\n")
sh <- stage_a_tbl[Subgroup == "SHH"]
if (nrow(sh) == 1) {
  cat(sprintf("  Key SHH numbers: n=%d, HR=%.3f, p=%.4f, CoxAge_p=%.4f\n",
              sh$n, sh$Cox_HR, sh$Cox_p, sh$CoxAge_p))
}

cat("\nSTAGE B: Xenium -> Cavalli subgroup matching\n")
cat("----------------------------------------------------------------\n")
fp_w <- dcast(fingerprint_match, MB_sample ~ Subgroup, value.var = "spearman_rho")
print(fp_w)
# best match per sample
best_match <- fingerprint_match[, .SD[which.max(spearman_rho)], by = MB_sample]
cat("\n  Best subgroup match per MB sample:\n")
print(best_match)

cat("\nSTAGE C: Leave-one-out stability\n")
cat("----------------------------------------------------------------\n")
cat("  ALL samples (full cohort):\n")
print(stage_c_tbl[scope == "ALL", .(excluded, n_genes, Cox_HR, Cox_p, CoxAge_p)])
cat("\n  Within SHH subgroup only:\n")
print(stage_c_tbl[scope == "SHH", .(excluded, n_genes, Cox_HR, Cox_p, B_logrank_p)])

cat("\nSTAGE D: 12-subtype profile\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  Overall 12-subtype ANOVA: p = %.3g\n", a_subtype))
cat("  Top 3 highest-score subtypes:\n")
print(head(sub_stats[order(-mean_score)], 3))
cat("  Bottom 3 lowest-score subtypes:\n")
print(tail(sub_stats[order(-mean_score)], 3))

cat("\n================================================================\n")
cat("  INTERPRETATION CHECKLIST\n")
cat("================================================================\n")
cat("\n[A] Subgroup-specificity:\n")
if (nrow(sh) == 1 && !is.na(sh$Cox_p)) {
  other_hrs <- stage_a_tbl[Subgroup %in% c("WNT", "Group3", "Group4") & !is.na(Cox_HR), Cox_HR]
  if (sh$Cox_HR < 0.5 && all(abs(log(other_hrs)) < abs(log(sh$Cox_HR)))) {
    cat("  >>> SHH-specific signal CONFIRMED (SHH has strongest HR)\n")
  } else {
    cat("  >>> Signal distribution mixed; report actual HRs per subgroup\n")
  }
}

cat("\n[B] Xenium-Cavalli concordance:\n")
shh_match_count <- sum(best_match$Subgroup == "SHH")
cat(sprintf("  %d / 4 MB samples best-match SHH subgroup\n", shh_match_count))
if (shh_match_count == 4) {
  cat("  >>> CONFIRMED: all MBEN samples fingerprint to SHH\n")
} else if (shh_match_count >= 2) {
  cat("  >>> Partial: majority MBEN samples match SHH\n")
} else {
  cat("  >>> Unexpected: MBEN samples do not match SHH; check fingerprint method\n")
}

cat("\n[C] Leave-one-out stability:\n")
loo_shh <- stage_c_tbl[scope == "SHH" & excluded != "none"]
loo_shh_pass <- sum(!is.na(loo_shh$B_logrank_p) & loo_shh$B_logrank_p < 0.05)
cat(sprintf("  %d / 4 LOO conditions retain SHH log-rank p<0.05\n", loo_shh_pass))
cat("  (original: SHH log-rank p = ", sh$Cox_p %||% NA, ")\n")

cat("\n[D] Subtype gradient:\n")
cat("  Top subtype:    ", as.character(sub_stats[order(-mean_score)][1, Subtype]), "\n", sep = "")
cat("  Bottom subtype: ", as.character(sub_stats[order(-mean_score)][.N, Subtype]), "\n", sep = "")
top_is_shh <- as.character(sub_stats[order(-mean_score)][1, Subgroup]) == "SHH"
if (top_is_shh) cat("  >>> Top-ranked subtype is SHH-derived (consistent with MBEN origin)\n")
sink()

cat("\n\nR13_SUMMARY.txt saved to: ", file.path(OUT_DIR, "R13_SUMMARY.txt"), "\n", sep = "")
cat("\nAll outputs:\n")
print(list.files(OUT_DIR, full.names = FALSE))
print(list.files(PLOT_DIR, full.names = FALSE))
cat("\nR13 complete.\n")
