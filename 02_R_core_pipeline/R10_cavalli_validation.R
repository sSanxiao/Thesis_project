#!/usr/bin/env Rscript
# ============================================================================
# R10_cavalli_validation.R  (主脚本 final 版)
#
# 功能: 用 Cavalli 2017 (GSE85217, 763 MB bulk + Table S1 临床+生存) 验证
#       R4 的 MB density signature
#
# 五个测试:
#   Test A1 - signature 能否区分 4 个 subgroup?                  ANOVA
#   Test A2 - signature 能否区分 12 个 subtype?                  ANOVA
#   Test B  - signature 能否分层生存 (log-rank)?
#   Test B-Cox - signature 连续变量 Cox (年龄调整)
#   Test C  - 单基因 Cox, 多少个方向与 Xenium 一致?
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
  library(survival)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ============================================================================
# 配置
# ============================================================================

# Configurable roots (see config/paths.R); override via environment variables.
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
R4_MB_DIR     <- file.path(RESULTS_DIR, "R4_Results", "Medulloblastoma_Human")
OUT_DIR       <- file.path(RESULTS_DIR, "R10_Cavalli")
EXT_DIR       <- file.path(EXTDATA_DIR, "Cavalli_GSE85217")
CLINICAL_CSV   <- file.path(EXT_DIR, "cavalli2017_mmc2_TableS1_clinical.csv")

MIN_SAMPLES_SHARED        <- 2
MIN_DIRECTION_CONSISTENCY <- 0.67

ANOVA_P_SUBGROUP_PASS    <- 0.001
ANOVA_P_SUBTYPE_PASS     <- 0.001
LOGRANK_P_PASS           <- 0.05
COX_P_PASS               <- 0.05
SINGLE_GENE_CONSISTENCY  <- 40

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat("R10: Cavalli 2017 bulk + clinical 验证\n")
cat("============================================\n")
cat("开始:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 阶段 1: 从 R4 构造 MB signature
# ============================================================================

cat("============================================\n")
cat("阶段 1: 构造 MB density signature\n")
cat("============================================\n")

mb_sample_dirs <- list.dirs(R4_MB_DIR, recursive = FALSE)
gene_records <- list()

for (d in mb_sample_dirs) {
  f <- file.path(d, "filtered_density_genes.csv")
  if (!file.exists(f)) next
  df <- fread(f)
  sname <- basename(d)
  t1 <- df[tier == "tier1_strict"]
  cat(sprintf("  %s: tier1_strict=%d (pos=%d, neg=%d)\n",
              sname, nrow(t1),
              sum(t1$direction == "positive"),
              sum(t1$direction == "negative")))
  if (nrow(t1) == 0) next
  for (i in seq_len(nrow(t1))) {
    gene_records[[length(gene_records) + 1]] <- list(
      gene = t1$gene[i], sample = sname,
      rho = t1$rho_knn_main[i], direction = t1$direction[i]
    )
  }
}

gi_df <- rbindlist(gene_records, fill = TRUE)
cat(sprintf("\n  总记录 %d, 唯一基因 %d, 样本 %d\n",
            nrow(gi_df), length(unique(gi_df$gene)), length(unique(gi_df$sample))))

gene_summary <- gi_df[, .(
  n_samples = .N,
  mean_rho = mean(rho, na.rm = TRUE),
  n_positive = sum(direction == "positive"),
  n_negative = sum(direction == "negative"),
  direction = {
    np <- sum(direction == "positive"); nn <- sum(direction == "negative")
    if (np > nn) "positive" else if (nn > np) "negative" else "mixed"
  },
  direction_consistency = {
    np <- sum(direction == "positive"); nn <- sum(direction == "negative")
    max(np, nn) / .N
  }
), by = gene]

sig_genes <- gene_summary[
  n_samples >= MIN_SAMPLES_SHARED &
  direction != "mixed" &
  direction_consistency >= MIN_DIRECTION_CONSISTENCY
][order(-n_samples, -abs(mean_rho))]

cat(sprintf("\n  最终 signature: %d 基因 (pos=%d, neg=%d)\n",
            nrow(sig_genes),
            sum(sig_genes$direction == "positive"),
            sum(sig_genes$direction == "negative")))
cat("  前 20:\n"); print(head(sig_genes, 20))
fwrite(sig_genes, file.path(OUT_DIR, "mb_density_signature_genes.csv"))

if (nrow(sig_genes) < 5) stop("signature 基因 < 5, 退出")

# ============================================================================
# 阶段 2: 读 GSE85217 + Table S1
# ============================================================================

cat("\n============================================\n")
cat("阶段 2: 读 GSE85217 + Table S1\n")
cat("============================================\n")

cat("  读 GSE85217 (本地缓存)...\n")
gse <- getGEO("GSE85217", destdir = EXT_DIR, getGPL = FALSE)
eset <- if (is.list(gse) && !is(gse, "ExpressionSet")) gse[[1]] else gse
cat(sprintf("    %d × %d\n", nrow(eset), ncol(eset)))

pd <- pData(eset)
pd$title_clean <- gsub('^"|"$', "", trimws(pd$title))

cat("\n  读 Table S1...\n")
# readxl: .name_repair='minimal' 保留原始列名含空格和括号
clin <- fread(CLINICAL_CSV)  # CSV 版本, 列名保留原始
cat(sprintf("    %d × %d\n", nrow(clin), ncol(clin)))
cat(sprintf("    原始列名前 10: %s\n",
            paste(head(names(clin), 10), collapse = " | ")))

# readxl 保留原名: "OS (years)", "Met status (1 Met, 0 M0)"
setnames(clin, old = "OS (years)",               new = "OS_years")
setnames(clin, old = "Met status (1 Met, 0 M0)",  new = "Met_status")

clin[, Age := as.numeric(Age)]
clin[, Dead := as.integer(Dead)]
clin[, OS_years := as.numeric(OS_years)]

cat(sprintf("\n    Dead 非NA: %d / %d\n", sum(!is.na(clin$Dead)), nrow(clin)))
cat(sprintf("    OS 非NA: %d / %d\n", sum(!is.na(clin$OS_years)), nrow(clin)))
cat(sprintf("    事件数 (Dead=1): %d\n", sum(clin$Dead == 1, na.rm = TRUE)))
cat("\n    Subgroup:\n"); print(table(clin$Subgroup, useNA = "ifany"))

# ============================================================================
# 阶段 3: 桥接
# ============================================================================

cat("\n============================================\n")
cat("阶段 3: GSM ID <-> Study_ID 桥接\n")
cat("============================================\n")

meta <- merge(
  data.table(gsm_id = rownames(pd),
             Study_ID = pd$title_clean,
             subgroup_gse = pd$`subgroup:ch1`,
             subtype_gse  = pd$`subtype:ch1`),
  clin, by = "Study_ID", all.x = TRUE
)

n_matched <- sum(!is.na(meta$Age))
cat(sprintf("  GSE85217 样本: %d, 匹配 Table S1: %d (%.1f%%)\n",
            nrow(pd), n_matched, 100 * n_matched / nrow(pd)))

cat("\n  一致性 (pData subgroup vs Table S1 Subgroup):\n")
print(table(GSE = meta$subgroup_gse, TableS1 = meta$Subgroup))

fwrite(meta, file.path(OUT_DIR, "cavalli_clinical_merged.csv"))

# ============================================================================
# 阶段 4: ENSG -> SYMBOL 映射
# ============================================================================

cat("\n============================================\n")
cat("阶段 4: ENSG -> SYMBOL 映射\n")
cat("============================================\n")

ensg_ids <- gsub("_at$", "", rownames(eset))
cat(sprintf("  ENSG IDs: %d\n", length(ensg_ids)))

ensg2sym <- suppressMessages(mapIds(
  org.Hs.eg.db, keys = ensg_ids, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"
))
cat(sprintf("    映射成功: %d / %d (%.1f%%)\n",
            sum(!is.na(ensg2sym)), length(ensg_ids),
            100 * sum(!is.na(ensg2sym)) / length(ensg_ids)))

sym_hits <- data.table(
  gene = sig_genes$gene,
  direction = sig_genes$direction,
  mean_rho = sig_genes$mean_rho,
  n_samples = sig_genes$n_samples,
  ensg = NA_character_,
  probe = NA_character_
)
for (i in seq_len(nrow(sym_hits))) {
  g <- sym_hits$gene[i]
  matched <- names(ensg2sym)[!is.na(ensg2sym) & ensg2sym == g]
  if (length(matched) == 0) next
  prb <- paste0(matched[1], "_at")
  if (prb %in% rownames(eset)) {
    sym_hits$ensg[i] <- matched[1]
    sym_hits$probe[i] <- prb
  }
}

n_sig_found <- sum(!is.na(sym_hits$probe))
cat(sprintf("\n  signature (%d) -> bulk probe: %d (%.1f%%)\n",
            nrow(sym_hits), n_sig_found, 100 * n_sig_found / nrow(sym_hits)))
cat("  未找到:", paste(sym_hits$gene[is.na(sym_hits$probe)], collapse = ", "), "\n")
fwrite(sym_hits, file.path(OUT_DIR, "ensg_symbol_map_used.csv"))

if (n_sig_found < 5) stop("匹配上的 signature 基因 < 5, 退出")

# ============================================================================
# 阶段 5: signature score
# ============================================================================

cat("\n============================================\n")
cat("阶段 5: signature score\n")
cat("============================================\n")

sig_hits <- sym_hits[!is.na(probe)]
expr_mat <- exprs(eset)
sub_expr <- expr_mat[sig_hits$probe, , drop = FALSE]

z_mat <- t(scale(t(sub_expr)))
z_mat[is.na(z_mat)] <- 0
neg_mask <- sig_hits$direction == "negative"
if (any(neg_mask)) z_mat[neg_mask, ] <- -z_mat[neg_mask, ]

sig_score <- colMeans(z_mat, na.rm = TRUE)
meta[, sig_score := sig_score[gsm_id]]
cat(sprintf("  score: n=%d, range=[%.3f, %.3f]\n",
            sum(!is.na(meta$sig_score)),
            min(meta$sig_score, na.rm = TRUE),
            max(meta$sig_score, na.rm = TRUE)))
fwrite(meta, file.path(OUT_DIR, "cavalli_clinical_merged.csv"))

# ============================================================================
# Test A: ANOVA (subgroup + subtype)
# ============================================================================

cat("\n============================================\n")
cat("Test A: ANOVA\n")
cat("============================================\n")

test_A_txt <- file(file.path(OUT_DIR, "test_A_anova.txt"), open = "w")

run_anova <- function(group_col, label) {
  m <- meta[!is.na(sig_score) & !is.na(get(group_col)) & get(group_col) != ""]
  if (nrow(m) < 20) {
    return(list(p = NA, n = nrow(m), means = NULL, data = NULL))
  }
  form <- as.formula(paste("sig_score ~", group_col))
  fit <- aov(form, data = m)
  p <- summary(fit)[[1]][["Pr(>F)"]][1]
  by_grp <- m[, .(mean = mean(sig_score), sd = sd(sig_score), n = .N),
              by = group_col]
  setorderv(by_grp, "mean")
  
  cat(sprintf("\n--- %s ---\n", label), file = test_A_txt)
  cat(sprintf("n=%d, ANOVA p=%.3e\n", nrow(m), p), file = test_A_txt)
  capture.output(print(by_grp), file = test_A_txt)
  
  cat(sprintf("  %s: p=%.3e, n=%d\n", label, p, nrow(m)))
  list(p = p, n = nrow(m), means = by_grp, data = m)
}

test_A1 <- run_anova("Subgroup", "Subgroup (k=4)")
test_A2 <- run_anova("Subtype", "Subtype (k=12)")
close(test_A_txt)

if (!is.null(test_A1$data)) {
  pA <- ggplot(test_A1$data, aes(x = Subgroup, y = sig_score, fill = Subgroup)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.2, alpha = 0.3) +
    theme_bw() +
    labs(title = sprintf("MB density signature across 4 subgroups\n(ANOVA p=%.2e, n=%d)",
                         test_A1$p, test_A1$n),
         y = "Signature score", x = "") +
    theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "test_A_subgroup_boxplot.pdf"),
         pA, width = 5, height = 4)
  
  pA2 <- ggplot(test_A2$data, aes(x = Subtype, y = sig_score, fill = Subgroup)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    theme_bw() +
    labs(title = sprintf("MB density signature across 12 subtypes\n(ANOVA p=%.2e, n=%d)",
                         test_A2$p, test_A2$n),
         y = "Signature score", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  ggsave(file.path(OUT_DIR, "test_A_subtype_boxplot.pdf"),
         pA2, width = 7, height = 5)
  cat("  ✓ boxplots saved\n")
}

# ============================================================================
# Test B: Survival
# ============================================================================

cat("\n============================================\n")
cat("Test B: Survival\n")
cat("============================================\n")

mB <- meta[!is.na(sig_score) & !is.na(OS_years) & !is.na(Dead) & OS_years > 0]
cat(sprintf("  n=%d, events=%d\n", nrow(mB), sum(mB$Dead == 1)))

test_B_result <- list(n = nrow(mB), n_events = sum(mB$Dead == 1))

if (nrow(mB) > 30 && sum(mB$Dead == 1) > 20) {
  mB[, sig_group := ifelse(sig_score > median(sig_score), "high", "low")]
  
  sv <- Surv(mB$OS_years, mB$Dead)
  sf <- survfit(sv ~ sig_group, data = mB)
  lr <- survdiff(sv ~ sig_group, data = mB)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  cox_c <- summary(coxph(sv ~ sig_score, data = mB))
  hr_c <- cox_c$coefficients[1, "exp(coef)"]
  p_c  <- cox_c$coefficients[1, "Pr(>|z|)"]
  
  mB_age <- mB[!is.na(Age)]
  cox_a <- summary(coxph(Surv(OS_years, Dead) ~ sig_score + Age, data = mB_age))
  hr_a <- cox_a$coefficients["sig_score", "exp(coef)"]
  p_a  <- cox_a$coefficients["sig_score", "Pr(>|z|)"]
  
  # per subgroup
  cox_grp <- list()
  for (sg in unique(mB$Subgroup)) {
    if (is.na(sg) || sg == "") next
    sub <- mB[Subgroup == sg]
    if (nrow(sub) < 30 || sum(sub$Dead == 1) < 10) next
    fit <- tryCatch(coxph(Surv(OS_years, Dead) ~ sig_score, data = sub),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      cox_grp[[sg]] <- data.frame(
        subgroup = sg, n = nrow(sub), events = sum(sub$Dead == 1),
        HR = round(s$coefficients[1, "exp(coef)"], 3),
        p = round(s$coefficients[1, "Pr(>|z|)"], 4)
      )
    }
  }
  cox_grp_df <- do.call(rbind, cox_grp)
  
  test_B_result$logrank_p <- lr_p
  test_B_result$hr_continuous <- hr_c
  test_B_result$p_continuous <- p_c
  test_B_result$hr_age_adj <- hr_a
  test_B_result$p_age_adj <- p_a
  test_B_result$per_subgroup <- cox_grp_df
  
  cat(sprintf("  log-rank p = %.3e\n", lr_p))
  cat(sprintf("  Cox continuous: HR=%.3f, p=%.3e\n", hr_c, p_c))
  cat(sprintf("  Cox age-adj:    HR=%.3f, p=%.3e\n", hr_a, p_a))
  cat("  per subgroup:\n"); print(cox_grp_df)
  
  sink(file.path(OUT_DIR, "test_B_cox.txt"))
  cat(sprintf("Test B: Survival\nn=%d, events=%d\n", nrow(mB), sum(mB$Dead == 1)))
  cat(sprintf("log-rank p = %.3e\n", lr_p))
  cat(sprintf("Cox continuous: HR=%.3f (p=%.3e)\n", hr_c, p_c))
  cat(sprintf("Cox age-adjusted: HR=%.3f (p=%.3e)\n", hr_a, p_a))
  cat("\nper subgroup:\n"); print(cox_grp_df)
  sink()
  
  pdf(file.path(OUT_DIR, "test_B_survival_KM.pdf"), width = 6, height = 5)
  par(mar = c(5, 5, 4, 2))
  plot(sf, col = c("red", "blue"), lwd = 2,
       xlab = "Time (years)", ylab = "Overall survival",
       main = sprintf("MB density signature vs OS\nlog-rank p=%.3e, n=%d (events=%d)",
                      lr_p, nrow(mB), sum(mB$Dead == 1)))
  legend("bottomleft",
         legend = c(sprintf("high (n=%d)", sum(mB$sig_group == "high")),
                    sprintf("low (n=%d)",  sum(mB$sig_group == "low"))),
         col = c("red", "blue"), lwd = 2, bty = "n")
  dev.off()
  cat("  ✓ KM saved\n")
}

# ============================================================================
# Test C: 单基因 Cox
# ============================================================================

cat("\n============================================\n")
cat("Test C: 单基因 Cox\n")
cat("============================================\n")

test_C_result <- list()
if (exists("mB") && nrow(mB) > 30 && sum(mB$Dead == 1) > 20) {
  gene_cox <- list()
  for (i in seq_len(nrow(sig_hits))) {
    g <- sig_hits$gene[i]; prb <- sig_hits$probe[i]; xdir <- sig_hits$direction[i]
    gvals <- as.numeric(expr_mat[prb, mB$gsm_id])
    if (all(is.na(gvals))) next
    fit <- tryCatch(coxph(Surv(mB$OS_years, mB$Dead) ~ gvals),
                    error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)$coefficients[1, ]
    hr <- s["exp(coef)"]; p <- s["Pr(>|z|)"]
    
    # 方向一致性:
    # Xenium +  (高密度高表达) & bulk HR>1 (高表达预后差) → 一致
    # Xenium -  (高密度低表达) & bulk HR<1 (高表达预后好) → 一致
    consistent <- (xdir == "positive" && hr > 1) ||
                  (xdir == "negative" && hr < 1)
    
    gene_cox[[g]] <- data.frame(
      gene = g, probe = prb, xenium_direction = xdir,
      bulk_HR = round(hr, 3), bulk_p = round(p, 4),
      bulk_direction = if (hr > 1) "hazard" else "protective",
      direction_consistent = consistent
    )
  }
  cox_df <- rbindlist(gene_cox, fill = TRUE)
  cox_df[, q := p.adjust(bulk_p, method = "BH")]
  setorder(cox_df, bulk_p)
  fwrite(cox_df, file.path(OUT_DIR, "test_C_single_gene_cox.csv"))
  
  test_C_result$n_total <- nrow(cox_df)
  test_C_result$n_p_sig <- sum(cox_df$bulk_p < 0.05, na.rm = TRUE)
  test_C_result$n_q_sig <- sum(cox_df$q < 0.05, na.rm = TRUE)
  test_C_result$n_consistent <- sum(cox_df$direction_consistent, na.rm = TRUE)
  test_C_result$pct_consistent <- 100 * test_C_result$n_consistent / nrow(cox_df)
  test_C_result$n_p_sig_consistent <- sum(
    cox_df$bulk_p < 0.05 & cox_df$direction_consistent, na.rm = TRUE)
  
  cat(sprintf("  基因: %d, p<0.05: %d, q<0.05: %d\n",
              test_C_result$n_total, test_C_result$n_p_sig, test_C_result$n_q_sig))
  cat(sprintf("  方向一致: %d (%.1f%%), p<0.05 且一致: %d\n",
              test_C_result$n_consistent, test_C_result$pct_consistent,
              test_C_result$n_p_sig_consistent))
  cat("  Top 10 by p:\n"); print(head(cox_df, 10))
}

# ============================================================================
# 阶段 6: 最终 Tier 判断
# ============================================================================

pass_A1 <- !is.na(test_A1$p) && test_A1$p < ANOVA_P_SUBGROUP_PASS
pass_A2 <- !is.na(test_A2$p) && test_A2$p < ANOVA_P_SUBTYPE_PASS
pass_B  <- !is.null(test_B_result$logrank_p) &&
           !is.na(test_B_result$logrank_p) &&
           test_B_result$logrank_p < LOGRANK_P_PASS
pass_B_cox <- !is.null(test_B_result$p_age_adj) &&
              !is.na(test_B_result$p_age_adj) &&
              test_B_result$p_age_adj < COX_P_PASS
pass_C  <- !is.null(test_C_result$pct_consistent) &&
           test_C_result$pct_consistent >= SINGLE_GENE_CONSISTENCY

n_pass <- sum(pass_A1, pass_A2, pass_B, pass_B_cox, pass_C)

verdict <- if (n_pass >= 4) {
  "TIER I LOCAL CHAPTER CONFIRMED"
} else if (n_pass == 3) {
  "STRONG TIER II / borderline TIER I"
} else if (n_pass == 2) {
  "TIER II CONFIRMED"
} else if (n_pass == 1) {
  "TIER II (weak)"
} else {
  "TIER III"
}

sink(file.path(OUT_DIR, "R10_SUMMARY.txt"))

cat("================================================================\n")
cat("  R10: MB density signature vs Cavalli 2017 bulk + clinical\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("================================================================\n\n")
cat("FINAL VERDICT:\n  ", verdict, "\n\n")
cat(sprintf("Tests passed: %d / 5\n\n", n_pass))
cat(sprintf("  Test A1 (Subgroup ANOVA, p<%.3f):         %s  [p=%.2e]\n",
            ANOVA_P_SUBGROUP_PASS,
            if (pass_A1) "PASS" else "FAIL",
            if (is.na(test_A1$p)) NaN else test_A1$p))
cat(sprintf("  Test A2 (Subtype ANOVA, p<%.3f):          %s  [p=%.2e]\n",
            ANOVA_P_SUBTYPE_PASS,
            if (pass_A2) "PASS" else "FAIL",
            if (is.na(test_A2$p)) NaN else test_A2$p))
cat(sprintf("  Test B (log-rank, p<%.2f):                %s  [p=%s]\n",
            LOGRANK_P_PASS,
            if (pass_B) "PASS" else "FAIL",
            if (is.null(test_B_result$logrank_p)) "NA"
            else format(test_B_result$logrank_p, digits = 3, scientific = TRUE)))
cat(sprintf("  Test B-Cox age-adj (p<%.2f):              %s  [p=%s]\n",
            COX_P_PASS,
            if (pass_B_cox) "PASS" else "FAIL",
            if (is.null(test_B_result$p_age_adj)) "NA"
            else format(test_B_result$p_age_adj, digits = 3, scientific = TRUE)))
cat(sprintf("  Test C (single-gene dir >=%d%%):           %s  [%.1f%%]\n",
            SINGLE_GENE_CONSISTENCY,
            if (pass_C) "PASS" else "FAIL",
            if (is.null(test_C_result$pct_consistent)) NA
            else test_C_result$pct_consistent))

cat("\n----------------------------------------------------------------\n")
cat("阶段 1 - MB signature\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  R4 tier1_strict 记录: %d, 唯一基因 %d\n",
            nrow(gi_df), length(unique(gi_df$gene))))
cat(sprintf("  通过筛: %d (pos=%d, neg=%d)\n",
            nrow(sig_genes),
            sum(sig_genes$direction == "positive"),
            sum(sig_genes$direction == "negative")))
cat(sprintf("  映射到 bulk: %d\n", n_sig_found))
cat("\n  Top 20:\n"); print(head(sig_genes, 20))

cat("\n----------------------------------------------------------------\n")
cat("阶段 2-3 - Cavalli 数据\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  GSE85217: %d, 匹配 S1: %d\n", nrow(pd), n_matched))
cat(sprintf("  有生存数据: %d (events %d)\n",
            sum(!is.na(meta$Dead) & !is.na(meta$OS_years)),
            sum(meta$Dead == 1, na.rm = TRUE)))
cat("\n  Subgroup 分布:\n"); print(table(meta$Subgroup, useNA = "ifany"))

cat("\n----------------------------------------------------------------\n")
cat("Test A - ANOVA\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  A1 Subgroup:  n=%d, p=%.3e\n", test_A1$n, test_A1$p))
if (!is.null(test_A1$means)) print(test_A1$means)
cat(sprintf("\n  A2 Subtype:   n=%d, p=%.3e\n", test_A2$n, test_A2$p))
if (!is.null(test_A2$means)) print(test_A2$means)

cat("\n----------------------------------------------------------------\n")
cat("Test B - Survival\n")
cat("----------------------------------------------------------------\n")
if (!is.null(test_B_result$logrank_p)) {
  cat(sprintf("  n=%d, events=%d\n", test_B_result$n, test_B_result$n_events))
  cat(sprintf("  log-rank p = %.3e\n", test_B_result$logrank_p))
  cat(sprintf("  Cox continuous: HR=%.3f (p=%.3e)\n",
              test_B_result$hr_continuous, test_B_result$p_continuous))
  cat(sprintf("  Cox age-adj: HR=%.3f (p=%.3e)\n",
              test_B_result$hr_age_adj, test_B_result$p_age_adj))
  if (!is.null(test_B_result$per_subgroup)) {
    cat("\n  per subgroup:\n"); print(test_B_result$per_subgroup)
  }
}

cat("\n----------------------------------------------------------------\n")
cat("Test C - 单基因 Cox\n")
cat("----------------------------------------------------------------\n")
if (!is.null(test_C_result$n_total)) {
  cat(sprintf("  基因: %d, p<0.05: %d, q<0.05: %d\n",
              test_C_result$n_total, test_C_result$n_p_sig, test_C_result$n_q_sig))
  cat(sprintf("  方向一致: %d (%.1f%%), p<0.05 且一致: %d\n",
              test_C_result$n_consistent, test_C_result$pct_consistent,
              test_C_result$n_p_sig_consistent))
  cat("\n  Top 20 by p:\n")
  top_df <- head(fread(file.path(OUT_DIR, "test_C_single_gene_cox.csv")), 20)
  print(top_df)
}

cat("\n================================================================\n")
cat("Tier 判断:\n")
cat("  5/5 或 4/5 -> TIER I, MB 独立章节\n")
cat("  3/5        -> 边缘 TIER I / 强 TIER II\n")
cat("  2/5        -> TIER II\n")
cat("  0-1/5      -> TIER III\n")
cat("================================================================\n")

sink()

cat("\n============================================\n")
cat("R10 完成, 结束:", format(Sys.time()), "\n")
cat("============================================\n")
cat("\n关键输出 (贴给 Claude):\n")
cat("  ", file.path(OUT_DIR, "R10_SUMMARY.txt"), "\n")
