# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ================================================================
#  R12_patch: 修 R12 阶段 4 的 setorder bug + SHH_alpha Cox 稳健化
#
#  R12 阶段 1/2/3 的结果已经写盘 (TestC / R10_reproduced / SHH_alpha KM),
#  只需要补跑:
#    a) 阶段 4: 99 基因样本来源分解
#    b) SHH_alpha 的稳健 Cox (解决 HR=0 的 separation 问题)
#    c) 最终 SUMMARY
# ================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

R4_DIR  <- "./results/R4_Results/Medulloblastoma_Human"
R11_DIR <- "./results/R11_SHH_Deepdive"
R12_DIR <- "./results/R12_Gaps"
GEO_DIR <- "./external_data/Cavalli_GSE85217"
CLINICAL <- file.path(GEO_DIR, "cavalli2017_mmc2_TableS1_clinical.csv")

MB_SAMPLES <- c("GSM8840046", "GSM8840047", "GSM8840048", "GSM8840049")

hline <- function(msg = "") {
  cat("\n"); cat(paste(rep("=", 64), collapse="")); cat("\n")
  if (nchar(msg) > 0) cat(msg, "\n")
  cat(paste(rep("=", 64), collapse="")); cat("\n")
}

# ================================================================
# 读已有材料
# ================================================================
sig_strict_p <- fread(file.path(R11_DIR, "sig_strict_genes.csv"))
cat(sprintf("Strict signature: %d 基因\n", nrow(sig_strict_p)))

# R4
r4_all <- list()
for (samp in MB_SAMPLES) {
  f <- file.path(R4_DIR, samp, "filtered_density_genes.csv")
  if (!file.exists(f)) next
  dt <- fread(f); dt[, sample := samp]
  r4_all[[samp]] <- dt
}
r4 <- rbindlist(r4_all, use.names = TRUE, fill = TRUE)

# ================================================================
# 阶段 4 (修复): 99 基因的样本来源分解
# ================================================================
hline("阶段 4: 99 基因样本来源分解")

r4_t1 <- r4[tier == "tier1_strict"]
r4_t1_wide <- dcast(r4_t1[, .(gene, sample, rho_knn_main, direction)],
                    gene ~ sample, value.var = c("rho_knn_main", "direction"))

provenance <- merge(sig_strict_p[, .(gene, direction_final, mean_rho, n_samples)],
                    r4_t1_wide, by = "gene", all.x = TRUE)

# 修 bug: 用已计算的列做排序
provenance[, abs_rho := abs(mean_rho)]
setorder(provenance, -n_samples, -abs_rho)
provenance[, abs_rho := NULL]

cat("\n  基因按出现样本数分布:\n")
print(table(provenance$n_samples))

# 各样本贡献: "在本样本 tier1" 的基因数 (一个基因可在多个样本同时出现)
cat("\n  各样本 tier1 贡献 (可重叠):\n")
for (samp in MB_SAMPLES) {
  rho_col <- paste0("rho_knn_main_", samp)
  n_total    <- sum(!is.na(provenance[[rho_col]]))
  n_unique   <- sum(!is.na(provenance[[rho_col]]) & provenance$n_samples == 1)
  cat(sprintf("    %s: %3d 基因 tier1 (其中 %3d 只出现在此样本)\n",
              samp, n_total, n_unique))
}

# MB266 独占占比
mb266_col  <- "rho_knn_main_GSM8840047"
mb266_only <- sum(provenance$n_samples == 1 & !is.na(provenance[[mb266_col]]))
multi      <- sum(provenance$n_samples >= 2)
other_one  <- sum(provenance$n_samples == 1) - mb266_only

cat(sprintf("\n  分解:\n"))
cat(sprintf("    仅 MB266 (GSM8840047):  %3d (%.0f%%)\n",
            mb266_only, 100*mb266_only/nrow(provenance)))
cat(sprintf("    其他单样本:              %3d (%.0f%%)\n",
            other_one, 100*other_one/nrow(provenance)))
cat(sprintf("    >=2 样本共现:            %3d (%.0f%%)\n",
            multi, 100*multi/nrow(provenance)))

# 保存
fwrite(provenance, file.path(R12_DIR, "sig_strict_99_provenance.csv"))

prov_slim <- provenance[, .(
  gene, direction_final, n_samples, mean_rho,
  in_GSM046 = !is.na(rho_knn_main_GSM8840046),
  in_MB266  = !is.na(rho_knn_main_GSM8840047),
  in_GSM048 = !is.na(rho_knn_main_GSM8840048),
  in_GSM049 = !is.na(rho_knn_main_GSM8840049)
)]
fwrite(prov_slim, file.path(R12_DIR, "sig_strict_99_provenance_slim.csv"))

# 柱图
contrib <- data.table(
  sample = c("MB266_only", "GSM046_only", "GSM048_only", "GSM049_only",
             "Multi (>=2)"),
  count  = c(
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840047)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840046)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840048)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840049)),
    sum(provenance$n_samples >= 2)
  )
)
contrib[, sample := factor(sample, levels = sample)]

g_c <- ggplot(contrib, aes(x = sample, y = count, fill = sample)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("MB266_only"   = "#d62728",
                                "GSM046_only" = "#aec7e8",
                                "GSM048_only" = "#98df8a",
                                "GSM049_only" = "#ffbb78",
                                "Multi (>=2)" = "#2ca02c"),
                    guide = "none") +
  labs(title = "Sample provenance of 99 strict-signature genes",
       subtitle = "How many genes come from each of the 4 MB Xenium samples?",
       x = "", y = "Number of genes") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(R12_DIR, "sig_strict_99_provenance_bar.png"),
       g_c, width = 7, height = 5, dpi = 150)

# ================================================================
# 补: SHH_alpha 的稳健 Cox
# ================================================================
hline("SHH_alpha 稳健 Cox (解决 HR=0.000 的 separation)")

# 重新算 score (不想再读 Cavalli, 直接读 R11 保存的 strict score)
score_dt <- fread(file.path(R11_DIR, "scores_strict.csv"))

# 读 clinical 做对接
clin <- fread(CLINICAL)
# 需要 GSM<->Study_ID 桥接; R11 保存的 scores 已经是 gsm keyed
# 从 clinical 再补桥接
pdat_csv <- NULL
# 偷懒: 直接从 R11 输出里拿 subgroup / subtype / OS / event
# 但 R11 scores 表只有 gsm 和 score, 所以需要重新 merge
# 用 R10 的 gsm<->Study_ID 方式 (GEO title = Study_ID)

# 先 GEO 的 pData 里读 title
library(GEOquery); library(Biobase)
gse <- getGEO(filename = file.path(GEO_DIR, "GSE85217_series_matrix.txt.gz"),
              getGPL = FALSE)
pdat <- pData(gse)
pdat_dt <- data.table(gsm = rownames(pdat), title = pdat$title)

d <- merge(score_dt, pdat_dt, by = "gsm")
d <- merge(d, clin, by.x = "title", by.y = "Study_ID", all.x = TRUE)
setnames(d, "OS (years)", "OS_years", skip_absent = TRUE)
setnames(d, "Dead", "event", skip_absent = TRUE)

d_alpha <- d[Subgroup == "SHH" & Subtype == "SHH_alpha" &
             !is.na(score) & !is.na(OS_years) & !is.na(event)]

# 事件分布 (诊断 separation)
d_alpha[, group2 := ifelse(score > median(score), "High", "Low")]
d_alpha[, group2 := factor(group2, levels = c("Low", "High"))]
cat(sprintf("  SHH_alpha n=%d, events=%d\n", nrow(d_alpha), sum(d_alpha$event)))
cat("\n  事件 × 组分布:\n")
print(table(group = d_alpha$group2, event = d_alpha$event))

# 标准 log-rank + Cox
sd2  <- survdiff(Surv(OS_years, event) ~ group2, data = d_alpha)
p_lr <- 1 - pchisq(sd2$chisq, df = 1)
cox0 <- coxph(Surv(OS_years, event) ~ score, data = d_alpha)
hr0  <- exp(coef(cox0)); p0 <- summary(cox0)$coefficients[1, "Pr(>|z|)"]
ci0  <- exp(confint(cox0))

cat(sprintf("\n  标准 Cox (可能有 separation):\n"))
cat(sprintf("    log-rank p = %.4f\n", p_lr))
cat(sprintf("    Cox HR     = %.4f  (95%% CI: %.4f - %.4f)\n", hr0, ci0[1], ci0[2]))
cat(sprintf("    Cox p      = %.4f\n", p0))

# Firth 校正 Cox (处理 separation) - 用 survival::coxph 不支持,
# 需要 coxphf 或 logistf 包. 看看装了没
has_coxphf <- requireNamespace("coxphf", quietly = TRUE)
if (has_coxphf) {
  cat("\n  Firth 校正 Cox (coxphf 包):\n")
  library(coxphf)
  fit_firth <- coxphf(Surv(OS_years, event) ~ score, data = d_alpha)
  cat(sprintf("    Firth HR = %.4f  (95%% CI: %.4f - %.4f)\n",
              exp(fit_firth$coefficients),
              exp(fit_firth$ci.lower), exp(fit_firth$ci.upper)))
  cat(sprintf("    Firth p  = %.4f\n", fit_firth$prob))
} else {
  cat("\n  coxphf 包未装, 改用其他方式报稳健效应量:\n")
  # 手工: 如果 High 组 events=0, 用 "median survival ratio" 和 log-rank 做报告
  evt_lo <- sum(d_alpha[group2=="Low"]$event)
  evt_hi <- sum(d_alpha[group2=="High"]$event)
  n_lo   <- sum(d_alpha$group2=="Low")
  n_hi   <- sum(d_alpha$group2=="High")
  cat(sprintf("    Low 组:  n=%d, events=%d (事件率 %.0f%%)\n",
              n_lo, evt_lo, 100*evt_lo/n_lo))
  cat(sprintf("    High 组: n=%d, events=%d (事件率 %.0f%%)\n",
              n_hi, evt_hi, 100*evt_hi/n_hi))
  if (evt_hi == 0) {
    cat("    High 组零事件 -> Cox HR 不可识别 (separation),\n")
    cat("    但 log-rank p=%.4f 仍有效; 论文报 log-rank + 事件率对比\n")
  }
  # 先 binning 连续 score 到三分位再做 Cox, 避免 separation
  q <- quantile(d_alpha$score, c(1/3, 2/3))
  d_alpha[, group3 := cut(score, breaks = c(-Inf, q[1], q[2], Inf),
                           labels = c("Low", "Mid", "High"))]
  cat("\n  3 分位替代 Cox (减少 separation):\n")
  cox3 <- coxph(Surv(OS_years, event) ~ group3, data = d_alpha)
  print(summary(cox3))
}

# Score 在两组的分布 (诊断)
cat("\n  Score 分布 (按组):\n")
print(d_alpha[, .(n = .N, mean = mean(score), sd = sd(score),
                   min = min(score), max = max(score)), by = group2])

# ================================================================
# 最终汇总
# ================================================================
hline("R12 Final Summary")
sink(file.path(R12_DIR, "R12_SUMMARY_final.txt"), split = TRUE)

cat("================================================================\n")
cat(sprintf("  R12 FINAL: strict signature full validation\n"))
cat(sprintf("  %s\n", Sys.time()))
cat("================================================================\n")

# Part 1 结果 (从 TestC CSV 读)
tc <- fread(file.path(R12_DIR, "TestC_strict_99genes.csv"))
consistency_all <- mean(tc$consistent)
consistency_sig <- ifelse(sum(tc$bulk_p<0.05)>0,
                           mean(tc[bulk_p<0.05]$consistent), NA)

cat("\nPART 1: Test C (单基因方向一致率)\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  全部 99 基因:          %d/%d = %.1f%%  [%s]\n",
            sum(tc$consistent), nrow(tc), 100*consistency_all,
            ifelse(consistency_all >= 0.40, "PASS", "FAIL")))
cat(sprintf("  bulk p<0.05 子集:      %d/%d = %.1f%%\n",
            sum(tc$consistent & tc$bulk_p<0.05),
            sum(tc$bulk_p<0.05), 100*consistency_sig))
cat(sprintf("  bulk q<0.05 基因数:    %d\n", sum(tc$q<0.05)))

# Part 2 (R10 复现)
cat("\nPART 2: 4 种 signature 对比\n")
cat("----------------------------------------------------------------\n")
if (file.exists(file.path(R12_DIR, "ALL_signatures_comparison.csv"))) {
  all_cmp <- fread(file.path(R12_DIR, "ALL_signatures_comparison.csv"))
  print(all_cmp)
}

# Part 3 (SHH_alpha)
cat("\nPART 3: SHH_alpha (driving subtype)\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  n=%d, events=%d\n", nrow(d_alpha), sum(d_alpha$event)))
cat(sprintf("  log-rank (2-part):  p = %.4f\n", p_lr))
cat(sprintf("  事件率: Low %d/%d (%.0f%%)  vs  High %d/%d (%.0f%%)\n",
            sum(d_alpha[group2=="Low"]$event),  sum(d_alpha$group2=="Low"),
            100*sum(d_alpha[group2=="Low"]$event)/sum(d_alpha$group2=="Low"),
            sum(d_alpha[group2=="High"]$event), sum(d_alpha$group2=="High"),
            100*sum(d_alpha[group2=="High"]$event)/sum(d_alpha$group2=="High")))
if (sum(d_alpha[group2=="High"]$event) == 0) {
  cat(sprintf("  >>> High 组零事件, Cox HR 不可识别 -> 报 log-rank + 事件率对比 <<<\n"))
}

# Part 4 (provenance)
cat("\nPART 4: 99 基因样本来源\n")
cat("----------------------------------------------------------------\n")
print(contrib)
cat(sprintf("\n  MB266 独占 %.0f%% | 其他单样本 %.0f%% | 多样本 %.0f%%\n",
            100*mb266_only/nrow(provenance),
            100*other_one/nrow(provenance),
            100*multi/nrow(provenance)))

# 最终 5 测试表
cat("\n================================================================\n")
cat("  最终 5 测试判定 (strict signature, 99 基因)\n")
cat("================================================================\n")
r11_cmp <- fread(file.path(R11_DIR, "signature_comparison.csv"))
s <- r11_cmp[signature == "strict"]
final_tests <- data.table(
  Test    = c("A1 Subgroup ANOVA", "A2 Subtype ANOVA", "B log-rank (all, n=612)",
              "B' Cox age-adj", "C dir consistency"),
  Value   = c(sprintf("p=%.2e", s$A1_Subgroup_p),
              sprintf("p=%.2e", s$A2_Subtype_p),
              sprintf("p=%.4f", s$B_logrank_p),
              sprintf("p=%.4f", s$Cox_ageadj_p),
              sprintf("%.1f%% (%d/%d)", 100*consistency_all,
                      sum(tc$consistent), nrow(tc))),
  Verdict = c(ifelse(s$A1_Subgroup_p  < 0.001, "PASS", "FAIL"),
              ifelse(s$A2_Subtype_p   < 0.001, "PASS", "FAIL"),
              ifelse(s$B_logrank_p    < 0.05,  "PASS", "FAIL"),
              ifelse(s$Cox_ageadj_p   < 0.05,  "PASS", "FAIL"),
              ifelse(consistency_all  >= 0.40, "PASS", "FAIL"))
)
print(final_tests)
n_pass <- sum(final_tests$Verdict == "PASS")
cat(sprintf("\n  ==>  %d / 5 PASS  ==>  ", n_pass))
if (n_pass == 5)      cat("TIER I (完全确认)\n")
else if (n_pass == 4) cat("TIER I\n")
else if (n_pass == 3) cat("边缘 TIER I / 强 TIER II\n")
else                  cat("TIER II\n")

cat("\n附加亮点:\n")
cat(sprintf("  - SHH subgroup HR=0.36 (R10 已报, n=172)\n"))
cat(sprintf("  - SHH_alpha log-rank p=%.4f (n=%d, R12 新增)\n",
            p_lr, nrow(d_alpha)))
cat("================================================================\n")
sink()

cat("\nR12 patch 完成.\n")
cat("关键输出:\n")
cat("  ", file.path(R12_DIR, "R12_SUMMARY_final.txt"), "\n")
cat("  ", file.path(R12_DIR, "sig_strict_99_provenance_bar.png"), "\n")
