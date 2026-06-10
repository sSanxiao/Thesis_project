#!/usr/bin/env Rscript
# ================================================================
#  R12: Close the 4 gaps on the R11 "strict" signature
#  1) Test C on strict (99 genes): bulk Cox direction vs Xenium direction
#  2) Reproduce the exact R10 signature formula as a baseline anchor
#  3) SHH_alpha standalone K-M (the driving subtype)
#  4) Sample-provenance of the 99 strict genes (MB266-dominated or not?)
# ================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(survival)
  library(ggplot2)
})

# -------- 路径 --------
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
R4_DIR      <- file.path(RESULTS_DIR, "R4_Results", "Medulloblastoma_Human")
GEO_DIR     <- file.path(EXTDATA_DIR, "Cavalli_GSE85217")
CLINICAL    <- file.path(GEO_DIR, "cavalli2017_mmc2_TableS1_clinical.csv")
R11_DIR     <- file.path(RESULTS_DIR, "R11_SHH_Deepdive")
OUT_DIR     <- file.path(RESULTS_DIR, "R12_Gaps")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "KM_plots"), showWarnings = FALSE)

MB_SAMPLES  <- c("GSM8840046", "GSM8840047", "GSM8840048", "GSM8840049")
MB266       <- "GSM8840047"

hline <- function(msg = "") {
  cat("\n"); cat(paste(rep("=", 64), collapse="")); cat("\n")
  if (nchar(msg) > 0) cat(msg, "\n")
  cat(paste(rep("=", 64), collapse="")); cat("\n")
}

# ================================================================
# 阶段 0: 读原材料 (R4 + Cavalli + R11 strict signature)
# ================================================================
hline("阶段 0: 读取材料")

# R4
r4_all <- list()
for (samp in MB_SAMPLES) {
  f <- file.path(R4_DIR, samp, "filtered_density_genes.csv")
  if (!file.exists(f)) { cat("  ! 缺", samp, "\n"); next }
  dt <- fread(f); dt[, sample := samp]
  r4_all[[samp]] <- dt
}
r4 <- rbindlist(r4_all, use.names = TRUE, fill = TRUE)
cat(sprintf("  R4: %d 样本, %d 行\n", length(r4_all), nrow(r4)))

# Cavalli bulk
cat("  读 GSE85217 (本地缓存)...\n")
gse <- getGEO(filename = file.path(GEO_DIR, "GSE85217_series_matrix.txt.gz"),
              getGPL = FALSE)
expr <- exprs(gse); pdat <- pData(gse)
cat(sprintf("    %d × %d\n", nrow(expr), ncol(expr)))

clin <- fread(CLINICAL)
pdat_dt <- data.table(gsm = rownames(pdat), title = pdat$title)
clin_merged <- merge(pdat_dt, clin, by.x = "title", by.y = "Study_ID", all.x = TRUE)
setnames(clin_merged, "OS (years)", "OS_years", skip_absent = TRUE)
setnames(clin_merged, "Dead", "event", skip_absent = TRUE)

# ENSG -> SYMBOL
probe_ids <- rownames(expr)
ensg_clean <- sub("_at$", "", probe_ids)
sym_map <- mapIds(org.Hs.eg.db, keys = ensg_clean, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
sym2probe <- data.table(probe = probe_ids, symbol = sym_map)
sym2probe <- sym2probe[!is.na(symbol)]; sym2probe <- sym2probe[!duplicated(symbol)]
cat(sprintf("  ENSG->SYMBOL: %d 可映射\n", nrow(sym2probe)))

# R11 strict signature
sig_strict_p <- fread(file.path(R11_DIR, "sig_strict_genes.csv"))
cat(sprintf("  R11 strict signature: %d 基因\n", nrow(sig_strict_p)))

# ================================================================
# 阶段 1: Test C on strict (99 基因的单基因 Cox 方向一致率)
# ================================================================
hline("阶段 1: Test C - 单基因 Cox 方向一致率 (strict, 99 基因)")

# 用有生存数据的全部样本 (n=612)
d_surv <- clin_merged[!is.na(OS_years) & !is.na(event)]
cat(sprintf("  生存数据样本: %d, events=%d\n", nrow(d_surv), sum(d_surv$event)))

# 每个 signature 基因单跑 Cox
test_c_results <- list()
for (i in seq_len(nrow(sig_strict_p))) {
  g     <- sig_strict_p$gene[i]
  pb    <- sig_strict_p$probe[i]
  xe_dir<- sig_strict_p$direction_final[i]
  e <- expr[pb, ]
  # 对齐到有生存的样本
  d_i <- copy(d_surv); d_i[, expr_val := e[match(d_i$gsm, names(e))]]
  d_i <- d_i[!is.na(expr_val)]

  cox <- tryCatch(coxph(Surv(OS_years, event) ~ expr_val, data = d_i),
                  error = function(e) NULL)
  if (is.null(cox)) next
  hr <- exp(coef(cox)); p <- summary(cox)$coefficients[1, "Pr(>|z|)"]
  # bulk 方向: HR>1 = hazard (高表达更差), HR<1 = protective (高表达更好)
  bulk_dir <- ifelse(hr > 1, "hazard", "protective")
  # Xenium 方向:
  #   positive = 高密度区该基因高表达; Xenium 里我们 signature 把 pos 基因 +z, neg 基因 -z
  #   所以 "signature 高 score 对应生存好" (HR<1) 时, 我们期望:
  #     positive 基因 -> bulk HR<1 (因为它"和高 score 同向", 高 score 是保护因子)
  #     negative 基因 -> bulk HR>1 (它"和高 score 反向", 所以它的"高表达"对应"低 score", 低 score 是危险)
  expected_bulk <- ifelse(xe_dir == "positive", "protective", "hazard")
  consistent <- (bulk_dir == expected_bulk)

  test_c_results[[i]] <- data.table(
    gene = g, probe = pb,
    xenium_direction = xe_dir,
    bulk_HR = hr, bulk_p = p,
    bulk_direction = bulk_dir,
    expected_bulk = expected_bulk,
    consistent = consistent
  )
}
tc <- rbindlist(test_c_results)
# FDR
tc[, q := p.adjust(bulk_p, method = "BH")]
setorder(tc, bulk_p)

cat(sprintf("\n  基因: %d\n", nrow(tc)))
cat(sprintf("  bulk p<0.05: %d (%.1f%%)\n",
            sum(tc$bulk_p < 0.05), 100*mean(tc$bulk_p < 0.05)))
cat(sprintf("  bulk q<0.05: %d (%.1f%%)\n",
            sum(tc$q < 0.05), 100*mean(tc$q < 0.05)))
cat(sprintf("  方向一致率 (所有基因): %d/%d = %.1f%%\n",
            sum(tc$consistent), nrow(tc), 100*mean(tc$consistent)))
cat(sprintf("  方向一致率 (p<0.05 子集): %d/%d = %.1f%%\n",
            sum(tc$consistent & tc$bulk_p < 0.05),
            sum(tc$bulk_p < 0.05),
            ifelse(sum(tc$bulk_p<0.05)>0,
                   100*mean(tc[bulk_p<0.05]$consistent), NA)))
cat(sprintf("  方向一致率 (q<0.05 子集): %d/%d = %.1f%%\n",
            sum(tc$consistent & tc$q < 0.05),
            sum(tc$q < 0.05),
            ifelse(sum(tc$q<0.05)>0,
                   100*mean(tc[q<0.05]$consistent), NA)))

# Test C 判定: 整体方向一致率 >=40% (原 R10 判定阈值)
test_c_pass_all <- mean(tc$consistent) >= 0.40
test_c_pass_sig <- ifelse(sum(tc$bulk_p<0.05)>=10,
                           mean(tc[bulk_p<0.05]$consistent) >= 0.40, NA)

cat(sprintf("\n  >>> Test C (全部基因): %s <<<\n",
            ifelse(test_c_pass_all, "PASS", "FAIL")))
if (!is.na(test_c_pass_sig)) {
  cat(sprintf("  >>> Test C (显著基因子集): %s <<<\n",
              ifelse(test_c_pass_sig, "PASS", "FAIL")))
}

cat("\n  Top 20 by bulk p:\n")
print(head(tc[, .(gene, xenium_direction, bulk_HR, bulk_p, q,
                   bulk_direction, consistent)], 20))

fwrite(tc, file.path(OUT_DIR, "TestC_strict_99genes.csv"))

# ================================================================
# 阶段 2: 复现 R10 原始 signature 作为 baseline anchor
# ================================================================
hline("阶段 2: 复现 R10 原始 signature (6 基因) 作为基准锚")

# R10 的筛选逻辑 (来自原脚本):
#   1. 取所有 MB 样本的 tier1_strict 记录
#   2. 按 gene 聚合
#   3. 要求 n_samples >= 2 且方向一致 >= 67%
# 这个严格比 R11 的任何一个都严格, 应该重现 6 基因

r10_cands <- r4[tier == "tier1_strict"]
agg10 <- r10_cands[, .(
  n_samples = .N,
  n_pos = sum(direction == "positive"),
  n_neg = sum(direction == "negative"),
  mean_rho = mean(rho_knn_main, na.rm = TRUE)
), by = gene]
agg10[, consistency := pmax(n_pos, n_neg) / n_samples]
agg10[, direction_final := ifelse(n_pos > n_neg, "positive",
                              ifelse(n_neg > n_pos, "negative", "tie"))]
sig_r10 <- agg10[n_samples >= 2 & consistency >= 0.67 & direction_final != "tie"]
cat(sprintf("  R10 原始配方: %d 基因 (R10 实际跑出 6 基因)\n", nrow(sig_r10)))
cat("  列表:\n")
print(sig_r10[, .(gene, n_samples, direction_final, mean_rho)])

# 跑 5 测试
sig_r10_p <- merge(sig_r10, sym2probe, by.x = "gene", by.y = "symbol")
cat(sprintf("  映射到 bulk: %d 基因\n", nrow(sig_r10_p)))

compute_score <- function(sig_p, expr) {
  probes <- sig_p$probe; dirs <- sig_p$direction_final
  M <- expr[probes, , drop = FALSE]
  Z <- t(scale(t(M)))
  sign_vec <- ifelse(dirs == "positive", 1, -1)
  colMeans(Z * sign_vec, na.rm = TRUE)
}
run_five_tests <- function(sig_p, clin_merged, label) {
  score <- compute_score(sig_p, expr)
  d <- copy(clin_merged); d$score <- score[match(d$gsm, names(score))]
  d <- d[!is.na(score) & !is.na(Subgroup)]

  a1 <- aov(score ~ Subgroup, data = d)
  A1_p <- summary(a1)[[1]][["Pr(>F)"]][1]
  d_sub <- d[!is.na(Subtype)]
  a2 <- aov(score ~ Subtype, data = d_sub)
  A2_p <- summary(a2)[[1]][["Pr(>F)"]][1]

  d_surv <- d[!is.na(OS_years) & !is.na(event)]
  d_surv[, score_hi := score > median(score, na.rm = TRUE)]
  sdiff <- survdiff(Surv(OS_years, event) ~ score_hi, data = d_surv)
  B_p <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)

  cox <- coxph(Surv(OS_years, event) ~ score, data = d_surv)
  Cox_p <- summary(cox)$coefficients[1, "Pr(>|z|)"]; Cox_HR <- exp(coef(cox))

  d_cox <- d_surv[!is.na(Age)]
  cox2 <- coxph(Surv(OS_years, event) ~ score + Age, data = d_cox)
  Cox_age_p <- summary(cox2)$coefficients["score", "Pr(>|z|)"]
  Cox_age_HR <- exp(coef(cox2))[1]

  data.table(label = label, n_genes = nrow(sig_p),
             A1_p = A1_p, A2_p = A2_p, B_p = B_p,
             Cox_p = Cox_p, Cox_HR = Cox_HR,
             Cox_age_p = Cox_age_p, Cox_age_HR = Cox_age_HR)
}

r10_tests <- run_five_tests(sig_r10_p, clin_merged, "R10_reproduced")
cat("\n  R10 复现 signature 的 5 测试:\n")
print(r10_tests)
fwrite(r10_tests, file.path(OUT_DIR, "R10_reproduced_tests.csv"))
fwrite(sig_r10_p, file.path(OUT_DIR, "R10_reproduced_genes.csv"))

# 汇总表: R10 复现 vs R11 strict/medium/loose
cat("\n  锚定对比表 (R10 复现 vs R11):\n")
r11_cmp <- fread(file.path(R11_DIR, "signature_comparison.csv"))
setnames(r10_tests, c("label","A1_p","A2_p","B_p","Cox_age_p","Cox_age_HR"),
         c("signature","A1_Subgroup_p","A2_Subtype_p","B_logrank_p",
           "Cox_ageadj_p","Cox_ageadj_HR"))
all_cmp <- rbind(r10_tests[, .(signature, n_genes, A1_Subgroup_p, A2_Subtype_p,
                                B_logrank_p, Cox_p, Cox_HR, Cox_ageadj_p, Cox_ageadj_HR)],
                 r11_cmp, fill = TRUE)
print(all_cmp)
fwrite(all_cmp, file.path(OUT_DIR, "ALL_signatures_comparison.csv"))

# ================================================================
# 阶段 3: SHH_alpha 单独 K-M (driving subtype)
# ================================================================
hline("阶段 3: SHH_alpha 单独 K-M")

# 用 strict signature
score_strict <- compute_score(sig_strict_p, expr)
d <- copy(clin_merged); d$score <- score_strict[match(d$gsm, names(score_strict))]

d_alpha <- d[Subgroup == "SHH" & Subtype == "SHH_alpha" &
             !is.na(score) & !is.na(OS_years) & !is.na(event)]
cat(sprintf("  SHH_alpha: n=%d, events=%d\n", nrow(d_alpha), sum(d_alpha$event)))

# 2 分位
d_alpha[, group2 := ifelse(score > median(score), "High", "Low")]
d_alpha[, group2 := factor(group2, levels = c("Low", "High"))]
sd2 <- survdiff(Surv(OS_years, event) ~ group2, data = d_alpha)
p2  <- 1 - pchisq(sd2$chisq, df = 1)
cox_a <- coxph(Surv(OS_years, event) ~ score, data = d_alpha)
hr_a  <- exp(coef(cox_a)); p_a <- summary(cox_a)$coefficients[1, "Pr(>|z|)"]
cat(sprintf("  2-part: log-rank p=%.4f, Cox HR=%.3f (p=%.4f)\n", p2, hr_a, p_a))

fit2 <- survfit(Surv(OS_years, event) ~ group2, data = d_alpha)
km_to_df <- function(fit) {
  s <- summary(fit, times = seq(0, max(fit$time, na.rm=TRUE), by = 0.1), extend = TRUE)
  data.frame(time = s$time, surv = s$surv, strata = s$strata)
}
df2 <- km_to_df(fit2)
g_a <- ggplot(df2, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c("group2=Low" = "#1f77b4", "group2=High" = "#d62728"),
                     labels = c("Low score", "High score")) +
  labs(title = "SHH_alpha — KM 2-partition (strict signature, 99 genes)",
       subtitle = sprintf("n=%d, events=%d | log-rank p=%.4f | Cox HR=%.3f (p=%.4f)",
                          nrow(d_alpha), sum(d_alpha$event), p2, hr_a, p_a),
       x = "Overall Survival (years)", y = "Survival probability", color = "") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"), legend.position = "top") +
  coord_cartesian(ylim = c(0, 1))
ggsave(file.path(OUT_DIR, "KM_plots", "SHH_alpha_KM_2partition_strict.png"),
       g_a, width = 7, height = 5, dpi = 150)

# 3 分位 (SHH_alpha n=51, 每组 ~17, 可能勉强)
q <- quantile(d_alpha$score, c(1/3, 2/3))
d_alpha[, group3 := cut(score, breaks = c(-Inf, q[1], q[2], Inf),
                         labels = c("Low", "Mid", "High"))]
sd3 <- survdiff(Surv(OS_years, event) ~ group3, data = d_alpha)
p3  <- 1 - pchisq(sd3$chisq, df = 2)
cat(sprintf("  3-part: log-rank p=%.4f\n", p3))

fit3 <- survfit(Surv(OS_years, event) ~ group3, data = d_alpha)
df3 <- km_to_df(fit3)
g_b <- ggplot(df3, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c("group3=Low" = "#1f77b4",
                                 "group3=Mid" = "#ff7f0e",
                                 "group3=High" = "#d62728"),
                     labels = c("Low", "Mid", "High")) +
  labs(title = "SHH_alpha — KM 3-partition (strict signature, 99 genes)",
       subtitle = sprintf("n=%d, events=%d | log-rank p=%.4f",
                          nrow(d_alpha), sum(d_alpha$event), p3),
       x = "Overall Survival (years)", y = "Survival probability", color = "Tertile") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"), legend.position = "top") +
  coord_cartesian(ylim = c(0, 1))
ggsave(file.path(OUT_DIR, "KM_plots", "SHH_alpha_KM_3partition_strict.png"),
       g_b, width = 7, height = 5, dpi = 150)

alpha_summary <- data.table(
  subtype = "SHH_alpha", n = nrow(d_alpha), events = sum(d_alpha$event),
  logrank_2p = p2, Cox_HR = hr_a, Cox_p = p_a, logrank_3p = p3
)
fwrite(alpha_summary, file.path(OUT_DIR, "SHH_alpha_summary.csv"))

# ================================================================
# 阶段 4: strict 99 基因的样本来源分解
# ================================================================
hline("阶段 4: 99 基因的样本来源分解 (是否 MB266 主导?)")

# 对每个 strict signature 基因, 查它在 4 个 MB 样本中哪些是 tier1_strict
r4_t1 <- r4[tier == "tier1_strict"]
r4_t1_wide <- dcast(r4_t1[, .(gene, sample, rho_knn_main, direction)],
                    gene ~ sample, value.var = c("rho_knn_main", "direction"))

provenance <- merge(sig_strict_p[, .(gene, direction_final, mean_rho, n_samples)],
                    r4_t1_wide, by = "gene", all.x = TRUE)
setorder(provenance, -n_samples, -abs(mean_rho))

# 统计: 每个基因在多少个 MB 样本里是 tier1
cat("\n  基因按出现样本数分布:\n")
print(table(provenance$n_samples))

# MB266 独占率: 基因只在 MB266 出现
mb266_only    <- sum(provenance$n_samples == 1 &
                     !is.na(provenance$rho_knn_main_GSM8840047))
other_single  <- sum(provenance$n_samples == 1) - mb266_only
multi_samples <- sum(provenance$n_samples >= 2)
cat(sprintf("\n  99 基因的样本来源:\n"))
cat(sprintf("    MB266 独占 (仅 GSM8840047): %d (%.0f%%)\n",
            mb266_only, 100*mb266_only/nrow(provenance)))
cat(sprintf("    其他单样本: %d (%.0f%%)\n",
            other_single, 100*other_single/nrow(provenance)))
cat(sprintf("    ≥2 样本共现: %d (%.0f%%)\n",
            multi_samples, 100*multi_samples/nrow(provenance)))

# 各样本贡献
for (samp in MB_SAMPLES) {
  rho_col <- paste0("rho_knn_main_", samp)
  n <- sum(!is.na(provenance[[rho_col]]))
  cat(sprintf("    %s tier1 贡献: %d 基因\n", samp, n))
}

fwrite(provenance, file.path(OUT_DIR, "sig_strict_99_provenance.csv"))

# 为论文做一个精简版: 只有关键列
prov_slim <- provenance[, .(
  gene, direction_final, n_samples, mean_rho,
  in_MB266     = !is.na(rho_knn_main_GSM8840047),
  in_GSM046    = !is.na(rho_knn_main_GSM8840046),
  in_GSM048    = !is.na(rho_knn_main_GSM8840048),
  in_GSM049    = !is.na(rho_knn_main_GSM8840049)
)]
fwrite(prov_slim, file.path(OUT_DIR, "sig_strict_99_provenance_slim.csv"))

# 画一个样本贡献的柱图
contrib <- data.table(
  sample = c("MB266_only", "GSM046_only", "GSM048_only", "GSM049_only",
             "≥2 samples"),
  count = c(
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840047)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840046)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840048)),
    sum(provenance$n_samples == 1 & !is.na(provenance$rho_knn_main_GSM8840049)),
    sum(provenance$n_samples >= 2)
  )
)
contrib[, sample := factor(sample, levels = sample)]  # 保持顺序

g_c <- ggplot(contrib, aes(x = sample, y = count, fill = sample)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("MB266_only" = "#d62728", "GSM046_only" = "#aec7e8",
                                "GSM048_only" = "#98df8a", "GSM049_only" = "#ffbb78",
                                "≥2 samples" = "#2ca02c"),
                    guide = "none") +
  labs(title = "Sample provenance of 99 strict-signature genes",
       subtitle = "How many genes come from each of the 4 MB Xenium samples?",
       x = "", y = "Number of genes") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT_DIR, "sig_strict_99_provenance_bar.png"),
       g_c, width = 7, height = 5, dpi = 150)

# ================================================================
# 阶段 5: 汇总报告
# ================================================================
hline("R12 Summary")
sink(file.path(OUT_DIR, "R12_SUMMARY.txt"), split = TRUE)

cat("================================================================\n")
cat(sprintf("  R12: Close the 4 gaps on the R11 strict signature\n"))
cat(sprintf("  %s\n", Sys.time()))
cat("================================================================\n")

cat("\nPART 1: Test C on strict (99 genes)\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  方向一致率 (全部 99 基因):      %d/%d = %.1f%%  [%s]\n",
            sum(tc$consistent), nrow(tc), 100*mean(tc$consistent),
            ifelse(test_c_pass_all, "PASS", "FAIL")))
if (sum(tc$bulk_p<0.05) >= 10) {
  cat(sprintf("  方向一致率 (bulk p<0.05 子集): %d/%d = %.1f%%  [%s]\n",
              sum(tc$consistent & tc$bulk_p < 0.05),
              sum(tc$bulk_p < 0.05),
              100*mean(tc[bulk_p<0.05]$consistent),
              ifelse(test_c_pass_sig, "PASS", "FAIL")))
}
cat(sprintf("  bulk p<0.05 基因数: %d (%.1f%%)\n",
            sum(tc$bulk_p<0.05), 100*mean(tc$bulk_p<0.05)))
cat(sprintf("  bulk q<0.05 基因数: %d (%.1f%%)\n",
            sum(tc$q<0.05), 100*mean(tc$q<0.05)))

cat("\nPART 2: R10 signature 复现作为 baseline anchor\n")
cat("----------------------------------------------------------------\n")
cat("4 种 signature 的全量对比:\n")
print(all_cmp)
cat("\n说明: R10_reproduced 应该与 R10_SUMMARY.txt 的结果接近\n")

cat("\nPART 3: SHH_alpha standalone K-M (driving subtype)\n")
cat("----------------------------------------------------------------\n")
print(alpha_summary)

cat("\nPART 4: 99 基因的样本来源\n")
cat("----------------------------------------------------------------\n")
print(contrib)
cat(sprintf("\n  MB266 独占占比: %.0f%%  ->  ", 100*mb266_only/nrow(provenance)))
if (mb266_only/nrow(provenance) >= 0.6) {
  cat("MB266 主导: 论文需明确声明\n")
} else if (mb266_only/nrow(provenance) >= 0.3) {
  cat("MB266 贡献大但非独占\n")
} else {
  cat("多样本共同支持, 不存在单样本主导\n")
}

cat("\n================================================================\n")
cat("  更新后的 5 测试判定 (strict signature):\n")
cat("----------------------------------------------------------------\n")
strict_row <- r11_cmp[signature == "strict"]
final_tests <- data.table(
  Test  = c("A1 Subgroup ANOVA", "A2 Subtype ANOVA", "B log-rank",
            "B' Cox age-adj", "C dir consistency"),
  Value = c(sprintf("p=%.2e", strict_row$A1_Subgroup_p),
            sprintf("p=%.2e", strict_row$A2_Subtype_p),
            sprintf("p=%.4f", strict_row$B_logrank_p),
            sprintf("p=%.4f", strict_row$Cox_ageadj_p),
            sprintf("%.1f%%", 100*mean(tc$consistent))),
  Verdict = c(
    ifelse(strict_row$A1_Subgroup_p < 0.001, "PASS", "FAIL"),
    ifelse(strict_row$A2_Subtype_p < 0.001, "PASS", "FAIL"),
    ifelse(strict_row$B_logrank_p < 0.05, "PASS", "FAIL"),
    ifelse(strict_row$Cox_ageadj_p < 0.05, "PASS", "FAIL"),
    ifelse(test_c_pass_all, "PASS", "FAIL")
  )
)
print(final_tests)
n_pass <- sum(final_tests$Verdict == "PASS")
cat(sprintf("\n  最终: %d / 5 PASS\n", n_pass))
cat(sprintf("  判定: "))
if (n_pass >= 4) cat("TIER I\n")
else if (n_pass == 3) cat("边缘 TIER I / 强 TIER II\n")
else cat("TIER II\n")

cat("\n输出文件 (",OUT_DIR,"):\n")
cat("  - TestC_strict_99genes.csv                 # Test C 详表\n")
cat("  - R10_reproduced_genes.csv                 # R10 原始 6 基因\n")
cat("  - R10_reproduced_tests.csv                 # R10 复现的 5 测试\n")
cat("  - ALL_signatures_comparison.csv            # 4 种 signature 对比总表\n")
cat("  - SHH_alpha_summary.csv                    # SHH_alpha 生存汇总\n")
cat("  - KM_plots/SHH_alpha_KM_{2,3}partition_strict.png\n")
cat("  - sig_strict_99_provenance.csv             # 99 基因样本来源完整\n")
cat("  - sig_strict_99_provenance_slim.csv        # 精简版 (布尔)\n")
cat("  - sig_strict_99_provenance_bar.png         # 贡献柱图\n")
cat("================================================================\n")
sink()

cat("\nR12 完成\n")
