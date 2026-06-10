#!/usr/bin/env Rscript
# ================================================================
#  R11: Signature Expansion + SHH Deep-dive
#  - R11.1: 三种 signature 策略, 在 Cavalli 上重跑 5 测试
#  - R11.2: SHH (n=172) 内部 2/3 分位 K-M, 每个 signature 一套
#  - R11.3: SHH 四个亚亚型 (alpha/beta/gamma/delta) 分层 Cox
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
OUT_DIR     <- file.path(RESULTS_DIR, "R11_SHH_Deepdive")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "KM_plots"), showWarnings = FALSE)

MB_SAMPLES  <- c("GSM8840046", "GSM8840047", "GSM8840048", "GSM8840049")

# 方便的分隔符
hline <- function(msg = "") {
  cat("\n"); cat(paste(rep("=", 64), collapse="")); cat("\n")
  if (nchar(msg) > 0) cat(msg, "\n")
  cat(paste(rep("=", 64), collapse="")); cat("\n")
}

# ================================================================
# 阶段 0: 读所有原材料 (R4 + Cavalli)
# ================================================================
hline("阶段 0: 读取 R4 + Cavalli")

# --- R4: 收集 tier1_strict 和 tier2_moderate 所有记录 ---
r4_all <- list()
for (samp in MB_SAMPLES) {
  f <- file.path(R4_DIR, samp, "filtered_density_genes.csv")
  if (!file.exists(f)) { cat("  ! 缺", samp, "\n"); next }
  dt <- fread(f)
  dt[, sample := samp]
  r4_all[[samp]] <- dt
}
r4 <- rbindlist(r4_all, use.names = TRUE, fill = TRUE)
cat(sprintf("  R4 读入: %d 样本, %d 行\n", length(r4_all), nrow(r4)))
cat("  tier 取值频次:\n")
print(table(r4$tier))

# --- Cavalli bulk (本地缓存, GEOquery 会认出) ---
cat("\n  读 GSE85217 (本地缓存)...\n")
gse <- getGEO(filename = file.path(GEO_DIR, "GSE85217_series_matrix.txt.gz"),
              getGPL = FALSE)
eset <- gse
expr <- exprs(eset)
pdat <- pData(eset)
cat(sprintf("    %d × %d\n", nrow(expr), ncol(expr)))

# --- Clinical (Table S1) ---
clin <- fread(CLINICAL)
cat(sprintf("  Clinical: %d × %d\n", nrow(clin), ncol(clin)))

# --- 对接 GSM ↔ Study_ID ---
# pdat 的 title 格式: MB_SubtypeStudy_XXXXX, clin 的 Study_ID 也是这个格式
pdat_dt <- data.table(
  gsm      = rownames(pdat),
  title    = pdat$title,
  subgroup_gse = pdat$`subgroup:ch1`,
  subtype_gse  = pdat$`subtype:ch1`
)
clin_merged <- merge(pdat_dt, clin, by.x = "title", by.y = "Study_ID", all.x = TRUE)
setnames(clin_merged, "OS (years)", "OS_years", skip_absent = TRUE)
setnames(clin_merged, "Dead", "event", skip_absent = TRUE)
cat(sprintf("  桥接后: %d 行 (有生存 %d, 有 SHH 亚型 %d)\n",
            nrow(clin_merged),
            sum(!is.na(clin_merged$OS_years) & !is.na(clin_merged$event)),
            sum(clin_merged$Subgroup == "SHH", na.rm = TRUE)))

# --- ENSG -> SYMBOL 映射 (去掉 _at 后缀) ---
probe_ids <- rownames(expr)
ensg_clean <- sub("_at$", "", probe_ids)
sym_map <- mapIds(org.Hs.eg.db, keys = ensg_clean, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
cat(sprintf("  ENSG -> SYMBOL: %d / %d 映射成功\n",
            sum(!is.na(sym_map)), length(sym_map)))

# 做一个 symbol -> probe 的反查表 (保留第一个出现)
sym2probe <- data.table(probe = probe_ids, symbol = sym_map)
sym2probe <- sym2probe[!is.na(symbol)]
sym2probe <- sym2probe[!duplicated(symbol)]
cat(sprintf("  unique symbols: %d\n", nrow(sym2probe)))

# ================================================================
# 阶段 1: 构造三种 signature (宽/中/严)
# ================================================================
hline("阶段 1: 构造三种 signature")

build_signature <- function(r4, strategy) {
  if (strategy == "strict") {
    # 严格: tier1_strict, 任何样本 (任何出现都算)
    cands <- r4[tier == "tier1_strict"]
  } else if (strategy == "medium") {
    # 中等: tier1_strict 或 tier2_moderate, 需 ≥2 样本 + 方向一致
    cands <- r4[tier %in% c("tier1_strict", "tier2_moderate")]
  } else if (strategy == "loose") {
    # 宽松: tier1 或 tier2, 单样本即可
    cands <- r4[tier %in% c("tier1_strict", "tier2_moderate")]
  }

  # 按 gene 聚合
  agg <- cands[, .(
    n_samples = .N,
    n_pos     = sum(direction == "positive"),
    n_neg     = sum(direction == "negative"),
    mean_rho  = mean(rho_knn_main, na.rm = TRUE)
  ), by = gene]
  agg[, consistency := pmax(n_pos, n_neg) / n_samples]
  agg[, direction_final := ifelse(n_pos > n_neg, "positive",
                            ifelse(n_neg > n_pos, "negative", "tie"))]

  if (strategy == "strict") {
    sig <- agg[n_samples >= 1]  # tier1 出现即可
  } else if (strategy == "medium") {
    sig <- agg[n_samples >= 2 & consistency >= 0.8]  # R10 同款
  } else if (strategy == "loose") {
    sig <- agg[n_samples >= 1]  # 单样本 tier2 也收
  }

  sig <- sig[direction_final != "tie"]
  sig[order(-abs(mean_rho))]
}

sig_strict <- build_signature(r4, "strict")
sig_medium <- build_signature(r4, "medium")
sig_loose  <- build_signature(r4, "loose")

cat(sprintf("  strict (tier1, ≥1 样本): %d 基因 (pos=%d, neg=%d)\n",
            nrow(sig_strict), sum(sig_strict$direction_final=="positive"),
            sum(sig_strict$direction_final=="negative")))
cat(sprintf("  medium (t1+t2, ≥2 样本 + 方向一致): %d 基因 (pos=%d, neg=%d)\n",
            nrow(sig_medium), sum(sig_medium$direction_final=="positive"),
            sum(sig_medium$direction_final=="negative")))
cat(sprintf("  loose  (t1+t2, ≥1 样本): %d 基因 (pos=%d, neg=%d)\n",
            nrow(sig_loose),  sum(sig_loose$direction_final=="positive"),
            sum(sig_loose$direction_final=="negative")))

# 映射到 probe
map_sig_to_probe <- function(sig) {
  sig2 <- merge(sig, sym2probe, by.x = "gene", by.y = "symbol")
  sig2
}
sig_strict_p <- map_sig_to_probe(sig_strict)
sig_medium_p <- map_sig_to_probe(sig_medium)
sig_loose_p  <- map_sig_to_probe(sig_loose)

cat(sprintf("\n  映射到 bulk probe: strict %d, medium %d, loose %d\n",
            nrow(sig_strict_p), nrow(sig_medium_p), nrow(sig_loose_p)))

# 保存 signature 基因列表
fwrite(sig_strict_p, file.path(OUT_DIR, "sig_strict_genes.csv"))
fwrite(sig_medium_p, file.path(OUT_DIR, "sig_medium_genes.csv"))
fwrite(sig_loose_p,  file.path(OUT_DIR, "sig_loose_genes.csv"))

# ================================================================
# 阶段 2: 为每个 signature 计算 score + 跑 5 测试
# ================================================================
hline("阶段 2: 三种 signature 各跑 5 测试")

compute_score <- function(sig_p, expr) {
  # 方向加权 z-score 平均
  probes <- sig_p$probe
  dirs   <- sig_p$direction_final
  M <- expr[probes, , drop = FALSE]
  # 每基因 z-score
  Z <- t(scale(t(M)))
  # 方向翻转: negative 基因的 z 取反, 使 "signature 高" = "高密度状态"
  sign_vec <- ifelse(dirs == "positive", 1, -1)
  Z_dir <- Z * sign_vec
  colMeans(Z_dir, na.rm = TRUE)
}

run_five_tests <- function(score, clin_merged, label) {
  # 对齐
  clin_merged$score <- score[match(clin_merged$gsm, names(score))]
  d <- clin_merged[!is.na(score) & !is.na(Subgroup)]

  out <- list(label = label, n_sig_genes = NA_integer_)

  # Test A1: Subgroup ANOVA (4)
  a1 <- aov(score ~ Subgroup, data = d)
  out$A1_p <- summary(a1)[[1]][["Pr(>F)"]][1]
  out$A1_n <- nrow(d)

  # Test A2: Subtype ANOVA (12)
  d_sub <- d[!is.na(Subtype)]
  a2 <- aov(score ~ Subtype, data = d_sub)
  out$A2_p <- summary(a2)[[1]][["Pr(>F)"]][1]
  out$A2_n <- nrow(d_sub)

  # Test B: log-rank (中位数分层), 整体
  d_surv <- d[!is.na(OS_years) & !is.na(event)]
  d_surv[, score_hi := score > median(score, na.rm = TRUE)]
  sdiff <- survdiff(Surv(OS_years, event) ~ score_hi, data = d_surv)
  out$B_p  <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
  out$B_n  <- nrow(d_surv)
  out$B_ev <- sum(d_surv$event == 1, na.rm = TRUE)

  # Cox 连续 + 年龄校正
  cox1 <- coxph(Surv(OS_years, event) ~ score, data = d_surv)
  out$Cox_HR     <- exp(coef(cox1))
  out$Cox_p      <- summary(cox1)$coefficients[1, "Pr(>|z|)"]

  if ("Age" %in% names(d_surv)) {
    d_cox <- d_surv[!is.na(Age)]
    cox2 <- coxph(Surv(OS_years, event) ~ score + Age, data = d_cox)
    out$Cox_age_HR <- exp(coef(cox2))[1]
    out$Cox_age_p  <- summary(cox2)$coefficients["score", "Pr(>|z|)"]
  } else {
    out$Cox_age_HR <- NA
    out$Cox_age_p  <- NA
  }

  out
}

results <- list()
for (strat in c("strict", "medium", "loose")) {
  sig_p <- get(paste0("sig_", strat, "_p"))
  cat(sprintf("\n--- [%s]  %d 基因 ---\n", strat, nrow(sig_p)))
  score <- compute_score(sig_p, expr)
  res <- run_five_tests(score, clin_merged, strat)
  res$n_sig_genes <- nrow(sig_p)
  results[[strat]] <- res
  # 保存每样本 score
  score_dt <- data.table(gsm = names(score), score = score)
  fwrite(score_dt, file.path(OUT_DIR, sprintf("scores_%s.csv", strat)))
}

# 对比表
cmp <- rbindlist(lapply(results, function(r) {
  data.table(
    signature    = r$label,
    n_genes      = r$n_sig_genes,
    A1_Subgroup_p = r$A1_p,
    A2_Subtype_p  = r$A2_p,
    B_logrank_p   = r$B_p,
    Cox_p         = r$Cox_p,
    Cox_HR        = r$Cox_HR,
    Cox_ageadj_p  = r$Cox_age_p,
    Cox_ageadj_HR = r$Cox_age_HR
  )
}))
cat("\n=== 三种 signature 对比 ===\n")
print(cmp)
fwrite(cmp, file.path(OUT_DIR, "signature_comparison.csv"))

# ================================================================
# 阶段 3: SHH 内部 K-M (2/3 分位) × 3 signatures = 6 张图
# ================================================================
hline("阶段 3: SHH 内 K-M (172 例)")

# R10 里 SHH n=172, events=38
d_all <- copy(clin_merged)

km_result_all <- list()

for (strat in c("strict", "medium", "loose")) {
  sig_p <- get(paste0("sig_", strat, "_p"))
  score <- compute_score(sig_p, expr)

  d <- copy(d_all)
  d$score <- score[match(d$gsm, names(score))]
  d_shh <- d[Subgroup == "SHH" & !is.na(score) & !is.na(OS_years) & !is.na(event)]
  cat(sprintf("\n  [%s] SHH n=%d, events=%d\n", strat, nrow(d_shh), sum(d_shh$event)))

  # --- 2 分位 ---
  d_shh[, group2 := ifelse(score > median(score), "High", "Low")]
  d_shh[, group2 := factor(group2, levels = c("Low", "High"))]
  sd2 <- survdiff(Surv(OS_years, event) ~ group2, data = d_shh)
  p2  <- 1 - pchisq(sd2$chisq, df = 1)
  cox2 <- coxph(Surv(OS_years, event) ~ score, data = d_shh)
  hr2  <- exp(coef(cox2)); p2_cox <- summary(cox2)$coefficients[1, "Pr(>|z|)"]

  # --- 3 分位 ---
  q <- quantile(d_shh$score, c(1/3, 2/3))
  d_shh[, group3 := cut(score, breaks = c(-Inf, q[1], q[2], Inf),
                         labels = c("Low", "Mid", "High"))]
  sd3 <- survdiff(Surv(OS_years, event) ~ group3, data = d_shh)
  p3  <- 1 - pchisq(sd3$chisq, df = 2)

  km_result_all[[strat]] <- list(
    n = nrow(d_shh), events = sum(d_shh$event),
    logrank2 = p2, Cox_HR = hr2, Cox_p = p2_cox,
    logrank3 = p3
  )

  # ---- 画 K-M (2 分位) ----
  fit2 <- survfit(Surv(OS_years, event) ~ group2, data = d_shh)
  # 手工转 ggplot 数据
  km_to_df <- function(fit) {
    s <- summary(fit, times = seq(0, max(fit$time, na.rm=TRUE), by = 0.1), extend = TRUE)
    data.frame(time = s$time, surv = s$surv, strata = s$strata)
  }
  df2 <- km_to_df(fit2)
  p_a <- ggplot(df2, aes(x = time, y = surv, color = strata)) +
    geom_step(size = 1) +
    scale_color_manual(values = c("group2=Low" = "#1f77b4", "group2=High" = "#d62728"),
                       labels = c("Low score", "High score")) +
    labs(title = sprintf("SHH KM, 2-partition  [%s signature, %d genes]",
                         strat, nrow(sig_p)),
         subtitle = sprintf("n=%d, events=%d | log-rank p=%.4f | Cox HR=%.3f (p=%.4f)",
                            nrow(d_shh), sum(d_shh$event), p2, hr2, p2_cox),
         x = "Overall Survival (years)", y = "Survival probability", color = "") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face="bold"),
          legend.position = "top") +
    coord_cartesian(ylim = c(0, 1))
  ggsave(file.path(OUT_DIR, "KM_plots", sprintf("SHH_KM_2partition_%s.png", strat)),
         p_a, width = 7, height = 5, dpi = 150)

  # ---- 画 K-M (3 分位) ----
  fit3 <- survfit(Surv(OS_years, event) ~ group3, data = d_shh)
  df3 <- km_to_df(fit3)
  p_b <- ggplot(df3, aes(x = time, y = surv, color = strata)) +
    geom_step(size = 1) +
    scale_color_manual(values = c("group3=Low" = "#1f77b4",
                                   "group3=Mid" = "#ff7f0e",
                                   "group3=High" = "#d62728"),
                       labels = c("Low", "Mid", "High")) +
    labs(title = sprintf("SHH KM, 3-partition  [%s signature, %d genes]",
                         strat, nrow(sig_p)),
         subtitle = sprintf("n=%d, events=%d | log-rank p=%.4f",
                            nrow(d_shh), sum(d_shh$event), p3),
         x = "Overall Survival (years)", y = "Survival probability", color = "Tertile") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face="bold"),
          legend.position = "top") +
    coord_cartesian(ylim = c(0, 1))
  ggsave(file.path(OUT_DIR, "KM_plots", sprintf("SHH_KM_3partition_%s.png", strat)),
         p_b, width = 7, height = 5, dpi = 150)

  cat(sprintf("    2-part: log-rank p=%.4f, Cox HR=%.3f (p=%.4f)\n", p2, hr2, p2_cox))
  cat(sprintf("    3-part: log-rank p=%.4f\n", p3))
}

# 汇总 K-M 结果
km_summary <- rbindlist(lapply(km_result_all, function(r) {
  data.table(n = r$n, events = r$events,
             logrank_2p = r$logrank2, Cox_HR = r$Cox_HR, Cox_p = r$Cox_p,
             logrank_3p = r$logrank3)
}), idcol = "signature")
cat("\n=== SHH K-M 对比 ===\n")
print(km_summary)
fwrite(km_summary, file.path(OUT_DIR, "SHH_KM_summary.csv"))

# ================================================================
# 阶段 4: SHH 四个亚亚型分层 Cox
# ================================================================
hline("阶段 4: SHH alpha/beta/gamma/delta 分层")

# Cavalli 的 SHH 亚型命名: SHH_alpha / SHH_beta / SHH_gamma / SHH_delta
# Subtype 列格式: "SHH alpha" 或 "SHH_alpha" - 检查一下
cat("  SHH 内 Subtype 取值:\n")
shh_sub_vals <- unique(clin_merged[Subgroup == "SHH", Subtype])
print(shh_sub_vals)

# 用 medium signature 做这个分析 (中等纳入标准最平衡)
score_m <- compute_score(sig_medium_p, expr)
d <- copy(clin_merged)
d$score <- score_m[match(d$gsm, names(score_m))]
d_shh_all <- d[Subgroup == "SHH" & !is.na(score) & !is.na(OS_years) & !is.na(event) & !is.na(Subtype)]

shh_subtype_results <- list()
for (st in sort(unique(d_shh_all$Subtype))) {
  d_st <- d_shh_all[Subtype == st]
  if (nrow(d_st) < 20 || sum(d_st$event) < 5) {
    cat(sprintf("  %-16s n=%d events=%d (skip, 样本/事件过少)\n",
                st, nrow(d_st), sum(d_st$event)))
    shh_subtype_results[[st]] <- data.table(subtype=st, n=nrow(d_st), events=sum(d_st$event),
                                             HR=NA, p=NA, note="skipped")
    next
  }
  cox <- tryCatch(coxph(Surv(OS_years, event) ~ score, data = d_st),
                  error = function(e) NULL)
  if (is.null(cox)) next
  hr <- exp(coef(cox)); p <- summary(cox)$coefficients[1, "Pr(>|z|)"]
  cat(sprintf("  %-16s n=%d events=%d  HR=%.3f  p=%.4f\n",
              st, nrow(d_st), sum(d_st$event), hr, p))
  shh_subtype_results[[st]] <- data.table(subtype=st, n=nrow(d_st), events=sum(d_st$event),
                                           HR=hr, p=p, note="ok")
}

shh_sub_tbl <- rbindlist(shh_subtype_results, fill = TRUE)
fwrite(shh_sub_tbl, file.path(OUT_DIR, "SHH_subtype_Cox.csv"))

# ================================================================
# 阶段 5: 汇总报告
# ================================================================
hline("R11 Summary")
sink(file.path(OUT_DIR, "R11_SUMMARY.txt"), split = TRUE)

cat("================================================================\n")
cat(sprintf("  R11: Signature Expansion + SHH Deep-dive\n"))
cat(sprintf("  %s\n", Sys.time()))
cat("================================================================\n\n")

cat("PART 1: Signature 扩展对比 (R10 结果对比)\n")
cat("----------------------------------------------------------------\n")
print(cmp)

# R10 原始结果放在这里做对比
cat("\n【R10 原始 signature = medium 策略】\n")
cat("  A1 Subgroup: p=3.06e-19\n")
cat("  A2 Subtype:  p=3.03e-25\n")
cat("  B  log-rank: p=0.022\n")
cat("  Cox age-adj: p=0.067  <- R10 的未过门槛项\n\n")

cat("Cox age-adj 最小 p (三策略中):\n")
min_p <- cmp[which.min(Cox_ageadj_p)]
cat(sprintf("  %s (%d 基因): p=%.4f, HR=%.3f\n",
            min_p$signature, min_p$n_genes, min_p$Cox_ageadj_p, min_p$Cox_ageadj_HR))
if (!is.na(min_p$Cox_ageadj_p) && min_p$Cox_ageadj_p < 0.05) {
  cat("  >>> Cox age-adj 推过 0.05, Test B' 升级为 PASS <<<\n")
}

cat("\nPART 2: SHH 内部 K-M (n=172)\n")
cat("----------------------------------------------------------------\n")
print(km_summary)

cat("\nPART 3: SHH 亚亚型分层 Cox (用 medium signature)\n")
cat("----------------------------------------------------------------\n")
print(shh_sub_tbl)

driving_subtype <- shh_sub_tbl[!is.na(p) & p < 0.05]
if (nrow(driving_subtype) > 0) {
  cat("\n  驱动 SHH 信号的亚亚型:\n")
  print(driving_subtype)
} else {
  cat("\n  亚亚型层面: 无单个亚亚型达 p<0.05 (样本量限制)\n")
}

cat("\nPART 4: 对论文写作的建议\n")
cat("----------------------------------------------------------------\n")
# 找"最好"的 signature - 综合 A1 + B + Cox age-adj
best_strat <- cmp[, .(score = -log10(A1_Subgroup_p) + -log10(B_logrank_p) +
                                 -log10(pmax(Cox_ageadj_p, 1e-10))),
                  by = .(signature, n_genes)]
best_strat <- best_strat[order(-score)]
cat("综合评分 (A1 + B + Cox age-adj 的 -log10(p) 总和):\n")
print(best_strat)
cat(sprintf("\n  论文推荐主 signature: %s (%d 基因)\n",
            best_strat$signature[1], best_strat$n_genes[1]))

cat("\n================================================================\n")
cat("输出文件 (",OUT_DIR,"):\n")
cat("  - signature_comparison.csv       # 三策略 × 5 测试对比\n")
cat("  - sig_strict_genes.csv / sig_medium_genes.csv / sig_loose_genes.csv\n")
cat("  - scores_strict.csv / scores_medium.csv / scores_loose.csv\n")
cat("  - SHH_KM_summary.csv             # SHH 172 例 K-M 汇总\n")
cat("  - KM_plots/SHH_KM_{2,3}partition_{strict,medium,loose}.png  # 6 张图\n")
cat("  - SHH_subtype_Cox.csv            # 4 亚亚型分层\n")
cat("================================================================\n")
sink()

cat("\nR11 全部完成\n")
