# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# install_R10_deps.R
# 功能: 安装 R10 需要的包, 带镜像 + 失败诊断 + 自动 fallback
# 运行: Rscript install_R10_deps.R 2>&1 | tee install_R10.log
# ============================================================================

cat("============================================\n")
cat("R10 依赖包安装\n")
cat("============================================\n")
cat("R version:", R.version.string, "\n")
cat("Lib path :", .libPaths()[1], "\n")
cat("Writable :", file.access(.libPaths()[1], 2) == 0, "\n\n")

# 显式设镜像 (以防 .Rprofile 没生效)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/")
options(install.packages.check.source = "no")
options(timeout = 600)  # 下载超时拉长到 10 分钟

# 需要检查的包
needed <- list(
  # CRAN
  survival       = list(source = "CRAN",        critical = TRUE),
  survminer      = list(source = "CRAN",        critical = TRUE),
  data.table     = list(source = "CRAN",        critical = TRUE),
  # Bioconductor
  GEOquery       = list(source = "Bioconductor", critical = TRUE),
  Biobase        = list(source = "Bioconductor", critical = TRUE)
)

# 检查已装
cat("检查已装状态:\n")
for (p in names(needed)) {
  have <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-15s [%s]  source=%s\n",
              p, if (have) "ok" else "MISSING", needed[[p]]$source))
}
cat("\n")

# 找出要装的
to_install <- names(needed)[!sapply(names(needed), requireNamespace, quietly = TRUE)]
if (length(to_install) == 0) {
  cat("所有包已安装, 无需操作\n")
  quit(save = "no", status = 0)
}

cat("需要安装:", paste(to_install, collapse = ", "), "\n\n")

# BiocManager 检查
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("BiocManager 不存在, 先装...\n")
  install.packages("BiocManager")
}
cat("BiocManager version:", as.character(BiocManager::version()), "\n\n")

# 逐个装, 失败了继续 (不要 all-or-nothing)
install_results <- list()

for (p in to_install) {
  meta <- needed[[p]]
  cat("========================================\n")
  cat(sprintf("  Installing: %s (from %s)\n", p, meta$source))
  cat("========================================\n")
  
  t0 <- Sys.time()
  result <- tryCatch({
    if (meta$source == "CRAN") {
      install.packages(p, dependencies = TRUE)
    } else if (meta$source == "Bioconductor") {
      BiocManager::install(p, update = FALSE, ask = FALSE, force = FALSE)
    }
    requireNamespace(p, quietly = TRUE)
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    FALSE
  }, warning = function(w) {
    cat("\nWARNING:", conditionMessage(w), "\n")
    requireNamespace(p, quietly = TRUE)  # warning 不一定是失败
  })
  
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  
  if (isTRUE(result)) {
    cat(sprintf("\n  ✓ %s OK (%.1fs)\n\n", p, elapsed))
    install_results[[p]] <- "OK"
  } else {
    cat(sprintf("\n  ✗ %s FAILED (%.1fs)\n\n", p, elapsed))
    install_results[[p]] <- "FAILED"
  }
}

# 总结
cat("============================================\n")
cat("安装结果总结:\n")
cat("============================================\n")
for (p in names(install_results)) {
  cat(sprintf("  %-15s  %s\n", p, install_results[[p]]))
}

n_fail <- sum(unlist(install_results) == "FAILED")
if (n_fail == 0) {
  cat("\n全部安装成功!\n")
  # 最终验证: library() 一遍
  cat("\n最终验证 (library load test):\n")
  for (p in names(needed)) {
    res <- tryCatch({
      suppressMessages(library(p, character.only = TRUE))
      "OK"
    }, error = function(e) paste("FAILED:", conditionMessage(e)))
    cat(sprintf("  %-15s  %s\n", p, res))
  }
} else {
  cat(sprintf("\n%d 个包安装失败, 把 install_R10.log 贴给 Claude 查看错误\n", n_fail))
}
