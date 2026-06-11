# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# install_R10_deps_v2.R
# 修复 v1 的两个问题:
#   1. survminer -> fastGHQuad 要 C++14, gcc 4.8.5 不支持 -> 放弃 survminer
#      (不影响生存分析, Cox 模型和 KM 曲线用 survival 包 + base plot 就够了)
#   2. BiocManager 镜像 -> 不覆盖 repos, 而是同时保留 CRAN 和 BioC 源
#
# 运行: Rscript install_R10_deps_v2.R 2>&1 | tee install_R10_v2.log
# ============================================================================

cat("============================================\n")
cat("R10 依赖包安装 v2\n")
cat("============================================\n")
cat("R version:", R.version.string, "\n")
cat("Lib path :", .libPaths()[1], "\n")
cat("Writable :", file.access(.libPaths()[1], 2) == 0, "\n\n")

options(timeout = 600)

# -------- 关键修复 1: 正确设置 CRAN + BioC 双仓库 --------
# BiocManager::repositories() 会返回 BioC 的所有仓库地址,
# 我们自己加上 CRAN 镜像, 让 install.packages / BiocManager 都能找到包

cran_mirror <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
# 也可以换: "https://mirrors.ustc.edu.cn/CRAN/" 或 "https://mirror.lzu.edu.cn/CRAN/"

# 中科大 BioC 镜像 (BiocManager 3.16 对应 BioC 3.16 版本)
biocm_base <- "https://mirrors.ustc.edu.cn/bioc/3.16"
bioc_repos <- c(
  BioCsoft       = paste0(biocm_base, "/bioc"),
  BioCann        = paste0(biocm_base, "/data/annotation"),
  BioCexp        = paste0(biocm_base, "/data/experiment"),
  BioCworkflows  = paste0(biocm_base, "/workflows")
)

# 合并: CRAN + BioC
all_repos <- c(CRAN = cran_mirror, bioc_repos)
options(repos = all_repos)

cat("已配置 repos:\n")
for (nm in names(all_repos)) {
  cat(sprintf("  %-15s %s\n", nm, all_repos[[nm]]))
}
cat("\n")

# 验证 BioC 镜像能连通 (试拉 PACKAGES 文件)
cat("验证 BioC 镜像连通性...\n")
test_url <- paste0(bioc_repos["BioCsoft"], "/src/contrib/PACKAGES.gz")
t <- tryCatch({
  con <- url(test_url)
  on.exit(close(con))
  readLines(con, n = 1)
  TRUE
}, error = function(e) FALSE, warning = function(w) FALSE)
cat(sprintf("  BioCsoft 连通: %s\n\n", t))

if (!t) {
  cat("⚠ 中科大 BioC 镜像连不上, 切回官方源\n")
  bioc_repos_official <- BiocManager::repositories()
  options(repos = c(CRAN = cran_mirror, bioc_repos_official))
  cat("已切到官方 BioC 源\n\n")
}

# -------- 关键修复 2: 放弃 survminer, 只装核心包 --------

needed <- list(
  survival       = list(source = "CRAN"),          # 你已装
  data.table     = list(source = "CRAN"),          # 你已装
  Biobase        = list(source = "Bioconductor"),  # 你已装
  GEOquery       = list(source = "Bioconductor"),  # 要装 (拉 GSE85217)
  ggplot2        = list(source = "CRAN"),          # 画 KM 曲线 (应该已装)
  ggpubr         = list(source = "CRAN")           # 组图 (可选, 无则跳过)
)

cat("检查已装状态:\n")
for (p in names(needed)) {
  have <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-15s [%s]\n", p, if (have) "ok" else "MISSING"))
}
cat("\n")

to_install <- names(needed)[!sapply(names(needed), requireNamespace, quietly = TRUE)]
if (length(to_install) == 0) {
  cat("所有包已安装\n")
  quit(save = "no", status = 0)
}
cat("需要安装:", paste(to_install, collapse = ", "), "\n\n")

# 逐个装
install_results <- list()
for (p in to_install) {
  cat("========================================\n")
  cat(sprintf("  Installing: %s (from %s)\n", p, needed[[p]]$source))
  cat("========================================\n")
  
  t0 <- Sys.time()
  result <- tryCatch({
    if (needed[[p]]$source == "CRAN") {
      install.packages(p, dependencies = TRUE)
    } else {
      # Bioconductor: 不再用 BiocManager::install,
      # 因为 options(repos) 里已经有 BioC 仓库, 直接 install.packages 就能找到
      install.packages(p, dependencies = TRUE)
    }
    requireNamespace(p, quietly = TRUE)
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    FALSE
  })
  
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  install_results[[p]] <- if (isTRUE(result)) "OK" else "FAILED"
  cat(sprintf("\n  %s %s (%.1fs)\n\n",
              if (isTRUE(result)) "✓" else "✗", p, elapsed))
}

# 总结
cat("============================================\n")
cat("安装结果:\n")
cat("============================================\n")
for (p in names(install_results)) {
  cat(sprintf("  %-15s  %s\n", p, install_results[[p]]))
}

# 最终 library 测试
cat("\n最终 library 测试:\n")
for (p in c("survival", "data.table", "GEOquery", "Biobase", "ggplot2")) {
  res <- tryCatch({
    suppressMessages(library(p, character.only = TRUE))
    "OK"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("  %-15s  %s\n", p, res))
}

# 功能性验证 (最重要): GEOquery 能不能 ping 到 GSE
if (requireNamespace("GEOquery", quietly = TRUE)) {
  cat("\nGEOquery 功能性验证 (不下载, 只查元数据):\n")
  t_geo <- tryCatch({
    suppressMessages({
      gse <- GEOquery::getGEO("GSE85217", GSEMatrix = FALSE, getGPL = FALSE,
                              destdir = tempdir())
    })
    cat("  ✓ GSE85217 metadata fetched OK\n")
    TRUE
  }, error = function(e) {
    cat("  ✗ FAILED:", conditionMessage(e), "\n")
    FALSE
  })
}

cat("\n============================================\n")
cat("完成\n")
cat("============================================\n")
