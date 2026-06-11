#!/usr/bin/env Rscript
# ============================================================================
# setup/install_deps.R
#
# Install all R packages required by the pipeline (CRAN + Bioconductor),
# one at a time with per-package error handling so a single failure does
# not abort the whole run. Consolidates the earlier install_R10_deps.R /
# install_R10_deps_v2.R helpers into one setup script covering R1-R21.
#
# Repositories default to the official CRAN / Bioconductor mirrors. To use
# a faster local mirror, set the CRAN_MIRROR / BIOC_MIRROR env vars, e.g.
#   export CRAN_MIRROR=https://mirrors.tuna.tsinghua.edu.cn/CRAN/
#   export BIOC_MIRROR=https://mirrors.ustc.edu.cn/bioc/
#
# Run: Rscript setup/install_deps.R 2>&1 | tee install_deps.log
#
# Note: Seurat 5.x is compatibility-sensitive; a frozen sessionInfo() is
#       recommended for exact reproduction (see README).
# ============================================================================

options(timeout = 600)  # lengthen download timeout to 10 min

cran_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_mirror))

bioc_mirror <- Sys.getenv("BIOC_MIRROR", unset = "")
if (nzchar(bioc_mirror)) options(BioC_mirror = bioc_mirror)

cat("============================================\n")
cat("Pipeline R dependency installation\n")
cat("============================================\n")
cat("R version:", R.version.string, "\n")
cat("Lib path :", .libPaths()[1], "\n")
cat("CRAN     :", cran_mirror, "\n\n")

# --- Package inventory (collected from R1-R21) -----------------------------
cran_pkgs <- c(
  "Matrix", "Seurat", "SeuratObject", "data.table", "jsonlite",
  "ggplot2", "ggrepel", "ggridges", "gridExtra", "patchwork",
  "pheatmap", "reshape2", "scales", "viridis", "matrixStats",
  "survival", "coxphf", "harmony"
)
bioc_pkgs <- c(
  "Biobase", "AnnotationDbi", "org.Hs.eg.db", "GEOquery",
  "SingleCellExperiment", "slingshot"
)

# --- Ensure BiocManager is present -----------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}
cat("BiocManager version:", as.character(BiocManager::version()), "\n\n")

install_one <- function(pkg, source) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-22s [already installed]\n", pkg))
    return("OK")
  }
  cat(sprintf("  Installing %-15s (%s)...\n", pkg, source))
  t0 <- Sys.time()
  ok <- tryCatch({
    if (source == "CRAN") {
      install.packages(pkg, dependencies = TRUE)
    } else {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n"); FALSE
  })
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  status <- if (isTRUE(ok)) "OK" else "FAILED"
  cat(sprintf("    %s (%.1fs)\n", status, elapsed))
  status
}

results <- list()
cat("--- CRAN packages ---\n")
for (p in cran_pkgs) results[[p]] <- install_one(p, "CRAN")
cat("\n--- Bioconductor packages ---\n")
for (p in bioc_pkgs) results[[p]] <- install_one(p, "Bioconductor")

# --- Summary ---------------------------------------------------------------
cat("\n============================================\n")
cat("Installation summary\n")
cat("============================================\n")
for (p in names(results)) cat(sprintf("  %-22s %s\n", p, results[[p]]))

n_fail <- sum(unlist(results) == "FAILED")
if (n_fail == 0) {
  cat("\nAll packages installed successfully.\n")
} else {
  cat(sprintf("\n%d package(s) failed; inspect install_deps.log for details.\n", n_fail))
}
