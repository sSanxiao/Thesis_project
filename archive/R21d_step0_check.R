# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R21d Step 0: Pre-flight check
# Verify SHH pathway genes are in Aldinger 2000 var genes
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

cat("================================================================\n")
cat("R21d Step 0: SHH pathway gene coverage in Aldinger\n")
cat("================================================================\n\n")

cat("Loading Aldinger Seurat object...\n")
aldinger <- readRDS("./results/R17_Aldinger2021/aldinger2021_scored.rds")

cat(sprintf("Aldinger cells: %d, genes: %d\n",
            ncol(aldinger), nrow(aldinger)))

# Check SHH pathway markers (Ghasemi 2024 quiescent CGNP markers)
shh_genes <- c("PTCH1", "SMO", "HHIP", "GLI1", "GLI2", "GLI3", "PTPRK",
               "BOC", "CDON", "GAS1", "SUFU")

# R21c top up genes (key sig_94-high markers)
r21c_top <- c("HHIP", "DCC", "BOC", "TMEM108", "KCNMB1", "LHX8", "EBF3",
              "TRPM3", "PLCH1", "CDCA7", "INSRR", "C3orf22", "BCL11A")

# sig_94 positive direction (need to load)
sig_prov <- read.csv("./results/R12_Gaps/sig_strict_99_provenance.csv")
sig_94_pos <- sig_prov$gene[sig_prov$direction_final == "positive" &
                              !sig_prov$gene %in% c("TENM1","SLC17A7","NRXN3","DCN","SV2B")]

cat("\n=== SHH pathway gene coverage ===\n")
present_shh <- intersect(shh_genes, rownames(aldinger))
absent_shh <- setdiff(shh_genes, rownames(aldinger))
cat(sprintf("Present (%d/%d): %s\n", length(present_shh), length(shh_genes),
            paste(present_shh, collapse = ", ")))
cat(sprintf("Absent: %s\n", paste(absent_shh, collapse = ", ")))

cat("\n=== R21c sig_94-high top markers coverage ===\n")
present_r21c <- intersect(r21c_top, rownames(aldinger))
absent_r21c <- setdiff(r21c_top, rownames(aldinger))
cat(sprintf("Present (%d/%d): %s\n", length(present_r21c), length(r21c_top),
            paste(present_r21c, collapse = ", ")))
cat(sprintf("Absent: %s\n", paste(absent_r21c, collapse = ", ")))

cat("\n=== sig_94 positive direction coverage ===\n")
present_pos <- intersect(sig_94_pos, rownames(aldinger))
cat(sprintf("Present (%d/%d): %.1f%%\n", length(present_pos), length(sig_94_pos),
            100 * length(present_pos) / length(sig_94_pos)))

cat("\n=== VERDICT ===\n")
if (length(present_shh) >= 3 && "PTCH1" %in% present_shh) {
  cat("✓ Sufficient SHH pathway coverage. Proceed with R21d full version.\n")
} else {
  cat("⚠ Limited SHH coverage. R21d will be a degraded version using available genes.\n")
}
