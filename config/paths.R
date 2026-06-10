# ============================================================
# config/paths.R
# ------------------------------------------------------------
# Central path configuration for the R pipeline.
#
# All scripts resolve input/output locations from environment
# variables instead of hardcoded absolute paths. Set these in your
# shell before running, or rely on the local-relative defaults below.
#
#   DATA_DIR     Xenium spatial input datasets   (was: Datasets_April)
#   EXTDATA_DIR  External validation datasets     (was: External_Data)
#   RESULTS_DIR  Pipeline outputs / results        (was: Results_New)
#
# Usage inside a script:
#   source("config/paths.R")
#   readRDS(file.path(RESULTS_DIR, "R1_seurat.rds"))
#
# Bash:        export DATA_DIR=/path/to/Xenium_datasets
# PowerShell:  $env:DATA_DIR = "D:\path\to\Xenium_datasets"
# ============================================================

DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = "./data")
EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
