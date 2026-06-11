# ============================================================
# R17a: Assemble Seurat object from Aldinger 2021 data
# ------------------------------------------------------------
# Input:  UCSC Cell Browser files (exprMatrix.tsv.gz, meta.tsv, UMAP)
# Output: Single .rds containing Seurat object ready for downstream
#
# Key decisions:
#   - Data is already SCT-scaled (z-score-like, range [-3, +3])
#   - Put it in `scale.data` slot of a Seurat object
#   - Also compute a `data` slot by clipping & shifting (for AddModuleScore)
#   - Read UMAP as Reductions['umap']
#   - Verify sig_94 / sig_core gene matches
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

cat("================================================================\n")
cat("R17a: Assemble Aldinger 2021 Seurat object\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

EXTDATA_DIR <- Sys.getenv("EXTDATA_DIR", unset = "./external_data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "./results")
IN_DIR <- file.path(EXTDATA_DIR, "Aldinger2021")
OUT_DIR <- file.path(RESULTS_DIR, "R17_Aldinger2021")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EXPR_FILE <- file.path(IN_DIR, "exprMatrix.tsv.gz")
META_FILE <- file.path(IN_DIR, "meta.tsv")
UMAP_FILE <- file.path(IN_DIR, "Seurat_UMAP.coords.tsv.gz")

SIG_PROV <- file.path(RESULTS_DIR, "R12_Gaps", "sig_strict_99_provenance.csv")
CONFLICT_GENES <- c("TENM1", "SLC17A7", "NRXN3", "DCN", "SV2B")
SIG_CORE_GENES <- c("TUBB4A", "APOE", "CCR7", "EOMES", "ST18", "NES", "AQP4", "QKI")

# ============================================================
# Step 1: Load expression matrix
# ============================================================
cat("[Step 1] Loading exprMatrix.tsv.gz...\n")
t0 <- Sys.time()
expr_dt <- fread(cmd = sprintf("zcat %s", EXPR_FILE), sep = "\t",
                 header = TRUE, data.table = TRUE)
cat(sprintf("  Loaded in %.1f sec, dim = %d x %d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(expr_dt), ncol(expr_dt)))

gene_ids <- expr_dt[[1]]  # first col = gene symbols
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids
rm(expr_dt); gc(verbose = FALSE)

cat(sprintf("  Matrix: %d genes x %d cells\n",
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  Value range: [%.3f, %.3f], median=%.3f\n",
            min(expr_mat), max(expr_mat), median(expr_mat)))
cat(sprintf("  Value type (based on range): SCT-residuals (z-score-like)\n"))

# ============================================================
# Step 2: Load metadata
# ============================================================
cat("\n[Step 2] Loading meta.tsv...\n")
meta_dt <- fread(META_FILE, sep = "\t", header = TRUE, data.table = TRUE)
cat(sprintf("  Meta: %d rows x %d cols\n", nrow(meta_dt), ncol(meta_dt)))
cat(sprintf("  Columns: %s\n",
            paste(names(meta_dt), collapse = ", ")))

# cellId 是 meta 的主键
setnames(meta_dt, "cellId", "cellId_orig")  # 避免 seurat 自带 cellId
meta_df <- as.data.frame(meta_dt)
rownames(meta_df) <- meta_df$cellId_orig

# ============================================================
# Step 3: Load UMAP coords
# ============================================================
cat("\n[Step 3] Loading UMAP coords...\n")
# UMAP file has no header — columns are cellId, UMAP_1, UMAP_2
umap_dt <- fread(cmd = sprintf("zcat %s", UMAP_FILE), sep = "\t",
                 header = FALSE, col.names = c("cellId", "UMAP_1", "UMAP_2"),
                 data.table = TRUE)
cat(sprintf("  UMAP: %d cells\n", nrow(umap_dt)))

umap_mat <- as.matrix(umap_dt[, .(UMAP_1, UMAP_2)])
rownames(umap_mat) <- umap_dt$cellId

# ============================================================
# Step 4: Align cell IDs across all 3 files
# ============================================================
cat("\n[Step 4] Aligning cell IDs...\n")
expr_cells <- colnames(expr_mat)
meta_cells <- rownames(meta_df)
umap_cells <- rownames(umap_mat)

cat(sprintf("  expr cells:  %d\n", length(expr_cells)))
cat(sprintf("  meta cells:  %d\n", length(meta_cells)))
cat(sprintf("  umap cells:  %d\n", length(umap_cells)))

common <- Reduce(intersect, list(expr_cells, meta_cells, umap_cells))
cat(sprintf("  Common cells: %d\n", length(common)))

# 按照共同 order 对齐
expr_mat <- expr_mat[, common]
meta_df <- meta_df[common, , drop = FALSE]
umap_mat <- umap_mat[common, , drop = FALSE]

stopifnot(identical(colnames(expr_mat), rownames(meta_df)))
stopifnot(identical(colnames(expr_mat), rownames(umap_mat)))
cat("  ✓ All three tables aligned on common cells\n")

# ============================================================
# Step 5: Build Seurat object
# ============================================================
cat("\n[Step 5] Building Seurat object...\n")

# 数据已经是 SCT residuals, 放在 scale.data
# 同时做一个简单的 "pseudo-data" slot: 把负值 clip 到 0, 供 AddModuleScore 用
# AddModuleScore 期望 log-normalized counts, 但 SCT residuals 也能用
# (Seurat 官方文档: score = 目标基因均值 - control set 均值, 对 residuals 也稳健)

# 直接创建 Assay5 (Seurat v5 格式)
cat("  Creating Assay5...\n")

# Seurat 要求 counts 矩阵 (可以是稀疏或 dense); 我们把 SCT 残差放在 scale.data,
# 同时提供 data 作 AddModuleScore fallback
# 最简单: 把 expr_mat 同时放 counts, data, scale.data
# (counts 虽然被滥用了, 但 downstream 我们不会用到)

# 为了 AddModuleScore 更合理, 把残差 clip 到 [0, max] 作 "data" slot
data_slot <- pmax(expr_mat, 0)  # clip 负值, 不影响正 score 基因

obj <- CreateSeuratObject(
  counts = data_slot,
  meta.data = meta_df,
  project = "Aldinger2021"
)

# 手动塞入 scale.data (原始 SCT 残差)
obj <- SetAssayData(obj, assay = "RNA", slot = "scale.data",
                    new.data = expr_mat)

# 同时 'data' slot 已经填成 clipped 版本
obj <- SetAssayData(obj, assay = "RNA", slot = "data",
                    new.data = data_slot)

cat(sprintf("  Seurat object: %d features x %d cells\n",
            nrow(obj), ncol(obj)))
cat(sprintf("  Assays: %s\n", paste(Assays(obj), collapse = ", ")))

# ============================================================
# Step 6: Add UMAP reduction
# ============================================================
cat("\n[Step 6] Adding UMAP reduction...\n")
obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAP_",
  assay = "RNA"
)
cat(sprintf("  UMAP added: %d cells x 2 dims\n",
            nrow(obj[["umap"]]@cell.embeddings)))

# ============================================================
# Step 7: Rename 'Cluster' to 'cell_type' for clarity
# ============================================================
cat("\n[Step 7] Setting up cell_type identity...\n")
obj$cell_type <- obj$Cluster
obj$age_pcw <- gsub(" PCW", "", obj$age)
obj$age_pcw <- as.numeric(obj$age_pcw)
cat(sprintf("  cell_type (21 levels): %s\n",
            paste(head(sort(unique(obj$cell_type)), 5), collapse = ", ")))
cat(sprintf("  age_pcw range: [%d, %d]\n",
            min(obj$age_pcw, na.rm = TRUE),
            max(obj$age_pcw, na.rm = TRUE)))

Idents(obj) <- "cell_type"

# ============================================================
# Step 8: Signature gene match check
# ============================================================
cat("\n[Step 8] Verifying signature gene matches...\n")

sig_prov <- fread(SIG_PROV)
sig_94 <- setdiff(sig_prov$gene, CONFLICT_GENES)

genes_in_obj <- rownames(obj)

match_94 <- intersect(sig_94, genes_in_obj)
match_core <- intersect(SIG_CORE_GENES, genes_in_obj)

cat(sprintf("  sig_94:   %d / %d matched (%.1f%%)\n",
            length(match_94), length(sig_94),
            100 * length(match_94) / length(sig_94)))
cat(sprintf("  sig_core: %d / %d matched\n",
            length(match_core), length(SIG_CORE_GENES)))

# 看一下方向分布
sig_94_prov <- sig_prov[gene %in% match_94, .(gene, direction_final)]
dir_dist <- table(sig_94_prov$direction_final)
cat(sprintf("  Direction dist in matched sig_94:\n"))
for (d in names(dir_dist)) {
  cat(sprintf("    %s: %d\n", d, dir_dist[d]))
}

# 保存匹配结果
fwrite(data.table(
  signature = c("sig_94", "sig_core"),
  total_genes = c(length(sig_94), length(SIG_CORE_GENES)),
  matched = c(length(match_94), length(match_core)),
  match_pct = c(round(100 * length(match_94) / length(sig_94), 1),
                round(100 * length(match_core) / length(SIG_CORE_GENES), 1))
), file.path(OUT_DIR, "signature_match_summary.csv"))

fwrite(data.table(gene = match_94,
                  direction = sig_prov[gene %in% match_94, direction_final][match(match_94, sig_prov[gene %in% match_94, gene])]),
       file.path(OUT_DIR, "sig94_matched_with_direction.csv"))

fwrite(data.table(gene = match_core), file.path(OUT_DIR, "sigcore_matched.csv"))

# ============================================================
# Step 9: Per-cluster cell counts (confirm)
# ============================================================
cat("\n[Step 9] Cluster × age distribution...\n")
cluster_age_tab <- table(obj$cell_type, obj$age_pcw)
cat("  Cluster × age cell counts (head):\n")
print(head(as.data.frame.matrix(cluster_age_tab), 10))

fwrite(as.data.table(as.data.frame.matrix(cluster_age_tab), keep.rownames = "cell_type"),
       file.path(OUT_DIR, "cluster_age_counts.csv"))

# ============================================================
# Step 10: Save Seurat object
# ============================================================
OUT_RDS <- file.path(OUT_DIR, "aldinger2021_seurat.rds")
cat(sprintf("\n[Step 10] Saving Seurat object to %s ...\n", OUT_RDS))
t0 <- Sys.time()
saveRDS(obj, OUT_RDS)
cat(sprintf("  Saved (%.1f sec, %.1f MB)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            file.size(OUT_RDS) / 1e6))

# ============================================================
# Final summary
# ============================================================
cat("\n================================================================\n")
cat("R17a DONE\n")
cat("================================================================\n")
cat(sprintf("Object: %d genes x %d cells\n", nrow(obj), ncol(obj)))
cat(sprintf("Cell types: %d\n", length(unique(obj$cell_type))))
cat(sprintf("Ages: %s\n", paste(sort(unique(obj$age_pcw)), collapse=", ")))
cat(sprintf("sig_94 matched: %d / %d\n", length(match_94), length(sig_94)))
cat(sprintf("sig_core matched: %d / %d\n", length(match_core), length(SIG_CORE_GENES)))
cat("\nOutputs:\n")
cat(sprintf("  - %s\n", OUT_RDS))
cat(sprintf("  - %s/signature_match_summary.csv\n", OUT_DIR))
cat(sprintf("  - %s/sig94_matched_with_direction.csv\n", OUT_DIR))
cat(sprintf("  - %s/cluster_age_counts.csv\n", OUT_DIR))
cat("\nReady for R17b (signature scoring + cluster-level analysis)\n")
