# [ARCHIVE] Debugging / version-iteration script — retained for provenance only.
# Not part of the clean pipeline. Server paths replaced by relative placeholders
# (./data, ./results, ./external_data); see config/paths.R for the main-script convention.
# ============================================================================
# ============================================================
# R21d UMAP fix - just redo the UMAP plot with correct column names
# Other 3 plots already saved successfully
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUT_DIR <- "./results/R21d_Aldinger_SHH"

cat("Loading Aldinger Seurat object...\n")
aldinger <- readRDS("./results/R17_Aldinger2021/aldinger2021_scored.rds")

aldinger$sig_94 <- aldinger$sig_94_zscore

# Get UMAP - try both possible column name conventions
umap_red <- Embeddings(aldinger, reduction = "umap")
cat(sprintf("UMAP column names: %s\n", paste(colnames(umap_red), collapse = ", ")))

umap_coords <- as.data.table(umap_red)
# Standardize column names
setnames(umap_coords, c("d1", "d2"))
umap_coords[, cell_id := colnames(aldinger)]
umap_coords[, sig_94 := as.numeric(aldinger$sig_94)]

# SHH gene expression
expr <- GetAssayData(aldinger, assay = "RNA", layer = "data")
shh_present <- intersect(c("PTCH1", "HHIP", "GLI1", "GLI2", "GLI3", "BOC", "PTPRK"),
                          rownames(aldinger))
for (g in shh_present) {
  umap_coords[[g]] <- as.numeric(expr[g, ])
}

# Read correlation table
cor_table <- fread(file.path(OUT_DIR, "shh_gene_correlation_table.csv"))

# UMAP for sig_94
umap_lim <- max(abs(quantile(umap_coords$sig_94, c(0.01, 0.99), na.rm=TRUE)))

p_umap_sig94 <- ggplot(umap_coords, aes(x = d1, y = d2, color = sig_94)) +
  geom_point(size = 0.15, alpha = 0.5) +
  scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0,
                         limits = c(-umap_lim, umap_lim),
                         oob = scales::squish, name = "sig_94") +
  coord_fixed() +
  labs(title = "sig_94 zscore", x = "UMAP_1", y = "UMAP_2") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

# Top 3 SHH genes
top3_shh <- cor_table$gene[1:3]
shh_umap_panels <- list(p_umap_sig94)

for (g in top3_shh) {
  rho <- cor_table[gene == g, rho_sig94]
  p <- ggplot(umap_coords, aes(x = d1, y = d2, color = .data[[g]])) +
    geom_point(size = 0.15, alpha = 0.5) +
    scale_color_gradient(low = "grey90", high = "#b2182b", name = g) +
    coord_fixed() +
    labs(title = sprintf("%s (ρ=%+.3f)", g, rho),
         x = "UMAP_1", y = "UMAP_2") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())
  shh_umap_panels[[length(shh_umap_panels) + 1]] <- p
}

p_umap_combined <- wrap_plots(shh_umap_panels, ncol = 2) +
  plot_annotation(title = "R21d: sig_94 vs SHH genes on Aldinger UMAP",
                  subtitle = "Visual co-localization check (cell-level Spearman shown per panel)",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(OUT_DIR, "umap_sig94_vs_SHH_genes.png"),
       p_umap_combined, width = 14, height = 12, dpi = 150, bg = "white")

cat("\nUMAP plot saved.\n")
cat(sprintf("Output: %s/umap_sig94_vs_SHH_genes.png\n", OUT_DIR))
