# =============================================================================
# make_test_data.R
# Generate a small synthetic Seurat object to demo the Seurat Viewer app.
#
# Usage (from this directory):
#   Rscript make_test_data.R
#
# It writes `test_data.rds` next to app.R. Launch the app with:
#   library(shiny); runApp(".")
# and pick `test_data.rds` in the sidebar.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

# --- 1. Define lineage-specific marker genes ---------------------------------
# Cluster labels follow the "<number>_<Lineage>" convention so the app's
# dynamic lineage-gradient palette is triggered automatically.
lineage_markers <- list(
  Epithelial          = c("EPCAM", "KRT19", "KRT5", "KRT14", "TP63",
                          "MUC5B", "CEACAM6", "CAV1"),
  Stromal_Endothelial = c("PECAM1", "VWF", "CDH5", "FLT1", "COL1A1", "COL3A1",
                          "PDGFRA", "THY1", "PRG4", "CD34", "ACTA2", "NOTCH3"),
  TNK_ILC             = c("CD3E", "CD4", "CD8A", "IL7R", "NKG7", "GNLY",
                          "GZMB", "GZMK", "FOXP3", "NCAM1", "CCL5"),
  B_Plasma            = c("CD19", "MS4A1", "CD79A", "MZB1", "XBP1", "IGHG1",
                          "PRDM1", "TCL1A"),
  Myeloid             = c("CD68", "C1QA", "LYZ", "CD14", "CD163", "FCGR3A",
                          "CLEC9A", "CD1C", "IRF7", "MARCO")
)
shared_genes <- c("MKI67", "TOP2A", "MALAT1", "ACTB", "HLA-DRA",
                  "PTPRC", "VIM", "S100A4")
all_genes <- unique(c(unlist(lineage_markers, use.names = FALSE), shared_genes))

# --- 2. Define clusters (number_lineage) -------------------------------------
clusters <- c(
  "0_Epithelial",
  "1_Epithelial",
  "2_Stromal_Endothelial",
  "3_TNK_ILC",
  "4_TNK_ILC",
  "5_B_Plasma",
  "6_Myeloid"
)
cluster_lineage <- sub("^\\d+_", "", clusters)
n_per_cluster <- 90
cell_clusters <- rep(clusters, each = n_per_cluster)
n_cells <- length(cell_clusters)
cell_ids <- sprintf("cell_%04d", seq_len(n_cells))

# --- 3. Simulate a counts matrix (genes x cells) -----------------------------
counts <- matrix(
  rpois(length(all_genes) * n_cells, lambda = 0.3),
  nrow = length(all_genes), ncol = n_cells,
  dimnames = list(all_genes, cell_ids)
)
for (i in seq_len(n_cells)) {
  lin <- cluster_lineage[match(cell_clusters[i], clusters)]
  mk <- intersect(lineage_markers[[lin]], all_genes)
  counts[mk, i] <- counts[mk, i] + rpois(length(mk), lambda = 8)
}
# Proliferation markers up in one epithelial and one T/NK cluster
prolif_cells <- cell_clusters %in% c("1_Epithelial", "4_TNK_ILC")
counts[c("MKI67", "TOP2A"), prolif_cells] <-
  counts[c("MKI67", "TOP2A"), prolif_cells] +
  rpois(2 * sum(prolif_cells), lambda = 6)
counts <- Matrix(counts, sparse = TRUE)

# --- 4. Build the Seurat object ----------------------------------------------
obj <- CreateSeuratObject(counts = counts, project = "demo")
obj <- NormalizeData(obj, verbose = FALSE)

# --- 5. Metadata -------------------------------------------------------------
# Conventional numeric cluster id + descriptive lineage label
obj$seurat_clusters <- factor(sub("_.*$", "", cell_clusters),
                              levels = as.character(0:6))
obj$cell_type <- factor(cell_clusters, levels = clusters)
obj$lineage   <- factor(cluster_lineage, levels = unique(cluster_lineage))

# Donors nested within site + study (for the Composition tab)
donor_table <- data.frame(
  donor = sprintf("D%02d", 1:8),
  study = rep(c("StudyA", "StudyB"), each = 4),
  site  = rep(c("Synovium", "Synovium", "Tonsil", "Tonsil"), times = 2),
  stringsAsFactors = FALSE
)
cell_donor_idx <- sample(seq_len(nrow(donor_table)), n_cells, replace = TRUE)
obj$donor <- donor_table$donor[cell_donor_idx]
obj$study <- donor_table$study[cell_donor_idx]
obj$site  <- donor_table$site[cell_donor_idx]

# Continuous score (numeric DEG mode: Top X% vs Bottom X%)
base_score <- c(Epithelial = 0.2, Stromal_Endothelial = 0.3,
                TNK_ILC = 0.7, B_Plasma = 0.5, Myeloid = 0.9)
obj$inflammation_score <- as.numeric(
  base_score[cluster_lineage] + rnorm(n_cells, sd = 0.15)
)

# --- 6. A fake UMAP embedding (clusters separated; no uwot needed) -----------
centers <- matrix(
  c( 6,  6,   # 0_Epithelial
     7,  4,   # 1_Epithelial
    -6,  5,   # 2_Stromal_Endothelial
    -5, -5,   # 3_TNK_ILC
    -7, -3,   # 4_TNK_ILC
     5, -6,   # 5_B_Plasma
     0, -1),  # 6_Myeloid
  ncol = 2, byrow = TRUE,
  dimnames = list(clusters, c("UMAP_1", "UMAP_2"))
)
umap_emb <- centers[cell_clusters, ] + matrix(rnorm(n_cells * 2, sd = 0.8),
                                              ncol = 2)
rownames(umap_emb) <- cell_ids
obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_emb, key = "UMAP_", assay = DefaultAssay(obj)
)

Idents(obj) <- obj$cell_type

# --- 7. Save -----------------------------------------------------------------
out <- "test_data.rds"
saveRDS(obj, out)
cat(sprintf("Saved %s : %d cells x %d genes\n", out, ncol(obj), nrow(obj)))
print(table(obj$cell_type))
