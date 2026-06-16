# =============================================================================
# make_spatial_test_data.R
# Generate a small synthetic *spatial* Seurat object to demo the Spatial tab.
#
# Usage (from this directory):
#   Rscript make_spatial_test_data.R
#
# Writes `spatial_test_data.rds` with x/y coordinates in meta.data, a `sample`
# column (2 tissue sections), plus cell_type / lineage labels — so you can pick a
# color variable and a sample and plot the spatial layout. It also includes
# **segmentation** (per-cell boundary polygons) in @images, so the Spatial tab's
# "Show segmentation" option can draw filled cell shapes.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(123)

# --- Lineage-specific marker genes (same convention as make_test_data.R) ------
lineage_markers <- list(
  Epithelial          = c("EPCAM", "KRT19", "KRT5", "MUC5B", "CAV1"),
  Stromal_Endothelial = c("PECAM1", "VWF", "COL1A1", "PDGFRA", "ACTA2"),
  TNK_ILC             = c("CD3E", "CD8A", "NKG7", "GZMB", "FOXP3"),
  B_Plasma            = c("CD19", "MS4A1", "MZB1", "XBP1", "IGHG1"),
  Myeloid             = c("CD68", "C1QA", "LYZ", "CD14", "MARCO")
)
shared_genes <- c("MKI67", "PTPRC", "VIM", "ACTB", "HLA-DRA")
all_genes <- unique(c(unlist(lineage_markers, use.names = FALSE), shared_genes))

clusters <- c("0_Epithelial", "1_Stromal_Endothelial", "2_TNK_ILC",
              "3_B_Plasma", "4_Myeloid")
cluster_lineage <- sub("^\\d+_", "", clusters)

# --- Two tissue sections; each cell gets a cluster + 2D position --------------
make_section <- function(section, n) {
  # Cluster centers placed in a 2D tissue (different layout per section)
  ang <- if (section == "sectionA") 0 else pi / 5
  base <- matrix(c(20, 20,  60, 25,  40, 55,  75, 65,  30, 80),
                 ncol = 2, byrow = TRUE)
  rot <- matrix(c(cos(ang), -sin(ang), sin(ang), cos(ang)), 2)
  centers <- base %*% rot + (if (section == "sectionB") 5 else 0)
  cl <- sample(seq_along(clusters), n, replace = TRUE,
               prob = c(0.30, 0.22, 0.20, 0.10, 0.18))
  xy <- centers[cl, ] + matrix(rnorm(n * 2, sd = 6), ncol = 2)
  data.frame(
    cluster = clusters[cl],
    lineage = cluster_lineage[cl],
    x = xy[, 1], y = xy[, 2],
    sample = section,
    stringsAsFactors = FALSE
  )
}

meta_df <- rbind(make_section("sectionA", 1200), make_section("sectionB", 1000))
n_cells <- nrow(meta_df)
cell_ids <- sprintf("cell_%04d", seq_len(n_cells))
rownames(meta_df) <- cell_ids

# --- Simulate counts (lineage markers up in their lineage) --------------------
counts <- matrix(rpois(length(all_genes) * n_cells, lambda = 0.3),
                 nrow = length(all_genes), ncol = n_cells,
                 dimnames = list(all_genes, cell_ids))
for (i in seq_len(n_cells)) {
  mk <- intersect(lineage_markers[[meta_df$lineage[i]]], all_genes)
  counts[mk, i] <- counts[mk, i] + rpois(length(mk), lambda = 8)
}
counts <- Matrix(counts, sparse = TRUE)

# --- Build the Seurat object --------------------------------------------------
obj <- CreateSeuratObject(counts = counts, project = "spatial_demo",
                          meta.data = meta_df)
obj <- NormalizeData(obj, verbose = FALSE)
obj$cell_type <- factor(meta_df$cluster, levels = clusters)
obj$lineage   <- factor(meta_df$lineage, levels = unique(cluster_lineage))
obj$seurat_clusters <- factor(sub("_.*$", "", meta_df$cluster), levels = as.character(0:4))
# continuous variable to demo numeric coloring
obj$density_score <- as.numeric(scale(meta_df$x + meta_df$y) + rnorm(n_cells, sd = 0.3))

# --- Segmentation: one irregular polygon (cell boundary) per cell -------------
# Build a FOV with both centroids and segmentation, one image per section.
make_polygon <- function(cx, cy, r, k = 8) {
  ang <- sort(runif(k, 0, 2 * pi))
  rad <- r * runif(k, 0.7, 1.1)        # irregular radius -> cell-like shape
  data.frame(x = cx + rad * cos(ang), y = cy + rad * sin(ang))
}
for (sec in c("sectionA", "sectionB")) {
  idx <- which(meta_df$sample == sec)
  cids <- rownames(meta_df)[idx]
  cent <- data.frame(x = meta_df$x[idx], y = meta_df$y[idx], cell = cids)
  poly <- do.call(rbind, lapply(seq_along(idx), function(j) {
    pg <- make_polygon(meta_df$x[idx[j]], meta_df$y[idx[j]], r = 2.2)
    pg$cell <- cids[j]; pg
  }))
  fov <- CreateFOV(
    coords = list(segmentation = CreateSegmentation(poly),
                  centroids    = CreateCentroids(cent)),
    type = c("segmentation", "centroids"), assay = "RNA")
  obj[[sec]] <- subset(fov, cells = cids)
}

saveRDS(obj, "spatial_test_data.rds")
cat(sprintf("Saved spatial_test_data.rds : %d cells x %d genes\n", ncol(obj), nrow(obj)))
cat("samples:\n"); print(table(obj$sample))
cat("images:", paste(names(obj@images), collapse = ", "), "\n")
cat("boundaries[1]:", paste(SeuratObject::Boundaries(obj@images[[1]]), collapse = ", "), "\n")
cat("x range:", round(range(meta_df$x), 1), " y range:", round(range(meta_df$y), 1), "\n")
