# shiny_visualization

A self-contained **Shiny app for exploring a Seurat object** (single-cell / spatial
single-cell data). Drop the app next to one or more `.rds` files, launch it, and
interactively browse gene expression, cluster composition, and differential
expression — no coding required. The interface is bilingual (English / 日本語).

A small synthetic dataset (`test_data.rds`) is included so you can try every
feature immediately.

---

## Features

| Tab | What it does |
|-----|--------------|
| 🎻 **Expression** | One tab with shared **gene** and **group-variable** controls and three inner sub-tabs: **Violin** (gene across groups), **Feature UMAP** (gene on the embedding), and **Group UMAP** (embedding colored by a metadata column; hover a point to see its cluster). |
| 📊 **DEG / GSEA** | Differential expression (`FindMarkers`) with a Volcano plot and a searchable table. Supports **categorical** comparisons (Group 1 vs one **or multiple** control groups, or "all others") and **numeric** comparisons (Top X% vs Bottom X% of a continuous score). The tested **genes** can be restricted to *all genes*, a *marker set*, a *set + your own genes*, or *only your own genes*. One click sends the significant **up- or down-DEGs to Enrichr** and opens the results. An inner **GSEA** sub-tab runs `fgsea` against **MSigDB** collections (`msigdbr`) using the **same DEG result** (genes ranked by avg_log2FC), with a top-pathways NES bar plot and a sortable table. |
| 🧱 **Composition** | Stacked proportion bar plot (e.g. cluster composition per donor), with optional **nested faceting** by additional variables (e.g. `site` + `study`). |
| 🗺️ **Spatial** | Plot cells at their tissue x/y coordinates, colored by any metadata variable (categorical → lineage palette, numeric → viridis) **or by a single gene's expression** (viridis), with a **sample** picker to focus on / facet specific tissue sections. If the object has **segmentation** (cell-boundary polygons in @images), a toggle draws filled cell shapes instead of points. You can **highlight selected clusters** (others greyed) and, in a **Neighborhood** sub-tab, plot the nearest-neighbor distance ECDF from the selected clusters to every other cluster (x = distance, log scale; y = cumulative fraction; mean+/-SD band across samples) — closer clusters' curves sit further left; hover a line to see which cluster. Two more sub-tabs reproduce squidpy-style analyses: **Co-occurrence** (P(j \| i,d)/P(j) vs distance, per-sample CI) and **Neighbors enrichment** (permutation z-score matrix on a spatial kNN graph). Pick which loaded dataset to plot from the sidebar. Works on objects that keep coordinates in meta.data (x/y) or in @images. |
| 🧭 **Azimuth** | Automated cell-type annotation with **Pan-Human Azimuth** via **CloudAzimuth** (`AzimuthAPI::CloudAzimuth`). One click sends the active dataset's expression to the satijalab cloud server and returns `azimuth_broad/medium/fine` labels. Shows a UMAP of predicted labels, a summary table (counts + mean confidence), a CSV download, and a **cluster × Azimuth correspondence** heatmap + CSV (cross-tabulate any metadata cluster column against the Azimuth labels, row/column-normalized). *Note: data is sent to an external cloud service.* |
| 🧩 **Niche** | Visualize spatial **niches identified by [tessera](https://github.com/JEFworks-Lab/tessera)**. This tab loads its **own** files (the niche data are plain `data.frame`s, not Seurat objects): a **tile_meta** `.rds` (one `sf` polygon per niche tile + a niche-cluster column such as `seurat_clusters`) and, optionally, a per-cell **annotation** `.rds`. Pick a **sample** and a **niche variable**, and the tiles are drawn as filled polygons colored by niche (interactive: **hover a tile to see its niche**; with optional **niche highlighting** and Y-flip). If a cell-annotation file is loaded, you choose which **cell-type column** to use and get a **niche × cell-type composition** heatmap (each niche row-normalized) + **CSV** of what cell types make up each niche, for the displayed sample or across all samples. The cell→tile link is `tile_id`; if the annotation file lacks it (e.g. a merged fine-cell-type `meta.rds`), it is **auto-recovered from a sibling `clusters_meta`** by matching cell IDs, and per-sample filtering is done by tile membership (so a different `sample_id` naming in the annotation file is fine). Samples without cell-type annotation still draw the map and show a clear message for the composition. |
| 🔥 **Heatmap / Dot** | One tab with shared gene/cluster controls and three inner sub-tabs: a **Heatmap** (mean-expression, optional per-gene Z-scoring and row/column clustering) and a **Dot Plot** (dot size = % expressing, color = scaled mean, optional marker-group facets), and a **Sub-cluster** sub-tab for iterative drill-down (re-cluster selected clusters into nested names, then re-plot + re-correlate vs the reference). With a reference, both plots have a **Combine** mode; the Dot Plot's combine also outputs a **correspondence table** (each active cluster → nearest reference by expression correlation, with a confidence *margin*) that you can **download as CSV**. |

### Highlights

- **Multi-dataset comparison** — load several `.rds` files at once, pick an
  **active** and a **reference** dataset, and every tab (Violin, Feature/Group
  UMAP, Composition, Heatmap, Dot Plot) renders the two **side by side**. The
  **group / cluster variable**, **clusters to draw**, **point size**, **plot
  height** and **plot width** can be set **independently** for the active and
  reference datasets. Set the reference to *(none)* for a single view. The
  Heatmap also has a **Combine with reference** option that plots the active and
  reference clusters in **one** hierarchically-clustered heatmap, so similar
  clusters from the two datasets sit next to each other, with a **dataset color
  bar** down the left so you can tell active vs reference rows at a glance (handy
  for transferring cell-type names across datasets).
- **Pre-registered marker sets** — curated panels (`Cell type`, `T cell`, `B cell`,
  `Myeloid`, `NK / ILC`, `Stromal / Tissue`, `All curated`) drive the Heatmap and
  Dot Plot tabs. On the Heatmap/Dot Plot/DEG tabs you can use a *set*, a
  *set + your own genes*, or *only your own genes*, and prune individual genes
  out of a set. The **Heatmap and Dot Plot live under one "Heatmap / Dot" tab and
  share one set of gene/cluster controls** (at the top of that tab, with inner
  sub-tabs for each plot), so you configure them once and both plots match.
- **Run on demand** — Composition, Heatmap, and Dot Plot render only when you
  click their **Plot** button (they don't recompute on every option change).
- **Interactive plots** — hover to inspect: a Group UMAP point shows its cluster;
  a DEG volcano point shows the gene and stats; a Composition bar segment shows
  the cluster, x value and proportion; a Dot Plot dot shows the gene, cluster,
  percent expressed and scaled mean; a Spatial cell shows its label (via plotly;
  all fall back to static plots if plotly is absent).
- **Dynamic lineage palette & ordering** — when cluster labels follow the
  `<number>_<Lineage>` convention (e.g. `0_B_Plasma`, `2_TNK_ILC`), the app
  automatically (a) assigns a light→dark color gradient within each lineage and
  (b) orders the cluster levels first by lineage, then by the leading number.
  Applied to **Group UMAP**, **Violin**, **Composition**, **Heatmap**,
  **Dot Plot**, and the **DEG** group selectors (so same-lineage clusters stay
  adjacent). Falls back to default colors / natural ordering for non-matching
  variables (e.g. `donor`, `site`). Supported lineage suffixes: `Epithelial`,
  `Stromal_Endothelial`, `TNK_ILC`, `B_Plasma`, `Myeloid`.
- **Cluster sub-selection** — the Composition, Heatmap, and Dot Plot tabs let you
  pick which clusters to draw. In Composition the proportions are **recomputed
  within the selected subset** (e.g. select only the T-cell clusters to see
  composition *within T cells*).
- **External database links** — each gene links out to NCBI Gene, ImmuNexUT, and
  AMP RA2.
- **Light / dark mode** and a **JP / EN** language toggle.
- Adjustable **plot height and width**, and point size.

---

## Requirements

- R (>= 4.2)
- R packages:

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "ggrepel",
                   "ggh4x", "plotly", "patchwork", "ggnewscale",
                   "FNN", "RANN", "msigdbr", "httr", "jsonlite", "sf"))
# Seurat (and its dependency Matrix):
install.packages("Seurat")
# fgsea (Bioconductor) for the GSEA tab:
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install("fgsea")
# AzimuthAPI (optional) for the Pan-Human Azimuth (CloudAzimuth) tab:
devtools::install_github("satijalab/AzimuthAPI")
```

Most packages degrade gracefully if missing: `ggh4x` enables nested facets in
the Composition tab (falls back to `facet_grid`); `plotly` powers the
interactive plots (falls back to static); `patchwork`/`ggnewscale` lay out the
side-by-side comparison and the dataset color bar; `FNN`/`RANN` are used by the
Spatial **Neighborhood / Co-occurrence / Neighbors-enrichment** analyses;
`msigdbr`+`fgsea` drive the **GSEA** sub-tab and `httr`+`jsonlite` the **Enrichr**
links; `sf` is required by the **Niche** tab to draw the tessera niche polygons.
The Heatmap is drawn with plain `ggplot2` (no `pheatmap` needed).

The heaviest analyses use multiple cores when available — the Neighbors-
enrichment permutation test runs via `parallel::mclapply` and `fgsea` via its
`nproc` — and stay **reproducible regardless of core count** (each iteration is
seeded independently). Set the environment variable `SHINY_VIZ_THREADS` to
override the thread count, or use the **Parallel threads** slider in the sidebar
(default = cores − 1, capped at 8; forking falls back to serial on Windows).

---

## Quick start (with the included test data)

```r
# from the repository directory
library(shiny)
runApp(".")
```

Then in the sidebar:

1. Select **`test_data.rds`** and click **Load**.
2. Try each tab. Suggested settings for the demo data:
   - **Group UMAP / Violin** → set *Group Variable* to **`cell_type`** to see the
     lineage-gradient palette.
   - **DEG** → categorical: *Group Variable* = `cell_type`, Group 1 = `5_B_Plasma`,
     Group 2 = one or several other clusters. Numeric: *Group Variable* =
     `inflammation_score` to compare Top vs Bottom percentiles.
   - **Composition** → *Cluster* = `cell_type`, *X-axis* = `donor`,
     *Facet* = `site`, `study`.
   - **Heatmap / Dot Plot** → *Cluster* = `cell_type`, *Marker Set* =
     `Cell type (grouped)`.

### Regenerating the test data

`test_data.rds` is committed to the repo, but you can rebuild it:

```r
Rscript make_test_data.R
```

It creates a 630-cell × 57-gene synthetic Seurat object with 7 clusters across
5 lineages, donor/site/study metadata, a continuous `inflammation_score`, and a
UMAP embedding.

---

## Using your own data

Place any Seurat `.rds` file in this directory and select it in the sidebar.
For the full experience, the object should ideally contain:

- A **UMAP reduction** (any reduction whose name contains `umap`) — required for
  the Feature/Group UMAP tabs.
- **Categorical metadata** columns (clusters, donor, condition, etc.).
- Optionally, **continuous metadata** (module scores, etc.) for numeric DEG.
- Cluster labels in `<number>_<Lineage>` form to enable the lineage palette
  (optional — any labels work, they just use default colors).

---

## Files

| File | Description |
|------|-------------|
| `app.R` | The Shiny application. |
| `make_test_data.R` | Script that generates `test_data.rds`. |
| `test_data.rds` | Synthetic demo dataset. |
| `test_data_reference.rds` | A subsampled variant of the demo data, for trying the reference comparison. |
| `make_spatial_test_data.R` | Script that generates `spatial_test_data.rds`. |
| `spatial_test_data.rds` | Synthetic spatial demo (2 sections, x/y coords) for the Spatial tab. |
