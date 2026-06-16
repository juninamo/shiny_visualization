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
| 🎻 **Violin** | Violin plot of a gene across a chosen grouping variable. |
| 🗺️ **Feature UMAP** | Gene expression overlaid on the UMAP embedding. |
| 🏷️ **Group UMAP** | UMAP colored by any categorical metadata column. |
| 📊 **DEG** | Differential expression (`FindMarkers`) with a Volcano plot and a searchable table. Supports **categorical** comparisons (Group 1 vs one **or multiple** control groups, or "all others") and **numeric** comparisons (Top X% vs Bottom X% of a continuous score). |
| 🧱 **Composition** | Stacked proportion bar plot (e.g. cluster composition per donor), with optional **nested faceting** by additional variables (e.g. `site` + `study`). |
| 🔥 **Heatmap** | Mean-expression heatmap of a pre-registered marker-gene set, with optional per-gene Z-scoring and row/column clustering. |
| 🔵 **Dot Plot** | Dot plot of a marker set: dot size = % expressing, color = scaled mean expression, optionally faceted by marker group. |

### Highlights

- **Pre-registered marker sets** — curated panels (`Cell type`, `T cell`, `B cell`,
  `Myeloid`, `NK / ILC`, `Stromal / Tissue`, `All curated`) drive the Heatmap and
  Dot Plot tabs out of the box.
- **Dynamic lineage palette** — when cluster labels follow the
  `<number>_<Lineage>` convention (e.g. `0_B_Plasma`, `2_TNK_ILC`), the app
  automatically assigns a light→dark color gradient within each lineage. This is
  applied to the **Composition**, **Group UMAP**, and **Violin** tabs, and falls
  back to default colors for non-matching variables (e.g. `donor`, `site`).
  Supported lineage suffixes: `Epithelial`, `Stromal_Endothelial`, `TNK_ILC`,
  `B_Plasma`, `Myeloid`.
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
                   "pheatmap", "ggh4x"))
# Seurat (and its dependency Matrix):
install.packages("Seurat")
```

`pheatmap` is required for the Heatmap tab; `ggh4x` enables nested facets in the
Composition tab (it gracefully falls back to `facet_grid` if absent).

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
