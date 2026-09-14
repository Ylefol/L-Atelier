# Marker gene heatmap

Visualises cluster-averaged expression of marker genes as a tiled
heatmap. When `markers` is a `keraunos_markers` object, genes are
grouped into facets by their source cluster so marker identity is
immediately clear. Columns are clusters; rows are genes; colour encodes
z-scored (or raw) mean expression.

## Usage

``` r
ASPIS_plot_marker_heatmap(
  sce,
  markers,
  cluster_col = "cluster",
  assay_name = "logcounts",
  top_n = 5L,
  scale = TRUE,
  palette = NULL,
  show_gene_labels = TRUE,
  title = "Marker gene heatmap"
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `logcounts` (or other) assay.

- markers:

  A `keraunos_markers` object returned by
  [`KERAUNOS_find_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_find_markers.md),
  or a character vector of gene names.

- cluster_col:

  Character. `colData` column containing cluster labels. Default
  `"cluster"`.

- assay_name:

  Character. Assay to pull expression from. Default `"logcounts"`.

- top_n:

  Integer. Maximum marker genes per cluster shown when `markers` is a
  `keraunos_markers` object. Default `5L`.

- scale:

  Logical. Z-score mean expression per gene across clusters. Default
  `TRUE`.

- palette:

  Character vector of colours for the fill gradient. `NULL` (default)
  uses a blue–white–red diverging palette.

- show_gene_labels:

  Logical. Show gene names on the y axis. Default `TRUE`.

- title:

  Character. Plot title. Default `"Marker gene heatmap"`.

## Value

A `ggplot` object.
