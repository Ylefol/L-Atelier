# Marker gene dot plot

Visualises cluster marker genes as a dot plot where dot **size** encodes
the fraction of cells expressing each gene and dot **colour** encodes
scaled (or raw) mean expression. Genes are ordered by source cluster and
dashed horizontal lines separate marker groups.

## Usage

``` r
ASPIS_plot_marker_dotplot(
  sce,
  markers,
  cluster_col = "cluster",
  assay_name = "logcounts",
  top_n = 5L,
  scale = TRUE,
  min_expr = 0,
  dot_max_size = 6,
  palette = NULL,
  title = "Marker genes"
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

  Logical. Z-score mean expression per gene across clusters before
  plotting. Default `TRUE`.

- min_expr:

  Numeric. Expression threshold used to count a cell as "expressing" the
  gene. Default `0`.

- dot_max_size:

  Numeric. Maximum dot diameter (ggplot size units). Default `6`.

- palette:

  Character vector of colours for the fill gradient. `NULL` (default)
  uses a blue–white–red diverging palette.

- title:

  Character. Plot title. Default `"Marker genes"`.

## Value

A `ggplot` object.
