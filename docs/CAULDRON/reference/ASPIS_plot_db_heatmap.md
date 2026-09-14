# Marker database heatmap

Visualises cluster-averaged expression of genes drawn from a
user-supplied marker database. Rows are genes grouped into facets by
cell type; columns are clusters. Genes absent from the SCE are retained
as grey `NA` tiles so coverage gaps are immediately visible. When a
`neg_markers` column is present, positive and negative markers form
nested sub-facets within each cell type.

## Usage

``` r
ASPIS_plot_db_heatmap(
  sce,
  markers,
  cluster_col = "cluster",
  assay_name = "logcounts",
  scale = TRUE,
  na_color = "grey80",
  palette = NULL,
  show_gene_labels = TRUE,
  title = "Marker database heatmap"
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `logcounts` (or other) assay.

- markers:

  A `data.frame` with columns:

  cell_type

  :   Cell type name used for row faceting.

  markers

  :   Positive marker genes, comma-separated string.

  neg_markers

  :   (Optional) Negative marker genes, comma-separated. When present,
      positive and negative genes are shown in nested sub-facets.

- cluster_col:

  Character. `colData` column containing cluster labels used for column
  grouping. Default `"cluster"`.

- assay_name:

  Character. Assay to pull expression from. Default `"logcounts"`.

- scale:

  Logical. Z-score mean expression per gene across clusters. Default
  `TRUE`.

- na_color:

  Character. Fill colour for genes absent from the SCE. Default
  `"grey80"`.

- palette:

  Character vector of colours for the fill gradient. `NULL` (default)
  uses a blue-white-red diverging palette.

- show_gene_labels:

  Logical. Show gene names on the y axis. Default `TRUE`.

- title:

  Character. Plot title. Default `"Marker database heatmap"`.

## Value

A `ggplot` object.
