# Plot cell type proportions as heatmap

Creates a heatmap of cell type proportions (cell types as rows, samples
as columns) using pheatmap.

## Usage

``` r
AETHER_plot_cibersort_heatmap(
  result,
  annotation_col = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color_palette = NULL,
  title = "Cell Type Proportions",
  legend_title = "Proportion",
  show_rownames = TRUE,
  show_colnames = TRUE,
  ...
)
```

## Arguments

- result:

  An artemis_cibersort object or a proportions matrix (samples x cell
  types).

- annotation_col:

  Data.frame of sample annotations for column annotation. Rownames must
  match sample names. Default: NULL.

- cluster_rows:

  Logical. Cluster cell types. Default: TRUE.

- cluster_cols:

  Logical. Cluster samples. Default: TRUE.

- scale:

  Character. Scale by "row", "column", or "none". Default: "none"
  (proportions are already comparable).

- color_palette:

  Character vector of colors for heatmap gradient. Default:
  white-to-dark-red sequential palette.

- title:

  Character. Plot title. Default: "Cell Type Proportions".

- legend_title:

  Character. Title for the color legend. Default: "Proportion".

- show_rownames:

  Logical. Show cell type names. Default: TRUE.

- show_colnames:

  Logical. Show sample names. Default: auto (TRUE if \<= 30 samples).

- ...:

  Additional arguments passed to pheatmap::pheatmap().

## Value

A pheatmap object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_cibersort_heatmap(cibersort_result)

# With sample annotations
annot <- data.frame(condition = groups, row.names = sample_names)
AETHER_plot_cibersort_heatmap(cibersort_result, annotation_col = annot)

} # }
```
