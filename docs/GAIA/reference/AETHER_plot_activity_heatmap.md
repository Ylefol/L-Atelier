# Plot activity heatmap

Creates a heatmap of activity scores (TFs or pathways vs samples).

## Usage

``` r
AETHER_plot_activity_heatmap(
  result,
  top_n = 50,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = NULL,
  color_palette = NULL,
  title = NULL,
  show_rownames = NULL,
  show_colnames = NULL,
  ...
)
```

## Arguments

- result:

  A decoupler_result or decoupler_comparison object, OR a matrix of
  activity scores (sources as rows, samples as columns).

- top_n:

  Integer. Number of top sources to display (by absolute mean activity).
  Default: 50. Use NULL to show all.

- cluster_rows:

  Logical. Cluster rows (sources). Default: TRUE.

- cluster_cols:

  Logical. Cluster columns (samples). Default: TRUE.

- scale:

  Character. Scale data by "row", "column", or "none". Default: "row".

- annotation_col:

  Data.frame of sample annotations (rownames = sample names). Default:
  NULL.

- color_palette:

  Character vector of colors for heatmap. Default: blue-white-red.

- title:

  Character. Plot title. Default: auto-generated.

- show_rownames:

  Logical. Show source names. Default: TRUE if \<= 50 sources.

- show_colnames:

  Logical. Show sample names. Default: TRUE if \<= 30 samples.

- ...:

  Additional arguments passed to pheatmap::pheatmap().

## Value

A pheatmap object (invisibly).

## Details

If a decoupler_comparison object is provided, uses the consensus scores
by default. To plot a specific method's results, extract it first:
comparison\$results\$ulm\$activities

## Examples

``` r
if (FALSE) { # \dontrun{
# From decoupler result
AETHER_plot_activity_heatmap(tf_result)

# With sample annotations
annot <- data.frame(
  condition = c("ctrl", "ctrl", "treat", "treat"),
  row.names = colnames(tf_result$activities)
)
AETHER_plot_activity_heatmap(tf_result, annotation_col = annot)

# From comparison (uses consensus)
AETHER_plot_activity_heatmap(comparison)

} # }
```
