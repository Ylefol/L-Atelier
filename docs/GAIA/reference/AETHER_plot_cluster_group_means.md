# Plot Cluster Mean Expression by Group

Companion plot to
[`AETHER_plot_part_heatmap()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_part_heatmap.md)
for cross-sectional (non-time-series) PART results. Shows the mean
expression per cluster per group as a dodged bar chart, offering a
compact summary view when the heatmap itself is too large to read
column-by-column.

## Usage

``` r
AETHER_plot_cluster_group_means(
  part_result,
  sample_info,
  sample_col = NULL,
  group_col = "group",
  clusters = NULL,
  group_colors = NULL,
  error_bars = c("se", "sd", "none"),
  title = NULL
)
```

## Arguments

- part_result:

  An `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).

- sample_info:

  Data.frame with sample metadata. Must have rownames matching colnames
  of the PART data matrix.

- group_col:

  Character. Group column in sample_info. Default: "group".

- clusters:

  Character vector of clusters to plot. NULL = all non-outlier clusters,
  in the same order as `part_result$cluster_colors` (matching the
  row-block order of the PART heatmap). Default: NULL.

- group_colors:

  Named character vector of colors for groups. NULL = generated with the
  same default palette/ordering logic used by
  [`AETHER_plot_part_heatmap()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_part_heatmap.md),
  so colors match between the two plots when both are called with the
  same `sample_info`/`group_col`. Default: NULL.

- error_bars:

  Character. "se" (standard error), "sd" (standard deviation), or
  "none". Default: "se".

- title:

  Character or NULL. Plot title. Default: NULL (auto-generated).

## Value

A ggplot object.

## Details

For each sample, the mean z-score across genes in a cluster is computed
first; bars then show the mean (± error) of those per-sample cluster
means within each group. Cluster x-axis labels are colored using
`part_result$cluster_colors` to visually match the heatmap's row
annotation blocks.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_cluster_group_means(part_result, sample_info, group_col = "Group")

} # }
```
