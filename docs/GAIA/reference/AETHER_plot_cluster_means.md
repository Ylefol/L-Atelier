# Plot Cluster Mean Trajectories (Overlay)

Creates a single panel showing mean expression trajectories for all
clusters, colored by cluster. Provides a compact overview of cluster
behavior across timepoints. When groups are present, lines are
distinguished by linetype.

## Usage

``` r
AETHER_plot_cluster_means(
  part_result,
  sample_info,
  time_col = "timepoint",
  group_col = "group",
  clusters = NULL,
  colors = NULL,
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

- time_col:

  Character. Column in sample_info for timepoint. Default: "timepoint".

- group_col:

  Character or NULL. Column in sample_info for group. If NULL, averages
  across all samples per timepoint. Default: "group".

- clusters:

  Character vector of clusters to plot. NULL = all non-outlier. Default:
  NULL.

- colors:

  Named character vector of colors for clusters. NULL = uses PART
  cluster_colors. Default: NULL.

- title:

  Character or NULL. Plot title. Default: NULL (auto-generated).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_cluster_means(part_result, sample_info)

} # }
```
