# Plot Cluster Expression Trajectories

Creates faceted line plots showing expression trajectories per PART
cluster. Matches TiSA's visualization style where each panel shows one
cluster for one group, with individual gene lines (thin, colored) and a
mean trajectory line (thick, gray).

## Usage

``` r
AETHER_plot_cluster_trajectories(
  part_result,
  sample_info,
  time_col = "timepoint",
  group_col = "group",
  norm_counts = NULL,
  clusters = NULL,
  scale_features = FALSE,
  ncol = 4,
  colors = NULL,
  alpha = 0.4,
  title_size = 10
)
```

## Arguments

- part_result:

  An `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).

- sample_info:

  Data.frame with sample metadata. Must have rownames matching colnames
  of the PART data matrix, and columns for timepoint and group.

- time_col:

  Character. Column in sample_info for timepoint. Default: "timepoint".

- group_col:

  Character. Column in sample_info for group. Default: "group".

- norm_counts:

  Optional matrix of normalized counts (genes x samples). If NULL, uses
  the z-scored data from part_result. Default: NULL.

- clusters:

  Character vector of clusters to plot (e.g., c("C1", "C3")). NULL = all
  non-outlier clusters. Default: NULL.

- scale_features:

  Logical. Apply scale feature sum transformation (value / row_sum) as
  in TiSA. Recommended when using norm_counts. Default: FALSE.

- ncol:

  Integer. Number of columns in facet grid. Default: 4.

- colors:

  Named character vector of colors for groups. NULL = auto. Default:
  NULL.

- alpha:

  Numeric. Transparency for individual gene lines. Default: 0.4.

- title_size:

  Numeric. Font size for facet titles. Default: 10.

## Value

A ggplot object.

## Details

This function is ported from TiSA's `plot_cluster_traj()` function. Key
features:

- Faceted by cluster AND group (e.g., "C1 - 50 genes \| IgM")

- Individual gene trajectories as thin colored lines

- Mean cluster trajectory as thick gray line

- X-axis shows timepoints, Y-axis shows expression

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_cluster_trajectories(part_result, sample_info)

} # }
```
