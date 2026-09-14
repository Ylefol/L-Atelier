# Plot Silhouette for Cluster Quality Assessment

Visualize silhouette widths for each sample, grouped by cluster.
Silhouette width measures how similar a sample is to its own cluster
compared to other clusters. Values range from -1 to 1, where higher is
better.

## Usage

``` r
AETHER_plot_silhouette(
  cluster_result,
  colors = NULL,
  show_avg = TRUE,
  title = "Silhouette Plot"
)
```

## Arguments

- cluster_result:

  An artemis_cluster object from ARTEMIS_cluster_mixed().

- colors:

  Named character vector. Custom colors for clusters. If NULL (default),
  uses a default palette.

- show_avg:

  Logical. Show average silhouette lines per cluster. Default = TRUE.

- title:

  Character. Plot title. Default = "Silhouette Plot".

## Value

A ggplot object

## Details

Interpretation of silhouette width:

- \> 0.7: Strong structure

- 0.5 - 0.7: Reasonable structure

- 0.25 - 0.5: Weak structure, could be artificial

- \< 0.25: No substantial structure

Negative values indicate samples that may be misclassified (closer to
another cluster than their assigned one).

## Examples

``` r
if (FALSE) { # \dontrun{
clust <- ARTEMIS_cluster_mixed(data)
AETHER_plot_silhouette(clust)

} # }
```
