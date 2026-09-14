# Plot Radar/Spider Chart for Cluster Profiles

Visualize how clusters differ across multiple variables using a radar
(spider) chart. Each axis represents a variable, and each polygon
represents a cluster's profile.

## Usage

``` r
AETHER_plot_cluster_radar(
  cluster_char,
  variables = NULL,
  data = NULL,
  clusters = NULL,
  normalize = TRUE,
  colors = NULL,
  alpha = 0.2,
  line_width = 1,
  title = "Cluster Profiles"
)
```

## Arguments

- cluster_char:

  An artemis_characterization object from
  ARTEMIS_characterize_clusters(), OR a named list of cluster means.

- variables:

  Character vector. Which variables to include in the plot. If NULL
  (default), uses top variables from characterization (max 10).

- data:

  Optional data frame. Required if cluster_char doesn't contain the
  original data or if you want to use different variables.

- clusters:

  Optional. Cluster assignments vector. Required if providing data
  instead of cluster_char.

- normalize:

  Logical. Normalize variables to 0-1 range for comparability. Default =
  TRUE.

- colors:

  Named character vector. Custom colors for clusters. If NULL (default),
  uses a default palette.

- alpha:

  Numeric. Fill transparency (0-1). Default = 0.2.

- line_width:

  Numeric. Width of polygon borders. Default = 1.

- title:

  Character. Plot title. Default = "Cluster Profiles".

## Value

A ggplot object

## Details

For quantitative variables, the plot shows the mean value per cluster.
Variables are normalized to 0-1 range by default so they're comparable
on the same axes.

Qualitative variables are excluded as they don't have a natural ordering
for radar plots. Use the cluster profiles table from characterization to
examine qualitative variable distributions.

## Examples

``` r
if (FALSE) { # \dontrun{
char <- ARTEMIS_characterize_clusters(clust)
AETHER_plot_cluster_radar(char)

# With specific variables
AETHER_plot_cluster_radar(char, variables = c("age", "score", "measure1"))

} # }
```
