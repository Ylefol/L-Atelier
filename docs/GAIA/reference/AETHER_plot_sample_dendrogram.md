# Plot sample dendrogram with trait colors

Creates a sample clustering dendrogram with trait values shown as
colored bars underneath.

## Usage

``` r
AETHER_plot_sample_dendrogram(
  cluster_result,
  traits = NULL,
  main = "Sample dendrogram and trait heatmap",
  ...
)
```

## Arguments

- cluster_result:

  A wgcna_cluster object from ARTEMIS_wgcna_cluster_samples(), or an
  hclust object.

- traits:

  Data.frame of traits (samples as rows) or wgcna_data object.

- main:

  Plot title.

- ...:

  Additional arguments passed to plotDendroAndColors().

## Value

Invisible NULL. Plots to current device.

## Examples

``` r
if (FALSE) { # \dontrun{
clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data)
AETHER_plot_sample_dendrogram(clust, wgcna_data)

} # }
```
