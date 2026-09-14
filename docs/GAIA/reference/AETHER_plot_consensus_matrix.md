# Plot Consensus Matrix Heatmap

Plots the averaged consensus matrix for a chosen k as a heatmap. Samples
are reordered by cluster assignment, making block structure visible.
Dark blue = always co-clustered; white = never co-clustered.

## Usage

``` r
AETHER_plot_consensus_matrix(cc_result, k, title = NULL)
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- k:

  Integer. Number of clusters to display.

- title:

  Character. Plot title. If NULL, auto-generated. Default = NULL.

## Value

A ggplot object.
