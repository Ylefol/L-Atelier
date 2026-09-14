# Multi-Panel Diagnostic Summary for Consensus Clustering

Combines CDF, WCSS elbow, and delta area into a single three-panel
diagnostic figure.

## Usage

``` r
AETHER_plot_consensus_summary(
  cc_result,
  title = "Consensus Clustering Diagnostics"
)
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- title:

  Character. Overall figure title. Default = "Consensus Clustering
  Diagnostics".

## Value

A ggplot object (patchwork layout).

## Details

Requires the `patchwork` package.
