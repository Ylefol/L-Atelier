# Plot Consensus CDF Curves

Plots empirical cumulative distribution functions (CDFs) of pairwise
consensus values for each k. A well-separated clustering shows a CDF
with values concentrated near 0 and 1 (flat/horizontal line). The
optimal k is where the CDF stops improving substantially.

## Usage

``` r
AETHER_plot_consensus_cdf(cc_result, k_range = NULL, title = "Consensus CDF")
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- k_range:

  Integer vector. Which k values to include. If NULL, uses all. Default
  = NULL.

- title:

  Character. Plot title. Default = "Consensus CDF".

## Value

A ggplot object.
