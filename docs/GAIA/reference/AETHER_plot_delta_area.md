# Plot Delta Area (Relative CDF Area Change)

Plots the relative change in CDF area between consecutive k values. A
large positive delta suggests that k explains substantially more
structure than k-1. Values plateau or drop near the optimal k.

## Usage

``` r
AETHER_plot_delta_area(cc_result, title = "Delta Area")
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- title:

  Character. Plot title. Default = "Delta Area".

## Value

A ggplot object.
