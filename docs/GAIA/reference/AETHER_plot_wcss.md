# Plot Within-Cluster Sum-of-Squares (WCSS) Elbow

Plots WCSS against k. The "elbow" — where adding another cluster gives
diminishing returns — suggests the optimal k.

## Usage

``` r
AETHER_plot_wcss(cc_result, title = "WCSS Elbow Plot")
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- title:

  Character. Plot title. Default = "WCSS Elbow Plot".

## Value

A ggplot object.
