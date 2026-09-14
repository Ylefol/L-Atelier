# Stability Histogram for Consensus Clustering

For each patient, plots how many of the `n_runs` independent runs
assigned them to the same (modal) cluster. A bar concentrated at
`n_runs` indicates high stability.

## Usage

``` r
AETHER_plot_stability_histogram(stability_result, title = NULL)
```

## Arguments

- stability_result:

  Output of
  [`ARTEMIS_consensus_stability()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_stability.md).

- title:

  Character. Plot title. If NULL, auto-generated. Default = NULL.

## Value

A ggplot object.
