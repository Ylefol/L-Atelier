# Alluvial Plot of Cluster Stability Across Runs

Shows how patients flow between cluster assignments across `n_runs`
independent runs. Stable patients remain within the same stratum;
unstable patients cross between strata.

## Usage

``` r
AETHER_plot_stability_alluvial(stability_result, title = NULL)
```

## Arguments

- stability_result:

  Output of
  [`ARTEMIS_consensus_stability()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_stability.md).

- title:

  Character. Plot title. If NULL, auto-generated. Default = NULL.

## Value

A ggplot object.

## Details

Requires the `ggalluvial` package. If not available, falls back to
[`AETHER_plot_stability_histogram()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_stability_histogram.md).
