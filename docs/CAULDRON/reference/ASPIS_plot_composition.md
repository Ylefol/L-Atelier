# Plot a 100% stacked bar chart of per-sample cell-type composition

Complements
[`ASPIS_plot_propeller`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_propeller.md):
that function shows one facet per cluster (a statistical test result in
isolation), which never makes the compositional nature of the data
visible. This shows the full picture instead – one bar per sample,
cluster proportions stacked to the full bar height, so it's immediately
clear that every sample's cell-type proportions sum to 100% (an increase
in one cluster necessarily comes at the expense of others).

## Usage

``` r
ASPIS_plot_composition(
  propeller_result,
  palette = NULL,
  facet_by_group = TRUE,
  title = NULL
)
```

## Arguments

- propeller_result:

  A `keraunos_propeller` object from
  [`KERAUNOS_propeller_proportions`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md)
  (uses its `$proportions` and `$results$baseline_prop`/`$FDR` fields).

- palette:

  Character vector of colours, one per cluster. Overrides any
  `cluster_colors` carried on `propeller_result`. Default `NULL`: uses
  `propeller_result$params$cluster_colors` if it was supplied to
  [`KERAUNOS_propeller_proportions`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md),
  otherwise falls back to the built-in 20-colour Tableau-inspired set
  (recycled if there are more clusters than colours).

- facet_by_group:

  Logical. Split samples into group-labelled facet panels so they're
  grouped visually by experimental group rather than in arbitrary order.
  Default `TRUE`.

- title:

  Character. Plot title. Default `NULL`.

## Value

A ggplot object.

## Details

Each legend entry is two lines: the cluster label, then `FDR = ...` from
the same propeller test – so the statistical context from
[`ASPIS_plot_propeller`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_propeller.md)
carries over without needing to cross-reference a second plot. Legend
key spacing is widened accordingly so the two-line entries don't crowd
each other.
