# Elbow plot of PCA variance explained

Plots the percentage of variance explained by each principal component,
with an optional cumulative variance overlay, a user-specified cut-off
line, and algorithmic suggestions to help choose the number of
components.

## Usage

``` r
ASPIS_plot_elbow(
  sce,
  n_show = 50L,
  n_pcs = NULL,
  show_cumulative = TRUE,
  suggest = NULL,
  cum_threshold = 80
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `reducedDims(sce)[["PCA"]]` populated by
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md).

- n_show:

  Integer. Number of components to display. Default `50`.

- n_pcs:

  Integer. If provided, draws a grey dashed line marking the user's
  chosen cut-off. Default `NULL`.

- show_cumulative:

  Logical. Overlay cumulative variance as a red dashed line. Default
  `TRUE`.

- suggest:

  Character or `NULL`. Algorithmic cut-off suggestions to display. One
  of `"elbow"`, `"cumulative"`, `"all"`, or `NULL` (default, no
  suggestion).

- cum_threshold:

  Numeric. Cumulative variance target (%) for the `"cumulative"` method.
  Default `80`.

## Value

A `ggplot` object.

## Details

Two suggestion methods are available and can be overlaid simultaneously:

- `"elbow"`:

  Finds the PC where the second derivative of the variance curve is
  maximised — i.e., where the rate of decline itself decelerates most
  sharply. Tends to suggest fewer PCs for data with a quick early
  drop-off.

- `"cumulative"`:

  Finds the minimum number of PCs needed to explain at least
  `cum_threshold`\\ naturally with how spread-out the variance is across
  components.

Suggestions are printed to the console and drawn as distinct vertical
lines (green = elbow, purple = cumulative) so they can be compared
against one another and against the user's chosen `n_pcs`.
