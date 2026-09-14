# Plot top-loading features for one view/factor

Creates a bar plot of the highest-magnitude feature weights for one view
on one MOFA factor.

## Usage

``` r
AETHER_plot_mofa_weights(mofa_result, view, factor, top_n = 20, title = NULL)
```

## Arguments

- mofa_result:

  A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
  HEPHAESTUS_load_mofa_model().

- view:

  Character. View name (must match a name in the X list passed to
  HEPHAESTUS_run_mofa()).

- factor:

  Numeric factor index (e.g. 1) or character factor name (e.g.
  "Factor1").

- top_n:

  Number of top-loading features to show, ranked by absolute weight.
  Default: 20.

- title:

  Plot title. If NULL, auto-generated.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
AETHER_plot_mofa_weights(result, view = "rna", factor = 1)

} # }
```
