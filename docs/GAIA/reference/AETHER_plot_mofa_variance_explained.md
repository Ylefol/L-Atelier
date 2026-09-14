# Plot variance explained per view and factor

Creates a heatmap showing how much variance each MOFA factor explains in
each omics view.

## Usage

``` r
AETHER_plot_mofa_variance_explained(
  mofa_result,
  show_values = TRUE,
  colors = NULL,
  title = "Variance Explained per View and Factor"
)
```

## Arguments

- mofa_result:

  A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
  HEPHAESTUS_load_mofa_model().

- show_values:

  Logical. Show variance explained (%) as cell text. Default: TRUE.

- colors:

  Two-color low/high gradient. Default: white to dark blue.

- title:

  Plot title. Default: "Variance Explained per View and Factor".

## Value

A ggplot object. Faceted by group if the model has more than one.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
AETHER_plot_mofa_variance_explained(result)

} # }
```
