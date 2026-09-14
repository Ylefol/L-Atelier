# Plot method comparison correlation heatmap

Visualizes the correlation between different methods' activity scores.

## Usage

``` r
AETHER_plot_method_correlation(
  comparison,
  title = "Method Correlation",
  color_palette = NULL,
  show_values = TRUE,
  ...
)
```

## Arguments

- comparison:

  A decoupler_comparison object.

- title:

  Character. Plot title. Default: "Method Correlation".

- color_palette:

  Character vector of colors. Default: white to blue.

- show_values:

  Logical. Display correlation values. Default: TRUE.

- ...:

  Additional arguments passed to pheatmap::pheatmap().

## Value

A pheatmap object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- ARTEMIS_decoupler_compare_methods(mat, network, methods = c("ulm", "mlm", "wsum"))
AETHER_plot_method_correlation(comparison)

} # }
```
