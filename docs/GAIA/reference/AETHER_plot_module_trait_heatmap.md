# Plot module-trait correlation heatmap

Creates a heatmap showing correlations between module eigengenes and
traits, with correlation values and p-values displayed in cells.

## Usage

``` r
AETHER_plot_module_trait_heatmap(
  trait_cor,
  use_padj = FALSE,
  show_values = TRUE,
  colors = NULL,
  main = "Module-trait relationships",
  cex.text = 0.7,
  ...
)
```

## Arguments

- trait_cor:

  A wgcna_trait_cor object from ARTEMIS_wgcna_module_traits().

- use_padj:

  Logical. Use adjusted p-values for display. Default: FALSE (uses raw
  p-values like original WGCNA).

- show_values:

  Logical. Show correlation and p-value in cells. Default: TRUE.

- colors:

  Color palette for heatmap. Default: blueWhiteRed(50).

- main:

  Plot title. Default: "Module-trait relationships".

- cex.text:

  Text size for cell values. Default: 0.7.

- ...:

  Additional arguments passed to labeledHeatmap().

## Value

Invisible NULL. Plots to current device.

## Examples

``` r
if (FALSE) { # \dontrun{
trait_cor <- ARTEMIS_wgcna_module_traits(modules, wgcna_data)
AETHER_plot_module_trait_heatmap(trait_cor)

} # }
```
