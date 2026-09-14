# Plot module-trait heatmap using ggplot2

Alternative ggplot2-based heatmap for module-trait correlations. More
customizable than the base R version.

## Usage

``` r
AETHER_plot_module_trait_heatmap_gg(
  trait_cor,
  use_padj = FALSE,
  sig_threshold = 0.05,
  show_values = "both",
  pval_format = "decimal",
  pval_digits = 3,
  colors = NULL,
  title = "Module-trait relationships"
)
```

## Arguments

- trait_cor:

  A wgcna_trait_cor object from ARTEMIS_wgcna_module_traits().

- use_padj:

  Logical. Use adjusted p-values. Default: FALSE.

- sig_threshold:

  P-value threshold for significance stars. Default: 0.05.

- show_values:

  What to show in cells: "cor", "pval", "both", or "stars". Default:
  "both".

- pval_format:

  Format for p-values: "scientific" (e.g., 1.2e-03) or "decimal" (e.g.,
  0.001). Default: "decimal".

- pval_digits:

  Integer. Decimal places for p-values (decimal format only). Default:
  3.

- colors:

  Color palette. Default: blue-white-red gradient.

- title:

  Plot title.

## Value

A ggplot object.
