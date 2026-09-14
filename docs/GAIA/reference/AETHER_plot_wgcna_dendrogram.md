# Plot module dendrogram with colors

Creates a dendrogram of genes colored by module assignment.

## Usage

``` r
AETHER_plot_wgcna_dendrogram(
  modules,
  block = 1,
  main = "Gene dendrogram and module colors",
  ...
)
```

## Arguments

- modules:

  A wgcna_modules object from ARTEMIS_wgcna_detect_modules().

- block:

  Which block's dendrogram to plot (for large datasets split into
  blocks). Default: 1.

- main:

  Plot title. Default: "Gene dendrogram and module colors".

- ...:

  Additional arguments passed to plotDendroAndColors().

## Value

Invisible NULL. Plots to current device.

## Examples

``` r
if (FALSE) { # \dontrun{
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
AETHER_plot_wgcna_dendrogram(modules)

} # }
```
