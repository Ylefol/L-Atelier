# Plot hub genes across modules

Bar plot showing the number of hub genes per module, with top hub gene
labeled.

## Usage

``` r
AETHER_plot_hub_summary(
  hubs,
  top_n = NULL,
  show_top_gene = TRUE,
  module_colors = NULL,
  title = "Hub genes per module"
)
```

## Arguments

- hubs:

  A wgcna_hubs object from ARTEMIS_wgcna_hub_genes().

- top_n:

  Number of modules to show. Default: NULL (all).

- show_top_gene:

  Logical. Label bars with top hub gene name. Default: TRUE.

- module_colors:

  Named character vector of module colors (module_name -\> hex). If
  NULL, generates colors automatically.

- title:

  Plot title.

## Value

A ggplot object.
