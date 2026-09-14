# Plot gene significance vs module membership

Creates a scatter plot of gene significance (GS) vs module membership
(MM) for genes in a specific module. Useful for identifying hub genes.

## Usage

``` r
AETHER_plot_gs_vs_mm(
  gene_sig,
  module,
  trait_name = NULL,
  mm_threshold = 0.8,
  gs_threshold = 0.2,
  n_label = 10,
  module_color = NULL,
  title = NULL
)
```

## Arguments

- gene_sig:

  A wgcna_gene_sig object from ARTEMIS_wgcna_gene_significance().

- module:

  Module color to plot.

- trait_name:

  Trait name for gene significance. If NULL, uses first trait.

- mm_threshold:

  Module membership threshold for highlighting. Default: 0.8.

- gs_threshold:

  Gene significance threshold for highlighting. Default: 0.2.

- n_label:

  Number of top genes to label. Default: 10.

- module_color:

  Hex color for hub gene points. If NULL, uses default green.

- title:

  Plot title. If NULL, auto-generated.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)
AETHER_plot_gs_vs_mm(gene_sig, module = "blue", trait_name = "severity")

} # }
```
