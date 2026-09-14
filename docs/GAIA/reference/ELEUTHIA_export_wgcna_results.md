# Export WGCNA results to files

Exports various WGCNA analysis results to CSV files, generates plots,
and optionally saves R objects for later use.

## Usage

``` r
ELEUTHIA_export_wgcna_results(
  output_dir,
  modules = NULL,
  trait_cor = NULL,
  gene_sig = NULL,
  hubs = NULL,
  enrichment = NULL,
  power_result = NULL,
  cluster_result = NULL,
  wgcna_data = NULL,
  save_plots = TRUE,
  plot_format = "png",
  prefix = "wgcna",
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- output_dir:

  Output directory path.

- modules:

  Optional. wgcna_modules object to export.

- trait_cor:

  Optional. wgcna_trait_cor object to export.

- gene_sig:

  Optional. wgcna_gene_sig object to export.

- hubs:

  Optional. wgcna_hubs object to export.

- enrichment:

  Optional. gost_enrichment object from APOLLO_enrich_gost(). If
  provided, exports enrichment results to \_enrichment/ subdirectory.

- power_result:

  Optional. wgcna_power object for power selection plot.

- cluster_result:

  Optional. wgcna_cluster object for sample dendrogram.

- wgcna_data:

  Optional. wgcna_data object (needed for sample dendrogram traits).

- save_plots:

  Logical. Generate and save plots. Default: TRUE.

- plot_format:

  Character. Plot format: "png", "pdf", or "both". Default: "png".

- prefix:

  Prefix for output filenames. Default: "wgcna".

- save_rds:

  Logical. Save R objects as RDS files. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible list of file paths created.

## Details

Creates the following files (depending on objects provided):

- \_gene_modules.csv: Gene-to-module assignments

- \_module_summary.csv: Module sizes

- \_module_eigengenes.csv: Module eigengene values

- \_trait_correlations.csv: Module-trait correlations

- \_trait_pvalues.csv: Correlation p-values

- \_significant_associations.csv: Significant module-trait pairs

- \_gene_significance.csv: Full gene info table

- \_hub_genes.csv: Hub gene list

- \_hub_summary.csv: Hub counts per module

- Module-specific gene lists in \_modules/ subdirectory

- Enrichment results in \_enrichment/ subdirectory (if enrichment
  provided)

- Plots in \_plots/ subdirectory (if save_plots = TRUE)

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_export_wgcna_results(
  output_dir = "results/wgcna",
  modules = modules,
  trait_cor = trait_cor,
  gene_sig = gene_sig,
  hubs = hubs,
  enrichment = gost_results,
  power_result = power_result,
  cluster_result = cluster_result
)

} # }
```
