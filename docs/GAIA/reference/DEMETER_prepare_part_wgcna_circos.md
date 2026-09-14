# Prepare data for PART-WGCNA circos plot

Combines PART clustering results, WGCNA module assignments, trait
correlations, and DEA log2 fold-change values into a structured object
for circos visualization.

## Usage

``` r
DEMETER_prepare_part_wgcna_circos(
  part_result,
  wgcna_modules,
  wgcna_trait_cor = NULL,
  dea_results = NULL,
  strip_version = TRUE,
  p_thresh = 0.05,
  sig_type = c("padj", "pvalue"),
  verbose = TRUE
)
```

## Arguments

- part_result:

  An `artemis_part` object from `ARTEMIS_part_clustering()`.

- wgcna_modules:

  A `wgcna_modules` object from WGCNA pipeline or
  [`DEMETER_load_wgcna()`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_wgcna.md).

- wgcna_trait_cor:

  A `wgcna_trait_cor` object from
  [`ARTEMIS_wgcna_module_traits()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_wgcna_module_traits.md)
  or
  [`DEMETER_load_wgcna()`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_wgcna.md).
  NULL if no trait correlations available.

- dea_results:

  Named list of DEA results. Each element should be either:

  - A result from
    [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
    (has `$results`)

  - A data.frame with feature_id and log2FoldChange columns Names become
    column labels in the circos (e.g., list(C5RO = dea_c5ro, CS1AN =
    dea_cs1an)). NULL if no DEA results to display.

- strip_version:

  Logical. Strip version suffix from Ensembl gene IDs (e.g.,
  ENSG00000103888.18 -\> ENSG00000103888) for matching between PART and
  WGCNA. Default: TRUE.

- p_thresh:

  Numeric. P-value threshold for trait significance in module_df.
  Default: 0.05.

- sig_type:

  Character. Type of significance to use for filtering: "padj" (adjusted
  p-value, default) or "pvalue" (raw p-value).

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `part_wgcna_circos` containing:

- PART_df:

  Data.frame with genes as rows, columns: cluster, module, and one L2FC
  column per DEA comparison

- module_df:

  Data.frame with modules as rows, columns: module, size, trait
  correlations, and significance flags

- association_matrix:

  Integer matrix (clusters x modules) of gene counts

- cluster_colors:

  Named character vector of cluster colors

- module_colors:

  Named character vector of module colors

- module_gene:

  Data.frame with gene and module columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare circos data from GAIA objects
circos_data <- DEMETER_prepare_part_wgcna_circos(
  part_result = part_res,
  wgcna_modules = wgcna$modules,
  wgcna_trait_cor = wgcna$trait_cor,
  dea_results = list(C5RO = dea_c5ro, CS1AN = dea_cs1an)
)

# Filter and plot
circos_data <- DEMETER_filter_circos_modules(circos_data, c("blue", "brown", "turquoise"))
AETHER_plot_part_wgcna_circos(circos_data, output_file = "circos.pdf")

} # }
```
