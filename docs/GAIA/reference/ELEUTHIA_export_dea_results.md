# Export DEA + PART + Enrichment Results

Comprehensive export function for the standard differential expression
analysis pipeline: DEA → gene selection → PART clustering → enrichment.

## Usage

``` r
ELEUTHIA_export_dea_results(
  output_dir,
  dea_result = NULL,
  part_result = NULL,
  gene_selection = NULL,
  enrichment = NULL,
  l2fc_thresh = 1,
  p_thresh = 0.05,
  sample_info = NULL,
  sample_col = NULL,
  group_col = "group",
  save_plots = TRUE,
  plot_format = "png",
  prefix = "dea",
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- output_dir:

  Character. Directory to save results. Created if needed.

- dea_result:

  DEA result from
  [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md).
  NULL to skip DEA export. Default: NULL.

- part_result:

  An `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).
  NULL to skip PART export. Default: NULL.

- gene_selection:

  Character vector of selected gene IDs (e.g., from
  [`ARTEMIS_select_de_genes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_select_de_genes.md)).
  NULL to skip. Default: NULL.

- enrichment:

  Enrichment result from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md).
  NULL to skip. Delegates to
  [`ELEUTHIA_export_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_enrichment.md).
  Default: NULL.

- l2fc_thresh:

  Numeric. log2 fold change threshold for significant genes export.
  Default: 1.0.

- p_thresh:

  Numeric. Adjusted p-value threshold for significant genes export.
  Default: 0.05.

- sample_info:

  Data.frame with sample metadata (for heatmap annotations). Must have
  rownames matching sample columns. Default: NULL.

- group_col:

  Character. Group column in sample_info. Default: "group".

- save_plots:

  Logical. Generate and save plots. Default: TRUE.

- plot_format:

  Character. "png", "pdf", or "both". Default: "png".

- prefix:

  Character. Prefix for output filenames. Default: "dea".

- save_rds:

  Logical. Save full R objects as RDS. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Depending on which arguments are provided, exports:

**DEA results** (when `dea_result` provided):

- All results CSV (full DESeq2 output)

- Significant results CSV (filtered by thresholds)

- Metadata text file (comparison, thresholds, counts)

**PART results** (when `part_result` provided):

- Cluster assignments CSV (gene → cluster)

- Cluster summary CSV (cluster, n_genes, color)

- Parameters text file

- Heatmap plot (if sample_info provided)

**Gene selection** (when `gene_selection` provided):

- Selected genes CSV

**Enrichment** (when `enrichment` provided):

- Delegates to
  [`ELEUTHIA_export_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_enrichment.md)
  in subdirectory

## Examples

``` r
if (FALSE) { # \dontrun{
# Full pipeline export
ELEUTHIA_export_dea_results(
  output_dir = "results/dea_analysis/",
  dea_result = dea_res,
  part_result = part_res,
  gene_selection = sig_genes,
  enrichment = enrich_res,
  sample_info = sample_sheet
)

# DEA only
ELEUTHIA_export_dea_results(
  output_dir = "results/",
  dea_result = dea_res,
  prefix = "wt_vs_ko"
)

} # }
```
