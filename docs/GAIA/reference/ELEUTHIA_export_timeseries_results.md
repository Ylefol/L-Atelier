# Export Time Series Analysis Results

Comprehensive export function for the time series analysis pipeline.
Accepts PART clustering results, time series DEA results, gene
selection, and enrichment results. Follows the same conventions as
[`ELEUTHIA_export_wgcna_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_wgcna_results.md)
and
[`ELEUTHIA_export_activity_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_activity_results.md).

## Usage

``` r
ELEUTHIA_export_timeseries_results(
  output_dir,
  part_result = NULL,
  ts_de = NULL,
  gene_selection = NULL,
  enrichment = NULL,
  sample_info = NULL,
  sample_col = NULL,
  time_col = "timepoint",
  group_col = "group",
  save_plots = TRUE,
  plot_format = "png",
  prefix = "timeseries",
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- output_dir:

  Character. Directory to save results. Created if needed.

- part_result:

  An `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).
  NULL to skip PART export. Default: NULL.

- ts_de:

  One or a list of `artemis_ts_de` objects from
  [`ARTEMIS_timeseries_conditional()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_conditional.md)
  or
  [`ARTEMIS_timeseries_temporal()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_temporal.md).
  NULL to skip. Default: NULL.

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

- sample_info:

  Data.frame with sample metadata (for trajectory plots). Must have
  rownames matching sample columns. Default: NULL.

- time_col:

  Character. Timepoint column in sample_info. Default: "timepoint".

- group_col:

  Character. Group column in sample_info. Default: "group".

- save_plots:

  Logical. Generate and save plots. Default: TRUE.

- plot_format:

  Character. "png", "pdf", or "both". Default: "png".

- prefix:

  Character. Prefix for output filenames. Default: "timeseries".

- save_rds:

  Logical. Save full R objects as RDS. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Depending on which arguments are provided, exports:

**PART results** (when `part_result` provided):

- Cluster assignments CSV (gene → cluster)

- Cluster map CSV (gene, cluster, color)

- Cluster summary CSV (cluster, size, color)

- Parameters text file

- Plots: heatmap, cluster trajectories (if sample_info provided),
  cluster means (if sample_info provided)

**Time series DEA** (when `ts_de` provided):

- Per-comparison full result tables in subdirectory

- Summary CSV (timepoints, n_up, n_down)

- Plots: DEA summary bar chart

**Gene selection** (when `gene_selection` provided):

- Selected genes CSV

- Gene summary CSV (per-comparison breakdown, if available)

**Enrichment** (when `enrichment` provided):

- Delegates to
  [`ELEUTHIA_export_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_enrichment.md)
  in subdirectory

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_export_timeseries_results(
  output_dir = "results/timeseries/",
  part_result = part_res,
  ts_de = list(cond_de, temp_de),
  gene_selection = selected_genes,
  sample_info = sample_info
)

} # }
```
