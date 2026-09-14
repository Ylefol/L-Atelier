# Export decoupleR activity results

Saves activity inference results (single method or multi-method
comparison) to CSV and RDS files. Handles both `decoupler_result` and
`decoupler_comparison` objects.

## Usage

``` r
ELEUTHIA_export_activity_results(
  result,
  output_dir,
  prefix = "activity",
  save_plots = TRUE,
  plot_format = "png",
  top_n = 30,
  save_rds = TRUE,
  save_long = TRUE,
  verbose = TRUE
)
```

## Arguments

- result:

  A `decoupler_result` or `decoupler_comparison` object from
  ARTEMIS_run_decoupler(), ARTEMIS_infer_tf_activity(),
  ARTEMIS_infer_pathway_activity(), or
  ARTEMIS_decoupler_compare_methods().

- output_dir:

  Character. Directory to save results. Created if needed.

- prefix:

  Character. Prefix for output filenames. Default: "activity".

- save_plots:

  Logical. Save diagnostic plots. Default: TRUE.

- plot_format:

  Character. Plot format: "png", "pdf", or "both". Default: "png".

- top_n:

  Integer. Number of top sources to show in heatmap and bar plots.
  Default: 30.

- save_rds:

  Logical. Save full R objects as RDS for exact reloading. Default:
  TRUE.

- save_long:

  Logical. Save full long-format results from decoupleR. Default: TRUE.
  Can produce large files for comparisons.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

For a `decoupler_result` (single method), exports:

- Activity matrix (sources x samples) as CSV

- P-value matrix as CSV (if available)

- Full long-format results as CSV (optional)

- Metadata text file (method, network info, gene overlap)

- Plots: activity heatmap, top activities bar chart

- RDS file of the full object (optional)

For a `decoupler_comparison` (multi-method), exports:

- Per-method activity matrices in subdirectory

- Consensus scores as CSV

- Summary table (per-source stats across methods) as CSV

- Method correlation matrix as CSV

- Method statistics as CSV

- Metadata text file

- Plots: activity heatmap, top activities, method correlation, method
  agreement

- RDS file of the full object (optional)

Saved RDS files can be reloaded with
[`DEMETER_load_activity()`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_activity.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Export single-method result
tf_result <- ARTEMIS_infer_tf_activity(expr_matrix, method = "ulm")
ELEUTHIA_export_activity_results(tf_result, "results/tf_activity/")

# Export multi-method comparison
comparison <- ARTEMIS_decoupler_compare_methods(
  expr_matrix, network, methods = c("ulm", "mlm", "wsum")
)
ELEUTHIA_export_activity_results(comparison, "results/tf_comparison/")

# Reload later
tf_result <- DEMETER_load_activity("results/tf_activity/activity_result.rds")

} # }
```
