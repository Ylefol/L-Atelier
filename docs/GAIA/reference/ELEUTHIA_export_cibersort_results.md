# Export CIBERSORT deconvolution results

Saves cell type proportions, quality diagnostics, metadata, and optional
plots from an `artemis_cibersort` object.

## Usage

``` r
ELEUTHIA_export_cibersort_results(
  result,
  output_dir,
  prefix = "cibersort",
  group_by = NULL,
  save_plots = TRUE,
  plot_format = "png",
  top_n = NULL,
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- result:

  An `artemis_cibersort` object from
  [`ARTEMIS_cibersort()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_cibersort.md).

- output_dir:

  Character. Directory to save results. Created if needed.

- prefix:

  Character. Prefix for output filenames. Default: "cibersort".

- group_by:

  Optional character/factor vector of condition labels per sample. Used
  for grouped boxplot. Supports named vectors (names = groups, values =
  colors). Default: NULL.

- save_plots:

  Logical. Save diagnostic plots. Default: TRUE.

- plot_format:

  Character. Plot format: "png", "pdf", or "both". Default: "png".

- top_n:

  Integer. Top N cell types for boxplot (by mean proportion). Default:
  NULL (show all non-zero).

- save_rds:

  Logical. Save full R object as RDS for reloading. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Exports the following files:

- Proportions matrix (samples x cell types) as CSV

- Diagnostics (correlation, RMSE, p-values) as CSV

- Metadata text file (parameters, quality summary)

- Plots: stacked bar, heatmap, boxplot (if group_by provided)

- RDS file of the full object (optional)

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_export_cibersort_results(cib_result, "results/cibersort/")

# With grouped boxplot
ELEUTHIA_export_cibersort_results(
  cib_result, "results/cibersort/",
  group_by = sample_conditions
)

} # }
```
