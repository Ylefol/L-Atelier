# Generate Stratification Pipeline Report

Generate a comprehensive report from the patient stratification
pipeline, including clustering results, stability assessment, variable
importance, and optionally save tables and plots.

## Usage

``` r
ELEUTHIA_stratification_report(
  cluster_result,
  stability_result = NULL,
  char_result = NULL,
  famd_result = NULL,
  comparison_result = NULL,
  qc_result = NULL,
  output_dir = NULL,
  prefix = "stratification",
  save_tables = TRUE,
  save_plots = TRUE,
  plot_width = 10,
  plot_height = 8,
  verbose = TRUE
)
```

## Arguments

- cluster_result:

  An artemis_cluster object from ARTEMIS_cluster_mixed().

- stability_result:

  Optional. An artemis_stability object from
  ARTEMIS_cluster_stability().

- char_result:

  Optional. An artemis_characterization object from
  ARTEMIS_characterize_clusters().

- famd_result:

  Optional. An artemis_famd object from ARTEMIS_famd().

- comparison_result:

  Optional. An artemis_variable_comparison object from
  ARTEMIS_compare_variable_importance().

- qc_result:

  Optional. A hades_qc object from HADES_qc_mixed_data().

- output_dir:

  Character. Directory to save output files. If NULL (default), only
  prints to console without saving files.

- prefix:

  Character. Prefix for output filenames. Default = "stratification".

- save_tables:

  Logical. Save key tables as CSV files. Default = TRUE (when output_dir
  is provided).

- save_plots:

  Logical. Save plots as PNG files. Default = TRUE (when output_dir is
  provided).

- plot_width:

  Numeric. Width of saved plots in inches. Default = 10.

- plot_height:

  Numeric. Height of saved plots in inches. Default = 8.

- verbose:

  Logical. Print report to console. Default = TRUE.

## Value

Invisibly returns a list containing:

- summary:

  Character vector of the text report

- tables:

  List of data frames that were/would be exported

- files_saved:

  Character vector of files saved (if output_dir provided)

## Details

This function consolidates results from the stratification pipeline into
a single report. It can:

- Print a formatted summary to console

- Save the summary as a text file

- Export key tables (cluster assignments, variable importance, etc.) as
  CSV

- Generate and save standard plots (silhouette, FAMD, radar)

At minimum, cluster_result is required. Other results are optional and
will be included if provided.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate console report only
ELEUTHIA_stratification_report(cluster_result)

# Full report with file export
ELEUTHIA_stratification_report(
  cluster_result = clust,
  stability_result = stab,
  char_result = char,
  famd_result = famd,
  output_dir = "results/stratification"
)

} # }
```
