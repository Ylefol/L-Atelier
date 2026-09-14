# Generate Discriminative Analysis Report

Generate a comprehensive report from the discriminative variable
analysis pipeline, including variable selection results, model
coefficients, and optionally save tables and plots.

## Usage

``` r
ELEUTHIA_discriminative_report(
  glasso_fit,
  cv_result,
  selected,
  final_model = NULL,
  coef_table = NULL,
  encoded = NULL,
  output_dir = NULL,
  prefix = "discriminative",
  save_tables = TRUE,
  save_plots = TRUE,
  plot_width = 10,
  plot_height = 8,
  verbose = TRUE
)
```

## Arguments

- glasso_fit:

  An artemis_group_lasso object from ARTEMIS_fit_group_lasso().

- cv_result:

  Output from ARTEMIS_select_lambda().

- selected:

  Output from ARTEMIS_extract_selected_variables().

- final_model:

  Optional. An artemis_final_model from ARTEMIS_fit_final_model().

- coef_table:

  Optional. Data.frame from ARTEMIS_extract_coefficients().

- encoded:

  Optional. Output from POSEIDON_encode_for_regression() for additional
  context.

- output_dir:

  Character. Directory to save output files. If NULL (default), only
  prints to console without saving files.

- prefix:

  Character. Prefix for output filenames. Default = "discriminative".

- save_tables:

  Logical. Save tables as CSV files. Default = TRUE.

- save_plots:

  Logical. Save plots as PNG files. Default = TRUE.

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

  List of data frames exported

- files_saved:

  Character vector of files saved

## Details

This function consolidates results from the discriminative analysis
pipeline:

- Group LASSO variable selection results

- Cross-validation lambda selection

- Selected variables and their effects

- Final model coefficients (odds ratios or coefficients)

At minimum, glasso_fit, cv_result, and selected are required. Final
model and coefficients are optional but recommended for complete
reporting.

## Examples

``` r
if (FALSE) { # \dontrun{
# Console report only
ELEUTHIA_discriminative_report(glasso_fit, cv_result, selected)

# Full report with file export
ELEUTHIA_discriminative_report(
  glasso_fit = fit,
  cv_result = cv_result,
  selected = selected,
  final_model = final,
  coef_table = coefs,
  output_dir = "results/discriminative"
)

} # }
```
