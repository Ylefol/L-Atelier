# Forest Plot of Regression Coefficients

Creates a forest plot displaying coefficients (or odds ratios) with
confidence intervals. Standard visualization for presenting regression
results in clinical research.

## Usage

``` r
AETHER_plot_coefficients(
  coef_table,
  null_line = 1,
  log_scale = TRUE,
  sort_by = "estimate",
  decreasing = TRUE,
  colors = NULL,
  point_size = 3,
  line_width = 0.8,
  text_size = 11,
  title = "Forest Plot",
  xlab = NULL,
  show_values = TRUE,
  sig_level = 0.05
)
```

## Arguments

- coef_table:

  A data.frame from
  [`ARTEMIS_extract_coefficients()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_extract_coefficients.md),
  or any data.frame with columns: variable, estimate, ci_lower,
  ci_upper.

- null_line:

  Numeric. Reference line value (default = 1 for odds ratios, use 0 for
  raw coefficients). Set to NULL to omit.

- log_scale:

  Logical. Use log scale for x-axis (default = TRUE for odds ratios,
  FALSE for coefficients).

- sort_by:

  Character. How to sort variables:

  - "estimate": By effect size (default)

  - "significance": By p-value

  - "name": Alphabetically

  - "none": Keep original order

- decreasing:

  Logical. Sort direction (default = TRUE, largest first).

- colors:

  Named vector or character. Colors for significant/non-significant:

  - Default: c(significant = "#E41A1C", nonsignificant = "#377EB8")

  - Or single color for all points

- point_size:

  Numeric. Size of point estimates (default = 3).

- line_width:

  Numeric. Width of confidence interval lines (default = 0.8).

- text_size:

  Numeric. Base text size (default = 11).

- title:

  Character. Plot title (default = "Forest Plot").

- xlab:

  Character. X-axis label (default based on log_scale).

- show_values:

  Logical. Show estimate values on right side (default = TRUE).

- sig_level:

  Numeric. Significance threshold for coloring (default = 0.05).

## Value

A ggplot object that can be further customized or saved.

## Details

Forest plots are the standard way to present regression coefficients in
medical/clinical research. Key features:

- Each row = one variable

- Point = effect estimate (odds ratio or coefficient)

- Horizontal line = confidence interval

- Vertical reference line at null value (1 for OR, 0 for coef)

Interpretation for odds ratios:

- OR \> 1 and CI doesn't cross 1: Significant positive association

- OR \< 1 and CI doesn't cross 1: Significant negative association

- CI crosses 1: Not statistically significant

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic forest plot of odds ratios
coef_table <- ARTEMIS_extract_coefficients(final_model, verbose = FALSE)
AETHER_plot_coefficients(coef_table)

# For raw coefficients
coef_table <- ARTEMIS_extract_coefficients(final_model, format = "coefficient",
                                            verbose = FALSE)
AETHER_plot_coefficients(coef_table, null_line = 0, log_scale = FALSE)

# Sort by p-value
AETHER_plot_coefficients(coef_table, sort_by = "significance")

} # }
```
