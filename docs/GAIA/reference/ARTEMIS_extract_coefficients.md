# Extract Coefficients with Confidence Intervals

Extracts model coefficients in an interpretable format with confidence
intervals and p-values. For logistic regression, converts coefficients
to odds ratios. Returns a publication-ready table.

## Usage

``` r
ARTEMIS_extract_coefficients(
  final_model,
  format = NULL,
  conf_level = 0.95,
  digits = 3,
  include_intercept = FALSE,
  verbose = TRUE
)
```

## Arguments

- final_model:

  An artemis_final_model object from
  [`ARTEMIS_fit_final_model()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_fit_final_model.md),
  or a standard glm object.

- format:

  Character string. Output format:

  - "odds_ratio" (default for binomial): Exponentiated coefficients

  - "coefficient": Raw coefficients (log-odds for binomial)

- conf_level:

  Numeric. Confidence level for intervals (default = 0.95).

- digits:

  Integer. Number of decimal places for rounding (default = 3).

- include_intercept:

  Logical. Include intercept in output (default = FALSE).

- verbose:

  Logical. Print formatted table (default = TRUE).

## Value

A data.frame with columns:

- variable:

  Variable name (and level for categorical)

- estimate:

  Coefficient or odds ratio

- ci_lower:

  Lower confidence bound

- ci_upper:

  Upper confidence bound

- std_error:

  Standard error

- z_value:

  Z-statistic (or t for gaussian)

- p_value:

  P-value

- significance:

  Significance stars (\* p\<0.05, \*\* p\<0.01, \*\*\* p\<0.001)

## Details

For logistic regression (family = binomial):

- Odds ratio \> 1: Higher values of predictor associated with higher
  probability of outcome = 1

- Odds ratio \< 1: Higher values associated with lower probability

- Odds ratio = 1: No association

For categorical variables, each level (except reference) gets its own
row. The odds ratio compares that level to the reference level.

Confidence intervals are Wald-based (coefficient ± z \* SE).

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting final model
final <- ARTEMIS_fit_final_model(data, "outcome", selected_vars, "binomial")

# Get odds ratios with 95% CI
coef_table <- ARTEMIS_extract_coefficients(final)

# Get raw coefficients instead
coef_table <- ARTEMIS_extract_coefficients(final, format = "coefficient")

# Use 99% confidence intervals
coef_table <- ARTEMIS_extract_coefficients(final, conf_level = 0.99)

} # }
```
