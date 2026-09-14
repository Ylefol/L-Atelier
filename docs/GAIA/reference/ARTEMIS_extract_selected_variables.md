# Extract Selected Variables from Group LASSO

Identifies which variable groups have non-zero coefficients at a given
lambda value. Returns the names of selected variables.

## Usage

``` r
ARTEMIS_extract_selected_variables(glasso_fit, lambda, verbose = TRUE)
```

## Arguments

- glasso_fit:

  An artemis_group_lasso object from
  [`ARTEMIS_fit_group_lasso()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_fit_group_lasso.md).

- lambda:

  Numeric. The lambda value at which to extract selected variables.
  Typically `cv_result$lambda.1se` or `cv_result$lambda.min`.

- verbose:

  Logical. Print selection summary (default = TRUE).

## Value

A list containing:

- selected_groups:

  Integer vector of selected group numbers

- selected_names:

  Character vector of selected variable names

- n_selected:

  Number of variables selected

- n_total:

  Total number of variables

- lambda:

  Lambda value used

- coefficients:

  Named vector of non-zero coefficients (all columns)

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting and CV
fit <- ARTEMIS_fit_group_lasso(X, y, groups, family = "binomial")
cv_result <- ARTEMIS_select_lambda(fit, X, y)

# Extract selected variables at lambda.1se
selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)
print(selected$selected_names)

} # }
```
