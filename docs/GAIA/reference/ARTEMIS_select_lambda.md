# Select Optimal Lambda via Cross-Validation

Performs cross-validation to select the optimal regularization parameter
(lambda) for a Group LASSO model. Returns both lambda.min (minimum CV
error) and lambda.1se (most parsimonious within 1 SE of minimum).

## Usage

``` r
ARTEMIS_select_lambda(
  glasso_fit,
  X,
  y,
  nfolds = 10,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- glasso_fit:

  An artemis_group_lasso object from
  [`ARTEMIS_fit_group_lasso()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_fit_group_lasso.md).

- X:

  The predictor matrix used to fit the model.

- y:

  The response vector used to fit the model.

- nfolds:

  Integer. Number of cross-validation folds (default = 10).

- seed:

  Integer. Random seed for reproducibility (default = NULL).

- verbose:

  Logical. Print CV summary (default = TRUE).

## Value

A list containing:

- lambda.min:

  Lambda with minimum CV error

- lambda.1se:

  Largest lambda within 1 SE of minimum (more parsimonious)

- cve:

  Cross-validation error for each lambda

- cvse:

  Standard error of CV error

- cv_fit:

  The cv.grpreg object for further inspection

- n_selected_min:

  Number of groups selected at lambda.min

- n_selected_1se:

  Number of groups selected at lambda.1se

## Details

Two lambda choices are provided:

- `lambda.min`: Minimizes prediction error. May include more variables
  than necessary.

- `lambda.1se`: The largest lambda (most regularized) whose error is
  within 1 standard error of the minimum. More parsimonious, recommended
  for variable selection.

For variable selection purposes, `lambda.1se` is typically preferred as
it provides a simpler model with similar predictive performance.

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting Group LASSO
fit <- ARTEMIS_fit_group_lasso(X, y, groups, family = "binomial")

# Select optimal lambda
cv_result <- ARTEMIS_select_lambda(fit, X, y, nfolds = 10)

# Use lambda.1se for variable selection (recommended)
selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)

} # }
```
