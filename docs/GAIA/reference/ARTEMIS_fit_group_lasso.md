# Fit Group LASSO for Variable Selection

Fits a Group LASSO model for variable selection in mixed data. Group
LASSO treats all dummy-coded columns from a single categorical variable
as a group, selecting entire variables rather than individual levels.

## Usage

``` r
ARTEMIS_fit_group_lasso(
  X,
  y,
  groups,
  family = "binomial",
  penalty = "grLasso",
  nlambda = 100,
  lambda = NULL,
  alpha = 1,
  group_names = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  Numeric matrix of predictors (samples x features), typically from
  [`POSEIDON_scale_predictors()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_scale_predictors.md).

- y:

  Response vector. For binary outcomes, should be 0/1 or a factor. For
  continuous outcomes, should be numeric.

- groups:

  Integer vector mapping each column in X to its variable group. Columns
  with the same group number belong to the same original variable.
  Typically from `encoded$var_mapping$groups`.

- family:

  Character string. Response type:

  - "binomial": Binary outcome (logistic regression)

  - "multinomial": Multi-class outcome (3+ categories)

  - "gaussian": Continuous outcome (linear regression)

  - "poisson": Count outcome (Poisson regression)

- penalty:

  Character string. Penalty type:

  - "grLasso" (default): Group LASSO - L1 penalty on group norms

  - "grMCP": Group MCP - non-convex, less biased estimates

  - "grSCAD": Group SCAD - non-convex alternative

- nlambda:

  Integer. Number of lambda values in the regularization path (default =
  100).

- lambda:

  Optional numeric vector. Specific lambda values to use. If NULL
  (default), automatically generated.

- alpha:

  Elastic net mixing parameter (default = 1 for pure Group LASSO).
  Values \< 1 add L2 penalty for stability with correlated predictors.

- group_names:

  Optional character vector. Names for each group (for nicer output). If
  NULL, uses "Group1", "Group2", etc.

- verbose:

  Logical. Print fitting summary (default = TRUE).

## Value

An object of class "artemis_group_lasso" containing:

- model:

  The fitted grpreg model object

- groups:

  Group assignments used

- group_names:

  Names for each group

- family:

  Response family used

- penalty:

  Penalty type used

- n_groups:

  Number of variable groups

- n_samples:

  Number of samples

- n_features:

  Number of features (columns in X)

- lambda:

  Vector of lambda values in path

## Details

Group LASSO is ideal for variable selection with categorical predictors
because:

- It treats all dummy columns from one variable as a unit

- A variable is either fully in or fully out of the model

- This avoids the awkward situation where only some levels are selected

- Results are interpretable: "Variable X is predictive of Y"

The regularization path fits models across a range of lambda values,
from high (few/no variables selected) to low (many variables selected).
Use
[`ARTEMIS_select_lambda()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_select_lambda.md)
to choose the optimal lambda via cross-validation.

## Examples

``` r
if (FALSE) { # \dontrun{
# After encoding and scaling
encoded <- POSEIDON_encode_for_regression(data, target_var = "outcome")
scaled <- POSEIDON_scale_predictors(encoded$X)

# Fit Group LASSO for binary outcome
fit <- ARTEMIS_fit_group_lasso(
  X = scaled$X_scaled,
  y = encoded$y,
  groups = encoded$var_mapping$groups,
  family = "binomial"
)

# With group names for better output
fit <- ARTEMIS_fit_group_lasso(
  X = scaled$X_scaled,
  y = encoded$y,
  groups = encoded$var_mapping$groups,
  group_names = encoded$var_mapping$group_names,
  family = "binomial"
)

} # }
```
