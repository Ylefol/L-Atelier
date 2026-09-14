# Scale Predictors for Penalized Regression

Scales predictor matrix for use in penalized regression methods.
Penalized regression (LASSO, Group LASSO, elastic net) requires scaling
because penalties are applied uniformly - without scaling, variables on
larger scales would be penalized less.

## Usage

``` r
POSEIDON_scale_predictors(
  X,
  method = "zscore",
  exclude_binary = TRUE,
  var_mapping = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  Numeric matrix of predictors (samples x features), typically from
  [`POSEIDON_encode_for_regression()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_encode_for_regression.md).

- method:

  Character string. Scaling method:

  - "zscore" (default): Standardize to mean=0, sd=1

  - "minmax": Scale to 0-1 range

  - "robust": Use median and IQR (robust to outliers)

- exclude_binary:

  Logical. If TRUE (default), columns that are binary (only 0 and 1
  values) are not scaled. Dummy-coded categoricals are already on a
  meaningful 0/1 scale.

- var_mapping:

  Optional. The var_mapping from
  [`POSEIDON_encode_for_regression()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_encode_for_regression.md).
  If provided, can use variable type information to guide scaling
  decisions.

- verbose:

  Logical. Print scaling summary (default = TRUE).

## Value

A list containing:

- X_scaled:

  Scaled predictor matrix

- scaling_params:

  List of parameters for each scaled column:

  - center: Value subtracted (mean, min, or median)

  - scale: Value divided by (sd, range, or IQR)

  These can be used to transform new data consistently.

- scaled_cols:

  Names of columns that were scaled

- unscaled_cols:

  Names of columns left unscaled (binary)

- method:

  Scaling method used

## Details

Why scaling matters for penalized regression:

- LASSO/Group LASSO apply the same penalty to all coefficients

- Without scaling, a variable measured in thousands would have a tiny
  coefficient (to compensate) and thus be penalized less

- Scaling puts all variables on equal footing for fair penalization

Binary columns (dummy variables) are typically excluded from scaling
because:

- They're already on a 0/1 scale with clear interpretation

- Scaling would change their interpretation (no longer
  "presence/absence")

- For Group LASSO, the group penalty handles their contribution

The scaling parameters are returned so you can apply the same
transformation to new/test data, ensuring consistency.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage after encoding
encoded <- POSEIDON_encode_for_regression(data, target_var = "outcome")
scaled <- POSEIDON_scale_predictors(encoded$X)

# Use robust scaling for data with outliers
scaled <- POSEIDON_scale_predictors(encoded$X, method = "robust")

# Scale everything including binary (not typical)
scaled <- POSEIDON_scale_predictors(encoded$X, exclude_binary = FALSE)

} # }
```
