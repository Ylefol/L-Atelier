# Extract Coefficients from Multinomial Model (Internal)

Internal helper function to extract and format coefficients from a
nnet::multinom model. Returns a table with one row per variable-class
combination.

## Usage

``` r
.extract_multinomial_coefficients(
  model,
  format,
  conf_level,
  digits,
  include_intercept,
  verbose
)
```

## Arguments

- model:

  A multinom model object

- format:

  "odds_ratio" or "coefficient"

- conf_level:

  Confidence level

- digits:

  Rounding digits

- include_intercept:

  Include intercepts?

- verbose:

  Print output?

## Value

data.frame of coefficients
