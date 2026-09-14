# Multiple Imputation by Chained Equations (MICE)

Imputes missing values in a numeric data frame using MICE with
predictive mean matching. Supports log1p-transformation of right-skewed
variables before imputation, with automatic back-transformation in the
output. Returns all completed datasets plus a mean-imputed dataset for
use with
[`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

## Usage

``` r
POSEIDON_mice_impute(
  data,
  n_imputations = 100,
  method = "pmm",
  max_iter = 20,
  log_vars = NULL,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- data:

  Numeric data frame. Should contain only the variables to be
  imputed/used for clustering — do not pass the full study dataset.

- n_imputations:

  Integer. Number of imputed datasets to generate. Oslo (Wick et al.)
  used 100; fewer imputations increase variability. Default = 100.

- method:

  Character. Imputation method passed to
  [`mice::mice()`](https://amices.org/mice/reference/mice.html). "pmm"
  (predictive mean matching) is appropriate for continuous variables.
  Default = "pmm".

- max_iter:

  Integer. Maximum MICE iterations per imputation. Default = 20.

- log_vars:

  Character vector. Columns to log1p-transform before imputation
  (back-transformed in all output datasets). Intended for right-skewed
  lab values (e.g., CRP, lactate, bilirubin). Default = NULL.

- seed:

  Integer. Random seed for reproducibility. Default = 42.

- verbose:

  Logical. Print progress. Default = TRUE.

## Value

A list containing:

- datasets:

  List of `n_imputations` completed data frames on original scale

- mean_dataset:

  Data frame with column means across all imputations (used for centroid
  computation)

- log_vars:

  Character vector of variables that were log1p-transformed

- n_imputations:

  Number of imputations generated

- n_obs:

  Number of observations

- n_vars:

  Number of variables

## Details

Requires the `mice` package (install with `install.packages("mice")`).
