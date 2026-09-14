# Filter Features by Variance

Removes features with zero or near-zero variance. Essential
preprocessing step before CCA to avoid numerical issues.

## Usage

``` r
POSEIDON_filter_zero_variance(X, variance_threshold = 0, verbose = TRUE)
```

## Arguments

- X:

  Data matrix or list of data matrices (samples x features)

- variance_threshold:

  Minimum variance to retain feature (default = 0)

- verbose:

  Print filtering summary (default = TRUE)

## Value

Filtered data matrix or list of matrices
