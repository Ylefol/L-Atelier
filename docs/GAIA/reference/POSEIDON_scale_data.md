# Scale and Center Data Matrix

Centers and scales columns of a data matrix. Wrapper around base::scale
with consistent handling.

## Usage

``` r
POSEIDON_scale_data(X, center = TRUE, scale = TRUE)
```

## Arguments

- X:

  Data matrix (samples x features)

- center:

  Logical, should columns be centered? (default = TRUE)

- scale:

  Logical, should columns be scaled to unit variance? (default = TRUE)

## Value

Scaled matrix with centering and scaling attributes
