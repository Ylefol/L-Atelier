# Compute Pseudoinverse of a Matrix

Computes the Moore-Penrose pseudoinverse using SVD. Handles
rank-deficient matrices by truncating near-zero singular values.

## Usage

``` r
DEMETER_gettingInverse(X, tolerance = 0.001)
```

## Arguments

- X:

  Input matrix

- tolerance:

  Singular values below this threshold are treated as zero (default:
  1e-3)

## Value

Pseudoinverse of X
