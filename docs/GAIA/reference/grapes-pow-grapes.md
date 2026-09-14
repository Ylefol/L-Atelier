# Matrix Power Operator

Raises a matrix to a given power using eigenvalue decomposition

## Usage

``` r
x %^% n
```

## Arguments

- x:

  A square matrix

- n:

  Power to raise the matrix to (can be negative for inverse powers)

## Value

Matrix x raised to the power n

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute inverse square root of a covariance matrix
Sigma <- cov(matrix(rnorm(100), 10, 10))
Sigma_inv_sqrt <- Sigma %^% (-0.5)

} # }
```
