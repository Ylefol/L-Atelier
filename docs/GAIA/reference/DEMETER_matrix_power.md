# Matrix Power Function (Function Form)

Raises a matrix to a given power using eigenvalue decomposition.
Function form of the %^% operator for explicit module referencing.

## Usage

``` r
DEMETER_matrix_power(x, n)
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
Sigma_inv_sqrt <- DEMETER_matrix_power(Sigma, -0.5)

} # }
```
