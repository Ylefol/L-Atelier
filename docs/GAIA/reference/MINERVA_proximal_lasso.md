# Apply Proximal Operator for LASSO (ADMM Framework)

Element-wise proximal operator for LASSO in ADMM. Used in RelPMDCCA
algorithm.

## Usage

``` r
MINERVA_proximal_lasso(omega, mu, c, tau, old)
```

## Arguments

- omega:

  Proximal gradient term

- mu:

  Step size parameter

- c:

  Cross-product term (gradient of correlation)

- tau:

  Sparsity tuning parameter

- old:

  Previous iteration value (fallback if update fails)

## Value

Updated value
