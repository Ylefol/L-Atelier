# Apply Proximal Operator for SCAD (ADMM Framework)

Element-wise proximal operator for SCAD in ADMM. Used in RelPMDCCA
algorithm.

## Usage

``` r
MINERVA_proximal_scad(omega, mu, c, tau, a = 3.7, old)
```

## Arguments

- omega:

  Proximal gradient term

- mu:

  Step size parameter

- c:

  Cross-product term

- tau:

  Sparsity tuning parameter

- a:

  SCAD shape parameter (default = 3.7)

- old:

  Previous iteration value (fallback)

## Value

Updated value
