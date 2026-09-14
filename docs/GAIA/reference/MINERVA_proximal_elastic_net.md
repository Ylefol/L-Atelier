# Apply Proximal Operator for Elastic-Net (ADMM Framework)

Element-wise proximal operator for Elastic-Net penalty. Combines L1
(LASSO) and L2 (Ridge) penalties.

## Usage

``` r
MINERVA_proximal_elastic_net(omega, mu, c, tau, tau_EN, old)
```

## Arguments

- omega:

  Proximal gradient term

- mu:

  Step size parameter

- c:

  Cross-product term

- tau:

  L2 penalty parameter (Ridge component)

- tau_EN:

  L1 penalty parameter (LASSO component)

- old:

  Previous iteration value (fallback)

## Value

Updated value
