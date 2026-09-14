# Update Canonical Vector for Two-Dataset ConvCCA

Updates one canonical vector while keeping the other fixed, using
soft-thresholding with LASSO or SCAD penalty

## Usage

``` r
HEPHAESTUS_updateW_convCCA(K, v, tau, penalty = "LASSO")
```

## Arguments

- K:

  Correlation matrix (transformed covariance structure)

- v:

  Current canonical vector for the other dataset

- tau:

  Sparsity tuning parameter (higher = more sparse)

- penalty:

  Penalty function to use: "LASSO" (default) or "SCAD"

## Value

Updated canonical vector (normalized)
