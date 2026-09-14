# Fit Sparse CCA Model (Unified Interface)

Wrapper function that fits any sCCA method with consistent interface

## Usage

``` r
MINERVA_scca_fit(
  X = NULL,
  X_1 = NULL,
  X_2 = NULL,
  method = c("ConvCCA", "RelPMDCCA", "PMDCCA"),
  tau = NULL,
  lambda = 10,
  nIter = 100,
  penalty = "LASSO",
  ...
)
```

## Arguments

- X:

  List of datasets (for multi-dataset) or use X_1, X_2 for two datasets

- X_1:

  First dataset (alternative to X for two-dataset methods)

- X_2:

  Second dataset (alternative to X for two-dataset methods)

- method:

  Method to use: "ConvCCA", "RelPMDCCA", or "PMDCCA"

- tau:

  List of tau values (one per dataset) or tauW_1, tauW_2 for two
  datasets

- lambda:

  Lambda parameter (RelPMDCCA only, default = 10)

- nIter:

  Maximum iterations

- penalty:

  Penalty function: "LASSO" (default) or "SCAD"

- ...:

  Additional arguments passed to the method

## Value

List with canonical vectors W (list of vectors, one per dataset)
