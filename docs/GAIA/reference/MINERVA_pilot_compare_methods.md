# Pilot Comparison of sCCA Methods

Quick comparison of all sCCA methods with default tau values to
determine which method works best for the data

## Usage

``` r
MINERVA_pilot_compare_methods(
  X = NULL,
  X_1 = NULL,
  X_2 = NULL,
  test_proportion = 0.2,
  default_tau_convCCA = NULL,
  default_tau_relPMDCCA = NULL,
  lambda = 10,
  nIter = 100,
  seed = 123,
  use_composite_score = TRUE,
  ideal_sparsity_min = 0.05,
  ideal_sparsity_max = 0.3,
  correlation_weight = 0.7,
  verbose = TRUE
)
```

## Arguments

- X:

  List of datasets (for multi-dataset) or use X_1, X_2

- X_1:

  First dataset (for two-dataset)

- X_2:

  Second dataset (for two-dataset)

- test_proportion:

  Proportion of data to hold out for testing (default = 0.2)

- default_tau_convCCA:

  Default tau for ConvCCA (default = 0.3 for all datasets)

- default_tau_relPMDCCA:

  Default tau for RelPMDCCA (default = 0.7 for all datasets)

- lambda:

  Lambda for RelPMDCCA (default = 10)

- nIter:

  Iterations for each method

- seed:

  Random seed

- use_composite_score:

  Use composite score combining correlation and sparsity quality
  (default = TRUE)

- ideal_sparsity_min:

  Minimum ideal sparsity (default = 0.05, or 5%)

- ideal_sparsity_max:

  Maximum ideal sparsity (default = 0.30, or 30%)

- correlation_weight:

  Weight for correlation in composite score (default = 0.7)

- verbose:

  Print results (default = TRUE)

## Value

List containing:

- best_method:

  Name of best-performing method

- results:

  Data frame with results for all methods

- models:

  List of fitted models
