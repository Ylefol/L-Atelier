# K-Fold Cross-Validation for Sparse CCA

Performs k-fold cross-validation to select optimal tau parameters

## Usage

``` r
MINERVA_cv_scca(
  X = NULL,
  X_1 = NULL,
  X_2 = NULL,
  method = c("ConvCCA", "RelPMDCCA", "PMDCCA"),
  tau_grid,
  lambda = 10,
  k = 5,
  nIter = 100,
  penalty = "LASSO",
  metric = "correlation",
  seed = 123,
  verbose = TRUE
)
```

## Arguments

- X:

  List of datasets (for multi-dataset) or use X_1, X_2 for two datasets

- X_1:

  First dataset (alternative for two-dataset methods)

- X_2:

  Second dataset (alternative for two-dataset methods)

- method:

  Method to use: "ConvCCA", "RelPMDCCA", or "PMDCCA"

- tau_grid:

  Data frame or list of tau combinations to test. For multi-dataset:
  data.frame(tau1 = ..., tau2 = ..., tau3 = ...) For two-dataset:
  data.frame(tau1 = ..., tau2 = ...)

- lambda:

  Lambda parameter (RelPMDCCA only, default = 10)

- k:

  Number of CV folds (default = 5)

- nIter:

  Maximum iterations for each fit

- penalty:

  Penalty function: "LASSO" (default) or "SCAD"

- metric:

  Metric to use: "correlation" (default, optimizes for correlation) or
  "both" (computes both correlation and sparsity, optimizes for
  correlation)

- seed:

  Random seed for reproducibility

- verbose:

  Print progress (default = TRUE)

## Value

List containing:

- best_tau:

  Optimal tau values (list)

- best_score:

  Best CV score

- cv_results:

  Data frame with all CV results

- best_model:

  Model fit on full data with best tau
