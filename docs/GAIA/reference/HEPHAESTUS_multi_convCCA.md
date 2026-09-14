# Multi-Dataset Convex Sparse Canonical Correlation Analysis

Performs sparse CCA on multiple (\>2) datasets simultaneously, finding
canonical vectors that maximize mutual correlation across all datasets

## Usage

``` r
HEPHAESTUS_multi_convCCA(
  X,
  tau,
  nIter = 1000,
  penalty = "LASSO",
  initW = NULL,
  R = 1
)
```

## Arguments

- X:

  List of M datasets, where each element is an (n x p_i) matrix. All
  datasets must have the same number of samples (rows)

- tau:

  List of M sparsity parameters, one per dataset (e.g., list(0.3, 0.5,
  0.2))

- nIter:

  Maximum number of iterations (default = 1000). Unlike convCCA,
  multi.convCCA requires a fixed nIter (no auto-convergence)

- penalty:

  Penalty function: "LASSO" (default) or "SCAD"

- initW:

  Optional list of M initial canonical vectors

- R:

  Number of canonical tuples to compute (default = 1)

## Value

List containing:

- W:

  List of M canonical vectors, one per dataset

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate 3 correlated datasets
set.seed(123)
n <- 50
u <- matrix(rnorm(n), ncol = 1)
X1 <- u %*% t(c(rep(1, 20), rep(0, 80))) + matrix(rnorm(n * 100), n, 100)
X2 <- u %*% t(c(rep(1, 15), rep(0, 85))) + matrix(rnorm(n * 100), n, 100)
X3 <- u %*% t(c(rep(1, 10), rep(0, 40))) + matrix(rnorm(n * 50), n, 50)

# Run multi.convCCA
result <- HEPHAESTUS_multi_convCCA(
  X = list(X1, X2, X3),
  tau = list(0.3, 0.3, 0.3),
  nIter = 100
)

# Check sparsity for each dataset
sapply(result$W, function(w) sum(abs(w) > 1e-6))

} # }
```
