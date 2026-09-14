# Multi-Dataset Relaxed Penalized Matrix Decomposition CCA

Performs sparse CCA on multiple (\>2) datasets simultaneously using ADMM
framework. Finds canonical vectors that maximize mutual correlation
across all datasets while allowing for dependent features.

## Usage

``` r
HEPHAESTUS_multi_relPMDCCA(
  X,
  lambda,
  tau,
  a = 3.7,
  sd.mu = sqrt(2),
  nIter = 1000,
  penalty = "LASSO",
  element_wise = TRUE,
  tau_EN = NULL,
  R = 1
)
```

## Arguments

- X:

  List of M datasets, where each element is an (n x p_i) matrix. All
  datasets must have the same number of samples (rows)

- lambda:

  Relaxation parameter controlling step size (typically fixed at 10)

- tau:

  List of M sparsity parameters, one per dataset (e.g., list(0.8, 0.8,
  0.7))

- a:

  SCAD shape parameter (default = 3.7)

- sd.mu:

  Standard deviation for step size initialization (default = sqrt(2))

- nIter:

  Maximum number of iterations (default = 1000). Stops early if
  converged.

- penalty:

  Penalty function: "LASSO" (default), "SCAD", or "ELASTIC-NET"

- element_wise:

  Use element-wise proximal updates (default = TRUE)

- tau_EN:

  Secondary tuning parameter for Elastic-Net penalty

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

# Run multi.relPMDCCA
result <- HEPHAESTUS_multi_relPMDCCA(
  X = list(X1, X2, X3),
  lambda = 10,
  tau = list(0.8, 0.8, 0.8),
  nIter = 100
)

# Check mutual correlations
scores <- lapply(1:3, function(i) result$X[[i]] %*% result$W[[i]])
cor(scores[[1]], scores[[2]])

} # }
```
