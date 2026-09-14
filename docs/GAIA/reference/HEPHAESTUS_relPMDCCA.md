# Relaxed Penalized Matrix Decomposition CCA for Two Datasets

Performs sparse CCA using ADMM framework with relaxed independence
assumption. Allows for dependent features within datasets. Based on Suo
et al. (2017).

## Usage

``` r
HEPHAESTUS_relPMDCCA(
  X_1,
  X_2,
  lambda,
  tauW_1,
  tauW_2,
  initW_1 = NULL,
  initW_2 = NULL,
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

- X_1:

  First dataset matrix (n x p1), samples in rows, features in columns

- X_2:

  Second dataset matrix (n x p2), same samples as X_1

- lambda:

  Relaxation parameter controlling step size (typically fixed at 10)

- tauW_1:

  Sparsity parameter for X_1 (higher = more sparse, typically 0.5-0.9)

- tauW_2:

  Sparsity parameter for X_2

- initW_1:

  Optional initial canonical vector for X_1 (default: rep(0.1, p1))

- initW_2:

  Optional initial canonical vector for X_2 (default: rep(0.1, p2))

- a:

  SCAD shape parameter (default = 3.7, from Fan & Li 2001)

- sd.mu:

  Standard deviation for step size initialization (default = sqrt(2))

- nIter:

  Maximum iterations (default = 1000). Algorithm stops early if
  converged.

- penalty:

  Penalty function: "LASSO" (default), "SCAD", or "ELASTIC-NET"

- element_wise:

  Use element-wise proximal updates (default = TRUE, recommended)

- tau_EN:

  Secondary tuning parameter for Elastic-Net penalty

- R:

  Number of canonical pairs to compute (default = 1)

## Value

List containing:

- W_1:

  Canonical vector for X_1 (length p1)

- W_2:

  Canonical vector for X_2 (length p2)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate two correlated datasets
set.seed(123)
n <- 50; p1 <- 100; p2 <- 80
u <- matrix(rnorm(n), ncol = 1)
v1 <- c(rep(1, 20), rep(0, 80))
v2 <- c(rep(1, 15), rep(0, 65))
X_1 <- u %*% t(v1) + matrix(rnorm(n * p1), n, p1)
X_2 <- u %*% t(v2) + matrix(rnorm(n * p2), n, p2)

# Run RelPMDCCA
result <- HEPHAESTUS_relPMDCCA(X_1, X_2, lambda = 10, tauW_1 = 0.8, tauW_2 = 0.8, nIter = 100)

# Check canonical correlation
cor(X_1 %*% result$W_1, X_2 %*% result$W_2)

} # }
```
