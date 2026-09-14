# Convex Sparse Canonical Correlation Analysis for Two Datasets

Performs sparse CCA on two datasets using convex optimization with LASSO
or SCAD penalization. Based on Parkhomenko et al. (2009).

## Usage

``` r
HEPHAESTUS_convCCA(
  X_1,
  X_2,
  tauW_1,
  tauW_2,
  initW_1 = NULL,
  initW_2 = NULL,
  nIter = NULL,
  penalty = "LASSO",
  R = 1
)
```

## Arguments

- X_1:

  First dataset matrix (n x p1), samples in rows, features in columns

- X_2:

  Second dataset matrix (n x p2), same samples as X_1

- tauW_1:

  Sparsity parameter for X_1 (higher = more sparse, typically 0.1-0.9)

- tauW_2:

  Sparsity parameter for X_2

- initW_1:

  Optional initial canonical vector for X_1 (default: rep(0.1, p1))

- initW_2:

  Optional initial canonical vector for X_2 (default: rep(0.1, p2))

- nIter:

  Maximum iterations (NULL = auto-converge until diff \< 1e-05,
  recommended)

- penalty:

  Penalty function: "LASSO" (default, L1) or "SCAD" (smoothly clipped
  absolute deviation)

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

# Run ConvCCA
result <- HEPHAESTUS_convCCA(X_1, X_2, tauW_1 = 0.3, tauW_2 = 0.3, nIter = NULL)

# Check sparsity
sum(abs(result$W_1) > 1e-6)  # Number of selected features in X_1

} # }
```
