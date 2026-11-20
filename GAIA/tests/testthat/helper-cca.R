# Shared fixture for the sparse-CCA tests and the AETHER weight-plot tests.

# n = 120 samples; datasets with 40 / 30 / 25 features; the first 8 / 6 / 5 carry signal.
simulate_cca <- function(seed = 1, n = 120, p = c(40, 30, 25), k = c(8, 6, 5), load = 1.5) {
  set.seed(seed)
  u <- rnorm(n)
  X <- lapply(seq_along(p), function(i) {
    v <- c(rep(load, k[i]), rep(0, p[i] - k[i]))
    u %*% t(v) + matrix(rnorm(n * p[i]), n)
  })
  list(X = X, k = k)
}

nonzero <- function(w) which(abs(w) > 1e-8)
