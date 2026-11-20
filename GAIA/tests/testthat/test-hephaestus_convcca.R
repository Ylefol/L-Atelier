# HEPHAESTUS: sparse CCA across two or more datasets. The simulated data has one
# shared latent variable; only the first k features of each dataset load on it,
# so a working method must put its non-zero weights on those features.

test_that("convCCA puts its non-zero weights on the features that carry the shared signal", {
  s <- simulate_cca()
  res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], tauW_1 = 0.15, tauW_2 = 0.15,
                                  nIter = 100))

  expect_named(res, c("W_1", "W_2"))
  expect_length(res$W_1, 40)
  expect_length(res$W_2, 30)

  nz1 <- nonzero(res$W_1); nz2 <- nonzero(res$W_2)
  expect_true(all(nz1 <= 8))          # no weight on pure-noise features
  expect_gte(length(nz1), 6)          # most of the signal features are found
  expect_true(all(nz2 <= 6))
  expect_gte(length(nz2), 4)
})

test_that("canonical vectors have unit length", {
  s <- simulate_cca()
  res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 100))
  expect_equal(sqrt(sum(res$W_1^2)), 1)
  expect_equal(sqrt(sum(res$W_2^2)), 1)
})

test_that("a larger tau gives sparser weights", {
  s <- simulate_cca()
  loose <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 100))
  tight <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.25, 0.25, nIter = 100))

  expect_lt(length(nonzero(tight$W_1)), length(nonzero(loose$W_1)))
})

test_that("the result is deterministic", {
  s <- simulate_cca()
  a <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 50))
  b <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 50))
  expect_identical(a, b)
})

test_that("a tau that is too high warns that every weight is zero", {
  s <- simulate_cca()
  expect_warning(
    res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.9, 0.9, nIter = 50)),
    "All canonical weights are zero"
  )
  expect_length(nonzero(res$W_1), 0)
})

test_that("the SCAD penalty can be used in place of the lasso", {
  s <- simulate_cca()
  res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 100,
                                  penalty = "SCAD"))
  expect_true(all(nonzero(res$W_1) <= 8))
  expect_gte(length(nonzero(res$W_1)), 6)
})

test_that("convCCA validates its arguments", {
  s <- simulate_cca()
  expect_error(quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 10, R = 0)),
               "positive integer")
  expect_error(quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = "many")),
               "numerical value")
})

test_that("multi_convCCA handles three datasets", {
  s <- simulate_cca()
  res <- quiet(HEPHAESTUS_multi_convCCA(s$X, tau = list(0.15, 0.15, 0.15), nIter = 100))

  expect_named(res, "W")
  expect_equal(vapply(res$W, length, integer(1)), c(40L, 30L, 25L))
  for (i in 1:3) {
    nz <- nonzero(res$W[[i]])
    expect_true(all(nz <= s$k[i]))
    expect_gte(length(nz), s$k[i] - 2)
  }
})

# ---- known issues (reported; fix pending approval) ---------------------------

test_that("convCCA with the default nIter = NULL terminates", {
  # KNOWN BUG: in the auto-converge loop the convergence measure is stored in
  # `diffU` but the loop tests `diff`, which stays at 1, so the loop only ends if
  # a value becomes NA. The documented default call never returns.
  # Remove this skip once HEPHAESTUS_convCCA_algo() updates `diff`.
  skip("known bug: nIter = NULL never converges (loop tests `diff`, assigns `diffU`)")
  s <- simulate_cca()
  setTimeLimit(elapsed = 30, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15))
  expect_named(res, c("W_1", "W_2"))
})

test_that("convCCA with R = 2 returns the second canonical pair", {
  # KNOWN BUG: `for (i in 1:R-1)` parses as (1:R)-1, so the deflation loop runs R
  # times instead of R-1 and R = 2 returns the third pair.
  # Remove this skip once the loop is `seq_len(R - 1)`.
  skip("known bug: R = 2 returns the third pair (for (i in 1:R-1) precedence)")
  s <- simulate_cca()
  X1 <- s$X[[1]]; X2 <- s$X[[2]]; tau <- 0.15
  first <- quiet(HEPHAESTUS_convCCA_algo(X1, X2, tau, tau, nIter = 100))
  deflate <- function(w1, w2) {
    list(X1 = rbind(X1, t(w1) %*% t(X1) %*% X1, t(w2) %*% t(X2) %*% X1),
         X2 = rbind(X2, t(w2) %*% t(X2) %*% X2, t(w1) %*% t(X1) %*% X2))
  }
  d <- deflate(first$w_1, first$w_2)
  second <- quiet(HEPHAESTUS_convCCA_algo(d$X1, d$X2, tau, tau, nIter = 100))

  res <- quiet(HEPHAESTUS_convCCA(X1, X2, tau, tau, nIter = 100, R = 2))
  expect_equal(as.numeric(res$W_1), as.numeric(second$w_1))
})
