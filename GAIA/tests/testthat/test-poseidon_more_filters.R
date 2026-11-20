# POSEIDON feature filters not covered in test-poseidon_filtering.R.

test_that("filter_low_variance removes the requested fraction of least-variable features", {
  set.seed(1)
  sds <- c(0.01, 0.02, 5, 6, 0.03, 7, 0.04, 8, 0.05, 9)   # 5 low, 5 high
  X <- sapply(sds, function(s) rnorm(30, sd = s))
  colnames(X) <- paste0("f", seq_along(sds))

  out <- quiet(POSEIDON_filter_low_variance(X, bottom_pct = 0.5, verbose = FALSE))
  expect_equal(ncol(out), 5)
  expect_equal(colnames(out), c("f3", "f4", "f6", "f8", "f10"))   # high-variance, original order
})

test_that("filter_low_variance validates bottom_pct and keeps at least two features", {
  X <- matrix(rnorm(60), 10)
  expect_error(POSEIDON_filter_low_variance(X, bottom_pct = 0), "between 0 and 1")
  expect_error(POSEIDON_filter_low_variance(X, bottom_pct = 1), "between 0 and 1")

  expect_warning(out <- quiet(POSEIDON_filter_low_variance(X, bottom_pct = 0.9, verbose = FALSE)),
                 "Keeping at least 2")
  expect_equal(ncol(out), 2)
})

test_that("filter_low_variance applies per dataset to a list", {
  set.seed(2)
  mk <- function(p) matrix(rnorm(20 * p, sd = rep(c(0.1, 5), length.out = p)), 20,
                           dimnames = list(NULL, paste0("f", 1:p)))
  out <- quiet(POSEIDON_filter_low_variance(list(a = mk(10), b = mk(6)),
                                            bottom_pct = 0.5, verbose = FALSE))
  expect_equal(vapply(out, ncol, integer(1)), c(a = 5L, b = 3L))
})

test_that("filter_pairwise_identical removes features where any two samples share a value", {
  X <- cbind(
    ok    = c(1, 2, 3, 4, 5),
    tie   = c(1, 2, 2, 4, 5),         # samples 2 and 3 identical
    zeros = c(0, 0, 3, 4, 5),         # a typical count-data tie
    ok2   = c(5.5, 1.1, 7.7, 3.3, 9.9)
  )
  out <- quiet(POSEIDON_filter_pairwise_identical(X, verbose = FALSE))
  expect_equal(colnames(out), c("ok", "ok2"))
})

test_that("filter_pairwise_identical honours the tolerance and accepts lists", {
  X <- cbind(near = c(1, 1.0000001, 3, 4), far = c(1, 2, 3, 4))
  expect_equal(colnames(quiet(POSEIDON_filter_pairwise_identical(X, tolerance = 1e-3,
                                                                 verbose = FALSE))), "far")
  expect_equal(colnames(quiet(POSEIDON_filter_pairwise_identical(X, tolerance = 1e-12,
                                                                 verbose = FALSE))),
               c("near", "far"))

  out <- quiet(POSEIDON_filter_pairwise_identical(list(a = X, b = X[, "far", drop = FALSE]),
                                                  tolerance = 1e-3, verbose = FALSE))
  expect_named(out, c("a", "b"))
  expect_equal(ncol(out$a), 1)
})

test_that("filter_pairwise_identical warns when there are fewer than two samples", {
  expect_warning(out <- POSEIDON_filter_pairwise_identical(matrix(1:3, 1), verbose = FALSE),
                 "at least 2 samples")
  expect_equal(dim(out), c(1L, 3L))
})
