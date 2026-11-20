# POSEIDON feature filters and dataset validation. These take samples x features
# matrices (or lists of them), so the same code serves every omic.

test_that("scale_data centres and scales each feature", {
  set.seed(1)
  X <- matrix(rnorm(60, mean = 5, sd = 3), nrow = 20)
  S <- POSEIDON_scale_data(X)

  expect_equal(unname(colMeans(S)), rep(0, 3), tolerance = 1e-12)
  expect_equal(unname(apply(S, 2, sd)), rep(1, 3), tolerance = 1e-12)

  # centre only: means are zero but spread is untouched
  C <- POSEIDON_scale_data(X, scale = FALSE)
  expect_equal(unname(apply(C, 2, sd)), unname(apply(X, 2, sd)))
})

test_that("filter_zero_variance drops constant features from a matrix", {
  X <- cbind(a = c(1, 2, 3, 4), const = rep(5, 4), b = c(2, 2, 3, 9))
  out <- POSEIDON_filter_zero_variance(X, verbose = FALSE)

  expect_equal(colnames(out), c("a", "b"))
  expect_equal(nrow(out), 4)
})

test_that("filter_zero_variance respects variance_threshold", {
  X <- cbind(low = c(1, 1.1, 1, 1.1), high = c(1, 5, 2, 9))
  expect_equal(colnames(POSEIDON_filter_zero_variance(X, verbose = FALSE)),
               c("low", "high"))
  expect_equal(
    colnames(POSEIDON_filter_zero_variance(X, variance_threshold = 1,
                                           verbose = FALSE)),
    "high"
  )
})

test_that("filter_zero_variance works on a list of datasets independently", {
  X1 <- cbind(a = 1:4, const = rep(1, 4))
  X2 <- cbind(x = rep(2, 4), y = c(1, 4, 2, 8), z = c(3, 1, 4, 1))
  out <- POSEIDON_filter_zero_variance(list(X1, X2), verbose = FALSE)

  expect_length(out, 2)
  expect_equal(colnames(out[[1]]), "a")
  expect_equal(colnames(out[[2]]), c("y", "z"))
})

test_that("filter_shared_variance keeps only features varying in every dataset", {
  # by_row = FALSE: features are columns
  X1 <- cbind(f1 = 1:4, f2 = rep(1, 4), f3 = c(4, 1, 3, 2))
  X2 <- cbind(f1 = c(2, 5, 1, 3), f2 = c(1, 2, 3, 4), f3 = rep(7, 4))
  out <- POSEIDON_filter_shared_variance(list(X1, X2), by_row = FALSE,
                                         verbose = FALSE)

  expect_equal(colnames(out[[1]]), "f1")
  expect_equal(colnames(out[[2]]), "f1")
})

test_that("filter_top_variable keeps the n most variable features, in original order", {
  set.seed(2)
  X <- cbind(low = rnorm(30, sd = 0.1), high = rnorm(30, sd = 10),
             mid = rnorm(30, sd = 2), tiny = rnorm(30, sd = 0.01))
  out <- POSEIDON_filter_top_variable(X, n_top = 2, verbose = FALSE)

  expect_equal(ncol(out), 2)
  expect_equal(colnames(out), c("high", "mid"))   # original column order kept
})

test_that("filter_top_variable keeps everything (with a warning) when n_top is too large", {
  X <- matrix(rnorm(40), nrow = 10)
  expect_warning(out <- POSEIDON_filter_top_variable(X, n_top = 100),
                 "Keeping all features")
  expect_equal(dim(out), dim(X))
})

test_that("validate_multi_dataset accepts consistent datasets", {
  X <- matrix(rnorm(40), nrow = 10, dimnames = list(paste0("s", 1:10), NULL))
  Y <- matrix(rnorm(50), nrow = 10, dimnames = list(paste0("s", 1:10), NULL))
  expect_true(POSEIDON_validate_multi_dataset(list(X, Y)))
})

test_that("validate_multi_dataset rejects missing values, mismatched samples and empty data", {
  X <- matrix(rnorm(40), nrow = 10)

  Xna <- X; Xna[3, 2] <- NA
  expect_error(POSEIDON_validate_multi_dataset(list(X, Xna)), "missing values")

  expect_error(POSEIDON_validate_multi_dataset(list(X, X[1:8, ])),
               "expected 10")

  expect_error(POSEIDON_validate_multi_dataset(list(X, X[, 0, drop = FALSE])),
               "zero features")

  expect_error(POSEIDON_validate_multi_dataset(X), "must be a list")
})

test_that("validate_multi_dataset warns when sample order differs", {
  X <- matrix(rnorm(20), nrow = 5, dimnames = list(letters[1:5], NULL))
  Y <- X[c(2, 1, 3, 4, 5), ]
  expect_warning(POSEIDON_validate_multi_dataset(list(X, Y),
                                                 sample_names = letters[1:5]),
                 "different sample names")
})

test_that("the simulator and the filters work together across datasets", {
  # Ties the pieces together: simulate three 'omics', filter, validate.
  out <- NULL
  capture.output(out <- DEMETER_generate_multi_datasets(
    n_datasets = 3, n_samples = 25, n_features = c(30, 40, 20),
    signal_samples = 10, signal_features = c(5, 5, 4), seed = 11))
  filtered <- POSEIDON_filter_zero_variance(out$datasets, verbose = FALSE)

  expect_true(POSEIDON_validate_multi_dataset(filtered))
  expect_equal(vapply(filtered, nrow, integer(1)), rep(25L, 3))
})
