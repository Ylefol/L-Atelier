# DEMETER_generate_multi_datasets() is the toolkit's own simulator for
# multi-omic data with a shared latent structure; other tests can build on it.

quiet_sim <- function(...) {
  out <- NULL
  capture.output(out <- DEMETER_generate_multi_datasets(...))
  out
}

test_that("simulated datasets have the requested dimensions", {
  sim <- quiet_sim(n_datasets = 2, n_samples = 30, n_features = c(40, 25),
                   signal_samples = 10, signal_features = c(5, 4), seed = 1)

  expect_length(sim$datasets, 2)
  expect_equal(dim(sim$datasets[[1]]), c(30L, 40L))
  expect_equal(dim(sim$datasets[[2]]), c(30L, 25L))
  expect_equal(dim(sim$shared_latent), c(30L, 1L))
  # latent is 1 for signal samples, 0 otherwise
  expect_equal(sum(sim$shared_latent), 10)
  # canonical vectors mark exactly the signal features
  expect_equal(sum(sim$true_canonical[[1]]), 5)
  expect_equal(sum(sim$true_canonical[[2]]), 4)
})

test_that("the same seed reproduces the same data, a different seed does not", {
  a <- quiet_sim(n_datasets = 2, n_samples = 20, n_features = c(10, 10),
                 signal_samples = 5, signal_features = c(3, 3), seed = 7)
  b <- quiet_sim(n_datasets = 2, n_samples = 20, n_features = c(10, 10),
                 signal_samples = 5, signal_features = c(3, 3), seed = 7)
  c2 <- quiet_sim(n_datasets = 2, n_samples = 20, n_features = c(10, 10),
                  signal_samples = 5, signal_features = c(3, 3), seed = 8)

  expect_identical(a$datasets, b$datasets)
  expect_false(identical(a$datasets, c2$datasets))
})

test_that("signal features carry the shared signal", {
  sim <- quiet_sim(n_datasets = 1, n_samples = 200, n_features = 50,
                   signal_samples = 100, signal_features = 10, seed = 3)
  X <- sim$datasets[[1]]
  signal_rows <- sim$shared_latent[, 1] == 1

  # signal features are shifted by ~1 in the signal samples, noise features are not
  shift <- colMeans(X[signal_rows, ]) - colMeans(X[!signal_rows, ])
  expect_gt(mean(shift[1:10]), 0.7)
  expect_lt(abs(mean(shift[11:50])), 0.2)
})

test_that("inconsistent argument lengths are rejected", {
  expect_error(quiet_sim(n_datasets = 3, n_features = c(10, 10)),
               "n_features must have length")
  expect_error(quiet_sim(n_datasets = 2, n_features = c(10, 10),
                         signal_features = 3),
               "signal_features must have length")
})
