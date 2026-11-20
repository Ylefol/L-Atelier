# POSEIDON batch correction and imputation. sva (ComBat) and mice are Suggests
# dependencies, so the tests skip where they are absent.

make_batch_counts <- function(seed = 6) {
  set.seed(seed)
  G <- 300; n <- 12
  group <- rep(c("A", "B"), 6)
  batch <- rep(c("b1", "b2"), each = 6)
  mu <- matrix(rgamma(G, 2, 0.02) + 20, G, n)
  mu[1:20,  group == "B"] <- mu[1:20, group == "B"] * 3      # a true group effect
  mu[1:150, batch == "b2"] <- mu[1:150, batch == "b2"] * 2.5 # a strong batch effect
  cnt <- matrix(rnbinom(G * n, mu = mu, size = 20), G)
  dimnames(cnt) <- list(sprintf("g%03d", 1:G), sprintf("s%02d", 1:n))
  list(counts = cnt,
       targets = data.frame(batch = batch, group = group, row.names = colnames(cnt),
                            stringsAsFactors = FALSE),
       batch = batch, group = group)
}

# share of the variance of a principal component explained by a factor
pc_r2 <- function(counts, f, pc = 1) {
  p <- stats::prcomp(t(log2(counts + 1)))
  summary(stats::lm(p$x[, pc] ~ f))$r.squared
}

test_that("ComBat removes the batch signal and keeps the group signal", {
  skip_if_not_installed("sva")
  d <- make_batch_counts()
  qr <- list(counts = d$counts, targets = d$targets)
  out <- quiet_msg(POSEIDON_correct_batch_counts(qr, "batch", "group"))

  expect_equal(dim(out$counts), dim(d$counts))
  expect_equal(dimnames(out$counts), dimnames(d$counts))
  expect_true(out$batch_corrected)
  expect_false(anyNA(out$counts))

  # before: PC1 is the batch; after: batch explains almost nothing, group dominates
  expect_gt(pc_r2(d$counts, factor(d$batch)), 0.8)
  expect_lt(pc_r2(out$counts, factor(d$batch)), 0.1)
  expect_gt(pc_r2(out$counts, factor(d$group)), 0.8)
})

test_that("batch correction leaves the sample table alone and rejects bad columns", {
  skip_if_not_installed("sva")
  d <- make_batch_counts()
  qr <- list(counts = d$counts, targets = d$targets)
  out <- quiet_msg(POSEIDON_correct_batch_counts(qr, "batch", "group"))
  expect_equal(out$targets, d$targets)

  expect_error(quiet_msg(POSEIDON_correct_batch_counts(qr, "zzz", "group")))
})

test_that("mice imputation fills the gaps with values close to the truth", {
  skip_if_not_installed("mice")
  set.seed(9)
  n <- 80
  x <- rnorm(n)
  # b follows a with realistic noise (R^2 about 0.8). Near-perfect collinearity
  # (e.g. sd = 0.2) makes mice itself impute poorly, independent of GAIA.
  full <- data.frame(a = x, b = x * 2 + rnorm(n, sd = 1), z = rnorm(n))
  d <- full
  gaps <- sample(n, 12)
  d$b[gaps] <- NA

  mi <- quiet_msg(POSEIDON_mice_impute(d, n_imputations = 5, max_iter = 5, seed = 1))

  expect_named(mi, c("datasets", "mean_dataset", "log_vars", "n_imputations", "n_obs", "n_vars"))
  expect_length(mi$datasets, 5)
  expect_equal(mi$n_imputations, 5)
  expect_false(anyNA(mi$mean_dataset))
  expect_equal(mi$mean_dataset$a, d$a)                                   # observed values untouched
  expect_equal(mi$mean_dataset$b[-gaps], d$b[-gaps])
  # b follows a, so imputed values must follow the truth
  expect_gt(cor(mi$mean_dataset$b[gaps], full$b[gaps]), 0.8)
})

test_that("mice imputation is reproducible for a given seed", {
  skip_if_not_installed("mice")
  set.seed(9)
  d <- data.frame(a = rnorm(40), b = rnorm(40)); d$b[1:5] <- NA
  f <- function() quiet_msg(POSEIDON_mice_impute(d, n_imputations = 3, max_iter = 3, seed = 4))$mean_dataset
  expect_identical(f(), f())
})
