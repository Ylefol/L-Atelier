# POSEIDON: preparing clinical / mixed tables for regression-type models.

make_clinical <- function(n = 30, seed = 1) {
  set.seed(seed)
  data.frame(
    outcome = rnorm(n),
    age = rnorm(n, 50, 10),
    sex = factor(rep(c("M", "F"), n / 2)),
    grp = factor(rep(c("a", "b", "c"), n / 3)),
    sev = factor(rep(c("low", "mid", "high"), n / 3),
                 levels = c("low", "mid", "high"), ordered = TRUE)
  )
}

encode <- function(d, ...) quiet(POSEIDON_encode_for_regression(d, target_var = "outcome", ...))

test_that("encode_for_regression dummy-codes categorical variables against a reference level", {
  e <- encode(make_clinical())

  # one column per non-reference level; the numeric variable is untouched
  expect_equal(colnames(e$X), c("age", "sex_M", "grp_b", "grp_c", "sev_mid", "sev_high"))
  expect_equal(nrow(e$X), 30)
  expect_equal(e$target_name, "outcome")
  expect_equal(e$var_mapping$reference_levels$sex, "F")
  expect_equal(e$var_mapping$reference_levels$grp, "a")
  expect_true(all(e$X[, c("sex_M", "grp_b", "grp_c")] %in% c(0, 1)))
})

test_that("encoded columns agree with the original data", {
  d <- make_clinical()
  e <- encode(d)

  expect_equal(unname(e$X[, "age"]), d$age)
  expect_equal(unname(e$X[, "sex_M"]), as.numeric(d$sex == "M"))
  expect_equal(unname(e$X[, "grp_c"]), as.numeric(d$grp == "c"))
  expect_equal(e$y, d$outcome)
  # exactly one dummy (or the reference) per categorical variable and sample
  expect_true(all(e$X[, "grp_b"] + e$X[, "grp_c"] <= 1))
})

test_that("the variable mapping links columns back to the original variables", {
  e <- encode(make_clinical())
  vm <- e$var_mapping

  expect_equal(vm$original_vars, c("age", "sex", "grp", "sev"))
  expect_equal(unname(vm$var_types[c("age", "sex")]), c("quantitative", "categorical"))
  expect_equal(vm$col_to_var, c("age", "sex", "grp", "grp", "sev", "sev"))
  expect_equal(unname(vm$groups), c(1, 2, 3, 3, 4, 4))
})

test_that("rows with a missing target are dropped, and counted", {
  d <- make_clinical()
  d$outcome[3] <- NA
  e <- encode(d)
  expect_equal(nrow(e$X), 29)
  expect_length(e$y, 29)
  expect_false("3" %in% e$sample_ids)
  expect_equal(e$n_removed_na, 1)
})

test_that("encode_for_regression rejects a missing target column", {
  expect_error(quiet(POSEIDON_encode_for_regression(make_clinical(), target_var = "nope")))
})

# ---- scale_predictors --------------------------------------------------------

test_that("zscore scaling centres and scales numeric columns but leaves dummies alone", {
  e <- encode(make_clinical())
  s <- quiet(POSEIDON_scale_predictors(e$X, method = "zscore", verbose = FALSE))

  expect_equal(s$scaled_cols, "age")
  expect_equal(s$unscaled_cols, c("sex_M", "grp_b", "grp_c", "sev_mid", "sev_high"))
  expect_equal(mean(s$X_scaled[, "age"]), 0, tolerance = 1e-12)
  expect_equal(sd(s$X_scaled[, "age"]), 1, tolerance = 1e-12)
  expect_equal(s$X_scaled[, "sex_M"], e$X[, "sex_M"])   # binary column unchanged
  expect_equal(s$method, "zscore")
})

test_that("minmax scaling maps numeric columns onto [0, 1]", {
  e <- encode(make_clinical())
  s <- quiet(POSEIDON_scale_predictors(e$X, method = "minmax", verbose = FALSE))
  expect_equal(range(s$X_scaled[, "age"]), c(0, 1))
})

test_that("robust scaling centres on the median", {
  e <- encode(make_clinical())
  s <- quiet(POSEIDON_scale_predictors(e$X, method = "robust", verbose = FALSE))
  expect_equal(stats::median(s$X_scaled[, "age"]), 0, tolerance = 1e-12)
})

test_that("exclude_binary = FALSE scales every column", {
  e <- encode(make_clinical())
  s <- quiet(POSEIDON_scale_predictors(e$X, method = "zscore", exclude_binary = FALSE,
                                       verbose = FALSE))
  expect_equal(length(s$scaled_cols), ncol(e$X))
  expect_equal(mean(s$X_scaled[, "sex_M"]), 0, tolerance = 1e-12)
})

test_that("scale_predictors validates its input", {
  expect_error(quiet(POSEIDON_scale_predictors(1:5)), "matrix or data.frame")
  expect_error(quiet(POSEIDON_scale_predictors(matrix(rnorm(10), 5), method = "bogus")),
               "must be 'zscore', 'minmax', or 'robust'")
})
