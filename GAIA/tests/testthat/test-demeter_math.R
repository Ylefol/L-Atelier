# DEMETER: matrix, correlation, formatting and grid helpers shared by the other
# modules.

make_pd <- function(p = 5, seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(p * 3 * p), ncol = p)
  crossprod(A) / nrow(A)   # symmetric positive definite
}

# ---- matrix utilities ---------------------------------------------------------

test_that("matrix_power reproduces integer powers and the inverse", {
  S <- make_pd()
  expect_equal(DEMETER_matrix_power(S, 2), S %*% S, ignore_attr = TRUE)
  expect_equal(DEMETER_matrix_power(S, 1), S, ignore_attr = TRUE)
  expect_equal(DEMETER_matrix_power(S, -1), solve(S), ignore_attr = TRUE)
})

test_that("the -1/2 power is the inverse square root, and %^% is the same operator", {
  S <- make_pd()
  R <- DEMETER_matrix_power(S, -0.5)
  expect_equal(R %*% R, solve(S), ignore_attr = TRUE)   # (S^-1/2)^2 = S^-1
  expect_equal(R %*% S %*% R, diag(5), ignore_attr = TRUE)  # whitening
  expect_equal(S %^% (-0.5), R)
})

test_that("gettingInverse inverts an invertible matrix", {
  S <- make_pd()
  expect_equal(DEMETER_gettingInverse(S) %*% S, diag(5), tolerance = 1e-8)
})

test_that("gettingInverse returns a generalised inverse for a rank-deficient matrix", {
  set.seed(2)
  B <- crossprod(matrix(rnorm(15), 3, 5))   # 5 x 5, rank 3
  Bp <- DEMETER_gettingInverse(B)
  expect_equal(B %*% Bp %*% B, B, tolerance = 1e-8)   # defining property
  expect_false(anyNA(Bp))
})

# ---- correlation matrix -------------------------------------------------------

test_that("correlation_matrix agrees with cor() and cor.test()", {
  set.seed(3)
  d <- data.frame(a = rnorm(40))
  d$b <- d$a + rnorm(40, sd = 0.5)
  d$c <- rnorm(40)
  res <- DEMETER_correlation_matrix(d, method = "pearson")

  expect_equal(res$correlation, cor(d), ignore_attr = FALSE)
  expect_equal(res$pvalue["a", "b"], cor.test(d$a, d$b)$p.value)
  expect_equal(res$pvalue["a", "c"], cor.test(d$a, d$c)$p.value)
  expect_true(isSymmetric(res$correlation))
  expect_true(isSymmetric(res$pvalue))
  expect_equal(unname(diag(res$correlation)), rep(1, 3))
  expect_equal(unname(diag(res$pvalue)), rep(0, 3))
  expect_equal(res$n["a", "b"], 40)
  expect_null(res$padj)
})

test_that("correlation_matrix handles missing values pairwise", {
  set.seed(4)
  d <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  d$a[1:5] <- NA
  res <- DEMETER_correlation_matrix(d, method = "pearson")

  expect_equal(res$n["a", "b"], 25)
  expect_equal(res$n["b", "c"], 30)
  expect_equal(res$correlation["a", "b"], cor(d$a, d$b, use = "complete.obs"))
})

test_that("correlation_matrix returns NA when fewer than 3 complete pairs exist", {
  d <- data.frame(a = c(1, 2, NA, NA, NA), b = c(NA, 2, 3, 4, 5), c = c(1, 2, 3, 4, 5))
  res <- DEMETER_correlation_matrix(d, method = "pearson")
  expect_true(is.na(res$correlation["a", "b"]))   # only 1 shared observation
  expect_false(is.na(res$correlation["b", "c"]))
})

test_that("correlation_matrix adjusts p-values over the upper triangle only", {
  set.seed(5)
  d <- as.data.frame(matrix(rnorm(40 * 4), ncol = 4))
  res <- DEMETER_correlation_matrix(d, method = "spearman", adjust = "BH")

  up <- upper.tri(res$pvalue)
  expect_equal(res$padj[up], p.adjust(res$pvalue[up], method = "BH"))
  expect_true(isSymmetric(res$padj))
  expect_true(all(res$padj[up] >= res$pvalue[up] - 1e-12))
  expect_equal(res$method, "spearman")
})

test_that("correlation_matrix names unnamed columns V1, V2, ...", {
  set.seed(6)
  res <- DEMETER_correlation_matrix(matrix(rnorm(60), ncol = 3))
  expect_equal(colnames(res$correlation), c("V1", "V2", "V3"))
})

# ---- p-value formatting -------------------------------------------------------

test_that("format_pvalues supports scientific, decimal, hybrid and star formats", {
  p <- c(0.00123, 0.5, NA)
  expect_equal(DEMETER_format_pvalues(p, "scientific"), c("1.23e-03", "5.00e-01", ""))

  expect_equal(DEMETER_format_pvalues(c(0.00123, 0.5, 0.0123456), "decimal"),
               c("0.0012", "0.5", "0.012"))
  expect_equal(DEMETER_format_pvalues(c(0.00123, 0.5, 0.0004), "hybrid"),
               c("0.0012", "0.5", "4.00e-04"))

  expect_equal(DEMETER_format_pvalues(c(0.0005, 0.005, 0.03, 0.07, 0.5), "stars"),
               c("***", "**", "*", ".", ""))
})

test_that("format_pvalues can hide non-significant values and keeps matrix shape", {
  p <- matrix(c(0.001, 0.2, 0.04, 0.9), 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  out <- DEMETER_format_pvalues(p, "scientific", hide_ns = TRUE, sig_threshold = 0.05)

  expect_equal(dim(out), c(2L, 2L))
  expect_equal(dimnames(out), dimnames(p))
  expect_equal(out["r1", "c1"], "1.00e-03")
  expect_equal(out["r2", "c1"], "")          # 0.2 hidden
  expect_equal(out["r1", "c2"], "4.00e-02")
  expect_equal(out["r2", "c2"], "")          # 0.9 hidden
})

test_that("format_pvalue_triangle blanks the diagonal and the opposite triangle", {
  p <- matrix(c(0, 0.01, 0.02,
                0.01, 0, 0.03,
                0.02, 0.03, 0), 3, byrow = TRUE)
  up <- DEMETER_format_pvalue_triangle(p, "upper", hide_ns = FALSE)
  lo <- DEMETER_format_pvalue_triangle(p, "lower", hide_ns = FALSE)

  expect_true(all(up[lower.tri(up, diag = TRUE)] == ""))
  expect_true(all(up[upper.tri(up)] != ""))
  expect_true(all(lo[upper.tri(lo, diag = TRUE)] == ""))
  expect_true(all(lo[lower.tri(lo)] != ""))
})

# ---- parameter grids ----------------------------------------------------------

test_that("coarse grid is the full factorial of the supplied values", {
  g <- quiet(DEMETER_create_coarse_grid(2, c(0.3, 0.5)))
  expect_equal(dim(g), c(4L, 2L))
  expect_equal(names(g), c("param1", "param2"))
  expect_equal(nrow(quiet(DEMETER_create_coarse_grid(3))), 27)   # default 3 values

  expect_error(DEMETER_create_coarse_grid(0), "positive integer")
  expect_error(DEMETER_create_coarse_grid(2, values = numeric(0)), "at least one element")
})

test_that("fine grid spaces points evenly across each range", {
  g <- quiet(DEMETER_create_fine_grid(list(c(0, 1), c(10, 20)), n_points = 3))
  expect_equal(dim(g), c(9L, 2L))
  expect_equal(sort(unique(g$param1)), c(0, 0.5, 1))
  expect_equal(sort(unique(g$param2)), c(10, 15, 20))
})

test_that("fine grid validates its ranges", {
  expect_error(DEMETER_create_fine_grid(list()), "non-empty list")
  expect_error(DEMETER_create_fine_grid(list(c(0, 1, 2))), "length 2")
  expect_error(DEMETER_create_fine_grid(list(c(1, 1))), "min must be less than max")
  expect_error(DEMETER_create_fine_grid(list(c(0, 1)), n_points = 1), "n_points")
})

test_that("ranges_around_values centres a window on each value and respects the bounds", {
  r <- quiet(DEMETER_create_ranges_around_values(c(0.5, 0.05), width = 0.2))
  expect_length(r, 2)
  expect_equal(r[[1]], c(0.4, 0.6))
  expect_equal(r[[2]], c(0, 0.15))   # clipped at the lower bound

  expect_error(DEMETER_create_ranges_around_values(numeric(0)), "at least one element")
  expect_error(DEMETER_create_ranges_around_values(0.5, width = -1), "positive number")
  expect_error(DEMETER_create_ranges_around_values(0.5, bounds = c(1, 0)), "bounds\\[1\\]")
})

test_that("ranges_around_values warns and falls back to the full bounds when a value is outside them", {
  expect_warning(r <- quiet(DEMETER_create_ranges_around_values(1.5, width = 0.2)),
                 "at boundary")
  expect_equal(r[[1]], c(0, 1))
})

# ---- colours ------------------------------------------------------------------

test_that("generate_colors returns the requested number of valid, bright-filtered colours", {
  cols <- DEMETER_generate_colors(12)
  expect_length(cols, 12)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
  expect_false(anyDuplicated(cols) > 0)

  # no near-white colours (mean RGB above the cutoff)
  expect_true(all(colMeans(grDevices::col2rgb(cols)) <= 200))
})
