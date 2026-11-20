# HADES: quality control for clinical / mixed-type tables. The fixture has one
# problem of each kind planted in it, so each check has something to find.

make_mixed <- function(n = 40, seed = 1) {
  set.seed(seed)
  d <- data.frame(
    id        = paste0("p", seq_len(n)),
    age       = round(rnorm(n, 50, 10)),
    bmi       = rnorm(n, 25, 4),
    const     = rep(1, n),                              # zero variance
    mostly_na = c(rep(NA, 25), rnorm(15)),              # 62.5% missing
    dominant  = rep("a", n),
    sex       = rep(c("M", "F"), n / 2),
    score_chr = as.character(round(rnorm(n, 10, 2))),   # numbers stored as text
    grade     = rep(c(0, 1, 2, 3), n / 4),              # looks categorical
    stringsAsFactors = FALSE
  )
  d$bmi2 <- d$bmi * 1.001 + rnorm(n, sd = 0.001)        # near-duplicate of bmi
  d$dominant[1] <- "b"
  d$score_chr[3] <- "unknown"                           # NA-like string
  d
}

qc_quiet <- function(d, ...) {
  suppressWarnings(quiet(HADES_qc_mixed_data(d, sample_id_col = "id", ...)))
}

test_that("qc_mixed_data returns a hades_qc report with the documented parts", {
  qc <- qc_quiet(make_mixed())
  expect_s3_class(qc, "hades_qc")
  expect_named(qc, c("variables", "samples", "correlations", "type_issues",
                     "id_issues", "summary", "thresholds"))
  expect_equal(nrow(qc$variables), 10)
  expect_equal(qc$summary$n_variables, 10)
  expect_equal(qc$summary$n_samples, 40)
})

test_that("variables are classified as quantitative or qualitative", {
  v <- qc_quiet(make_mixed())$variables
  type <- setNames(v$type, v$variable)
  expect_equal(unname(type[c("age", "bmi", "grade")]), rep("quantitative", 3))
  expect_equal(unname(type[c("id", "sex", "dominant")]), rep("qualitative", 3))
})

test_that("each planted problem is flagged, and clean variables are not", {
  v <- qc_quiet(make_mixed())$variables
  rownames(v) <- v$variable

  expect_true(v["mostly_na", "flagged_na"])
  expect_equal(v["mostly_na", "pct_missing"], 25 / 40)
  expect_true(v["const", "flagged_variance"])
  expect_true(v["score_chr", "has_na_strings"])       # "unknown" is an NA-like string
  expect_true(v["score_chr", "potential_numeric"])    # numbers stored as text

  expect_false(any(v[c("age", "bmi", "sex"), "flagged_any"]))
})

test_that("highly correlated variable pairs are reported", {
  qc <- qc_quiet(make_mixed())
  expect_equal(nrow(qc$correlations), 1)
  expect_setequal(c(qc$correlations$var1, qc$correlations$var2), c("bmi", "bmi2"))
  expect_gt(qc$correlations$correlation, 0.99)
  expect_equal(qc$summary$n_high_correlation_pairs, 1)
})

test_that("thresholds are applied as given", {
  d <- make_mixed()
  strict <- qc_quiet(d, na_threshold = 0.7)$variables
  expect_false(strict$flagged_na[strict$variable == "mostly_na"])   # 62.5% < 70%

  qc <- qc_quiet(d, check_correlations = FALSE, check_samples = FALSE)
  expect_null(qc$correlations)   # not computed when switched off
})

test_that("type issues are collected per variable", {
  ti <- qc_quiet(make_mixed())$type_issues
  expect_equal(names(ti$potential_numeric), "score_chr")
  expect_equal(names(ti$na_like_strings), "score_chr")
  expect_setequal(names(ti$potential_categorical), c("const", "grade"))
})

test_that("sample IDs are checked", {
  d <- make_mixed()
  ok <- qc_quiet(d)
  expect_false(ok$summary$id_issues_any)

  d$id[2] <- d$id[1]
  dup <- qc_quiet(d)
  expect_true(dup$summary$id_issues_any)
  expect_gt(dup$summary$id_duplicates, 0)
})

# ---- convert_to_factor -------------------------------------------------------

test_that("convert_to_factor converts the columns the QC report flagged as categorical", {
  d <- make_mixed()
  qc <- qc_quiet(d)
  out <- suppressWarnings(quiet(HADES_convert_to_factor(d, qc_report = qc)))

  expect_s3_class(out$grade, "factor")
  expect_s3_class(out$const, "factor")
  expect_type(out$age, "double")      # continuous columns are left alone
  expect_equal(nrow(out), nrow(d))
})

test_that("convert_to_factor honours explicit columns, ordering and exclusions", {
  d <- make_mixed()
  out <- quiet(HADES_convert_to_factor(d, columns = c("grade", "sex"),
                                       ordered = list(grade = TRUE)))
  expect_s3_class(out$grade, "ordered")
  expect_equal(levels(out$grade), c("0", "1", "2", "3"))   # numeric order
  expect_s3_class(out$sex, "factor")
  expect_false(is.ordered(out$sex))

  qc <- qc_quiet(d)
  out2 <- quiet(HADES_convert_to_factor(d, qc_report = qc, exclude = "grade"))
  expect_type(out2$grade, "double")
})

# ---- drop_high_na_variables --------------------------------------------------

test_that("drop_high_na_variables removes columns above the NA threshold", {
  d <- make_mixed()
  out <- quiet(HADES_drop_high_na_variables(d, na_threshold = 0.5))
  expect_equal(setdiff(names(d), names(out)), "mostly_na")

  # a looser threshold keeps everything
  keep <- quiet(HADES_drop_high_na_variables(d, na_threshold = 0.9))
  expect_equal(names(keep), names(d))
})

test_that("drop_high_na_variables never drops protected columns", {
  d <- make_mixed()
  out <- quiet(HADES_drop_high_na_variables(d, na_threshold = 0.5, protect = "mostly_na"))
  expect_true("mostly_na" %in% names(out))
})

test_that("drop_high_na_variables treats NA-like strings as missing", {
  d <- data.frame(a = c("N/A", "missing", "NA", "1", "2"),
                  b = c("1", "2", "3", "4", "5"), stringsAsFactors = FALSE)
  out <- quiet(HADES_drop_high_na_variables(d, na_threshold = 0.5))
  expect_equal(names(out), "b")                      # a is 60% NA-like

  d2 <- data.frame(x = c("N/A", rep("1.5", 9)), stringsAsFactors = FALSE)
  out2 <- quiet(HADES_drop_high_na_variables(d2, na_threshold = 0.9))
  expect_true(is.na(out2$x[1]))                      # the string was converted to NA
})
