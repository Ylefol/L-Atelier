# POSEIDON operations on 'quant_result' lists (counts = features x samples,
# targets = one row per sample) that are not covered in test-poseidon_samples.R.

make_quant2 <- function() {
  # 3 biological replicates x 2 technical replicates, 2 groups, 5 features
  targets <- data.frame(
    sample_id = paste0("s", 1:6),
    bio_rep   = rep(1:3, each = 2),
    tech_rep  = rep(1:2, times = 3),
    group     = rep(c("ctrl", "ctrl", "treat"), each = 2),
    batch     = "b1",
    stringsAsFactors = FALSE
  )
  counts <- matrix(c(10, 20,  30, 50,  5, 5,
                      0,  0,   1,  1,  0, 2,
                     100, 100, 90, 110, 200, 220,
                      3,  4,   2,  3,  8, 9,
                      0,  0,   0,  0,  0, 0),
                   nrow = 5, byrow = TRUE,
                   dimnames = list(paste0("f", 1:5), targets$sample_id))
  list(counts = counts, targets = targets,
       annotation = data.frame(id = paste0("f", 1:5), stringsAsFactors = FALSE))
}

# ---- filter_low_counts -------------------------------------------------------

test_that("filter_low_counts keeps features with enough samples above the count threshold", {
  q <- make_quant2()
  out <- quiet(POSEIDON_filter_low_counts(q, min_count = 10, min_samples = 2, verbose = FALSE))

  # f1 reaches 10 in 4 samples, f3 in all 6; f2, f4 and f5 never do
  expect_equal(rownames(out$counts), c("f1", "f3"))
  expect_equal(ncol(out$counts), 6)           # samples untouched
  # the annotation table is filtered along with the counts
  expect_equal(out$annotation$id, c("f1", "f3"))
})

test_that("filter_low_counts thresholds are inclusive and adjustable", {
  q <- make_quant2()
  # >= 3 in at least 2 samples additionally keeps f4 (3, 4, 2, 3, 8, 9)
  out <- quiet(POSEIDON_filter_low_counts(q, min_count = 3, min_samples = 2, verbose = FALSE))
  expect_equal(rownames(out$counts), c("f1", "f3", "f4"))

  strict <- quiet(POSEIDON_filter_low_counts(q, min_count = 10, min_samples = 6, verbose = FALSE))
  expect_equal(rownames(strict$counts), "f3")
})

test_that("filter_low_counts works when there is no annotation", {
  q <- make_quant2(); q$annotation <- NULL
  out <- quiet(POSEIDON_filter_low_counts(q, verbose = FALSE))
  expect_null(out$annotation)
  expect_gt(nrow(out$counts), 0)
})

# ---- filter_by_group ---------------------------------------------------------

test_that("filter_by_group keeps only the requested groups", {
  q <- make_quant2()
  out <- quiet(POSEIDON_filter_by_group(q, groups = "treat", verbose = FALSE))
  expect_equal(colnames(out$counts), c("s5", "s6"))
  expect_equal(out$targets$group, c("treat", "treat"))

  both <- quiet(POSEIDON_filter_by_group(q, groups = c("ctrl", "treat"), verbose = FALSE))
  expect_equal(ncol(both$counts), 6)
})

test_that("filter_by_group reports unknown groups and columns", {
  q <- make_quant2()
  expect_error(quiet(POSEIDON_filter_by_group(q, groups = "nope", verbose = FALSE)))
  expect_error(quiet(POSEIDON_filter_by_group(q, group_col = "nope", groups = "ctrl",
                                              verbose = FALSE)))
})

# ---- average_tech_reps -------------------------------------------------------

test_that("average_tech_reps collapses technical replicates to one column per biological sample", {
  q <- make_quant2()
  out <- quiet(POSEIDON_average_tech_reps(q, verbose = FALSE))

  expect_equal(ncol(out$counts), 3)
  expect_equal(nrow(out$counts), 5)
  expect_true(out$averaged)
  expect_equal(out$averaging_method, "mean")
  expect_equal(out$targets$n_tech_reps, c(2, 2, 2))
  expect_equal(out$targets$bio_rep, 1:3)
})

test_that("average_tech_reps averages the counts of the replicates", {
  q <- make_quant2()
  out <- quiet(POSEIDON_average_tech_reps(q, verbose = FALSE))

  expect_equal(unname(out$counts["f1", ]), c(15, 40, 5))
  expect_equal(unname(out$counts["f3", ]), c(100, 100, 210))
  # the column names match the sample IDs recorded in the targets
  expect_equal(colnames(out$counts), out$targets$sample_id)
})

test_that("average_tech_reps can use the median instead of the mean", {
  q <- make_quant2()
  # three technical replicates for one biological sample, one of them an outlier
  q$counts <- cbind(q$counts, s7 = c(1000, 0, 0, 0, 0))
  q$targets <- rbind(q$targets,
                     data.frame(sample_id = "s7", bio_rep = 1, tech_rep = 3,
                                group = "ctrl", batch = "b1", stringsAsFactors = FALSE))
  mean_out <- quiet(POSEIDON_average_tech_reps(q, method = "mean", verbose = FALSE))
  med_out  <- quiet(POSEIDON_average_tech_reps(q, method = "median", verbose = FALSE))

  expect_equal(unname(med_out$counts["f1", 1]), 20)          # median of 10, 20, 1000
  expect_gt(unname(mean_out$counts["f1", 1]), 300)           # pulled up by the outlier
  expect_equal(med_out$targets$n_tech_reps[1], 3)
})

test_that("average_tech_reps validates its column arguments", {
  q <- make_quant2()
  expect_error(quiet(POSEIDON_average_tech_reps(q, bio_rep_col = "nope", verbose = FALSE)),
               "bio_rep_col 'nope' not found")
  expect_error(quiet(POSEIDON_average_tech_reps(q, tech_rep_col = "nope", verbose = FALSE)),
               "tech_rep_col 'nope' not found")
})
