# ARTEMIS global shift test: do signals at a set of anchor regions move
# systematically between two groups?

make_anchor <- function(shift = TRUE, seed = 4, n_regions = 100) {
  set.seed(seed)
  counts <- matrix(rpois(n_regions * 6, 50), n_regions, 6,
                   dimnames = list(paste0("r", seq_len(n_regions)), paste0("s", 1:6)))
  if (shift) counts[, 4:6] <- matrix(rpois(n_regions * 3, 65), n_regions, 3)  # group B higher everywhere
  list(counts = counts,
       targets = data.frame(group = rep(c("WT", "KO"), each = 3),
                            row.names = colnames(counts), stringsAsFactors = FALSE))
}

shift_test <- function(q, ...) {
  quiet_msg(ARTEMIS_global_shift_test(q, group_a = "WT", group_b = "KO", ...))
}

test_that("a genuine global shift is detected, in the right direction", {
  res <- shift_test(make_anchor(shift = TRUE))

  expect_named(res, c("per_region", "global_test", "summary", "groups"))
  expect_equal(res$groups, c("WT", "KO"))
  expect_equal(nrow(res$per_region), 100)
  expect_lt(res$global_test$p.value, 1e-6)
  expect_match(res$global_test$test, "Wilcoxon")
  expect_gt(res$summary$median_log2FC, 0.25)      # true shift is log2(65/50) = 0.38
  expect_lt(res$summary$median_log2FC, 0.5)
  expect_gt(res$summary$n_up, 90)
})

test_that("no shift gives no significant result", {
  res <- shift_test(make_anchor(shift = FALSE))
  expect_gt(res$global_test$p.value, 0.01)
  expect_lt(abs(res$summary$median_log2FC), 0.1)
})

test_that("per-region values are consistent with the group means", {
  q <- make_anchor()
  res <- shift_test(q)
  r <- res$per_region[res$per_region$region_id == "r1", ]

  expect_equal(r$mean_WT, mean(q$counts["r1", 1:3]))
  expect_equal(r$mean_KO, mean(q$counts["r1", 4:6]))
  expect_equal(r$raw_diff, r$mean_KO - r$mean_WT)
  expect_equal(r$log2FC, log2((r$mean_KO + 1) / (r$mean_WT + 1)))
  expect_equal(res$summary$n_regions,
               res$summary$n_up + res$summary$n_down + res$summary$n_unchanged)
})

test_that("the median summary and the t-test option can be selected", {
  q <- make_anchor()
  med <- shift_test(q, summary_method = "median")
  tt  <- shift_test(q, test = "t.test")

  expect_equal(med$per_region$median_WT[1], median(q$counts["r1", 1:3]))
  expect_false(grepl("Wilcoxon", tt$global_test$test))
  expect_lt(tt$global_test$p.value, 1e-6)
})

test_that("global shift validates its input", {
  q <- make_anchor()
  expect_error(quiet_msg(ARTEMIS_global_shift_test(list(counts = q$counts))),
               "must be a list with 'counts' and 'targets'")
  expect_error(shift_test(q, group_col = "zzz"), "group_col 'zzz' not found")
  expect_error(quiet_msg(ARTEMIS_global_shift_test(q, group_a = "WT", group_b = "nope")),
               "group_b 'nope' not found")
  expect_error(quiet_msg(ARTEMIS_global_shift_test(q, group_a = "nope", group_b = "KO")),
               "group_a 'nope' not found")
})

test_that("the global shift plot is built from the result", {
  res <- shift_test(make_anchor())
  p <- quiet_msg(ARTEMIS_plot_global_shift(res))
  expect_named(p, c("density", "paired"))
  expect_s3_class(p$density, "ggplot")
  expect_s3_class(p$paired, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p$density))
})
