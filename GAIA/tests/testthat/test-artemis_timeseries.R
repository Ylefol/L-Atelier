# ARTEMIS time-series differential analysis (DESeq2), on simulated data with
# two planted effects; see helper-timeseries.R.

norm_ts <- function() {
  d <- make_timeseries()
  list(d = d, norm = quiet_msg(ARTEMIS_normalize_timeseries(d$counts, d$targets, verbose = FALSE)))
}

test_that("normalize_timeseries returns an artemis_ts_norm object with size factors", {
  x <- norm_ts()
  expect_s3_class(x$norm, "artemis_ts_norm")
  expect_named(x$norm, c("dds", "norm_counts", "targets", "size_factors", "parameters"))
  expect_equal(dim(x$norm$norm_counts), dim(x$d$counts))
  expect_true(all(x$norm$size_factors > 0))
  expect_equal(length(x$norm$size_factors), 24)

  out <- utils::capture.output(print(x$norm))
  expect_true(any(grepl("Samples: 24", out)))
  expect_true(any(grepl("Genes: 150", out)))
})

test_that("normalize_timeseries validates its inputs", {
  d <- make_timeseries()
  expect_error(quiet_msg(ARTEMIS_normalize_timeseries(d$counts, d$targets,
                                                      group_col = "zzz", verbose = FALSE)),
               "'zzz' column not found in targets")
  expect_error(quiet_msg(ARTEMIS_normalize_timeseries(d$counts, d$targets,
                                                      time_col = "zzz", verbose = FALSE)),
               "'zzz' column not found")
})

test_that("temporal comparisons find genes that change over time, pooling both groups", {
  x <- norm_ts()
  te <- quiet_msg(ARTEMIS_timeseries_temporal(x$norm, verbose = FALSE))

  expect_s3_class(te, "artemis_ts_de")
  expect_named(te$results, c("TP_T2_vs_TP_T1", "TP_T3_vs_TP_T2"))   # consecutive by default

  # T2 vs T1 pools both groups, so only the group-specific genes (which rise in
  # 'trt' at T2) stand out; T3 vs T2 picks up the genes that rise in both
  s21 <- sig_genes(te$results$TP_T2_vs_TP_T1)
  s32 <- sig_genes(te$results$TP_T3_vs_TP_T2)
  expect_gte(sum(x$d$group_time %in% s21), 8)
  expect_gte(sum(x$d$temporal %in% s32), 8)
  expect_lte(sum(x$d$null %in% s21), 4)
  expect_lte(sum(x$d$null %in% s32), 4)
})

test_that("comparisons = 'all' adds the non-consecutive pair", {
  x <- norm_ts()
  te <- quiet_msg(ARTEMIS_timeseries_temporal(x$norm, comparisons = "all", verbose = FALSE))
  expect_named(te$results, c("TP_T2_vs_TP_T1", "TP_T3_vs_TP_T1", "TP_T3_vs_TP_T2"))
  expect_error(quiet_msg(ARTEMIS_timeseries_temporal(x$norm, comparisons = "zzz",
                                                     verbose = FALSE)),
               "should be one of")
})

test_that("conditional comparisons find the group difference at each timepoint", {
  x <- norm_ts()
  co <- quiet_msg(ARTEMIS_timeseries_conditional(x$norm, reference = "ctrl",
                                                 experiment = "trt", verbose = FALSE))

  expect_named(co$results, c("trt_vs_ctrl_TPT1", "trt_vs_ctrl_TPT2", "trt_vs_ctrl_TPT3"))
  # no difference at T1, the planted group-specific genes at T2 and T3
  expect_lte(length(sig_genes(co$results$trt_vs_ctrl_TPT1)), 2)
  for (tp in c("trt_vs_ctrl_TPT2", "trt_vs_ctrl_TPT3")) {
    s <- sig_genes(co$results[[tp]])
    expect_gte(sum(x$d$group_time %in% s), 8)
    expect_lte(sum(!(s %in% x$d$group_time)), 4)
  }
  # and the genes that rise in both groups are NOT group differences
  expect_lte(sum(x$d$temporal %in% sig_genes(co$results$trt_vs_ctrl_TPT3)), 2)
})

test_that("the two time-series analyses reject bad input clearly", {
  x <- norm_ts()
  expect_error(ARTEMIS_timeseries_temporal(list(), verbose = FALSE),
               "must be an artemis_ts_norm object")
  expect_error(quiet_msg(ARTEMIS_timeseries_conditional(x$norm, "nope", "trt", verbose = FALSE)),
               "reference 'nope' not found in 'group' column")
})

test_that("time-series DE results feed the same gene-selection step as a plain comparison", {
  x <- norm_ts()
  te <- quiet_msg(ARTEMIS_timeseries_temporal(x$norm, comparisons = "all", verbose = FALSE))
  sel <- ARTEMIS_select_de_genes(te, l2fc_thresh = 1, verbose = FALSE)

  expect_gte(length(intersect(sel, c(x$d$temporal, x$d$group_time))), 16)
  expect_named(attr(sel, "per_comparison"), names(te$results))
})
