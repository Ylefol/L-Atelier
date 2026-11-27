# HADES Olink QC: every step records what it did in $qc_log (arguments actually
# used, counts, and the IDs removed or flagged with a reason), and
# HADES_get_olink_qc_log() reads it back. Data are simulated with planted problems
# (see helper-olink.R), so each removal has a known answer.

# ---- the loader starts the log ------------------------------------------------

test_that("the loader records itself as the first log entry, including what it dropped", {
  ol <- load_sim_olink()

  expect_length(ol$qc_log, 1)
  e <- ol$qc_log[[1]]
  expect_equal(e$step, "ELEUTHIA_load_olink")
  expect_equal(e$action, "loaded")
  expect_true(e$params$restrict_to_metadata)
  expect_equal(e$params$npx_file, "npx.csv")           # file name only, no local path
  expect_equal(e$samples_before, 25L)                  # 24 in metadata + X1
  expect_equal(e$samples_after, 24L)
  expect_equal(e$proteins_after, 30L)
  expect_equal(e$removed_samples$sample_id, "X1")
  expect_match(e$removed_samples$reason, "not in metadata")
  expect_equal(e$details$instrument_rows_dropped, 25L) # one ext_ctrl row per sample
  expect_equal(e$details$blank_sample_id_rows_dropped, 3L)
  expect_equal(dim(ol$wide), c(30L, 24L))
})

# ---- sample filter ------------------------------------------------------------

test_that("HADES_filter_olink logs which samples it removed and why", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))

  expect_equal(ncol(ol$wide), 20L)
  e <- ol$qc_log[[2]]
  expect_equal(e$step, "HADES_filter_olink")
  expect_equal(e$action, "removed")
  expect_equal(e$params, list(keep_controls = FALSE, filter_qc = TRUE))
  expect_equal(c(e$samples_before, e$samples_after), c(24L, 20L))

  rs <- e$removed_samples
  expect_setequal(rs$sample_id, c("N1", "N2", "S21", "S22"))
  expect_true(all(grepl("SampleType", rs$reason[rs$sample_id %in% c("N1", "N2")])))
  expect_true(all(rs$reason[rs$sample_id %in% c("S21", "S22")] == "SampleQC == FAIL"))
  expect_false(any(c("N1", "N2", "S21", "S22") %in% colnames(ol$wide)))
})

test_that("the arguments recorded are the ones the step ran with, defaults included", {
  raw <- load_sim_olink()
  keep <- quiet(HADES_filter_olink(raw, keep_controls = TRUE, filter_qc = FALSE))

  expect_equal(keep$qc_log[[2]]$params, list(keep_controls = TRUE, filter_qc = FALSE))
  expect_equal(nrow(keep$qc_log[[2]]$removed_samples), 0L)
  expect_equal(ncol(keep$wide), 24L)

  # a call that relies on defaults still records them
  d <- quiet(HADES_filter_olink(raw))
  expect_false(d$qc_log[[2]]$params$keep_controls)
  expect_true(d$qc_log[[2]]$params$filter_qc)
})

test_that("samples with a missing SampleType or SampleQC do not produce NA rows", {
  raw <- load_sim_olink()
  raw$sample_meta$SampleQC[raw$sample_meta$SampleID == "S03"]   <- NA   # not a failure
  raw$sample_meta$SampleType[raw$sample_meta$SampleID == "S04"] <- NA   # not a SAMPLE
  ol <- quiet(HADES_filter_olink(raw))

  expect_false(anyNA(ol$sample_meta$SampleID))
  expect_true("S03" %in% ol$sample_meta$SampleID)
  expect_false("S04" %in% ol$sample_meta$SampleID)
  rs <- ol$qc_log[[2]]$removed_samples
  expect_match(rs$reason[rs$sample_id == "S04"], "NA")
})

# ---- protein filters ----------------------------------------------------------

test_that("the missingness filter logs the proteins removed with their NA fractions", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_filter_olink_proteins(ol, max_na_fraction = 0.2))

  e <- ol$qc_log[[3]]
  expect_equal(e$step, "HADES_filter_olink_proteins")
  expect_equal(e$params, list(max_na_fraction = 0.2))
  expect_equal(c(e$proteins_before, e$proteins_after), c(30L, 28L))
  expect_equal(c(e$samples_before, e$samples_after), c(20L, 20L))

  rp <- e$removed_proteins
  expect_setequal(rp$OlinkID, c("OID001", "OID002"))
  expect_equal(rp$value, c(0.5, 0.5))
  expect_setequal(rp$Assay, c("GENE001", "GENE002"))
  expect_true(all(grepl("NA fraction > 0.2", rp$reason)))
  expect_false(any(c("OID001", "OID002") %in% rownames(ol$wide)))
})

test_that("the AssayQC warn filter logs the proteins removed with their warn fractions", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_filter_olink_assay_warn(ol, max_warn_fraction = 0.2))

  e <- ol$qc_log[[3]]
  expect_equal(e$step, "HADES_filter_olink_assay_warn")
  expect_equal(e$params, list(max_warn_fraction = 0.2))
  expect_equal(e$removed_proteins$OlinkID, "OID003")
  expect_equal(e$removed_proteins$value, 0.5)
  expect_true("OID004" %in% rownames(ol$wide))     # 1 warn in 24 samples: kept
  expect_equal(c(e$proteins_before, e$proteins_after), c(30L, 29L))
})

test_that("a threshold that removes nothing is still logged, with empty removed tables", {
  ol <- quiet(HADES_filter_olink_proteins(load_sim_olink(), max_na_fraction = 1))

  e <- ol$qc_log[[2]]
  expect_equal(nrow(e$removed_proteins), 0L)
  expect_equal(e$proteins_before, e$proteins_after)
})

# ---- outliers: flagged vs removed --------------------------------------------

test_that("outlier detection without filter logs the sample as flagged, and keeps it", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_detect_outliers_olink(ol, n_pcs = 3, threshold = 3))

  e <- ol$qc_log[[3]]
  expect_equal(e$step, "HADES_detect_outliers_olink")
  expect_equal(e$action, "flagged")
  expect_equal(e$flagged_samples$sample_id, "S05")
  expect_equal(nrow(e$removed_samples), 0L)
  expect_equal(c(e$samples_before, e$samples_after), c(20L, 20L))
  expect_true("S05" %in% colnames(ol$wide))              # still in the data
  expect_true(ol$sample_meta$outlier[ol$sample_meta$SampleID == "S05"])
  expect_equal(e$params, list(n_pcs = 3, threshold = 3, method = "sd", filter = FALSE))
  expect_true(is.numeric(e$details$cutoff))
  expect_gte(e$flagged_samples$value, e$details$cutoff)
})

test_that("outlier detection with filter logs the sample as removed", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_detect_outliers_olink(ol, n_pcs = 3, threshold = 3, filter = TRUE))

  e <- ol$qc_log[[3]]
  expect_equal(e$action, "removed")
  expect_equal(e$removed_samples$sample_id, "S05")
  expect_equal(nrow(e$flagged_samples), 0L)
  expect_equal(c(e$samples_before, e$samples_after), c(20L, 19L))
  expect_false("S05" %in% colnames(ol$wide))
})

test_that("each outlier call records its own arguments, so a mismatch between calls is visible", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_detect_outliers_olink(ol, n_pcs = 10, threshold = 3))
  ol <- quiet(HADES_detect_outliers_olink(ol, threshold = 3, filter = TRUE))  # default n_pcs

  p1 <- ol$qc_log[[3]]$params
  p2 <- ol$qc_log[[4]]$params
  expect_equal(p1$n_pcs, 10)
  expect_equal(p2$n_pcs, 5)                 # what the removal really used
  expect_false(identical(p1$n_pcs, p2$n_pcs))
})

test_that("a requested n_pcs larger than the data allows is recorded as requested and used", {
  ol <- quiet(HADES_filter_olink(load_sim_olink()))
  ol <- quiet(HADES_detect_outliers_olink(ol, n_pcs = 500))

  p <- ol$qc_log[[3]]$params
  expect_equal(p$n_pcs_requested, 500L)
  expect_lt(p$n_pcs, 500)
})

# ---- the accessor -------------------------------------------------------------

test_that("HADES_get_olink_qc_log returns the filter chain as a tidy table", {
  ol <- run_olink_chain(load_sim_olink())
  st <- HADES_get_olink_qc_log(ol)

  expect_s3_class(st, "data.frame")
  expect_equal(st$step, c("ELEUTHIA_load_olink", "HADES_filter_olink",
                          "HADES_filter_olink_proteins", "HADES_filter_olink_assay_warn",
                          "HADES_detect_outliers_olink", "HADES_detect_outliers_olink"))
  expect_equal(st$action, c("loaded", "removed", "removed", "removed", "flagged", "removed"))
  expect_equal(st$step_id, 1:6)
  expect_true(all(c("parameters", "samples_before", "samples_after", "proteins_before",
                    "proteins_after", "n_samples_removed", "n_proteins_removed",
                    "n_samples_flagged", "timestamp") %in% names(st)))
  expect_match(st$parameters[2], "keep_controls=FALSE")
  expect_equal(st$n_samples_removed, c(1L, 4L, 0L, 0L, 0L, 1L))
  expect_equal(st$n_proteins_removed, c(0L, 0L, 2L, 1L, 0L, 0L))
  expect_equal(st$n_samples_flagged, c(0L, 0L, 0L, 0L, 1L, 0L))
})

test_that("the log is internally consistent and ends at the object's real dimensions", {
  ol <- run_olink_chain(load_sim_olink())
  st <- HADES_get_olink_qc_log(ol)

  # each step starts where the previous one ended
  # (the load entry records the pre-restriction sample count, so start at step 2)
  n <- nrow(st)
  expect_equal(st$samples_before[3:n],  st$samples_after[2:(n - 1)])
  expect_equal(st$proteins_before[3:n], st$proteins_after[2:(n - 1)])
  expect_equal(st$samples_after[1], st$samples_before[2])
  # and the last step matches the data
  expect_equal(st$samples_after[n],  ncol(ol$wide))
  expect_equal(st$proteins_after[n], nrow(ol$wide))
  expect_equal(c(ncol(ol$wide), nrow(ol$wide)), c(19L, 27L))
})

test_that("the removed and flagged tables list every ID with the step that removed it", {
  ol <- run_olink_chain(load_sim_olink())

  rs <- HADES_get_olink_qc_log(ol, "removed_samples")
  expect_named(rs, c("step_id", "step", "sample_id", "reason", "value"))
  expect_setequal(rs$sample_id, c("X1", "N1", "N2", "S21", "S22", "S05"))
  expect_equal(rs$step[rs$sample_id == "X1"], "ELEUTHIA_load_olink")
  expect_equal(rs$step[rs$sample_id == "S05"], "HADES_detect_outliers_olink")

  rp <- HADES_get_olink_qc_log(ol, "removed_proteins")
  expect_named(rp, c("step_id", "step", "OlinkID", "Assay", "reason", "value"))
  expect_setequal(rp$OlinkID, c("OID001", "OID002", "OID003"))
  expect_equal(rp$step[rp$OlinkID == "OID003"], "HADES_filter_olink_assay_warn")

  fl <- HADES_get_olink_qc_log(ol, "flagged_samples")
  expect_equal(fl$sample_id, "S05")
  expect_equal(fl$step_id, 5L)
})

test_that("an empty removed table keeps its columns", {
  ol <- quiet(HADES_filter_olink(load_sim_olink(), keep_controls = TRUE, filter_qc = FALSE))
  rp <- HADES_get_olink_qc_log(ol, "removed_proteins")

  expect_equal(nrow(rp), 0L)
  expect_named(rp, c("step_id", "step", "OlinkID", "Assay", "reason", "value"))
})

test_that("the accessor fails clearly on a non-olink object or one without a log", {
  expect_error(HADES_get_olink_qc_log(list()), "olink_data")
  ol <- load_sim_olink(); ol$qc_log <- NULL
  expect_error(HADES_get_olink_qc_log(ol), "no \\$qc_log")
  expect_error(HADES_get_olink_qc_log(load_sim_olink(), what = "nope"), "should be one of")
})

test_that("the log does not change the data: filtering results are the same as before", {
  ol <- run_olink_chain(load_sim_olink())

  expect_setequal(rownames(ol$wide), setdiff(sprintf("OID%03d", 1:30),
                                             c("OID001", "OID002", "OID003")))
  expect_false(anyNA(ol$sample_meta$SampleID))
  expect_equal(nrow(ol$data), nrow(ol$wide) * ncol(ol$wide))
  expect_equal(nrow(ol$assay_meta), nrow(ol$wide))
})

test_that("print.olink_data reports the number of recorded steps", {
  ol <- run_olink_chain(load_sim_olink())
  out <- utils::capture.output(print(ol))

  expect_true(any(grepl("QC log\\s*: 6 step", out)))
})

test_that("a step on an object without a log starts one rather than failing", {
  ol <- load_sim_olink(); ol$qc_log <- NULL
  ol <- quiet(HADES_filter_olink_proteins(ol, max_na_fraction = 0.2))

  expect_length(ol$qc_log, 1)
  expect_equal(ol$qc_log[[1]]$step, "HADES_filter_olink_proteins")
})
