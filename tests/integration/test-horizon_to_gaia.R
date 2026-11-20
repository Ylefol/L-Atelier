# HORIZON -> GAIA: what HORIZON writes must be what GAIA reads.
#
# Two routes lead from the same per-sample count files to a differential
# analysis, and both are tested:
#   A. HORIZON_aggregate_counts()  ->  count matrix + sample metadata  ->  ARTEMIS
#   B. GAIA's own sample sheet     ->  ELEUTHIA_load_rnaseq_from_sheet() -> ARTEMIS
# and they must agree.

new_run <- function(skip = character(0)) {
  sim  <- simulate_counts()
  root <- file.path(tempfile("horizon_run_"))
  dir.create(root)
  ss <- write_horizon_run(sim, root, skip = skip)
  list(sim = sim, root = root, sheet = ss)
}

aggregate <- function(run) {
  suppressWarnings(quiet(HORIZON_aggregate_counts(run$sheet, output_root = run$root)))
}

agg_dir <- function(run) file.path(run$root, "aggregated", "counts")

test_that("the HORIZON sample sheet used here is valid", {
  run <- new_run()
  ss <- quiet(HORIZON_validate_sample_sheet(run$sheet, check_files = FALSE))
  expect_equal(nrow(ss), 8)
})

# ---- Route A -----------------------------------------------------------------

test_that("HORIZON_aggregate_counts writes the matrix and the sample metadata", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  mat <- aggregate(run)

  expect_true(file.exists(file.path(agg_dir(run), "count_matrix.csv")))
  expect_true(file.exists(file.path(agg_dir(run), "count_matrix.rds")))
  expect_true(file.exists(file.path(agg_dir(run), "sample_metadata.csv")))

  # genes x samples, identical to what was written per sample
  expect_equal(dim(mat), dim(run$sim$counts))
  expect_identical(rownames(mat), rownames(run$sim$counts))
  expect_identical(colnames(mat), colnames(run$sim$counts))
  expect_equal(unname(mat), unname(run$sim$counts))
})

test_that("the saved .rds and .csv matrices agree with each other", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  mat <- aggregate(run)

  from_rds <- readRDS(file.path(agg_dir(run), "count_matrix.rds"))
  from_csv <- as.matrix(utils::read.csv(file.path(agg_dir(run), "count_matrix.csv"),
                                        row.names = 1, check.names = FALSE))

  expect_equal(from_rds, mat)
  expect_equal(from_csv, mat)
})

test_that("route A: aggregated counts + metadata go straight into ARTEMIS and recover the planted genes", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  aggregate(run)

  counts  <- readRDS(file.path(agg_dir(run), "count_matrix.rds"))
  targets <- targets_from_metadata(file.path(agg_dir(run), "sample_metadata.csv"))
  dea <- run_dea(counts, targets, group_col = "condition")

  sig <- dea$results$feature_id[!is.na(dea$results$padj) & dea$results$padj < 0.05]
  expect_gte(sum(run$sim$up %in% sig), 8)
  expect_gte(sum(run$sim$down %in% sig), 8)
  expect_lte(sum(run$sim$null %in% sig), 5)
  expect_true(all(dea$results$log2FoldChange[dea$results$feature_id %in% run$sim$up] > 0))
})

test_that("the metadata CSV carries the sample sheet columns the analysis needs", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  aggregate(run)

  meta <- targets_from_metadata(file.path(agg_dir(run), "sample_metadata.csv"))
  expect_true(all(c("sample_id", "condition", "batch") %in% colnames(meta)))
  expect_identical(rownames(meta), colnames(run$sim$counts))
  expect_identical(meta$condition, run$sim$condition)
})

test_that("a batch column from the metadata can be used in the design", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  aggregate(run)

  counts  <- readRDS(file.path(agg_dir(run), "count_matrix.rds"))
  targets <- targets_from_metadata(file.path(agg_dir(run), "sample_metadata.csv"))
  norm <- quiet(ARTEMIS_normalize_counts(counts, targets, group_col = "condition",
                                         batch_col = "batch", verbose = FALSE))
  dea <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))

  expect_equal(norm$parameters$batch_col, "batch")
  expect_equal(nrow(dea$results), nrow(counts))
})

test_that("a sample whose count file is missing is skipped, and the rest still flows into GAIA", {
  run <- new_run(skip = "sample_8")
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)

  expect_warning(
    quiet(HORIZON_aggregate_counts(run$sheet, output_root = run$root)),
    "Count file not found for sample 'sample_8'"
  )

  counts  <- readRDS(file.path(agg_dir(run), "count_matrix.rds"))
  targets <- targets_from_metadata(file.path(agg_dir(run), "sample_metadata.csv"))

  # metadata and matrix stay consistent: the skipped sample is in neither
  expect_equal(ncol(counts), 7)
  expect_false("sample_8" %in% rownames(targets))
  expect_setequal(rownames(targets), colnames(counts))

  dea <- run_dea(counts, targets, group_col = "condition")
  expect_equal(nrow(dea$results), nrow(counts))
})

test_that("aggregation stops with a clear message when no sample has counts", {
  run <- new_run(skip = colnames(simulate_counts()$counts))
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)
  expect_error(
    suppressWarnings(quiet(HORIZON_aggregate_counts(run$sheet, output_root = run$root))),
    "No count files found"
  )
})

# ---- Route B -----------------------------------------------------------------

gaia_sheet <- function(run) {
  data.frame(
    sample_id = run$sheet$sample_id,
    file_loc  = file.path(run$sheet$output_dir, run$sheet$sample_id, "counts"),
    file_name = paste0(run$sheet$sample_id, "_counts.txt"),
    group     = run$sheet$condition,
    batch     = run$sheet$batch,
    omics     = "RNAseq",
    format    = "counts",
    stringsAsFactors = FALSE
  )
}

test_that("GAIA's sample sheet finds HORIZON's count files", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)

  v <- quiet(ELEUTHIA_validate_sample_sheet(gaia_sheet(run), check_files = TRUE,
                                            check_beds = FALSE, verbose = FALSE))
  expect_true(v$valid)
  expect_length(v$errors, 0)
})

test_that("route B: GAIA loads HORIZON's per-sample count files directly", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)

  v   <- quiet(ELEUTHIA_validate_sample_sheet(gaia_sheet(run), check_files = TRUE,
                                              check_beds = FALSE, verbose = FALSE))
  rna <- quiet(ELEUTHIA_load_rnaseq_from_sheet(v$sample_sheet, verbose = FALSE))

  expect_equal(dim(rna$counts), dim(run$sim$counts))
  expect_equal(unname(rna$counts), unname(run$sim$counts))
  expect_identical(rownames(rna$counts), rownames(run$sim$counts))
  expect_identical(rownames(rna$targets), colnames(rna$counts))
  expect_true("group" %in% colnames(rna$targets))
})

test_that("both routes produce the same matrix and the same differential results", {
  run <- new_run()
  on.exit(unlink(run$root, recursive = TRUE), add = TRUE)

  # route A
  mat_a   <- aggregate(run)
  targets <- targets_from_metadata(file.path(agg_dir(run), "sample_metadata.csv"))
  dea_a   <- run_dea(mat_a, targets, group_col = "condition")

  # route B
  v     <- quiet(ELEUTHIA_validate_sample_sheet(gaia_sheet(run), check_files = TRUE,
                                                check_beds = FALSE, verbose = FALSE))
  rna_b <- quiet(ELEUTHIA_load_rnaseq_from_sheet(v$sample_sheet, verbose = FALSE))
  dea_b <- run_dea(rna_b$counts, rna_b$targets, group_col = "group")

  expect_equal(unname(mat_a), unname(rna_b$counts))
  expect_equal(dea_a$results, dea_b$results)
})
