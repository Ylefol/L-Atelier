# PARSE Biosciences sample sheet: template creation and validation. File-existence
# checks use temporary empty files, so no sequencing data is needed.

parse_fixture <- function(dir, ids = c("run_1", "run_2")) {
  gd <- file.path(dir, "genome"); dir.create(gd, showWarnings = FALSE)
  r1 <- file.path(dir, paste0(ids, "_R1.fastq.gz")); r2 <- file.path(dir, paste0(ids, "_R2.fastq.gz"))
  file.create(r1, r2)
  data.frame(run_id = ids, fastq_r1 = r1, fastq_r2 = r2, output_dir = dir, genome_dir = gd,
             chemistry = "v3", kit = "WT", stringsAsFactors = FALSE)
}

test_that("create_parse_sheet returns the template, and writes it when asked", {
  t <- quiet(HORIZON_create_parse_sheet(n_runs = 3))
  expect_equal(nrow(t), 3L)
  expect_equal(names(t), c("run_id", "fastq_r1", "fastq_r2", "output_dir", "genome_dir",
                           "chemistry", "kit"))
  expect_equal(t$run_id, c("run_1", "run_2", "run_3"))
  expect_true(all(t$chemistry == "v3") && all(t$kit == "WT"))

  out <- tempfile(fileext = ".csv"); on.exit(unlink(out), add = TRUE)
  quiet(HORIZON_create_parse_sheet(out, n_runs = 2))
  expect_equal(nrow(utils::read.csv(out)), 2L)
})

test_that("validate_parse_sheet accepts a valid sheet, as a data.frame or a CSV path", {
  d <- tempfile("ps_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ss <- parse_fixture(d)

  expect_output(res <- HORIZON_validate_parse_sheet(ss), "2 run\\(s\\) validated")
  expect_equal(res, ss)

  csv <- file.path(d, "sheet.csv"); utils::write.csv(ss, csv, row.names = FALSE)
  expect_equal(quiet(HORIZON_validate_parse_sheet(csv)), ss)
})

test_that("validate_parse_sheet rejects missing columns, duplicate runs and bad chemistry or kit", {
  d <- tempfile("ps_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ss <- parse_fixture(d)

  expect_error(HORIZON_validate_parse_sheet(ss[, setdiff(names(ss), "kit")]), "Missing required columns: kit")
  dup <- ss; dup$run_id[2] <- dup$run_id[1]
  expect_error(HORIZON_validate_parse_sheet(dup), "Duplicate run_id")
  bad_chem <- ss; bad_chem$chemistry[1] <- "v9"
  expect_error(HORIZON_validate_parse_sheet(bad_chem), "Invalid chemistry")
  bad_kit <- ss; bad_kit$kit[2] <- "XL"
  expect_error(HORIZON_validate_parse_sheet(bad_kit), "Invalid kit")
})

test_that("validate_parse_sheet checks that files and the genome directory exist, unless told not to", {
  d <- tempfile("ps_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ss <- parse_fixture(d)

  no_r1 <- ss; no_r1$fastq_r1[1] <- file.path(d, "gone_R1.fastq.gz")
  expect_error(HORIZON_validate_parse_sheet(no_r1), "R1 FASTQ")
  no_r2 <- ss; no_r2$fastq_r2[2] <- file.path(d, "gone_R2.fastq.gz")
  expect_error(HORIZON_validate_parse_sheet(no_r2), "R2 FASTQ")
  no_genome <- ss; no_genome$genome_dir[1] <- file.path(d, "no_genome")
  expect_error(HORIZON_validate_parse_sheet(no_genome), "Genome directory")

  expect_output(HORIZON_validate_parse_sheet(no_r1, check_files = FALSE), "validated")
})

test_that("validate_parse_sheet rejects a missing file path and a wrong input type", {
  expect_error(HORIZON_validate_parse_sheet(file.path(tempdir(), "nope.csv")), "not found")
  expect_error(HORIZON_validate_parse_sheet(42), "file path string or a data.frame")
})
