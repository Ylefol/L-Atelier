# Tests for sample sheet creation and validation. These need no external tools
# and no real sequencing data: file-existence checks use temporary empty files.

make_sheet <- function(...) {
  ss <- data.frame(
    sample_id    = c("s1", "s2"),
    fastq_r1     = c("a_R1.fq", "b_R1.fq"),
    fastq_r2     = c("a_R2.fq", "b_R2.fq"),
    paired_end   = c(TRUE, TRUE),
    output_dir   = c("out", "out"),
    strandedness = c("unstranded", "reverse"),
    condition    = c("ctrl", "treat"),
    stringsAsFactors = FALSE
  )
  overrides <- list(...)
  for (nm in names(overrides)) ss[[nm]] <- overrides[[nm]]
  ss
}

test_that("create_sample_sheet returns a template with the required columns", {
  ss <- HORIZON_create_sample_sheet(n_samples = 3)
  expect_s3_class(ss, "data.frame")
  expect_equal(nrow(ss), 3)
  expect_true(all(c("sample_id", "fastq_r1", "fastq_r2", "paired_end",
                    "output_dir", "strandedness", "condition") %in% colnames(ss)))
  expect_equal(ss$sample_id, c("sample_1", "sample_2", "sample_3"))
})

test_that("create_sample_sheet writes a CSV that can be read back", {
  out <- tempfile(fileext = ".csv")
  on.exit(unlink(out), add = TRUE)
  capture.output(HORIZON_create_sample_sheet(output_path = out, n_samples = 2))
  expect_true(file.exists(out))
  back <- read.csv(out, stringsAsFactors = FALSE)
  expect_equal(nrow(back), 2)
})

test_that("validate_sample_sheet accepts a valid sheet (check_files = FALSE)", {
  expect_output(res <- HORIZON_validate_sample_sheet(make_sheet(), check_files = FALSE),
                "validated")
  expect_equal(nrow(res), 2)
})

test_that("validate_sample_sheet rejects missing required columns", {
  ss <- make_sheet()
  ss$strandedness <- NULL
  expect_error(HORIZON_validate_sample_sheet(ss, check_files = FALSE),
               "Missing required columns")
})

test_that("validate_sample_sheet rejects duplicate sample ids", {
  expect_error(
    HORIZON_validate_sample_sheet(make_sheet(sample_id = c("s1", "s1")),
                                  check_files = FALSE),
    "Duplicate sample_id"
  )
})

test_that("validate_sample_sheet rejects invalid strandedness", {
  expect_error(
    HORIZON_validate_sample_sheet(make_sheet(strandedness = c("unstranded", "sideways")),
                                  check_files = FALSE),
    "Invalid strandedness"
  )
})

test_that("validate_sample_sheet rejects bad input types and missing files", {
  expect_error(HORIZON_validate_sample_sheet(42), "file path string or a data.frame")
  expect_error(HORIZON_validate_sample_sheet("does_not_exist.csv"), "not found")
})

test_that("validate_sample_sheet checks that FASTQ files exist when asked", {
  # Non-existent paths must fail when check_files = TRUE
  expect_error(HORIZON_validate_sample_sheet(make_sheet(), check_files = TRUE),
               "R1 FASTQ file")

  # Existing (empty) files must pass
  d  <- tempfile("fq_")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  r1 <- file.path(d, c("a_R1.fq", "b_R1.fq"))
  r2 <- file.path(d, c("a_R2.fq", "b_R2.fq"))
  file.create(c(r1, r2))
  ss <- make_sheet(fastq_r1 = r1, fastq_r2 = r2)
  expect_output(HORIZON_validate_sample_sheet(ss, check_files = TRUE), "validated")
})
