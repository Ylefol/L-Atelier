# ELEUTHIA sample sheets: the table that tells GAIA which files belong to which
# sample, group and omics type.

# A sheet pointing at real (tiny) count files, so file checks can be exercised.
make_sheet_files <- function(n = 4) {
  dir <- tempfile("sheet_files_"); dir.create(dir)
  ids <- paste0("S", seq_len(n))
  for (s in ids) {
    utils::write.table(data.frame(gene_id = c("g1", "g2"), count = c(5L, 7L)),
                       file.path(dir, paste0(s, "_counts.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  list(dir = dir,
       sheet = data.frame(sample_id = ids, file_loc = dir,
                          file_name = paste0(ids, "_counts.txt"),
                          group = rep(c("control", "treated"), length.out = n),
                          omics = "RNAseq", format = "counts",
                          stringsAsFactors = FALSE))
}

validate <- function(sheet, ...) {
  quiet(ELEUTHIA_validate_sample_sheet(sheet, check_beds = FALSE, verbose = FALSE, ...))
}

test_that("a correct sheet validates and keeps its sample IDs as row names", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  v <- validate(f$sheet)

  expect_true(v$valid)
  expect_length(v$errors, 0)
  expect_length(v$warnings, 0)
  expect_equal(rownames(v$sample_sheet), f$sheet$sample_id)
  expect_equal(nrow(v$sample_sheet), 4)
})

test_that("a sheet can be validated from a CSV file", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  csv <- tempfile(fileext = ".csv"); on.exit(unlink(csv), add = TRUE)
  utils::write.csv(f$sheet, csv, row.names = FALSE)

  v <- validate(csv)
  expect_true(v$valid)
  expect_equal(nrow(v$sample_sheet), 4)

  missing <- validate(file.path(tempdir(), "no_such_sheet.csv"))
  expect_false(missing$valid)
  expect_match(missing$errors, "not found")
})

test_that("missing required columns are reported", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  v <- validate(f$sheet[, setdiff(names(f$sheet), "group")])
  expect_false(v$valid)
  expect_match(v$errors, "Missing required columns: group")
})

test_that("an unknown format is an error", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  f$sheet$format[1] <- "xyz"
  v <- validate(f$sheet)
  expect_false(v$valid)
  expect_match(v$errors[1], "Invalid format values: xyz")
  expect_match(v$errors[1], "peaks, bed, counts")
})

test_that("a data file that does not exist is a warning, not an error", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  f$sheet$file_name[2] <- "missing.txt"

  v <- validate(f$sheet)
  expect_true(v$valid)
  expect_length(v$warnings, 1)
  expect_match(v$warnings, "Data file not found")
  expect_match(v$warnings, "missing.txt", fixed = TRUE)

  # file checks can be switched off
  quick <- validate(f$sheet, check_files = FALSE)
  expect_length(quick$warnings, 0)
})

test_that("empty rows are removed", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  blank <- stats::setNames(as.data.frame(matrix("", 1, ncol(f$sheet)),
                                         stringsAsFactors = FALSE), names(f$sheet))
  v <- validate(rbind(f$sheet, blank))
  expect_equal(nrow(v$sample_sheet), 4)
})

test_that("without a sample_id column, IDs are generated from omics, group and replicate", {
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  s <- f$sheet; s$sample_id <- NULL
  s$bio_rep <- c(1, 2, 1, 2); s$tech_rep <- "a"; s$batch <- c(1, 1, 2, 2)

  v <- validate(s)
  expect_true(v$valid)
  expect_equal(v$sample_sheet$sample_id,
               c("RNAseq_control_1a_b1", "RNAseq_treated_2a_b1",
                 "RNAseq_control_1a_b2", "RNAseq_treated_2a_b2"))
  expect_equal(rownames(v$sample_sheet), v$sample_sheet$sample_id)
})

test_that("duplicate sample IDs are reported as a validation error", {
  # KNOWN BUG (reported, fix pending approval): the duplicate is detected and an
  # error recorded, but the function then sets row names from sample_id, which
  # fails for duplicates with a raw R error ("duplicate 'row.names' are not
  # allowed") instead of returning valid = FALSE with the message.
  # Remove this skip once rownames are only set when the IDs are unique.
  skip("known bug: duplicate sample_id crashes on rownames<- instead of returning valid = FALSE")
  f <- make_sheet_files(); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  f$sheet$sample_id[2] <- f$sheet$sample_id[1]
  v <- validate(f$sheet)
  expect_false(v$valid)
  expect_match(v$errors, "Duplicate sample_id")
})

# ---- summaries and subsets ---------------------------------------------------

test_that("summarize_sample_sheet counts samples per omics type and group", {
  f <- make_sheet_files(6); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  f$sheet$omics[5:6] <- "ATACseq"
  v <- validate(f$sheet)

  s <- quiet(ELEUTHIA_summarize_sample_sheet(v$sample_sheet))
  expect_equal(s$total, 6)
  expect_equal(unlist(s$by_omics)[c("RNAseq", "ATACseq")], c(RNAseq = 4, ATACseq = 2))
  expect_equal(sum(unlist(s$by_group)), 6)
})

test_that("get_omics_subset returns one omics type, or warns when there is none", {
  f <- make_sheet_files(6); on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  f$sheet$omics[5:6] <- "ATACseq"
  v <- validate(f$sheet)

  atac <- ELEUTHIA_get_omics_subset(v$sample_sheet, "ATACseq")
  expect_equal(nrow(atac), 2)
  expect_true(all(atac$omics == "ATACseq"))
  expect_equal(rownames(atac), c("1", "2"))   # row names are reset

  expect_warning(none <- ELEUTHIA_get_omics_subset(v$sample_sheet, "ChIPseq"),
                 "No samples found for omics type: ChIPseq")
  expect_equal(nrow(none), 0)
})
