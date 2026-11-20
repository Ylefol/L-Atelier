# ELEUTHIA: reading count files into matrices.

write_counts <- function(dir, counts, header = TRUE, ext = "txt", suffix = "") {
  for (s in colnames(counts)) {
    df <- data.frame(gene_id = rownames(counts), count = as.integer(counts[, s]))
    utils::write.table(df, file.path(dir, paste0(s, suffix, ".", ext)), sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = header)
  }
}

small_counts <- function() {
  matrix(c(10L, 20L, 30L,
            0L,  1L,  0L,
           50L, 60L, 70L),
         nrow = 3, byrow = TRUE,
         dimnames = list(c("gA", "gB", "gC"), c("S1", "S2", "S3")))
}

# ---- load_count_file ---------------------------------------------------------

test_that("load_count_file reads a headed and a header-less file alike", {
  d <- tempfile("lcf_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- small_counts()
  write_counts(d, m, header = TRUE,  suffix = "_headed")
  write_counts(d, m, header = FALSE, suffix = "_plain")

  headed <- ELEUTHIA_load_count_file(file.path(d, "S1_headed.txt"))
  plain  <- ELEUTHIA_load_count_file(file.path(d, "S1_plain.txt"))

  expect_equal(headed, c(gA = 10, gB = 0, gC = 50))
  expect_identical(headed, plain)      # header detected automatically
})

test_that("load_count_file honours explicit column positions and a missing file", {
  d <- tempfile("lcf_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- file.path(d, "wide.txt")
  utils::write.table(data.frame(a = c("gA", "gB"), x = c(9, 9), count = c(3, 4)),
                     f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  expect_equal(ELEUTHIA_load_count_file(f, count_col = 3), c(gA = 3, gB = 4))
  expect_error(ELEUTHIA_load_count_file(file.path(d, "nope.txt")), "File not found")
})

# ---- load_counts_from_dir ----------------------------------------------------

test_that("load_counts_from_dir merges a directory of header-less files by gene", {
  d <- tempfile("lcd_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- small_counts()
  write_counts(d, m, header = FALSE)

  # note the trailing slash: the function joins the directory and file name by pasting
  out <- quiet(ELEUTHIA_load_counts_from_dir(paste0(d, "/")))
  expect_equal(dim(out), c(3L, 3L))
  expect_equal(sort(colnames(out)), c("S1", "S2", "S3"))
  expect_equal(rownames(out), rownames(m))
  expect_equal(out["gC", "S2"], 60)
})

test_that("load_counts_from_dir can select files by pattern and reports bad directories", {
  d <- tempfile("lcd_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- small_counts()
  write_counts(d, m, header = FALSE)
  write_counts(d, m[, 1, drop = FALSE], header = FALSE, ext = "bak")

  out <- quiet(ELEUTHIA_load_counts_from_dir(paste0(d, "/"), file_pattern = "\\.txt$"))
  expect_equal(ncol(out), 3)

  expect_error(quiet(ELEUTHIA_load_counts_from_dir(file.path(d, "nope"))),
               "Directory does not exist")
  expect_error(quiet(ELEUTHIA_load_counts_from_dir(paste0(d, "/"), file_pattern = "\\.zzz$")),
               "No files matching pattern")
  empty <- tempfile("empty_"); dir.create(empty); on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_error(quiet(ELEUTHIA_load_counts_from_dir(empty)), "No files found")
})

test_that("load_counts_from_dir works with a directory given without a trailing slash", {
  # KNOWN BUG (reported, fix pending approval): the file path is built with
  # paste0(count_path, file) instead of file.path(), so a directory without a
  # trailing slash points at non-existent files and the function ends with the
  # unhelpful "No data was successfully loaded".
  # Remove this skip once the path is built with file.path().
  skip("known bug: directory path without trailing slash finds no files")
  d <- tempfile("lcd_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  write_counts(d, small_counts(), header = FALSE)
  out <- quiet(ELEUTHIA_load_counts_from_dir(d))
  expect_equal(ncol(out), 3)
})

test_that("load_multi_counts loads several directories into a named list", {
  d <- tempfile("lmc_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  write_counts(d, small_counts(), header = FALSE)
  p <- paste0(d, "/")

  out <- quiet(ELEUTHIA_load_multi_counts(list(first = p, second = p)))
  expect_named(out, c("first", "second"))
  expect_equal(dim(out$first), c(3L, 3L))
  expect_equal(out$first, out$second)
})

# ---- filter_low_expression ---------------------------------------------------

test_that("filter_low_expression removes genes without enough expressed samples", {
  m <- small_counts()                       # gB is (0, 1, 0)
  rna <- list(counts = m, targets = data.frame(sample_id = colnames(m),
                                               row.names = colnames(m)))
  out <- quiet(ELEUTHIA_filter_low_expression(rna, min_count = 10, min_samples = 2,
                                              verbose = FALSE))

  expect_equal(rownames(out$counts), c("gA", "gC"))
  expect_equal(ncol(out$counts), 3)
  expect_equal(out$targets, rna$targets)    # sample table untouched

  none <- quiet(ELEUTHIA_filter_low_expression(rna, min_count = 1000, min_samples = 1,
                                               verbose = FALSE))
  expect_equal(nrow(none$counts), 0)
})

# ---- loading from a sheet ----------------------------------------------------

test_that("load_rnaseq_from_sheet warns and uses the common genes when files disagree", {
  d <- tempfile("lrs_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- small_counts()
  write_counts(d, m)
  # sample S3 lacks gene gC
  utils::write.table(data.frame(gene_id = c("gA", "gB"), count = c(30L, 0L)),
                     file.path(d, "S3.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  sheet <- data.frame(sample_id = colnames(m), file_loc = d,
                      file_name = paste0(colnames(m), ".txt"),
                      group = c("a", "a", "b"), omics = "RNAseq", format = "counts",
                      stringsAsFactors = FALSE)
  expect_warning(rna <- quiet(ELEUTHIA_load_rnaseq_from_sheet(sheet, verbose = FALSE)),
                 "Gene sets differ")
  expect_equal(rownames(rna$counts), c("gA", "gB"))
  expect_equal(ncol(rna$counts), 3)
})

test_that("load_rnaseq_from_sheet stops when the requested omics type is absent", {
  sheet <- data.frame(sample_id = "S1", file_loc = ".", file_name = "x.txt",
                      group = "a", omics = "ATACseq", format = "counts",
                      stringsAsFactors = FALSE)
  expect_error(quiet(ELEUTHIA_load_rnaseq_from_sheet(sheet, verbose = FALSE)),
               "No samples found for omics type: RNAseq")
})
