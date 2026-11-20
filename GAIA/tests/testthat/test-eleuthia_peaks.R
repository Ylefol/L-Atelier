# ELEUTHIA: reading peak and BED files and building a consensus peak set, the
# basis of every ChIP / ATAC / CUT&Tag analysis. Files are written in the same
# formats HORIZON's peak callers produce.

write_tsv_rows <- function(path, rows) {
  writeLines(vapply(rows, function(r) paste(r, collapse = "\t"), ""), path)
  path
}

# narrowPeak: chrom start end name score strand signalValue pValue qValue peak
make_narrowpeak <- function(dir, name, rows) {
  write_tsv_rows(file.path(dir, paste0(name, ".narrowPeak")), rows)
}

peak_files <- function() {
  d <- tempfile("peaks_"); dir.create(d)
  a <- make_narrowpeak(d, "A", list(c("chr1", 100, 200, "a1", 500, ".", 5.0, 10, 8, 50),
                                    c("chr1", 1000, 1100, "a2", 300, ".", 3.0, 6, 4, 50),
                                    c("chr2", 500, 700, "a3", 800, ".", 9.0, 20, 15, 100)))
  b <- make_narrowpeak(d, "B", list(c("chr1", 150, 250, "b1", 400, ".", 4.0, 8, 6, 50),
                                    c("chr1", 5000, 5100, "b2", 200, ".", 2.0, 4, 3, 50),
                                    c("chr2", 520, 640, "b3", 700, ".", 8.0, 18, 12, 60)))
  list(dir = d, A = a, B = b)
}

load_peaks <- function(pf) {
  list(A = quiet(ELEUTHIA_load_narrowpeak(pf$A)), B = quiet(ELEUTHIA_load_narrowpeak(pf$B)))
}

# ---- narrowPeak and BED files ------------------------------------------------

test_that("load_narrowpeak reads the ten narrowPeak columns and adds the summit position", {
  pf <- peak_files(); on.exit(unlink(pf$dir, recursive = TRUE), add = TRUE)
  p <- quiet(ELEUTHIA_load_narrowpeak(pf$A))

  expect_s3_class(p, "data.frame")
  expect_equal(nrow(p), 3)
  expect_equal(names(p)[1:4], c("chr", "start", "end", "name"))
  expect_true(all(c("signal_value", "pvalue", "qvalue", "peak_summit", "summit_pos") %in% names(p)))
  expect_equal(p$chr, c("chr1", "chr1", "chr2"))
  expect_equal(p$start, c(100, 1000, 500))
  expect_equal(p$summit_pos, p$start + p$peak_summit)      # summit is an offset from the start
  expect_equal(p$qvalue, c(8, 4, 15))
})

test_that("load_narrowpeak reports a missing file", {
  expect_error(quiet(ELEUTHIA_load_narrowpeak(file.path(tempdir(), "no_such.narrowPeak"))),
               "not found", ignore.case = TRUE)
})

test_that("load_bed names the columns by position and rejects a missing file", {
  d <- tempfile("bed_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  bed3 <- write_tsv_rows(file.path(d, "three.bed"), list(c("chr1", 100, 200), c("chr2", 5, 50)))
  bed6 <- write_tsv_rows(file.path(d, "six.bed"),
                         list(c("chr1", 100, 200, "r1", 7, "+"), c("chr2", 5, 50, "r2", 9, "-")))

  b3 <- ELEUTHIA_load_bed(bed3)
  expect_equal(names(b3), c("chr", "start", "end"))
  expect_equal(nrow(b3), 2)

  b6 <- ELEUTHIA_load_bed(bed6)
  expect_equal(names(b6), c("chr", "start", "end", "name", "score", "strand"))
  expect_equal(b6$strand, c("+", "-"))

  expect_error(ELEUTHIA_load_bed(file.path(d, "none.bed")), "File not found")
})

test_that("load_bed ignores comment lines", {
  d <- tempfile("bed_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- file.path(d, "c.bed")
  writeLines(c("# header comment", "chr1\t1\t10", "chr1\t20\t30"), f)
  expect_equal(nrow(ELEUTHIA_load_bed(f)), 2)
})

# ---- consensus peaks ---------------------------------------------------------

test_that("overlapping peaks from different samples are merged into one consensus peak", {
  pf <- peak_files(); on.exit(unlink(pf$dir, recursive = TRUE), add = TRUE)
  cons <- quiet(ELEUTHIA_create_consensus_peaks(load_peaks(pf)))

  # chr1 100-200 + 150-250 -> 100-250 ; chr1 1000-1100 ; chr1 5000-5100 ; chr2 500-700 + 520-640 -> 500-700
  expect_equal(nrow(cons), 4)
  expect_true(all(c("chr", "start", "end", "peak_id", "n_samples") %in% names(cons)))

  m <- cons[cons$chr == "chr1" & cons$start == 100, ]
  expect_equal(m$end, 250)
  expect_equal(m$n_samples, 2)

  chr2 <- cons[cons$chr == "chr2", ]
  expect_equal(c(chr2$start, chr2$end), c(500, 700))     # the nested peak does not extend it
  expect_equal(chr2$n_samples, 2)

  expect_equal(sort(cons$n_samples), c(1, 1, 2, 2))
  expect_false(anyDuplicated(cons$peak_id) > 0)
})

test_that("min_overlap keeps only peaks seen in enough samples", {
  pf <- peak_files(); on.exit(unlink(pf$dir, recursive = TRUE), add = TRUE)
  both <- quiet(ELEUTHIA_create_consensus_peaks(load_peaks(pf), min_overlap = 2))

  expect_equal(nrow(both), 2)
  expect_true(all(both$n_samples == 2))
})

test_that("merge_distance joins peaks that are close but do not overlap", {
  pf <- peak_files(); on.exit(unlink(pf$dir, recursive = TRUE), add = TRUE)
  far  <- quiet(ELEUTHIA_create_consensus_peaks(load_peaks(pf), merge_distance = 0))
  near <- quiet(ELEUTHIA_create_consensus_peaks(load_peaks(pf), merge_distance = 4000))

  expect_equal(nrow(far), 4)
  expect_equal(nrow(near), 2)                    # on chr1 everything within 4 kb becomes one peak
  chr1 <- near[near$chr == "chr1", ]
  expect_equal(c(chr1$start, chr1$end), c(100, 5100))
})

test_that("the same peak twice in one sample counts as one sample", {
  p <- data.frame(chr = c("chr1", "chr1"), start = c(100, 120), end = c(200, 220))
  cons <- quiet(ELEUTHIA_create_consensus_peaks(list(only = p)))
  expect_equal(nrow(cons), 1)
  expect_equal(cons$n_samples, 1)
})

test_that("an empty peak list is rejected", {
  expect_error(quiet(ELEUTHIA_create_consensus_peaks(list())), "peak_list is empty")
})

test_that("consensus peaks convert to the 1-based SAF format used for counting", {
  pf <- peak_files(); on.exit(unlink(pf$dir, recursive = TRUE), add = TRUE)
  cons <- quiet(ELEUTHIA_create_consensus_peaks(load_peaks(pf)))
  saf <- ELEUTHIA_peaks_to_saf(cons)

  expect_named(saf, c("GeneID", "Chr", "Start", "End", "Strand"))
  expect_equal(saf$GeneID, cons$peak_id)
  expect_equal(saf$Start, cons$start + 1)         # BED is 0-based, SAF is 1-based
  expect_equal(saf$End, cons$end)
  expect_true(all(saf$Strand == "*"))
})

# ---- expanding regions -------------------------------------------------------

test_that("expand_regions widens each region symmetrically and clips at zero", {
  r <- data.frame(chr = c("chr1", "chr1"), start = c(1000, 50), end = c(1100, 80))

  same <- quiet(ELEUTHIA_expand_regions(r, window_size = 0, verbose = FALSE))
  expect_equal(same$start, r$start)

  out <- quiet(ELEUTHIA_expand_regions(r, window_size = 100, verbose = FALSE))
  expect_equal(out$start[1], 900)
  expect_equal(out$end[1], 1200)
  expect_equal(out$start[2], 0)                   # 50 - 100 would be negative
  expect_equal(out$end[2], 180)
})
