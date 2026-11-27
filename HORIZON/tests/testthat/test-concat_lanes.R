# Concatenating multi-lane FASTQs. gzip members can be appended byte for byte, so
# the merged file must decompress to the lanes in order; that is what is checked.

write_gz_fastq <- function(path, reads) {
  con <- gzfile(path, "wb"); on.exit(close(con))
  writeLines(unlist(lapply(reads, function(r) c(paste0("@", r), "ACGT", "+", "IIII"))), con)
  path
}
read_gz <- function(path) readLines(gzfile(path))

lane_fixture <- function(dir) {
  list(
    r1 = c(write_gz_fastq(file.path(dir, "L1_R1.fastq.gz"), c("a1", "a2")),
           write_gz_fastq(file.path(dir, "L2_R1.fastq.gz"), c("b1"))),
    r2 = c(write_gz_fastq(file.path(dir, "L1_R2.fastq.gz"), c("a1", "a2")),
           write_gz_fastq(file.path(dir, "L2_R2.fastq.gz"), c("b1")))
  )
}

test_that("concat_lanes merges paired-end lanes in order, per read direction", {
  d <- tempfile("lanes_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- lane_fixture(d)
  res <- quiet(HORIZON_concat_lanes(f$r1, f$r2, output_dir = file.path(d, "out"), sample_id = "S1"))

  expect_equal(basename(res$fastq_r1), "S1_R1.fastq.gz")
  expect_equal(basename(res$fastq_r2), "S1_R2.fastq.gz")
  expect_equal(read_gz(res$fastq_r1)[c(1, 5, 9)], c("@a1", "@a2", "@b1"))   # lane order kept
  expect_equal(length(read_gz(res$fastq_r1)), 12L)                          # 3 reads x 4 lines
  expect_equal(read_gz(res$fastq_r2), read_gz(res$fastq_r1))
})

test_that("concat_lanes handles single-end data and returns NULL for R2", {
  d <- tempfile("lanes_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- lane_fixture(d)
  res <- quiet(HORIZON_concat_lanes(f$r1, output_dir = file.path(d, "out"), sample_id = "S1"))

  expect_null(res$fastq_r2)
  expect_false(file.exists(file.path(d, "out", "S1_R2.fastq.gz")))
  expect_equal(length(read_gz(res$fastq_r1)), 12L)
})

test_that("concat_lanes does not overwrite existing output unless force = TRUE", {
  d <- tempfile("lanes_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- lane_fixture(d)
  out <- file.path(d, "out")
  res <- quiet(HORIZON_concat_lanes(f$r1, output_dir = out, sample_id = "S1"))
  writeLines("sentinel", res$fastq_r1)

  quiet(HORIZON_concat_lanes(f$r1, output_dir = out, sample_id = "S1"))
  expect_equal(readLines(res$fastq_r1), "sentinel")               # skipped, not overwritten

  quiet(HORIZON_concat_lanes(f$r1, output_dir = out, sample_id = "S1", force = TRUE))
  expect_equal(length(read_gz(res$fastq_r1)), 12L)                # rebuilt
})

test_that("concat_lanes validates its inputs", {
  d <- tempfile("lanes_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- lane_fixture(d)
  out <- file.path(d, "out")

  expect_error(HORIZON_concat_lanes(character(0), output_dir = out, sample_id = "S"), "at least one")
  expect_error(HORIZON_concat_lanes(c(f$r1, file.path(d, "nope.fastq.gz")), output_dir = out,
                                    sample_id = "S"), "not found")
  expect_error(HORIZON_concat_lanes(f$r1, f$r2[1], output_dir = out, sample_id = "S"), "same length")
  expect_error(HORIZON_concat_lanes(f$r1, c(f$r2[1], file.path(d, "nope_R2.fastq.gz")),
                                    output_dir = out, sample_id = "S"), "R2 FASTQ")
})
