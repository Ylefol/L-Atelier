# Chromosome sizes from a samtools .fai index. The FASTA route calls samtools via
# the conda environment and is not tested here; the .fai route is plain R.

write_fai <- function(path) {
  utils::write.table(data.frame(c("chr1", "chr2", "chrX"), c(1000L, 800L, 500L),
                                c(6L, 1020L, 1840L), c(60L, 60L, 60L), c(61L, 61L, 61L)),
                     path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  path
}

test_that("get_chrom_sizes reads a .fai into chr / size columns", {
  d <- tempfile("fai_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fai <- write_fai(file.path(d, "genome.fa.fai"))
  s <- quiet(HORIZON_get_chrom_sizes(fai, cache_dir = file.path(d, "cache")))

  expect_equal(names(s), c("chr", "size"))
  expect_equal(s$chr, c("chr1", "chr2", "chrX"))
  expect_equal(s$size, c(1000L, 800L, 500L))
  expect_type(s$size, "integer")
})

test_that("get_chrom_sizes can restrict to chosen chromosomes", {
  d <- tempfile("fai_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fai <- write_fai(file.path(d, "genome.fa.fai"))
  s <- quiet(HORIZON_get_chrom_sizes(fai, chromosomes = c("chrX", "chr1"),
                                     cache_dir = file.path(d, "cache")))

  expect_equal(s$chr, c("chr1", "chrX"))        # file order, not request order
  expect_equal(s$size, c(1000L, 500L))
})

test_that("get_chrom_sizes caches the full table and reuses it", {
  d <- tempfile("fai_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fai   <- write_fai(file.path(d, "genome.fa.fai"))
  cache <- file.path(d, "cache")
  quiet(HORIZON_get_chrom_sizes(fai, chromosomes = "chr1", cache_dir = cache))

  rds <- file.path(cache, "genome.fa_chrom_sizes.rds")
  expect_true(file.exists(rds))
  expect_equal(nrow(readRDS(rds)), 3L)           # the cache holds all chromosomes, not the subset

  # change the source: the cached table is what comes back until force = TRUE
  utils::write.table(data.frame("chrNew", 5L, 6L, 60L, 61L), fai, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  expect_equal(nrow(quiet(HORIZON_get_chrom_sizes(fai, cache_dir = cache))), 3L)
  expect_equal(quiet(HORIZON_get_chrom_sizes(fai, cache_dir = cache, force = TRUE))$chr, "chrNew")
})

test_that("get_chrom_sizes rebuilds the cache with force = TRUE and honours cache_name", {
  d <- tempfile("fai_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fai   <- write_fai(file.path(d, "genome.fa.fai"))
  cache <- file.path(d, "cache")
  quiet(HORIZON_get_chrom_sizes(fai, cache_dir = cache, cache_name = "mine"))
  expect_true(file.exists(file.path(cache, "mine.rds")))

  saveRDS(data.frame(chr = "stale", size = 1L), file.path(cache, "mine.rds"))
  expect_equal(quiet(HORIZON_get_chrom_sizes(fai, cache_dir = cache, cache_name = "mine"))$chr, "stale")
  fresh <- quiet(HORIZON_get_chrom_sizes(fai, cache_dir = cache, cache_name = "mine", force = TRUE))
  expect_equal(fresh$chr, c("chr1", "chr2", "chrX"))
})

test_that("get_chrom_sizes errors on a missing file", {
  expect_error(HORIZON_get_chrom_sizes(file.path(tempdir(), "nope.fa.fai")), "File not found")
})
