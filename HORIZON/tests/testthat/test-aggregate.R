# Aggregating per-sample HORIZON outputs into matrices. Inputs are simulated files
# laid out like HORIZON's own output (see helper-outputs.R), so each cell of the
# result has a known value.

# ---- gene counts --------------------------------------------------------------

test_that("aggregate_counts combines per-sample files into a genes x samples matrix", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out)
  m <- quiet(HORIZON_aggregate_counts(sim_sheet(out), output_root = out))

  expect_equal(dim(m), c(5L, 3L))
  expect_equal(colnames(m), c("s1", "s2", "s3"))
  expect_equal(rownames(m), paste0("g", 1:5))
  expect_equal(unname(m["g1", ]), c(11, 21, 31))      # i * 10 + gene index
  expect_equal(unname(m["g5", "s3"]), 35)
})

test_that("aggregate_counts writes the matrix, an RDS copy and the sample metadata", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out)
  m <- quiet(HORIZON_aggregate_counts(sim_sheet(out), output_root = out))
  agg <- file.path(out, "aggregated", "counts")

  expect_true(all(file.exists(file.path(agg, c("count_matrix.csv", "count_matrix.rds",
                                               "sample_metadata.csv")))))
  csv <- utils::read.csv(file.path(agg, "count_matrix.csv"), row.names = 1)
  expect_equal(as.matrix(csv), m, ignore_attr = TRUE)
  expect_equal(readRDS(file.path(agg, "count_matrix.rds")), m)
  meta <- utils::read.csv(file.path(agg, "sample_metadata.csv"))
  expect_equal(meta$sample_id, c("s1", "s2", "s3"))
})

test_that("aggregate_counts can skip the RDS and metadata outputs", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out)
  quiet(HORIZON_aggregate_counts(sim_sheet(out), output_root = out,
                                 save_rds = FALSE, save_metadata = FALSE))
  agg <- file.path(out, "aggregated", "counts")

  expect_true(file.exists(file.path(agg, "count_matrix.csv")))
  expect_false(file.exists(file.path(agg, "count_matrix.rds")))
  expect_false(file.exists(file.path(agg, "sample_metadata.csv")))
})

test_that("aggregate_counts infers output_root from the first sample and says so", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out)
  txt <- utils::capture.output(HORIZON_aggregate_counts(sim_sheet(out)))

  expect_true(file.exists(file.path(out, "aggregated", "counts", "count_matrix.csv")))
  expect_true(any(grepl("output_root inferred from first sample", txt)))
})

test_that("aggregate_counts output is a [HORIZON] headline with indented detail lines", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out)
  txt <- utils::capture.output(HORIZON_aggregate_counts(sim_sheet(out)))

  expect_true(grepl("^\\[HORIZON\\]", txt[1]))
  expect_true(any(grepl("^    Aggregated count matrix: 5 genes x 3 samples$", txt)))
  expect_true(any(grepl("^    Count matrix written to: ", txt)))
  expect_true(all(grepl("^(\\[HORIZON\\]|    )", txt)))    # every line is a headline or indented
  expect_gte(length(txt), 4)                                # each message ends in a newline
})

test_that("a sample without a count file is skipped with a warning", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out, ids = c("s1", "s3"))
  expect_warning(m <- quiet(HORIZON_aggregate_counts(sim_sheet(out), output_root = out)),
                 "s2")

  expect_equal(colnames(m), c("s1", "s3"))
  meta <- utils::read.csv(file.path(out, "aggregated", "counts", "sample_metadata.csv"))
  expect_equal(meta$sample_id, c("s1", "s3"))      # metadata only for aggregated samples
})

test_that("aggregate_counts errors when no sample has a count file", {
  out <- tempfile("agg_"); dir.create(out); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_error(suppressWarnings(quiet(HORIZON_aggregate_counts(sim_sheet(out), output_root = out))),
               "No count files found")
})

test_that("samples with different gene sets are aligned on the union, with NA for the gaps", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_counts(out, ids = "s1", genes = paste0("g", 1:3))
  write_sim_counts(out, ids = "s2", genes = paste0("g", 2:4))
  m <- quiet(HORIZON_aggregate_counts(sim_sheet(out, c("s1", "s2")), output_root = out))

  expect_setequal(rownames(m), paste0("g", 1:4))
  expect_true(is.na(m["g1", "s2"]))
  expect_true(is.na(m["g4", "s1"]))
})

# ---- salmon -------------------------------------------------------------------

test_that("aggregate_salmon builds TPM and NumReads matrices from quant.sf files", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_salmon(out)
  res <- quiet(HORIZON_aggregate_salmon(sim_sheet(out), output_root = out))

  expect_named(res, c("tpm", "counts"))
  expect_equal(dim(res$tpm), c(4L, 3L))
  expect_equal(rownames(res$tpm), paste0("tx", 1:4))
  expect_equal(unname(res$tpm["tx2", ]), c(102, 202, 302))
  expect_equal(unname(res$counts["tx2", ]), c(1002, 2002, 3002))

  agg <- file.path(out, "aggregated", "salmon")
  expect_true(all(file.exists(file.path(agg, c("tpm_matrix.csv", "counts_matrix.csv",
                                               "tpm_matrix.rds", "counts_matrix.rds",
                                               "sample_metadata.csv")))))
})

test_that("aggregate_salmon skips a missing sample with a warning and errors if none exist", {
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_salmon(out, ids = c("s1", "s2"))
  expect_warning(res <- quiet(HORIZON_aggregate_salmon(sim_sheet(out), output_root = out)), "s3")
  expect_equal(colnames(res$tpm), c("s1", "s2"))

  empty <- tempfile("agg_"); dir.create(empty); on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_error(suppressWarnings(quiet(HORIZON_aggregate_salmon(sim_sheet(empty), output_root = empty))),
               "No quant.sf files found")
})

# ---- TEcount ------------------------------------------------------------------

test_that("aggregate_tecount keeps only the TE rows, identified from the TE GTF", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("S4Vectors")
  out <- tempfile("agg_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_tecount(out)
  gtf <- write_te_gtf(file.path(out, "te.gtf"))
  m <- quiet(HORIZON_aggregate_tecount(sim_sheet(out), te_gtf = gtf, output_root = out))

  expect_equal(dim(m), c(2L, 3L))                              # the two TE subfamilies
  expect_setequal(rownames(m), c("AluY:Alu:SINE", "L1PA2:L1:LINE"))
  expect_false(any(grepl("^ENSG", rownames(m))))               # gene rows dropped
  expect_equal(unname(m["AluY:Alu:SINE", ]), c(10, 20, 30))
  expect_equal(unname(m["L1PA2:L1:LINE", ]), c(20, 40, 60))
  expect_true(file.exists(file.path(out, "aggregated", "tecount", "te_count_matrix.csv")))
})

test_that("aggregate_tecount reports a TE GTF without the expected attributes", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("S4Vectors")
  out <- tempfile("agg_"); dir.create(out); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_sim_tecount(out)
  gtf <- file.path(out, "bad.gtf")
  writeLines(paste("chr1", "x", "exon", 1, 100, ".", "+", ".", 'gene_id "a";', sep = "\t"), gtf)

  expect_error(quiet(HORIZON_aggregate_tecount(sim_sheet(out), te_gtf = gtf, output_root = out)),
               "family_id")
})

# ---- merging gene and TE counts ----------------------------------------------

gene_te_fixture <- function() {
  g <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  t <- matrix(11:16, nrow = 2, dimnames = list(c("AluY:Alu:SINE", "L1PA2:L1:LINE"), c("s1", "s2", "s3")))
  list(g = g, t = t)
}

test_that("merge_gene_te_counts stacks genes and TEs on the shared samples", {
  f <- gene_te_fixture()
  m <- quiet(HORIZON_merge_gene_te_counts(f$g, f$t))

  expect_equal(dim(m), c(4L, 3L))
  expect_equal(rownames(m), c("g1", "g2", "AluY:Alu:SINE", "L1PA2:L1:LINE"))
  expect_equal(m["g2", "s3"], 6L)
  expect_equal(m["L1PA2:L1:LINE", "s1"], 12L)
})

test_that("merge_gene_te_counts keeps only samples present in both and warns about the rest", {
  f <- gene_te_fixture()
  f$t <- f$t[, c("s1", "s2")]
  expect_warning(m <- quiet(HORIZON_merge_gene_te_counts(f$g, f$t)), "s3")

  expect_equal(colnames(m), c("s1", "s2"))
})

test_that("merge_gene_te_counts refuses ambiguous or incomplete input", {
  f <- gene_te_fixture()

  both <- f$t; rownames(both)[1] <- "g1"
  expect_error(HORIZON_merge_gene_te_counts(f$g, both), "collision")

  expect_error(HORIZON_merge_gene_te_counts(unname(f$g), f$t), "rownames")
  expect_error(HORIZON_merge_gene_te_counts(f$g, unname(f$t)), "rownames")

  g_nocols <- f$g; colnames(g_nocols) <- NULL
  expect_error(HORIZON_merge_gene_te_counts(g_nocols, f$t), "colnames")

  other <- f$t; colnames(other) <- c("x1", "x2", "x3")
  expect_error(HORIZON_merge_gene_te_counts(f$g, other), "in common")
})
