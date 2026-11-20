# POSEIDON_match_samples() / filter_by_metadata() operate on 'quant_result'
# lists (counts = features x samples, targets = one row per sample). The same
# structure is used for every assay, which is what lets datasets be aligned.

make_quant <- function(sample_ids, bio_reps, n_features = 6) {
  counts <- matrix(seq_len(n_features * length(sample_ids)),
                   nrow = n_features,
                   dimnames = list(paste0("f", seq_len(n_features)), sample_ids))
  targets <- data.frame(sample_id = sample_ids, bio_rep = bio_reps,
                        stringsAsFactors = FALSE)
  list(counts = counts, targets = targets)
}

test_that("filter_by_metadata subsets counts and targets together", {
  q <- make_quant(c("a", "b", "c", "d"), c(1, 1, 2, 3))
  out <- POSEIDON_filter_by_metadata(q, column = "bio_rep", values = c(1, 3),
                                     verbose = FALSE)

  expect_equal(colnames(out$counts), c("a", "b", "d"))
  expect_equal(out$targets$sample_id, c("a", "b", "d"))
  expect_equal(ncol(out$counts), nrow(out$targets))
})

test_that("filter_by_metadata errors clearly on bad input", {
  q <- make_quant(c("a", "b"), c(1, 2))
  expect_error(POSEIDON_filter_by_metadata(q, "nope", 1, verbose = FALSE),
               "not found in targets")
  expect_error(POSEIDON_filter_by_metadata(q, "bio_rep", 99, verbose = FALSE),
               "No samples found")
})

test_that("match_samples keeps only the biological replicates present in every dataset", {
  atac <- make_quant(c("a1", "a2", "a3", "a4"), c(1, 2, 3, 4))
  rna  <- make_quant(c("r1", "r2", "r3"),       c(2, 3, 5))
  chip <- make_quant(c("c1", "c2", "c3", "c4"), c(1, 2, 2, 3))   # replicate 2 twice

  res <- POSEIDON_match_samples(list(ATAC = atac, RNA = rna, ChIP = chip),
                                match_col = "bio_rep", verbose = FALSE)

  expect_setequal(res$matched_values, c(2, 3))
  expect_named(res$datasets, c("ATAC", "RNA", "ChIP"))
  expect_equal(res$datasets$ATAC$targets$bio_rep, c(2, 3))
  expect_equal(res$datasets$RNA$targets$bio_rep, c(2, 3))
  expect_equal(sort(res$datasets$ChIP$targets$bio_rep), c(2, 2, 3))

  # the sample map reports how many samples each dataset contributes per value
  expect_equal(res$sample_map$ChIP_n[res$sample_map$bio_rep == 2], 2)
  expect_equal(res$sample_map$ATAC_n[res$sample_map$bio_rep == 2], 1)
})

test_that("match_samples validates its input", {
  q <- make_quant(c("a", "b"), c(1, 2))
  expect_error(POSEIDON_match_samples(list(q), verbose = FALSE),
               "at least 2 datasets")
  expect_error(POSEIDON_match_samples(list(q, q), match_col = "nope",
                                      verbose = FALSE),
               "not found in targets")

  other <- make_quant(c("x", "y"), c(8, 9))
  expect_error(POSEIDON_match_samples(list(q, other), verbose = FALSE),
               "No matching values")
})
