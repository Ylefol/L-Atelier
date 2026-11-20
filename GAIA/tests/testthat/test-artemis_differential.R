# ARTEMIS differential analysis (DESeq2) on simulated counts with a known
# answer: 10 features up, 10 down, 180 unchanged, 4 vs 4 samples.

# ---- normalisation -----------------------------------------------------------

test_that("normalize_counts returns an artemis_norm object with the expected parts", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))

  expect_s3_class(norm, "artemis_norm")
  expect_equal(dim(norm$norm_counts), dim(d$counts))
  expect_equal(norm$parameters$group_col, "group")
  expect_equal(norm$parameters$size_factors_source, "deseq2_estimated")
  expect_equal(names(norm$size_factors), colnames(d$counts))
  expect_true(all(norm$size_factors > 0))
})

test_that("estimated size factors recover the simulated sequencing depths", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))

  expect_gt(cor(log(norm$size_factors), log(d$depth)), 0.95)
  # normalised counts are the raw counts divided by the size factors
  expect_equal(norm$norm_counts,
               sweep(d$counts, 2, norm$size_factors, "/"),
               ignore_attr = TRUE, tolerance = 1e-8)
})

test_that("normalize_counts validates its inputs", {
  d <- make_test_counts()

  bad <- d$counts; bad[1, 1] <- 1.5
  expect_error(quiet(ARTEMIS_normalize_counts(bad, d$targets, verbose = FALSE)),
               "must be integers")

  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              group_col = "nope", verbose = FALSE)),
               "group_col 'nope' not found")
  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              batch_col = "nope", verbose = FALSE)),
               "batch_col 'nope' not found")

  nm <- d$counts; rownames(nm) <- NULL
  expect_error(quiet(ARTEMIS_normalize_counts(nm, d$targets, verbose = FALSE)),
               "must have rownames")

  t2 <- d$targets; rownames(t2) <- paste0("other", 1:8)
  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, t2, verbose = FALSE)),
               "No matching samples")
})

test_that("samples are matched by name, so column order does not matter", {
  d <- make_test_counts()
  shuffled <- d$counts[, c(5, 1, 8, 2, 3, 7, 4, 6)]
  a <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  b <- quiet(ARTEMIS_normalize_counts(shuffled, d$targets, verbose = FALSE))

  expect_equal(unname(b$size_factors[colnames(d$counts)]),
               unname(a$size_factors), tolerance = 1e-8)
})

test_that("external size factors are used, validated and recorded", {
  d <- make_test_counts()
  sf <- setNames(rep(1, 8), colnames(d$counts))

  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                         size_factors = sf, verbose = FALSE))
  expect_equal(norm$parameters$size_factors_source, "external")
  expect_equal(unname(norm$size_factors), rep(1, 8))
  expect_equal(norm$norm_counts, d$counts, ignore_attr = TRUE)

  # a single number applies to every sample
  one <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                        size_factors = 2, verbose = FALSE))
  expect_equal(unname(one$size_factors), rep(2, 8))

  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              size_factors = rep(1, 8),
                                              verbose = FALSE)),
               "must be named")
  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              size_factors = sf[1:5], verbose = FALSE)),
               "missing values for")
  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              size_factors = -1, verbose = FALSE)),
               "strictly positive")
  expect_error(quiet(ARTEMIS_normalize_counts(d$counts, d$targets,
                                              size_factors = "a", verbose = FALSE)),
               "must be NULL, a named numeric")
})

test_that("a single-level group still normalises (size factors need no design)", {
  d <- make_test_counts()
  d$targets$all <- "everyone"
  expect_no_error(
    norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, group_col = "all",
                                           verbose = FALSE))
  )
  expect_equal(dim(norm$norm_counts), dim(d$counts))
})

# ---- differential testing ----------------------------------------------------

test_that("differential_counts finds the planted up and down features", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  dea  <- quiet(ARTEMIS_differential_counts(norm, reference = "control",
                                            experiment = "treated", verbose = FALSE))
  res <- dea$results

  expect_s3_class(dea, "artemis_dea")
  expect_equal(dea$comparison, "treated vs control")
  expect_named(res, c("feature_id", "baseMean", "log2FoldChange", "lfcSE",
                      "stat", "pvalue", "padj"))
  expect_equal(nrow(res), 200)
  expect_false(is.unsorted(res$padj, na.rm = TRUE))   # most significant first

  sig <- res$feature_id[!is.na(res$padj) & res$padj < 0.05]
  up   <- res[res$feature_id %in% d$up, ]
  down <- res[res$feature_id %in% d$down, ]

  # direction and size of the effect (true log2FC is +2 / -2)
  expect_true(all(up$log2FoldChange > 0))
  expect_true(all(down$log2FoldChange < 0))
  expect_equal(median(up$log2FoldChange), 2, tolerance = 0.4)
  expect_equal(median(down$log2FoldChange), -2, tolerance = 0.4)

  # most planted features are detected; few unchanged features are called
  expect_gte(sum(d$up %in% sig), 8)
  expect_gte(sum(d$down %in% sig), 8)
  expect_lte(sum(d$null %in% sig), 5)
})

test_that("the summary agrees with the results table", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  dea  <- quiet(ARTEMIS_differential_counts(norm, "control", "treated",
                                            alpha = 0.05, verbose = FALSE))
  res <- dea$results
  s   <- dea$summary

  expect_equal(s$n_features, nrow(res))
  expect_equal(s$n_tested, sum(!is.na(res$padj)))
  expect_equal(s$n_significant, sum(res$padj < 0.05, na.rm = TRUE))
  expect_equal(s$n_significant, s$n_sig_up + s$n_sig_down)
})

test_that("swapping reference and experiment flips the sign of the fold changes", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  ab <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))$results
  ba <- quiet(ARTEMIS_differential_counts(norm, "treated", "control", verbose = FALSE))$results

  m <- merge(ab, ba, by = "feature_id", suffixes = c(".ab", ".ba"))
  # equal up to the numerical noise of DESeq2's dispersion fitting
  expect_lt(max(abs(m$log2FoldChange.ab + m$log2FoldChange.ba)), 1e-4)
  expect_lt(cor(m$log2FoldChange.ab, m$log2FoldChange.ba), -0.99999)
})

test_that("a batch term can be included in the design", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, batch_col = "batch",
                                         verbose = FALSE))
  dea <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))

  expect_equal(norm$parameters$batch_col, "batch")
  expect_equal(nrow(dea$results), 200)
  expect_true(all(dea$results$log2FoldChange[dea$results$feature_id %in% d$up] > 0))
})

test_that("differential_counts rejects bad input with helpful messages", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))

  expect_error(ARTEMIS_differential_counts(list(), "control", "treated"),
               "must be an artemis_norm object")
  expect_error(quiet(ARTEMIS_differential_counts(norm, "nope", "treated", verbose = FALSE)),
               "reference 'nope' not found.*control, treated")
  expect_error(quiet(ARTEMIS_differential_counts(norm, "control", "nope", verbose = FALSE)),
               "experiment 'nope' not found")
})

test_that("a design with no residual degrees of freedom is refused up front", {
  d <- make_test_counts()
  # one control vs one treated sample: nothing left to estimate dispersion from
  keep <- c("S1", "S5")
  norm <- quiet(ARTEMIS_normalize_counts(d$counts[, keep], d$targets[keep, ],
                                         verbose = FALSE))
  expect_error(quiet(ARTEMIS_differential_counts(norm, "control", "treated",
                                                 verbose = FALSE)),
               "Not enough samples for this design")
})

# ---- selecting DE genes ------------------------------------------------------

make_de_table <- function(genes, lfc, padj) {
  data.frame(feature_id = genes, log2FoldChange = lfc, padj = padj,
             stringsAsFactors = FALSE)
}

test_that("select_de_genes applies the log2FC and adjusted p-value thresholds", {
  tab <- make_de_table(c("a", "b", "c", "d", "e"),
                       lfc  = c(2, -2, 0.5, 3, 2),
                       padj = c(0.01, 0.04, 0.001, 0.2, NA))
  sel <- ARTEMIS_select_de_genes(tab, l2fc_thresh = 1, p_thresh = 0.05, verbose = FALSE)

  # a, b pass; c fails on effect size, d on p-value, e has no p-value
  expect_setequal(as.character(sel), c("a", "b"))
})

test_that("select_de_genes combines comparisons by union or intersection", {
  t1 <- make_de_table(c("a", "b", "c"), c(2, 2, 2), c(0.01, 0.01, 0.01))
  t2 <- make_de_table(c("b", "c", "d"), c(-2, 2, 2), c(0.01, 0.01, 0.01))

  u <- ARTEMIS_select_de_genes(list(first = t1, second = t2), union = TRUE, verbose = FALSE)
  i <- ARTEMIS_select_de_genes(list(first = t1, second = t2), union = FALSE, verbose = FALSE)

  expect_setequal(as.character(u), c("a", "b", "c", "d"))
  expect_setequal(as.character(i), c("b", "c"))
  # how many comparisons each gene appears in is kept for later use
  gs <- attr(u, "gene_summary")
  expect_equal(gs$n_comparisons[gs$gene == "b"], 2)
  expect_equal(gs$n_comparisons[gs$gene == "a"], 1)
  expect_named(attr(u, "per_comparison"), c("first", "second"))
})

test_that("select_de_genes warns and returns an empty set when nothing passes", {
  # Regression test: this used to fail inside rowSums() when every comparison
  # was read correctly but none had a gene passing the thresholds.
  tab <- make_de_table(c("a", "b"), c(0.1, 0.2), c(0.9, 0.8))
  expect_warning(sel <- ARTEMIS_select_de_genes(tab, verbose = FALSE),
                 "No significant genes")
  expect_length(sel, 0)
  expect_named(attr(sel, "per_comparison"), "comparison_1")

  # same when several comparisons are all empty
  expect_warning(sel2 <- ARTEMIS_select_de_genes(list(x = tab, y = tab), verbose = FALSE),
                 "No significant genes")
  expect_length(sel2, 0)
  expect_named(attr(sel2, "per_comparison"), c("x", "y"))

  # one empty comparison among non-empty ones is not an error and does not warn
  good <- make_de_table(c("a", "b"), c(2, 2), c(0.01, 0.01))
  expect_no_warning(
    sel3 <- ARTEMIS_select_de_genes(list(x = tab, y = good), verbose = FALSE)
  )
  expect_setequal(as.character(sel3), c("a", "b"))
})

test_that("select_de_genes accepts a DESeq2 result object directly", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  dea  <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))
  sel  <- ARTEMIS_select_de_genes(dea, l2fc_thresh = 1, verbose = FALSE)

  expect_gte(length(intersect(sel, c(d$up, d$down))), 16)
  expect_lte(length(intersect(sel, d$null)), 5)
})
