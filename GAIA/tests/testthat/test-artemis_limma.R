# ARTEMIS limma differential analysis for continuous data (proteomics, arrays,
# log-scale measurements). This is the non-count counterpart of the DESeq2
# analysis; both feed the same downstream steps.
#
# limma is a Suggests dependency, so these tests are skipped where it is absent.

make_limma_data <- function(seed = 5) {
  set.seed(seed)
  G <- 200; n <- 16
  group <- rep(c("ctrl", "trt"), each = 8)
  batch <- rep(c("b1", "b2"), 8)
  M <- matrix(rnorm(G * n, 8, 1), G, n)
  M[1:15,  group == "trt"] <- M[1:15,  group == "trt"] + 2     # 15 up
  M[16:30, group == "trt"] <- M[16:30, group == "trt"] - 2     # 15 down
  M[, batch == "b2"] <- M[, batch == "b2"] + 0.5               # a batch offset
  dimnames(M) <- list(sprintf("p%03d", 1:G), sprintf("s%02d", 1:n))
  list(M = M,
       meta = data.frame(sample = colnames(M), group = group, batch = batch,
                         stringsAsFactors = FALSE),
       up = sprintf("p%03d", 1:15), down = sprintf("p%03d", 16:30),
       null = sprintf("p%03d", 31:G))
}

run_limma <- function(d, ...) {
  quiet_msg(ARTEMIS_limma_de(d$M, d$meta, group_col = "group", reference = "ctrl",
                             experiment = "trt", sample_col = "sample", ...))
}

test_that("limma recovers planted up- and down-regulated features", {
  skip_if_not_installed("limma")
  d <- make_limma_data()
  res <- run_limma(d)

  expect_s3_class(res, "artemis_limma")
  expect_equal(nrow(res$results), 200)
  expect_true(all(c("feature_id", "log2FoldChange", "AveExpr", "t", "pvalue", "padj", "sig")
                  %in% names(res$results)))

  r <- res$results
  sig <- r$feature_id[r$sig]
  expect_gte(sum(d$up %in% sig), 12)
  expect_gte(sum(d$down %in% sig), 12)
  expect_lte(sum(d$null %in% sig), 5)

  # direction and size of the effect (true difference is +2 / -2)
  expect_true(all(r$log2FoldChange[r$feature_id %in% d$up] > 0))
  expect_true(all(r$log2FoldChange[r$feature_id %in% d$down] < 0))
  expect_equal(median(r$log2FoldChange[r$feature_id %in% d$up]), 2, tolerance = 0.3)
})

test_that("the summary agrees with the results table", {
  skip_if_not_installed("limma")
  res <- run_limma(make_limma_data())
  r <- res$results; s <- res$summary

  expect_equal(s$n_features, 200)
  expect_equal(s$n_sig_up, sum(r$sig & r$log2FoldChange > 0))
  expect_equal(s$n_sig_down, sum(r$sig & r$log2FoldChange < 0))
  expect_equal(s$fdr, 0.05)
  expect_match(res$comparison, "trt")
})

test_that("a stricter fdr or a fold-change cutoff reduces the number of significant features", {
  skip_if_not_installed("limma")
  d <- make_limma_data()
  loose  <- run_limma(d)
  strict <- run_limma(d, fdr = 1e-6)
  big    <- run_limma(d, lfc = 3)

  expect_lt(sum(strict$results$sig), sum(loose$results$sig))
  expect_lt(sum(big$results$sig), sum(loose$results$sig))
})

test_that("a covariate can be included in the model", {
  skip_if_not_installed("limma")
  d <- make_limma_data()
  adj <- run_limma(d, covariates = "batch")

  expect_equal(nrow(adj$results), 200)
  expect_gte(sum(d$up %in% adj$results$feature_id[adj$results$sig]), 12)
})

test_that("limma rejects unknown groups and mismatched inputs", {
  skip_if_not_installed("limma")
  d <- make_limma_data()
  expect_error(quiet_msg(ARTEMIS_limma_de(d$M, d$meta, group_col = "group",
                                          reference = "nope", experiment = "trt",
                                          sample_col = "sample")))
  expect_error(quiet_msg(ARTEMIS_limma_de(d$M, d$meta, group_col = "zzz",
                                          reference = "ctrl", experiment = "trt",
                                          sample_col = "sample")))
})

test_that("limma results go through the same ranking step as DESeq2 results", {
  skip_if_not_installed("limma")
  d <- make_limma_data()
  res <- run_limma(d)
  ranked <- APOLLO_rank_from_de(res, rank_by = "t", verbose = FALSE)   # artemis_limma input

  expect_named(ranked)
  expect_false(is.unsorted(rev(ranked)))
  expect_gte(sum(names(ranked)[1:20] %in% d$up), 12)
  expect_gte(sum(tail(names(ranked), 20) %in% d$down), 12)
})
