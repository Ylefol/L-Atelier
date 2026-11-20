# APOLLO: preparing ranked lists for enrichment and filtering enrichment tables.
# The point of these tests is the shared interface: DE results from different
# assays (DESeq2 counts, limma-style proteomics) go through the same code and
# come out in the same form, ready for the same enrichment functions.

make_de_df <- function() {
  data.frame(
    feature_id     = c("g1", "g2", "g3", "g4", "g5"),
    log2FoldChange = c(2.0, -1.5, 0.2, -3.0, NA),
    stat           = c(5.0, -3.0, 0.4, -6.0, NA),
    pvalue         = c(1e-6, 1e-3, 0.5, 1e-9, NA),
    stringsAsFactors = FALSE
  )
}

test_that("rank_from_de returns a named, descending numeric vector without NAs", {
  r <- APOLLO_rank_from_de(make_de_df(), rank_by = "t", verbose = FALSE)

  expect_type(r, "double")
  expect_named(r)
  expect_false(anyNA(r))
  expect_false(is.unsorted(rev(r)))                 # descending
  expect_equal(names(r), c("g1", "g3", "g2", "g4"))  # g5 (NA) dropped
  expect_equal(unname(r["g4"]), -6)
})

test_that("each ranking statistic uses the right column", {
  df <- make_de_df()

  lfc <- APOLLO_rank_from_de(df, rank_by = "log2FoldChange", verbose = FALSE)
  expect_equal(unname(lfc["g1"]), 2)

  sl <- APOLLO_rank_from_de(df, rank_by = "signed_neg_log10p", verbose = FALSE)
  expect_equal(unname(sl["g1"]), 6)                  # +sign, -log10(1e-6)
  expect_equal(unname(sl["g4"]), -9)                 # -sign, -log10(1e-9)
})

test_that("limma-style results (column 't') rank the same way as DESeq2 ('stat')", {
  deseq <- make_de_df()
  limma <- data.frame(feature_id = deseq$feature_id, logFC = deseq$log2FoldChange,
                      t = deseq$stat, stringsAsFactors = FALSE)
  limma_obj <- structure(list(results = limma), class = c("artemis_limma", "list"))

  r_deseq <- APOLLO_rank_from_de(deseq, rank_by = "t", verbose = FALSE)
  r_limma <- APOLLO_rank_from_de(limma_obj, rank_by = "t", verbose = FALSE)

  # same interface in, same form out, whichever assay produced the statistics
  expect_identical(r_deseq, r_limma)
})

test_that("rank_from_de removes duplicated gene names, keeping the first", {
  df <- rbind(make_de_df(),
              data.frame(feature_id = "g1", log2FoldChange = -9, stat = -99,
                         pvalue = 0.1, stringsAsFactors = FALSE))
  r <- APOLLO_rank_from_de(df, rank_by = "t", verbose = FALSE)

  expect_equal(sum(names(r) == "g1"), 1)
  expect_equal(unname(r["g1"]), 5)   # the first occurrence, not -99
})

test_that("rank_from_de handles time-series results by comparison name or position", {
  ts <- structure(
    list(results = list(
      early = list(results = make_de_df()),
      late  = list(results = transform(make_de_df(), stat = -stat))
    )),
    class = c("artemis_ts_de", "list")
  )

  first <- APOLLO_rank_from_de(ts, rank_by = "t", verbose = FALSE)
  late  <- APOLLO_rank_from_de(ts, rank_by = "t", comparison = "late", verbose = FALSE)
  byidx <- APOLLO_rank_from_de(ts, rank_by = "t", comparison = 2, verbose = FALSE)

  expect_equal(unname(first["g1"]), 5)      # default: first comparison
  expect_equal(unname(late["g1"]), -5)
  expect_identical(late, byidx)
  expect_error(APOLLO_rank_from_de(ts, comparison = "nope", verbose = FALSE),
               "comparison 'nope' not found.*early, late")
})

test_that("rank_from_de explains missing columns and unsupported input", {
  df <- make_de_df()

  expect_error(APOLLO_rank_from_de(df, gene_col = "symbol", verbose = FALSE),
               "gene_col 'symbol' not found")
  expect_error(APOLLO_rank_from_de(df[, c("feature_id", "log2FoldChange")],
                                   rank_by = "t", verbose = FALSE),
               "requires a 't'")
  expect_error(APOLLO_rank_from_de(df[, c("feature_id", "pvalue")],
                                   rank_by = "log2FoldChange", verbose = FALSE),
               "requires a 'log2FoldChange'")
  expect_error(APOLLO_rank_from_de(df[, c("feature_id", "log2FoldChange")],
                                   rank_by = "signed_neg_log10p", verbose = FALSE),
               "requires 'log2FoldChange' and 'pvalue'")
  expect_error(APOLLO_rank_from_de(1:3, verbose = FALSE),
               "must be an artemis_limma")
})

# ---- filtering enrichment tables ---------------------------------------------

make_enrichment <- function() {
  combined <- data.frame(
    module    = rep(c("m1", "m2"), each = 6),
    source    = rep(rep(c("GO:BP", "KEGG"), each = 3), times = 2),
    term_name = paste0("term", 1:12),
    p_value   = c(0.001, 0.01, 0.02,  0.03, 0.04, 0.05,
                  0.002, 0.003, 0.2,  0.5, 0.6, 0.7),
    stringsAsFactors = FALSE
  )
  structure(list(combined = combined), class = c("gost_enrichment", "list"))
}

test_that("filter_enrichment returns everything when no filter is given", {
  enr <- make_enrichment()
  expect_equal(nrow(APOLLO_filter_enrichment(enr)), 12)
})

test_that("filter_enrichment filters by source and by module", {
  enr <- make_enrichment()

  go <- APOLLO_filter_enrichment(enr, sources = "GO:BP")
  expect_equal(nrow(go), 6)
  expect_true(all(go$source == "GO:BP"))

  m2 <- APOLLO_filter_enrichment(enr, modules = "m2")
  expect_true(all(m2$module == "m2"))

  both <- APOLLO_filter_enrichment(enr, sources = "KEGG", modules = "m1")
  expect_equal(both$term_name, c("term4", "term5", "term6"))
})

test_that("top_n keeps the most significant terms per module and source", {
  enr <- make_enrichment()
  top <- APOLLO_filter_enrichment(enr, top_n = 2)

  # 2 modules x 2 sources x top 2 = 8 rows
  expect_equal(nrow(top), 8)
  counts <- table(top$module, top$source)
  expect_true(all(counts == 2))
  # within a group the smallest p-values were kept
  m1_go <- top[top$module == "m1" & top$source == "GO:BP", ]
  expect_equal(sort(m1_go$p_value), c(0.001, 0.01))
  expect_equal(rownames(top), as.character(seq_len(nrow(top))))
})

test_that("filter_enrichment copes with empty results and wrong input", {
  enr <- make_enrichment()
  expect_equal(nrow(APOLLO_filter_enrichment(enr, sources = "REAC")), 0)

  empty <- structure(list(combined = enr$combined[0, ]),
                     class = c("gost_enrichment", "list"))
  expect_equal(nrow(APOLLO_filter_enrichment(empty, top_n = 3)), 0)

  expect_error(APOLLO_filter_enrichment(list(combined = data.frame())),
               "must be a gost_enrichment object")
})

# ---- the chain: counts -> DE -> ranked list ----------------------------------

test_that("counts -> normalise -> DE -> ranked list runs end to end", {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  dea  <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))

  ranked <- quiet(APOLLO_rank_from_de(dea$results, rank_by = "t", verbose = FALSE))

  expect_length(ranked, sum(!is.na(dea$results$stat)))
  # the planted up-regulated features sit at the top, down-regulated at the bottom
  expect_gte(sum(names(ranked)[1:15] %in% d$up), 8)
  expect_gte(sum(tail(names(ranked), 15) %in% d$down), 8)

  # the same genes selected as DE are the ones at the extremes of the ranking
  sel <- ARTEMIS_select_de_genes(dea, l2fc_thresh = 1, verbose = FALSE)
  extremes <- c(head(names(ranked), 15), tail(names(ranked), 15))
  expect_gte(length(intersect(sel, extremes)), 12)
})
