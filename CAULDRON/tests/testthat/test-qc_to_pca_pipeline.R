# The core preprocessing chain on a simulated SingleCellExperiment:
#   AEGIS (QC) -> PYRI (normalise, HVGs) -> TALOS (PCA)
# Each step is tested on its own, then the chain is run end to end.
# The fixture (make_test_sce) has 3 cell types with their own marker genes.

quiet <- function(expr) {
  # these functions print progress with cat(); keep test output readable
  utils::capture.output(res <- force(expr))
  res
}

# ---- AEGIS: QC metrics -------------------------------------------------------

test_that("compute_qc_metrics adds the expected per-cell metrics", {
  sce <- quiet(AEGIS_compute_qc_metrics(make_test_sce(), mito_prefix = "mt-",
                                        verbose = FALSE))
  cd <- SummarizedExperiment::colData(sce)

  expect_true(all(c("sum", "detected", "subsets_mt_percent",
                    "subsets_ribo_percent", "log10_genes_per_umi") %in% names(cd)))
  expect_equal(ncol(sce), 120)   # no cells added or removed
  expect_true(all(cd$subsets_mt_percent >= 0 & cd$subsets_mt_percent <= 100))
})

test_that("mitochondrial percentage matches an independent calculation", {
  base <- make_test_sce()
  sce  <- quiet(AEGIS_compute_qc_metrics(base, mito_prefix = "mt-", verbose = FALSE))

  counts  <- SummarizedExperiment::assay(base, "counts")
  is_mito <- grepl("^mt-", rownames(counts))
  expected <- 100 * colSums(counts[is_mito, ]) / colSums(counts)

  expect_equal(unname(sce$subsets_mt_percent), unname(expected))
  expect_equal(unname(sce$sum), unname(colSums(counts)))
})

test_that("compute_qc_metrics warns when the mito prefix matches nothing", {
  expect_warning(
    quiet(AEGIS_compute_qc_metrics(make_test_sce(), mito_prefix = "MT-",
                                   verbose = FALSE)),
    "No mitochondrial genes matched"
  )
})

# ---- AEGIS: filtering --------------------------------------------------------

test_that("filter_cells needs QC metrics first", {
  expect_error(AEGIS_filter_cells(make_test_sce(), remove_doublets = FALSE,
                                  verbose = FALSE),
               "Run AEGIS_compute_qc_metrics")
})

test_that("fixed-mode filtering removes exactly the cells outside the thresholds", {
  sce <- quiet(AEGIS_compute_qc_metrics(make_test_sce(), verbose = FALSE))
  thr <- stats::median(sce$sum)

  out <- quiet(AEGIS_filter_cells(sce, mode = "fixed", min_counts = thr,
                                  remove_doublets = FALSE, verbose = FALSE))

  expect_equal(ncol(out), sum(sce$sum >= thr))
  expect_true(all(out$sum >= thr))
  expect_equal(nrow(out), nrow(sce))   # genes untouched
})

test_that("adaptive filtering removes a planted low-quality cell", {
  base <- make_test_sce()
  # one near-empty cell: library size cut to ~2% of its original
  counts <- SummarizedExperiment::assay(base, "counts")
  counts[, 1] <- rbinom(nrow(counts), counts[, 1], 0.02)
  SummarizedExperiment::assay(base, "counts") <- counts

  sce <- quiet(AEGIS_compute_qc_metrics(base, verbose = FALSE))
  out <- quiet(AEGIS_filter_cells(sce, mode = "adaptive", n_mads = 3,
                                  remove_doublets = FALSE, verbose = FALSE))

  expect_false("cell1" %in% colnames(out))
  expect_gt(ncol(out), 100)   # but the rest of the cells are kept
})

test_that("filter_cells warns and carries on when doublet calls are missing", {
  sce <- quiet(AEGIS_compute_qc_metrics(make_test_sce(), verbose = FALSE))
  expect_warning(
    out <- quiet(AEGIS_filter_cells(sce, mode = "fixed", min_counts = 0,
                                    remove_doublets = TRUE, verbose = FALSE)),
    "scDblFinder.class not found"
  )
  expect_equal(ncol(out), ncol(sce))
})

# ---- PYRI: normalisation and HVGs --------------------------------------------

test_that("normalize adds log-normalised counts without altering the raw counts", {
  base <- make_test_sce()
  sce  <- quiet(PYRI_normalize(base, method = "lognorm", verbose = FALSE))

  expect_true(all(c("counts", "logcounts") %in% SummarizedExperiment::assayNames(sce)))
  expect_identical(SummarizedExperiment::assay(sce, "counts"),
                   SummarizedExperiment::assay(base, "counts"))
  expect_equal(dim(SummarizedExperiment::assay(sce, "logcounts")), dim(base))
  expect_false(anyNA(SummarizedExperiment::assay(sce, "logcounts")))

  # library-size factors are centred on 1 by construction
  expect_equal(mean(scuttle::librarySizeFactors(base)), 1, tolerance = 1e-8)
})

test_that("normalize reports a missing assay", {
  expect_error(PYRI_normalize(make_test_sce(), assay_name = "nope", verbose = FALSE),
               "not found.*Available: counts")
})

test_that("select_hvg stores the selection in rowData and metadata", {
  sce <- quiet(PYRI_normalize(make_test_sce(), verbose = FALSE))
  sce <- quiet(PYRI_select_hvg(sce, n_hvgs = 40, verbose = FALSE))

  expect_equal(sum(SummarizedExperiment::rowData(sce)$is_hvg), 40)
  expect_length(S4Vectors::metadata(sce)$hvg, 40)
  expect_setequal(S4Vectors::metadata(sce)$hvg,
                  rownames(sce)[SummarizedExperiment::rowData(sce)$is_hvg])
})

test_that("select_hvg picks the planted marker genes", {
  sce <- quiet(PYRI_normalize(make_test_sce(), verbose = FALSE))
  sce <- quiet(PYRI_select_hvg(sce, n_hvgs = 40, verbose = FALSE))

  n_marker_hits <- sum(grepl("^marker", S4Vectors::metadata(sce)$hvg))
  # 45 of the 255 genes are markers (18%); a working HVG selection is far above that
  expect_gte(n_marker_hits, 30)
})

test_that("select_hvg validates its inputs and caps n_hvgs at the gene count", {
  sce <- quiet(PYRI_normalize(make_test_sce(), verbose = FALSE))

  expect_error(PYRI_select_hvg(make_test_sce(), verbose = FALSE),
               "Run PYRI_normalize")
  expect_error(PYRI_select_hvg(sce, block_col = "nope", verbose = FALSE),
               "block_col 'nope' not found")

  # n_hvgs is an upper bound. scran::getTopHVGs() only returns genes with
  # positive biological variance, so asking for more than exist returns fewer.
  out <- quiet(PYRI_select_hvg(sce, n_hvgs = 100000, verbose = FALSE))
  n_sel <- length(S4Vectors::metadata(out)$hvg)
  expect_lte(n_sel, nrow(sce))
  expect_gt(n_sel, 0)
  expect_equal(sum(SummarizedExperiment::rowData(out)$is_hvg), n_sel)
  expect_lt(n_sel, nrow(sce))   # some genes have non-positive biological variance

  # the printed summary reports the number actually selected, not the request
  # (the request is capped at the number of genes before it is reported)
  txt <- paste(utils::capture.output(
    invisible(PYRI_select_hvg(sce, n_hvgs = 100000, verbose = TRUE))), collapse = "\n")
  expect_match(txt, sprintf("Selected   : %d / %d genes", n_sel, nrow(sce)),
               fixed = TRUE)
  expect_match(txt, sprintf("requested %d", nrow(sce)), fixed = TRUE)
})

# ---- TALOS: PCA --------------------------------------------------------------

test_that("run_pca stores a cells x PCs embedding", {
  skip_if_not_installed("irlba")
  sce <- quiet(PYRI_normalize(make_test_sce(), verbose = FALSE))
  sce <- quiet(PYRI_select_hvg(sce, n_hvgs = 60, verbose = FALSE))
  sce <- quiet(TALOS_run_pca(sce, n_pcs = 10, verbose = FALSE))

  pca <- SingleCellExperiment::reducedDim(sce, "PCA")
  expect_equal(dim(pca), c(ncol(sce), 10L))
  expect_equal(rownames(pca), colnames(sce))
  expect_false(anyNA(pca))
  pct <- attr(pca, "percentVar")
  expect_length(pct, 10)
  expect_true(all(diff(pct) <= 1e-8))   # variance explained decreases
})

test_that("run_pca warns without HVGs and validates the assay", {
  skip_if_not_installed("irlba")
  sce <- quiet(PYRI_normalize(make_test_sce(), verbose = FALSE))

  expect_warning(quiet(TALOS_run_pca(sce, n_pcs = 5, verbose = FALSE)),
                 "no HVGs found")
  expect_error(TALOS_run_pca(make_test_sce(), verbose = FALSE),
               "Run PYRI_normalize")
})

# ---- End to end --------------------------------------------------------------

test_that("QC -> normalise -> HVG -> PCA runs end to end and separates the cell types", {
  skip_if_not_installed("irlba")
  sce <- make_test_sce()

  sce <- quiet(AEGIS_compute_qc_metrics(sce, verbose = FALSE))
  sce <- quiet(AEGIS_filter_cells(sce, mode = "fixed", min_counts = 1,
                                  remove_doublets = FALSE, verbose = FALSE))
  sce <- quiet(PYRI_normalize(sce, verbose = FALSE))
  sce <- quiet(PYRI_select_hvg(sce, n_hvgs = 60, verbose = FALSE))
  sce <- quiet(TALOS_run_pca(sce, n_pcs = 10, verbose = FALSE))

  expect_true(KHALKOS_validate_sce(
    sce, require_assays = c("counts", "logcounts"),
    require_reduceddims = "PCA", require_metadata = "hvg", verbose = FALSE))

  # The three simulated cell types should separate in the first two PCs
  pcs <- SingleCellExperiment::reducedDim(sce, "PCA")[, 1:2]
  group <- factor(sce$cell_type)
  r2 <- vapply(1:2, function(i) summary(lm(pcs[, i] ~ group))$r.squared, numeric(1))
  expect_gt(max(r2), 0.6)
})
