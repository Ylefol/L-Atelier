# ==============================================================================
# PYRI - Normalization & Preprocessing
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The sacred forge fire of Hephaestus. Fire purifies raw ore into refined
# metal — PYRI refines raw count data into analysis-ready form, removing
# technical noise and batch effects.
#
# Core responsibilities:
#   - Normalization (library-size / scran / SCTransform)
#   - Highly variable gene (HVG) selection
#   - Feature scaling
#   - Cell cycle scoring and regression
#   - Batch effect correction (Harmony, fastMNN, etc.)
#
# All functions prefixed: PYRI_
# ==============================================================================


#' Normalise a SingleCellExperiment
#'
#' Computes normalised log-counts and stores them as the \code{"logcounts"}
#' assay.  Two methods are available:
#'
#' \describe{
#'   \item{\code{"lognorm"}}{Library-size normalisation followed by a log1p
#'     transform (\code{scater::logNormCounts}).  Fast and appropriate for
#'     most datasets.  Each cell is scaled to a common size factor derived
#'     from its total count.}
#'   \item{\code{"scran"}}{Pooling-based size factor estimation
#'     (\code{scran::quickCluster} + \code{scran::computeSumFactors}) followed
#'     by \code{scater::logNormCounts}.  More robust for datasets with large
#'     differences in cell type composition or sequencing depth, at the cost
#'     of additional compute time.}
#' }
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"counts"} assay.
#' @param method Character. \code{"lognorm"} (default) or \code{"scran"}.
#' @param assay_name Character. Input count assay. Default \code{"counts"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return The input SCE with a \code{"logcounts"} assay added and, for
#'   \code{method = "scran"}, size factors stored in
#'   \code{sizeFactors(sce)}.
#' @export
PYRI_normalize <- function(sce,
                            method     = c("lognorm", "scran"),
                            assay_name = "counts",
                            verbose    = TRUE) {

  method <- match.arg(method)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(assayNames(sce), collapse = ", "), call. = FALSE)

  if (method == "lognorm") {

    sce <- scater::logNormCounts(sce, assay.type = assay_name)

    if (isTRUE(verbose)) {
      sf <- scater::librarySizeFactors(sce)
      .pyri_print_norm_summary("lognorm", sf)
    }

  } else {

    message("Running scran: quick clustering for size factor estimation...")
    clusters <- scran::quickCluster(sce, assay.type = assay_name)
    sce      <- scran::computeSumFactors(sce, clusters = clusters,
                                          assay.type = assay_name)
    sce      <- scater::logNormCounts(sce, assay.type = assay_name)

    if (isTRUE(verbose)) {
      sf <- BiocGenerics::sizeFactors(sce)
      .pyri_print_norm_summary("scran", sf)
    }
  }

  sce
}


#' Select highly variable genes
#'
#' Models the mean-variance relationship across genes using
#' \code{scran::modelGeneVar} and selects the top \code{n_hvgs} genes by
#' biological variance component.  Results are stored in two places so
#' downstream functions can find them easily:
#'
#' \itemize{
#'   \item \code{rowData(sce)$is_hvg}  — logical vector flagging selected genes
#'   \item \code{metadata(sce)$hvg}    — character vector of selected gene names
#' }
#'
#' For multi-sample datasets it is strongly recommended to set
#' \code{block_col} to a \code{colData} column identifying sample of origin
#' (e.g. \code{"animal_id"}).  This models the mean-variance trend per sample
#' and combines the results, preventing any single sample with an unusual depth
#' from dominating the selection.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"logcounts"} assay
#'   (i.e. \code{\link{PYRI_normalize}} must have been run first).
#' @param n_hvgs Integer. Number of top HVGs to select. Default \code{2000}.
#' @param block_col Character. Column in \code{colData} to use as a blocking
#'   factor (typically sample or animal ID).  \code{NULL} (default) fits one
#'   trend across all cells.
#' @param assay_name Character. Assay to model variance from. Default
#'   \code{"logcounts"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return The input SCE with \code{rowData(sce)$is_hvg} and
#'   \code{metadata(sce)$hvg} populated.
#' @export
PYRI_select_hvg <- function(sce,
                             n_hvgs     = 2000L,
                             block_col  = NULL,
                             assay_name = "logcounts",
                             verbose    = TRUE) {

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  if (!is.null(block_col) && !block_col %in% names(colData(sce)))
    stop("block_col '", block_col, "' not found in colData. ",
         "Available columns: ", paste(names(colData(sce)), collapse = ", "),
         call. = FALSE)

  n_hvgs <- min(as.integer(n_hvgs), nrow(sce))

  block <- if (!is.null(block_col)) colData(sce)[[block_col]] else NULL

  # BPCells IterableMatrix is not compatible with beachmat::rowBlockApply()
  # used internally by scran::modelGeneVar.  Use a BPCells-native path instead.
  if (inherits(assay(sce, assay_name), "IterableMatrix")) {
    gene_var <- .pyri_bpcells_model_gene_var(assay(sce, assay_name), block)
    hvg_list <- rownames(gene_var)[order(gene_var$bio, decreasing = TRUE)][seq_len(n_hvgs)]
  } else {
    gene_var <- scran::modelGeneVar(sce, assay.type = assay_name, block = block)
    hvg_list <- scran::getTopHVGs(gene_var, n = n_hvgs)
  }

  rowData(sce)$is_hvg   <- rownames(sce) %in% hvg_list
  metadata(sce)$hvg     <- hvg_list

  if (isTRUE(verbose)) {
    bio_var <- gene_var$bio
    message(sprintf(
      "\u2500\u2500 PYRI: HVG selection %s\n  Selected   : %s / %s genes\n  Block      : %s\n  Bio var    : %.4f \u2013 %.4f (median %.4f)\n%s",
      strrep("\u2500", 34),
      format(n_hvgs, big.mark = ","),
      format(nrow(sce), big.mark = ","),
      if (!is.null(block_col)) block_col else "none",
      min(bio_var, na.rm = TRUE),
      max(bio_var, na.rm = TRUE),
      median(bio_var, na.rm = TRUE),
      strrep("\u2500", 56)
    ))
  }

  sce
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# BPCells-native mean-variance modelling for HVG selection.
#
# Replaces scran::modelGeneVar when the assay is a BPCells IterableMatrix.
# Uses BPCells::matrix_stats() (streaming, no full matrix load) for per-gene
# mean and variance, then scran::fitTrendVar() to model the technical trend and
# extract the biological component.
#
# For blocked data: fits a trend per block and averages the biological
# component across blocks, matching scran's blocked workflow conceptually.
#
# Returns a data.frame with a $bio column (biological variance per gene),
# row-named by gene, compatible with the verbose block in PYRI_select_hvg.
.pyri_bpcells_model_gene_var <- function(mat, block) {
  gene_names <- rownames(mat)

  .get_bio <- function(m) {
    means <- BPCells::matrix_stats(m, row_stats = "mean")$row_stats["mean", ]
    vars  <- BPCells::matrix_stats(m, row_stats = "variance")$row_stats["variance", ]
    names(means) <- gene_names
    names(vars)  <- gene_names
    trend <- scran::fitTrendVar(means, vars)
    vars - trend$trend(means)
  }

  if (is.null(block)) {
    bio <- .get_bio(mat)
  } else {
    lvls    <- levels(factor(block))
    bio_mat <- matrix(NA_real_, nrow = length(lvls), ncol = length(gene_names))
    for (i in seq_along(lvls))
      bio_mat[i, ] <- .get_bio(mat[, block == lvls[i], drop = FALSE])
    bio <- colMeans(bio_mat, na.rm = TRUE)
    names(bio) <- gene_names
  }

  data.frame(bio = bio, row.names = gene_names)
}


.pyri_print_norm_summary <- function(method, size_factors) {
  bar <- strrep("\u2500", 56)
  message(sprintf(
    "\u2500\u2500 PYRI: Normalisation (%s) %s\n  Size factors : %s \u2013 %s (median %.3f)\n  Output assay : logcounts\n%s",
    method, bar,
    format(round(min(size_factors),  3), nsmall = 3),
    format(round(max(size_factors),  3), nsmall = 3),
    median(size_factors),
    bar
  ))
}
