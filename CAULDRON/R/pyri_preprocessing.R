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

    cat("Running scran: quick clustering for size factor estimation...\n")
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
    cat(sprintf(
      "\u2500\u2500 PYRI: HVG selection %s\n  Selected   : %s / %s genes\n  Block      : %s\n  Bio var    : %.4f \u2013 %.4f (median %.4f)\n%s\n",
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


#' Sweep HVG count to find a data-driven optimum
#'
#' Tests a range of \code{n_hvgs} values and, for each, scores the resulting
#' gene set on two complementary metrics:
#'
#' \enumerate{
#'   \item \strong{Cumulative variance explained} — the percentage of total
#'     gene-expression variance captured by the top \code{n_pcs} principal
#'     components computed from the selected HVGs.  Increasing \code{n_hvgs}
#'     initially raises this sharply, then flattens as noise genes dilute the
#'     signal.  The elbow of this curve (detected via the kneedle method) is
#'     used as the suggested optimum.
#'   \item \strong{Mean biological variance} — the mean of the scran biological
#'     variance component across selected genes.  This falls monotonically as
#'     lower-variance genes enter the set and serves as a sanity check: the
#'     suggested \code{n_hvgs} should lie in the region where mean bio variance
#'     is still meaningfully above zero.
#' }
#'
#' The mean-variance model is fitted \strong{once} across all genes; the sweep
#' only changes the cutoff on the pre-ranked list, so runtime is dominated by
#' the repeated lightweight PCA steps, not by modelling.  PCA is performed with
#' \code{irlba::irlba} directly (BPCells-aware — no full matrix materialisation
#' required).
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"logcounts"} assay
#'   (\code{\link{PYRI_normalize}} must have been run first).
#' @param n_hvgs_range Integer vector of candidate HVG counts to test.
#'   Values exceeding \code{nrow(sce)} are silently dropped.
#'   Default \code{c(500, 1000, 2000, 3000, 5000, 7500)}.
#' @param n_pcs Integer.  Number of PCA components used to measure variance
#'   explained at each step.  Smaller values are faster; \code{20} is usually
#'   sufficient for a sweep.  Default \code{20}.
#' @param scale Logical.  Whether to scale genes to unit variance before PCA,
#'   matching \code{\link{TALOS_run_pca}}'s \code{scale} argument.
#'   Default \code{FALSE}.
#' @param block_col Character.  Column in \code{colData} for per-sample
#'   blocking in the mean-variance model (matches \code{\link{PYRI_select_hvg}}).
#'   \code{NULL} (default) fits one trend across all cells.
#' @param assay_name Character.  Assay to model.  Default \code{"logcounts"}.
#' @param seed Integer.  Random seed for irlba reproducibility.  Default
#'   \code{42L}.
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return A \code{pyri_hvg_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame: \code{n_hvgs}, \code{cum_var},
#'       \code{mean_bio_var}.}
#'     \item{\code{plot}}{Two-panel ggplot (cum var + mean bio var vs n_hvgs).}
#'     \item{\code{best_n_hvgs}}{Suggested HVG count (elbow in cum_var).}
#'     \item{\code{params}}{Sweep parameters.}
#'   }
#' @export
PYRI_tune_hvg <- function(sce,
                           n_hvgs_range = c(500L, 1000L, 2000L, 3000L, 5000L, 7500L),
                           n_pcs        = 20L,
                           scale        = FALSE,
                           block_col    = NULL,
                           assay_name   = "logcounts",
                           seed         = 42L,
                           verbose      = TRUE) {

  if (!requireNamespace("irlba", quietly = TRUE))
    stop("Package 'irlba' is required. Install via: install.packages(\"irlba\")",
         call. = FALSE)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  if (!is.null(block_col) && !block_col %in% names(colData(sce)))
    stop("block_col '", block_col, "' not found in colData. ",
         "Available columns: ", paste(names(colData(sce)), collapse = ", "),
         call. = FALSE)

  n_hvgs_range <- sort(unique(as.integer(n_hvgs_range)))
  n_hvgs_range <- n_hvgs_range[n_hvgs_range >= 2L & n_hvgs_range <= nrow(sce)]
  if (length(n_hvgs_range) == 0L)
    stop("No valid n_hvgs values after clipping to [2, nrow(sce)].", call. = FALSE)

  n_pcs   <- as.integer(n_pcs)
  n_cells <- ncol(sce)
  n_genes <- nrow(sce)
  mat     <- assay(sce, assay_name)
  block   <- if (!is.null(block_col)) colData(sce)[[block_col]] else NULL

  # ── Step 1: Fit mean-variance model ONCE over all genes ──────────────────────
  if (verbose) cat("  Modelling mean-variance relationship (all genes) ...\n")

  if (inherits(mat, "IterableMatrix")) {
    gene_var_df <- .pyri_bpcells_model_gene_var(mat, block)
    gene_bio    <- setNames(gene_var_df$bio, rownames(gene_var_df))
  } else {
    gene_var <- scran::modelGeneVar(sce, assay.type = assay_name, block = block)
    gene_bio <- setNames(as.numeric(gene_var$bio), rownames(gene_var))
  }

  # Rank genes by biological variance (descending) — done once
  ranked_genes <- names(sort(gene_bio, decreasing = TRUE))

  if (verbose)
    cat(sprintf(
      "\u2500\u2500 PYRI: HVG sweep %s\n  n_hvgs range : %s\n  n_pcs        : %d  |  scale: %s\n  block_col    : %s\n  Total genes  : %s\n%s\n",
      strrep("\u2500", 37),
      paste(format(n_hvgs_range, big.mark = ","), collapse = ", "),
      n_pcs, scale,
      if (!is.null(block_col)) block_col else "none",
      format(n_genes, big.mark = ","),
      strrep("\u2500", 56)))

  # ── Step 2: Sweep ────────────────────────────────────────────────────────────
  results <- vector("list", length(n_hvgs_range))

  for (i in seq_along(n_hvgs_range)) {
    n <- n_hvgs_range[i]
    if (verbose)
      cat(sprintf("  (%d/%d) n_hvgs = %s ...\n",
                      i, length(n_hvgs_range), format(n, big.mark = ",")))

    genes_i  <- ranked_genes[seq_len(n)]
    mat_i    <- mat[genes_i, , drop = FALSE]
    cum_var  <- .pyri_quick_pca_var(mat_i, n_pcs, scale, seed)
    mean_bio <- mean(gene_bio[genes_i], na.rm = TRUE)

    results[[i]] <- data.frame(
      n_hvgs       = n,
      cum_var      = round(cum_var,  2),
      mean_bio_var = round(mean_bio, 4)
    )
  }

  res_df <- do.call(rbind, results)

  # ── Step 3: Elbow detection on cum_var ───────────────────────────────────────
  best_n_hvgs <- .pyri_elbow(res_df$n_hvgs, res_df$cum_var)

  p <- .pyri_tune_hvg_plot(res_df, best_n_hvgs)

  if (verbose)
    cat(sprintf(
      "%s\n  Suggested n_hvgs: %s (elbow in variance explained)\n%s\n",
      strrep("\u2500", 56),
      format(best_n_hvgs, big.mark = ","),
      strrep("\u2500", 56)))

  structure(
    list(results     = res_df,
         plot        = p,
         best_n_hvgs = best_n_hvgs,
         params      = list(n_hvgs_range = n_hvgs_range,
                            n_pcs        = n_pcs,
                            scale        = scale,
                            block_col    = block_col,
                            assay_name   = assay_name,
                            seed         = seed)),
    class = "pyri_hvg_sweep"
  )
}


#' @export
print.pyri_hvg_sweep <- function(x, ...) {
  cat("\u2500\u2500 PYRI HVG sweep \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  n_hvgs tested   : %s\n",
              paste(format(x$params$n_hvgs_range, big.mark = ","), collapse = ", ")))
  cat(sprintf("  Suggested n_hvgs: %s  (elbow in variance explained)\n\n",
              format(x$best_n_hvgs, big.mark = ",")))
  print(x$results, row.names = FALSE)
  invisible(x)
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


# Run a lightweight irlba PCA on a gene-subset matrix and return the cumulative
# percentage of total variance explained by the top n_pcs components.
#
# mat   : genes × cells (dense, sparse Matrix, or BPCells IterableMatrix)
# n_pcs : number of PCs to compute
# scale : logical — scale genes to unit variance before PCA
# seed  : integer random seed for irlba reproducibility
#
# Returns a single numeric: sum(d^2 / (n_cells-1)) / total_var * 100
.pyri_quick_pca_var <- function(mat, n_pcs, scale, seed) {
  n_genes <- nrow(mat)
  n_cells <- ncol(mat)
  n_pcs   <- min(n_pcs, n_genes - 1L, n_cells - 1L)

  is_bpcells <- inherits(mat, "IterableMatrix")

  # Per-gene means and variances (streaming for BPCells, vectorised otherwise)
  if (is_bpcells) {
    gene_means <- BPCells::matrix_stats(mat, row_stats = "mean")$row_stats["mean", ]
    gene_vars  <- BPCells::matrix_stats(mat, row_stats = "variance")$row_stats["variance", ]
  } else {
    gene_means <- rowMeans(mat)
    # E[X^2] - E[X]^2, then Bessel correction — works for dense and sparse Matrix
    gene_vars  <- (rowMeans(mat * mat) - gene_means^2) * n_cells / (n_cells - 1L)
  }

  scale_vec <- if (isTRUE(scale)) {
    pmax(sqrt(gene_vars), .Machine$double.eps)
  } else {
    FALSE
  }

  # Transpose for cells × genes layout required by irlba.
  # BiocGenerics::t() dispatches S4 methods — handles dgCMatrix and BPCells
  # IterableMatrix alike; base::t.default() fails on both.
  tmat <- BiocGenerics::t(mat)

  set.seed(seed)
  result <- irlba::irlba(tmat, nv = n_pcs, center = gene_means, scale = scale_vec)

  total_var <- if (isTRUE(scale)) n_genes else sum(gene_vars)
  sum(result$d^2 / (n_cells - 1L)) / total_var * 100
}


# Kneedle elbow detection on a monotone curve.
# Finds the point of maximum perpendicular distance from the line connecting
# the first and last (x, y) point — the "elbow" where returns diminish.
# x : numeric vector (n_hvgs values, sorted ascending)
# y : numeric vector (cum_var values, same length)
# Returns the x value at the detected elbow.
.pyri_elbow <- function(x, y) {
  n <- length(x)
  if (n <= 2L) return(x[which.max(y)])

  # Normalise both axes to [0,1]
  x_n <- (x - min(x)) / (max(x) - min(x))
  y_n <- (y - min(y)) / (max(y) - min(y) + .Machine$double.eps)

  # Direction vector of the line from first to last point
  dx <- x_n[n] - x_n[1L]
  dy <- y_n[n] - y_n[1L]

  # Perpendicular distance of each point from that line
  dists <- abs(dy * x_n - dx * y_n + x_n[n] * y_n[1L] - y_n[n] * x_n[1L]) /
           sqrt(dx^2 + dy^2)

  x[which.max(dists)]
}


# Two-panel line plot for the HVG sweep.
# Panel 1: cumulative variance explained (%) vs n_hvgs
# Panel 2: mean biological variance vs n_hvgs
# Dashed red vertical line at best_n_hvgs.
.pyri_tune_hvg_plot <- function(res_df, best_n_hvgs) {

  metrics <- c("cum_var", "mean_bio_var")
  lab_map <- c(
    cum_var      = "Cumulative variance\nexplained (%) in top PCs",
    mean_bio_var = "Mean biological\nvariance of selected HVGs"
  )

  df_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(n_hvgs = res_df$n_hvgs,
               metric = lab_map[m],
               value  = res_df[[m]])
  }))
  df_long$metric <- factor(df_long$metric, levels = lab_map[metrics])

  ggplot(df_long, aes(x = n_hvgs, y = value)) +
    geom_line(colour = "#4E79A7", linewidth = 0.8) +
    geom_point(colour = "#4E79A7", size = 2.5) +
    geom_vline(xintercept = best_n_hvgs,
               linetype = "dashed", colour = "#E15759", linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x       = "Number of HVGs",
         y       = NULL,
         title   = "HVG count sweep",
         caption = sprintf("Red dashed: suggested n_hvgs = %s",
                           format(best_n_hvgs, big.mark = ","))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


.pyri_print_norm_summary <- function(method, size_factors) {
  bar <- strrep("\u2500", 56)
  cat(sprintf(
    "\u2500\u2500 PYRI: Normalisation (%s) %s\n  Size factors : %s \u2013 %s (median %.3f)\n  Output assay : logcounts\n%s\n",
    method, bar,
    format(round(min(size_factors),  3), nsmall = 3),
    format(round(max(size_factors),  3), nsmall = 3),
    median(size_factors),
    bar
  ))
}
