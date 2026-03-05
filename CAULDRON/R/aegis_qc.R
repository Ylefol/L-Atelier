# ==============================================================================
# AEGIS - Quality Control & Filtering
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The divine shield forged by Hephaestus. AEGIS protects the dataset from
# bad data — filtering low-quality cells, doublets, and ambient RNA
# contamination before any analysis proceeds.
#
# Core responsibilities:
#   - Per-cell QC metric computation (nUMI, nGenes, MT%, ribo%)
#   - Adaptive threshold filtering (MAD-based)
#   - Doublet detection (scDblFinder)
#   - Ambient RNA removal (SoupX / DecontX)
#   - QC visualisation (violin plots, scatter plots, knee plots)
#   - Sample-level QC summaries
#
# All functions prefixed: AEGIS_
# ==============================================================================


#' Compute per-cell QC metrics
#'
#' Adds standard per-cell quality metrics to the \code{colData} of a
#' \code{SingleCellExperiment}:
#'
#' \itemize{
#'   \item \code{sum}          — total UMI / read counts per cell
#'   \item \code{detected}     — number of genes with non-zero counts
#'   \item \code{subsets_mt_percent}   — \% counts from mitochondrial genes
#'   \item \code{subsets_ribo_percent} — \% counts from ribosomal genes
#'     (if ribosomal genes are found)
#'   \item \code{log10_genes_per_umi}  — complexity score:
#'     \code{log10(detected) / log10(sum)}; values near 0.8 indicate
#'     high-complexity libraries
#' }
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_name Character. Assay to compute metrics from.  If \code{NULL}
#'   (default), tries \code{"counts"} then \code{"X"} then the first assay.
#' @param mito_prefix Character. Prefix identifying mitochondrial genes.
#'   Default \code{"mt-"} (mouse).  Use \code{"MT-"} for human.
#' @param ribo_prefixes Character vector. Prefixes identifying ribosomal genes.
#'   Default \code{c("Rps", "Rpl")} (mouse).  Use \code{c("RPS", "RPL")} for
#'   human.  Pass \code{NULL} to skip ribosomal metrics.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return The input \code{SingleCellExperiment} with QC columns added to
#'   \code{colData}.
#' @export
AEGIS_compute_qc_metrics <- function(sce,
                                      assay_name    = NULL,
                                      mito_prefix   = "mt-",
                                      ribo_prefixes = c("Rps", "Rpl"),
                                      verbose       = TRUE) {

  assay_name <- .aegis_resolve_assay(sce, assay_name)
  genes      <- rownames(sce)

  # ── Identify mitochondrial genes ────────────────────────────────────────────
  is_mito <- logical(nrow(sce))
  if (!is.null(mito_prefix) && nzchar(mito_prefix))
    is_mito <- grepl(paste0("^", mito_prefix), genes)

  if (!any(is_mito) && !is.null(mito_prefix))
    warning("No mitochondrial genes matched prefix '", mito_prefix,
            "'. Check rownames(sce) and adjust mito_prefix.", call. = FALSE)

  # ── Identify ribosomal genes ─────────────────────────────────────────────────
  is_ribo <- logical(nrow(sce))
  if (!is.null(ribo_prefixes) && length(ribo_prefixes) > 0) {
    pat     <- paste0("^(", paste(ribo_prefixes, collapse = "|"), ")")
    is_ribo <- grepl(pat, genes)
  }

  # ── Build subset list and compute metrics ────────────────────────────────────
  subsets <- list()
  if (any(is_mito)) subsets$mt   <- is_mito
  if (any(is_ribo)) subsets$ribo <- is_ribo

  # BPCells IterableMatrix does not implement S4Arrays::extract_array(), so
  # beachmat::colBlockApply() (used internally by scuttle) fails.  Use
  # BPCells-native streaming operations instead.
  if (inherits(assay(sce, assay_name), "IterableMatrix")) {
    sce <- .aegis_bpcells_qc_metrics(sce, assay_name, is_mito, is_ribo)
  } else {
    sce <- scuttle::addPerCellQCMetrics(sce, subsets = subsets,
                                         assay.type = assay_name)
  }

  # ── Complexity score ─────────────────────────────────────────────────────────
  s <- colData(sce)$sum
  d <- colData(sce)$detected
  colData(sce)$log10_genes_per_umi <- ifelse(
    s > 1 & d > 0, log10(d) / log10(s), NA_real_
  )

  if (isTRUE(verbose))
    .aegis_print_metrics_summary(sce, is_mito, is_ribo, assay_name)

  sce
}


#' Detect doublets with scDblFinder
#'
#' Identifies likely doublets (droplets containing two cells) using
#' \code{\link[scDblFinder]{scDblFinder}}.  Doublet scores and classifications
#' are appended to \code{colData} as \code{scDblFinder.score} and
#' \code{scDblFinder.class}.
#'
#' When the dataset contains multiple samples, pass the \code{sample_col}
#' argument so that doublet simulation is performed per sample, which
#' substantially improves accuracy.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param sample_col Character. Column in \code{colData} identifying sample of
#'   origin (e.g. \code{"mouse_id"}).  Doublet detection is then run per
#'   sample.  \code{NULL} (default) treats the entire dataset as one sample.
#' @param clusters_col Character. Column in \code{colData} with cluster labels.
#'   Providing pre-computed clusters improves sensitivity.  \code{NULL}
#'   (default) lets scDblFinder cluster internally.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallelisation.
#'   Default \code{SerialParam()} (single-threaded).
#' @param seed Integer. Random seed for reproducibility. Default \code{42L}.
#' @param force_bpcells Logical. When the counts assay is a
#'   \pkg{BPCells}-backed \code{IterableMatrix}, \code{scDblFinder} will
#'   internally coerce the full matrix into RAM, which may exhaust memory on
#'   large datasets.  By default (\code{FALSE}), doublet detection is skipped
#'   and a message is printed.  Set \code{TRUE} to run regardless — only do
#'   this if you are confident the dataset fits in available RAM.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return The input \code{SingleCellExperiment} with doublet annotations added
#'   to \code{colData}, or the unchanged SCE if the step was skipped due to a
#'   BPCells backend.
#' @export
AEGIS_detect_doublets <- function(sce,
                                   sample_col    = NULL,
                                   clusters_col  = NULL,
                                   BPPARAM       = BiocParallel::SerialParam(),
                                   seed          = 42L,
                                   force_bpcells = FALSE,
                                   verbose       = TRUE) {

  # ── BPCells safeguard ────────────────────────────────────────────────────────
  if (inherits(counts(sce), "IterableMatrix") && !isTRUE(force_bpcells)) {
    message(
      "\u2500\u2500 AEGIS: Doublet detection skipped \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n",
      "  The counts assay is BPCells-backed (on-disk IterableMatrix).\n",
      "  scDblFinder would load the full count matrix into RAM to run, which\n",
      "  defeats the purpose of using BPCells and risks an out-of-memory crash.\n",
      "\n",
      "  To run doublet detection anyway (if RAM permits):\n",
      "    AEGIS_detect_doublets(sce, force_bpcells = TRUE)\n",
      strrep("\u2500", 56)
    )
    return(invisible(sce))
  }

  if (!is.null(sample_col) && !sample_col %in% names(colData(sce)))
    stop("sample_col '", sample_col, "' not found in colData. ",
         "Available columns: ", paste(names(colData(sce)), collapse = ", "),
         call. = FALSE)

  if (!is.null(clusters_col) && !clusters_col %in% names(colData(sce)))
    stop("clusters_col '", clusters_col, "' not found in colData.",
         call. = FALSE)

  samples  <- if (!is.null(sample_col))   colData(sce)[[sample_col]]   else NULL
  clusters <- if (!is.null(clusters_col)) colData(sce)[[clusters_col]] else NULL

  # scDblFinder uses counts(sce) internally; TALARIA_load_h5ad guarantees the
  # assay is named "counts" so no aliasing is needed here.
  set.seed(seed)
  sce <- scDblFinder::scDblFinder(sce,
                                   samples  = samples,
                                   clusters = clusters,
                                   BPPARAM  = BPPARAM)

  if (isTRUE(verbose)) {
    cls   <- colData(sce)$scDblFinder.class
    n_dbl <- sum(cls == "doublet", na.rm = TRUE)
    n_tot <- length(cls)
    message(sprintf(
      "\u2500\u2500 AEGIS: Doublet detection \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n  Singlets : %s (%s%%)\n  Doublets : %s (%s%%)\n%s",
      format(n_tot - n_dbl, big.mark = ","),
      round(100 * (n_tot - n_dbl) / n_tot, 1),
      format(n_dbl, big.mark = ","),
      round(100 * n_dbl / n_tot, 1),
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' Filter low-quality cells
#'
#' Removes cells that fail QC thresholds, returning a filtered
#' \code{SingleCellExperiment}.
#'
#' Two modes are supported:
#' \describe{
#'   \item{\code{"adaptive"}}{Outliers are identified per metric using
#'     median-absolute-deviation (MAD) via
#'     \code{\link[scuttle]{perCellQCFilters}}.  Cells more than
#'     \code{n_mads} MADs below the median for \code{sum} or \code{detected},
#'     or above the median for MT\%, are flagged.  This approach automatically
#'     adapts to dataset depth and is the recommended default.}
#'   \item{\code{"fixed"}}{Hard thresholds supplied by the user.  At least one
#'     of \code{min_counts}, \code{max_counts}, \code{min_features},
#'     \code{max_features}, or \code{max_pct_mt} must be provided.}
#' }
#'
#' In either mode, \code{max_pct_mt} acts as an additional hard ceiling if
#' supplied alongside \code{mode = "adaptive"}.
#'
#' Requires \code{\link{AEGIS_compute_qc_metrics}} to have been run first.
#' If \code{remove_doublets = TRUE}, also requires
#' \code{\link{AEGIS_detect_doublets}}.
#'
#' @param sce A \code{SingleCellExperiment} with QC columns in
#'   \code{colData}.
#' @param mode Character. \code{"adaptive"} (default) or \code{"fixed"}.
#' @param n_mads Numeric. Number of MADs for adaptive outlier detection.
#'   Default \code{3}.
#' @param min_counts,max_counts Numeric. Hard bounds on total counts
#'   (\code{sum}).
#' @param min_features,max_features Numeric. Hard bounds on detected genes
#'   (\code{detected}).
#' @param max_pct_mt Numeric. Hard ceiling on \% mitochondrial counts.
#'   Applied in both adaptive and fixed modes when provided.
#' @param remove_doublets Logical. Also remove cells classified as doublets by
#'   \code{scDblFinder.class}.  Default \code{TRUE}.
#' @param verbose Logical. Print a removal summary. Default \code{TRUE}.
#'
#' @return Filtered \code{SingleCellExperiment}.
#' @export
AEGIS_filter_cells <- function(sce,
                                mode            = c("adaptive", "fixed"),
                                n_mads          = 3,
                                min_counts      = NULL,
                                max_counts      = NULL,
                                min_features    = NULL,
                                max_features    = NULL,
                                max_pct_mt      = NULL,
                                remove_doublets = TRUE,
                                verbose         = TRUE) {

  mode <- match.arg(mode)

  cd <- colData(sce)

  # ── Prerequisite checks ──────────────────────────────────────────────────────
  if (!"sum" %in% names(cd))
    stop("QC metrics not found. Run AEGIS_compute_qc_metrics() first.",
         call. = FALSE)

  if (isTRUE(remove_doublets) && !"scDblFinder.class" %in% names(cd)) {
    warning("scDblFinder.class not found in colData; skipping doublet removal. ",
            "Run AEGIS_detect_doublets() first, or set remove_doublets = FALSE.",
            call. = FALSE)
    remove_doublets <- FALSE
  }

  n_before <- ncol(sce)
  discard  <- rep(FALSE, n_before)
  reasons  <- list()

  # ── Adaptive mode ────────────────────────────────────────────────────────────
  if (mode == "adaptive") {
    sub_fields <- intersect("subsets_mt_percent", names(cd))

    filters <- scuttle::perCellQCFilters(
      as.data.frame(cd),
      sum.field      = "sum",
      detected.field = "detected",
      sub.fields     = if (length(sub_fields)) sub_fields else NULL,
      nmads          = n_mads
    )

    for (col in setdiff(names(filters), "discard")) {
      mask          <- as.logical(filters[[col]])
      reasons[[col]] <- sum(mask, na.rm = TRUE)
      discard       <- discard | mask
    }
  }

  # ── Fixed mode ───────────────────────────────────────────────────────────────
  if (mode == "fixed") {
    if (!is.null(min_counts)) {
      mask                      <- cd$sum < min_counts
      reasons[["low_counts"]]   <- sum(mask, na.rm = TRUE)
      discard                   <- discard | mask
    }
    if (!is.null(max_counts)) {
      mask                      <- cd$sum > max_counts
      reasons[["high_counts"]]  <- sum(mask, na.rm = TRUE)
      discard                   <- discard | mask
    }
    if (!is.null(min_features)) {
      mask                      <- cd$detected < min_features
      reasons[["low_features"]] <- sum(mask, na.rm = TRUE)
      discard                   <- discard | mask
    }
    if (!is.null(max_features)) {
      mask                       <- cd$detected > max_features
      reasons[["high_features"]] <- sum(mask, na.rm = TRUE)
      discard                    <- discard | mask
    }
  }

  # ── Hard MT ceiling (both modes) ─────────────────────────────────────────────
  if (!is.null(max_pct_mt) && "subsets_mt_percent" %in% names(cd)) {
    mask                   <- cd$subsets_mt_percent > max_pct_mt
    reasons[["hard_mt_cap"]] <- sum(mask, na.rm = TRUE)
    discard                <- discard | mask
  }

  # ── Doublet removal ───────────────────────────────────────────────────────────
  if (isTRUE(remove_doublets)) {
    mask                 <- cd$scDblFinder.class == "doublet"
    reasons[["doublets"]] <- sum(mask, na.rm = TRUE)
    discard              <- discard | mask
  }

  sce_filtered <- sce[, !discard]

  if (isTRUE(verbose))
    .aegis_print_filter_summary(n_before, ncol(sce_filtered), reasons, mode,
                                 n_mads)

  sce_filtered
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.aegis_resolve_assay <- function(sce, assay_name) {
  an <- assayNames(sce)
  if (!is.null(assay_name)) {
    if (!assay_name %in% an)
      stop("Assay '", assay_name, "' not found. Available: ",
           paste(an, collapse = ", "), call. = FALSE)
    return(assay_name)
  }
  for (p in c("counts", "X", "logcounts")) if (p %in% an) return(p)
  message("Note: using assay '", an[1], "' for QC metrics.")
  an[1]
}


.aegis_print_metrics_summary <- function(sce, is_mito, is_ribo, assay_name) {
  bar <- strrep("\u2500", 56)
  cd  <- colData(sce)

  fmt_range <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return("N/A")
    sprintf("%s \u2013 %s", format(round(min(x), 1), big.mark = ","),
            format(round(max(x), 1), big.mark = ","))
  }

  message("\u2500\u2500 AEGIS: QC metrics computed (assay: ", assay_name, ") ", bar)
  message(sprintf("  Cells              : %s", format(ncol(sce), big.mark = ",")))
  message(sprintf("  Mito genes found   : %s", sum(is_mito)))
  message(sprintf("  Ribo genes found   : %s", sum(is_ribo)))
  message(sprintf("  Total counts (sum) : %s", fmt_range(cd$sum)))
  message(sprintf("  Genes detected     : %s", fmt_range(cd$detected)))
  if ("subsets_mt_percent" %in% names(cd))
    message(sprintf("  MT %%               : %s", fmt_range(cd$subsets_mt_percent)))
  if ("subsets_ribo_percent" %in% names(cd))
    message(sprintf("  Ribo %%             : %s", fmt_range(cd$subsets_ribo_percent)))
  message(bar)
}


# Compute per-cell QC metrics when the counts assay is a BPCells IterableMatrix.
# Adds the same colData columns that scuttle::addPerCellQCMetrics() would add:
#   sum, detected, subsets_mt_sum, subsets_mt_detected, subsets_mt_percent,
#   subsets_ribo_sum, subsets_ribo_detected, subsets_ribo_percent
.aegis_bpcells_qc_metrics <- function(sce, assay_name, is_mito, is_ribo) {
  mat      <- assay(sce, assay_name)
  col_sums <- colSums(mat)
  detected <- BPCells::matrix_stats(mat, col_stats = "nonzero")$col_stats["nonzero", ]

  colData(sce)$sum      <- col_sums
  colData(sce)$detected <- detected

  if (any(is_mito)) {
    mito_mat <- mat[is_mito, , drop = FALSE]
    mito_sum <- colSums(mito_mat)
    mito_det <- BPCells::matrix_stats(mito_mat, col_stats = "nonzero")$col_stats["nonzero", ]
    colData(sce)$subsets_mt_sum      <- mito_sum
    colData(sce)$subsets_mt_detected <- mito_det
    colData(sce)$subsets_mt_percent  <- ifelse(col_sums > 0,
                                               mito_sum / col_sums * 100, 0)
  }

  if (any(is_ribo)) {
    ribo_mat <- mat[is_ribo, , drop = FALSE]
    ribo_sum <- colSums(ribo_mat)
    ribo_det <- BPCells::matrix_stats(ribo_mat, col_stats = "nonzero")$col_stats["nonzero", ]
    colData(sce)$subsets_ribo_sum      <- ribo_sum
    colData(sce)$subsets_ribo_detected <- ribo_det
    colData(sce)$subsets_ribo_percent  <- ifelse(col_sums > 0,
                                                  ribo_sum / col_sums * 100, 0)
  }

  sce
}


.aegis_print_filter_summary <- function(n_before, n_after, reasons, mode,
                                         n_mads) {
  bar      <- strrep("\u2500", 56)
  n_rm     <- n_before - n_after
  pct_rm   <- round(100 * n_rm / n_before, 1)
  mode_str <- if (mode == "adaptive") paste0("adaptive (nmads = ", n_mads, ")")
              else "fixed"

  message("\u2500\u2500 AEGIS: Cell filtering (", mode_str, ") ", bar)
  message(sprintf("  Before   : %s cells", format(n_before, big.mark = ",")))
  message(sprintf("  Removed  : %s cells (%s%%)",
                  format(n_rm, big.mark = ","), pct_rm))
  message(sprintf("  Retained : %s cells", format(n_after, big.mark = ",")))
  if (length(reasons)) {
    message("  Reasons (non-exclusive):")
    for (nm in names(reasons))
      message(sprintf("    %-32s: %s", nm, format(reasons[[nm]], big.mark = ",")))
  }
  message(bar)
}
