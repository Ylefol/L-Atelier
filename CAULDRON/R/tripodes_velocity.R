# ==============================================================================
# TRIPODES - RNA Velocity Functions
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Functions for spliced/unspliced QC and RNA velocity estimation via scVelo,
# running in a CAULDRON-managed basilisk environment (Python 3.10,
# numpy==1.23.5, scvelo==0.2.5).  All Python logic lives in
# inst/python/scvelo_run.py — no velociraptor dependency.
#
# Functions:
#   TRIPODES_check_velocity_ready() — fast SCE validation
#   TRIPODES_assess_splicing()      — per-cell/gene splicing QC (Step 2)
#   TRIPODES_run_velocity()         — scVelo runner (Step 3)
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

.tripodes_msg <- function(..., type = "info") {
  prefix <- switch(type,
    info    = "[TRIPODES]",
    warn    = "[TRIPODES] Warning:",
    ok      = "[TRIPODES] \u2713",
    .       = "[TRIPODES]"
  )
  message(prefix, " ", paste0(...))
}

# Known assay naming conventions for spliced / unspliced matrices.
# Each entry: list(spliced = <name>, unspliced = <name>, note = <explanation>)
.tripodes_known_conventions <- list(
  canonical = list(
    spliced   = "spliced",
    unspliced = "unspliced",
    note      = NULL   # this is the required convention — no rename needed
  ),
  scvelo_moments = list(
    spliced   = "Ms",
    unspliced = "Mu",
    note      = paste0(
      "'Ms' and 'Mu' are scVelo moment (smoothed) matrices, not raw counts. ",
      "RNA velocity in CAULDRON requires raw spliced/unspliced count matrices. ",
      "Load the raw counts (e.g. from a STARsolo velocyto loom file or alevin-fry ",
      "output) and name them 'spliced' and 'unspliced'."
    )
  ),
  counts_prefixed = list(
    spliced   = "counts_spliced",
    unspliced = "counts_unspliced",
    note      = paste0(
      "Found 'counts_spliced'/'counts_unspliced'. CAULDRON expects 'spliced'/'unspliced'. ",
      "Rename with: assay(sce, 'spliced') <- assay(sce, 'counts_spliced'); ",
      "assay(sce, 'unspliced') <- assay(sce, 'counts_unspliced')"
    )
  )
)


# ------------------------------------------------------------------------------
#' Check whether an SCE object is ready for RNA velocity analysis
#'
#' Validates that a SingleCellExperiment contains the spliced and unspliced
#' count assays required by \code{TRIPODES_run_velocity}. Detects alternate or
#' incompatible naming conventions and reports actionable guidance. Also checks
#' data type (must be raw integer counts) and reports summary splicing metrics.
#'
#' This function is non-destructive — it does \strong{not} modify the SCE.
#' Run it before \code{TRIPODES_assess_splicing} or
#' \code{TRIPODES_run_velocity} to confirm your object is usable.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param spliced_assay   Name of the spliced counts assay (default
#'   \code{"spliced"}).
#' @param unspliced_assay Name of the unspliced counts assay (default
#'   \code{"unspliced"}).
#' @param verbose Logical; print diagnostic report (default \code{TRUE}).
#'
#' @return Invisibly returns a named list:
#'   \describe{
#'     \item{ready}{Logical — \code{TRUE} only if both assays are present, have
#'       matching dimensions, contain integer-like counts, and no threshold
#'       warnings were triggered.}
#'     \item{spliced_assay, unspliced_assay}{Confirmed assay names, or
#'       \code{NA} if not found.}
#'     \item{n_cells, n_genes}{Dimensions of the SCE.}
#'     \item{median_spliced, median_unspliced}{Median total UMI per cell for
#'       each modality (\code{NA} if assays absent).}
#'     \item{splicing_ratio}{Global splicing ratio: sum(spliced) /
#'       (sum(spliced) + sum(unspliced)); \code{NA} if assays absent.}
#'     \item{issues}{Character vector of all warnings/errors found (empty if
#'       \code{ready = TRUE}).}
#'   }
#'
#' @export
TRIPODES_check_velocity_ready <- function(sce,
                                           spliced_assay   = "spliced",
                                           unspliced_assay = "unspliced",
                                           verbose         = TRUE) {

  avail    <- assayNames(sce)
  issues   <- character(0L)   # hard errors — block ready=TRUE
  warnings <- character(0L)   # soft warnings — inform only

  # ------------------------------------------------------------------
  # 1. Detect alternate / incompatible naming conventions
  # ------------------------------------------------------------------
  for (conv_name in names(.tripodes_known_conventions)) {
    conv <- .tripodes_known_conventions[[conv_name]]
    if (is.null(conv$note)) next   # skip canonical — handled below

    found_s <- conv$spliced   %in% avail
    found_u <- conv$unspliced %in% avail

    if (found_s || found_u) {
      found_names <- c(
        if (found_s) conv$spliced   else character(0L),
        if (found_u) conv$unspliced else character(0L)
      )
      issues <- c(issues, paste0(
        "Detected '", paste(found_names, collapse = "' / '"), "': ", conv$note
      ))
    }
  }

  # ------------------------------------------------------------------
  # 2. Check for canonical assays
  # ------------------------------------------------------------------
  has_spliced   <- spliced_assay   %in% avail
  has_unspliced <- unspliced_assay %in% avail

  if (!has_spliced || !has_unspliced) {
    missing_desc <- c(
      if (!has_spliced)   paste0("'", spliced_assay,   "' (spliced counts)")   else NULL,
      if (!has_unspliced) paste0("'", unspliced_assay, "' (unspliced counts)") else NULL
    )
    issues <- c(issues, paste0(
      "Required assay(s) not found: ",
      paste(missing_desc, collapse = ", "), ". ",
      "Available assays: ", paste(avail, collapse = ", "), "."
    ))

    if (verbose) {
      .tripodes_msg("RNA velocity readiness check: FAILED", type = "warn")
      for (iss in issues) message("  [!] ", iss)
    }

    return(invisible(list(
      ready            = FALSE,
      spliced_assay    = NA_character_,
      unspliced_assay  = NA_character_,
      n_cells          = ncol(sce),
      n_genes          = nrow(sce),
      median_spliced   = NA_real_,
      median_unspliced = NA_real_,
      splicing_ratio   = NA_real_,
      issues           = issues
    )))
  }

  # ------------------------------------------------------------------
  # 3. Canonical assays found — run metric checks
  # ------------------------------------------------------------------
  mat_s <- assay(sce, spliced_assay)
  mat_u <- assay(sce, unspliced_assay)
  n_cells <- ncol(sce)
  n_genes <- nrow(sce)

  # Dimension match
  if (!identical(dim(mat_s), dim(mat_u))) {
    issues <- c(issues, paste0(
      "Spliced and unspliced matrices have different dimensions: ",
      paste(dim(mat_s), collapse = "x"), " vs ",
      paste(dim(mat_u), collapse = "x"), "."
    ))
  }

  # Data type — should be raw integer counts, not log-normalised
  .is_continuous <- function(mat) {
    if (inherits(mat, "IterableMatrix")) return(FALSE)  # BPCells: trust raw load
    samp <- if (is(mat, "sparseMatrix")) {
      mat@x[seq_len(min(500L, length(mat@x)))]
    } else {
      as.vector(mat[seq_len(min(10L, nrow(mat))), seq_len(min(10L, ncol(mat)))])
    }
    samp <- samp[samp != 0]
    length(samp) > 0L && any(samp != floor(samp))
  }

  if (.is_continuous(mat_s) || .is_continuous(mat_u)) {
    issues <- c(issues, paste0(
      "Spliced or unspliced matrix appears to contain non-integer (continuous) ",
      "values. RNA velocity requires raw count matrices, not log-normalised data."
    ))
  }


  # ------------------------------------------------------------------
  # 4. Compute summary metrics (BPCells-aware)
  # ------------------------------------------------------------------
  if (inherits(mat_s, "IterableMatrix")) {
    cs_s <- BPCells::matrix_stats(mat_s, col_stats = "sum")$col_stats["sum", ]
    cs_u <- BPCells::matrix_stats(mat_u, col_stats = "sum")$col_stats["sum", ]
  } else {
    cs_s <- Matrix::colSums(mat_s)
    cs_u <- Matrix::colSums(mat_u)
  }

  median_spliced   <- stats::median(cs_s)
  median_unspliced <- stats::median(cs_u)
  total_s          <- sum(cs_s)
  total_u          <- sum(cs_u)
  splicing_ratio   <- if ((total_s + total_u) > 0) total_s / (total_s + total_u) else NA_real_

  # Soft threshold warnings — inform only, do not block ready
  if (!is.na(median_unspliced) && median_unspliced < 0.5) {
    warnings <- c(warnings, paste0(
      "Median unspliced UMI per cell is very low (",
      round(median_unspliced, 2),
      "). Velocity estimates may be unreliable. Verify that unspliced ",
      "counts were correctly generated at alignment (STARsolo/alevin-fry with ",
      "velocity output enabled)."
    ))
  }

  if (!is.na(splicing_ratio)) {
    if (splicing_ratio < 0.30) {
      warnings <- c(warnings, paste0(
        "Global splicing ratio is low (", round(splicing_ratio * 100L, 1L), "% spliced / ",
        round((1 - splicing_ratio) * 100L, 1L), "% unspliced). ",
        "This can be expected for certain capture chemistries (e.g. Parse Biosciences ",
        "random hexamer priming) but is unusual for oligo-dT-based protocols. ",
        "Verify this is consistent with your library preparation method before proceeding."
      ))
    } else if (splicing_ratio > 0.95) {
      warnings <- c(warnings, paste0(
        "Global splicing ratio is high (", round(splicing_ratio * 100L, 1L), "% spliced). ",
        "Very few unspliced reads were detected. This may be expected depending on your ",
        "protocol, but could also indicate that unspliced counts were not correctly ",
        "captured at alignment. Verify before proceeding."
      ))
    }
  }

  # ------------------------------------------------------------------
  # 5. Report
  # ------------------------------------------------------------------
  # ready = TRUE as long as there are no hard errors (missing assays,
  # dim mismatch, continuous data). Soft threshold warnings are surfaced
  # but do not block the pipeline — the user decides how to proceed.
  ready <- length(issues) == 0L

  if (verbose) {
    status_label <- if (ready) "PASSED" else "FAILED"
    .tripodes_msg(paste0("RNA velocity readiness check: ", status_label))
    .tripodes_msg(paste0("  Cells:                 ", n_cells))
    .tripodes_msg(paste0("  Genes:                 ", n_genes))
    .tripodes_msg(paste0("  Median spliced UMI:    ", round(median_spliced,   1L)))
    .tripodes_msg(paste0("  Median unspliced UMI:  ", round(median_unspliced, 1L)))
    if (!is.na(splicing_ratio)) {
      .tripodes_msg(paste0(
        "  Global splicing ratio: ",
        round(splicing_ratio * 100L, 1L), "% spliced / ",
        round((1 - splicing_ratio) * 100L, 1L), "% unspliced"
      ))
    }
    if (length(issues) > 0L) {
      message("")
      for (iss in issues) message("  [!] ", iss)
    }
    if (length(warnings) > 0L) {
      message("")
      for (w in warnings) message("  [~] ", w)
    }
  }

  invisible(list(
    ready            = ready,
    spliced_assay    = spliced_assay,
    unspliced_assay  = unspliced_assay,
    n_cells          = n_cells,
    n_genes          = n_genes,
    median_spliced   = median_spliced,
    median_unspliced = median_unspliced,
    splicing_ratio   = splicing_ratio,
    issues           = issues,
    warnings         = warnings
  ))
}


# ------------------------------------------------------------------------------
#' Assess spliced/unspliced read distribution across cells and genes
#'
#' Computes per-cell and per-gene splicing metrics and stores them in the SCE.
#' This is the detailed QC step to run after
#' \code{\link{TRIPODES_check_velocity_ready}} and before
#' \code{\link{TRIPODES_run_velocity}}. Results are consumed by
#' \code{ASPIS_plot_splicing} for visualisation.
#'
#' Per-cell metrics added to \code{colData}:
#' \describe{
#'   \item{splicing_ratio}{Spliced UMI / (spliced + unspliced) UMI per cell.}
#'   \item{total_spliced}{Total spliced UMI count per cell.}
#'   \item{total_unspliced}{Total unspliced UMI count per cell.}
#' }
#'
#' Per-gene metrics added to \code{rowData}:
#' \describe{
#'   \item{mean_spliced}{Mean spliced counts across cells.}
#'   \item{mean_unspliced}{Mean unspliced counts across cells.}
#'   \item{gene_splicing_ratio}{mean_spliced / (mean_spliced + mean_unspliced).}
#'   \item{velocity_reliable}{Logical — gene has sufficient unspliced signal for
#'     reliable velocity estimation (see \code{min_unspliced_mean} and
#'     \code{ratio_range}).}
#' }
#'
#' @param sce A \code{SingleCellExperiment}; must pass
#'   \code{TRIPODES_check_velocity_ready}.
#' @param spliced_assay   Name of the spliced counts assay (default
#'   \code{"spliced"}).
#' @param unspliced_assay Name of the unspliced counts assay (default
#'   \code{"unspliced"}).
#' @param min_unspliced_mean Minimum mean unspliced counts for a gene to be
#'   flagged \code{velocity_reliable = TRUE} (default \code{0.01}).
#' @param ratio_range Numeric vector of length 2; genes whose
#'   \code{gene_splicing_ratio} falls outside this range are flagged as
#'   unreliable (default \code{c(0.05, 0.99)}). Genes with near-zero or
#'   near-total unspliced fractions have poorly constrained dynamics.
#' @param verbose Logical; print summary report (default \code{TRUE}).
#'
#' @return The input SCE with per-cell metrics added to \code{colData} and
#'   per-gene metrics added to \code{rowData}.
#'   \code{metadata(sce)$splicing_assessed} is set to \code{TRUE}.
#'
#' @export
TRIPODES_assess_splicing <- function(sce,
                                      spliced_assay      = "spliced",
                                      unspliced_assay    = "unspliced",
                                      min_unspliced_mean = 0.01,
                                      ratio_range        = c(0.05, 0.99),
                                      verbose            = TRUE) {

  # ------------------------------------------------------------------
  # 1. Validate — must pass check_velocity_ready first
  # ------------------------------------------------------------------
  chk <- TRIPODES_check_velocity_ready(sce,
                                        spliced_assay   = spliced_assay,
                                        unspliced_assay = unspliced_assay,
                                        verbose         = FALSE)
  if (!chk$ready) {
    stop(
      "[TRIPODES] SCE did not pass velocity readiness check. ",
      "Run TRIPODES_check_velocity_ready(sce, verbose=TRUE) for details."
    )
  }

  mat_s <- assay(sce, spliced_assay)
  mat_u <- assay(sce, unspliced_assay)
  is_bpcells <- inherits(mat_s, "IterableMatrix")

  if (verbose) .tripodes_msg("Assessing spliced/unspliced distribution...")

  # ------------------------------------------------------------------
  # 2. Per-cell metrics
  # ------------------------------------------------------------------
  if (is_bpcells) {
    cs_s <- BPCells::matrix_stats(mat_s, col_stats = "sum")$col_stats["sum", ]
    cs_u <- BPCells::matrix_stats(mat_u, col_stats = "sum")$col_stats["sum", ]
  } else {
    cs_s <- Matrix::colSums(mat_s)
    cs_u <- Matrix::colSums(mat_u)
  }

  cell_total      <- cs_s + cs_u
  cell_ratio      <- ifelse(cell_total > 0, cs_s / cell_total, NA_real_)

  colData(sce)$total_spliced   <- as.numeric(cs_s)
  colData(sce)$total_unspliced <- as.numeric(cs_u)
  colData(sce)$splicing_ratio  <- as.numeric(cell_ratio)

  # ------------------------------------------------------------------
  # 3. Per-gene metrics
  # ------------------------------------------------------------------
  n_cells <- ncol(sce)

  if (is_bpcells) {
    rs_s <- BPCells::matrix_stats(mat_s, row_stats = "mean")$row_stats["mean", ]
    rs_u <- BPCells::matrix_stats(mat_u, row_stats = "mean")$row_stats["mean", ]
  } else {
    rs_s <- Matrix::rowMeans(mat_s)
    rs_u <- Matrix::rowMeans(mat_u)
  }

  gene_total <- rs_s + rs_u
  gene_ratio <- ifelse(gene_total > 0, rs_s / gene_total, NA_real_)

  velocity_reliable <- (
    !is.na(gene_ratio) &
    rs_u >= min_unspliced_mean &
    gene_ratio >= ratio_range[1] &
    gene_ratio <= ratio_range[2]
  )

  rowData(sce)$mean_spliced        <- as.numeric(rs_s)
  rowData(sce)$mean_unspliced      <- as.numeric(rs_u)
  rowData(sce)$gene_splicing_ratio <- as.numeric(gene_ratio)
  rowData(sce)$velocity_reliable   <- velocity_reliable

  # ------------------------------------------------------------------
  # 4. Summary report
  # ------------------------------------------------------------------
  n_reliable     <- sum(velocity_reliable, na.rm = TRUE)
  n_genes        <- nrow(sce)
  pct_reliable   <- round(n_reliable / n_genes * 100, 1)
  med_cell_ratio <- stats::median(cell_ratio, na.rm = TRUE)
  pct_low_ratio  <- round(mean(cell_ratio < 0.10, na.rm = TRUE) * 100, 1)

  if (verbose) {
    .tripodes_msg("Splicing assessment complete.")
    .tripodes_msg(paste0("  Cells assessed:          ", ncol(sce)))
    .tripodes_msg(paste0("  Genes assessed:          ", n_genes))
    .tripodes_msg(paste0(
      "  Median cell splicing ratio: ",
      round(med_cell_ratio * 100, 1), "% spliced"
    ))
    .tripodes_msg(paste0(
      "  Cells with <10% spliced:    ", pct_low_ratio, "%",
      if (pct_low_ratio > 20) " [notable — review splicing_ratio distribution]" else ""
    ))
    .tripodes_msg(paste0(
      "  Velocity-reliable genes:    ",
      n_reliable, " / ", n_genes, " (", pct_reliable, "%)"
    ))
    if (pct_reliable < 20) {
      message("  [!] Fewer than 20% of genes are flagged as velocity-reliable. ",
              "Consider reviewing alignment parameters or unspliced count generation.")
    }
  }

  metadata(sce)$splicing_assessed <- TRUE
  sce
}


# ------------------------------------------------------------------------------
#' Estimate RNA velocity via scVelo
#'
#' Runs RNA velocity estimation using the Python \pkg{scVelo} library,
#' called directly through a \pkg{basilisk}-managed environment and a bundled
#' script (\code{inst/python/scvelo_run.py}) — not via the \pkg{velociraptor}
#' Bioconductor wrapper (CAULDRON has no dependency on \pkg{velociraptor}).
#' Spliced and unspliced count assays must be present in the SCE. Results are
#' stored back into the input SCE and returned.
#'
#' If \code{TRIPODES_assess_splicing} has not been run, a warning is issued but
#' execution continues. If the SCE fails the hard readiness checks (missing
#' assays, dimension mismatch, continuous data), execution stops with a message
#' directing the user to \code{TRIPODES_check_velocity_ready}.
#'
#' Stored results:
#' \describe{
#'   \item{\code{assay(sce, "velocity")}}{Gene × cell velocity matrix.}
#'   \item{\code{assay(sce, "Ms")}, \code{assay(sce, "Mu")}}{Moments-smoothed
#'     spliced/unspliced counts (gene × cell), as fit against by
#'     \code{scv.pp.moments()}. Only stored when \code{store_moments = TRUE}
#'     (default \code{FALSE}) — these are full dense gene × cell matrices and
#'     substantially increase the SCE's size.}
#'   \item{\code{colData(sce)$velocity_confidence}}{Per-cell velocity
#'     confidence score from scVelo.}
#'   \item{\code{colData(sce)$velocity_length}}{Per-cell velocity vector
#'     magnitude.}
#'   \item{\code{reducedDim(sce, "velocity_<use_dimred>")}}{Projected velocity
#'     coordinates on the chosen embedding (if \code{use_dimred} is present in
#'     the SCE).}
#' }
#'
#' @param sce A \code{SingleCellExperiment}; must contain spliced and unspliced
#'   count assays.
#' @param spliced_assay   Name of the spliced counts assay (default
#'   \code{"spliced"}).
#' @param unspliced_assay Name of the unspliced counts assay (default
#'   \code{"unspliced"}).
#' @param mode scVelo fitting mode. One of \code{"deterministic"} (default,
#'   fast first-order OLS fit), \code{"stochastic"} (accounts for
#'   transcriptional noise via second-order moments; requires more memory), or
#'   \code{"dynamical"} (most accurate, substantially slower — fits full
#'   kinetic model per gene via \code{scv.tl.recover_dynamics()}, run
#'   automatically first when this mode is selected).
#' @param use_dimred Name of the PCA-like dimensionality reduction in the SCE
#'   to pass to scVelo for neighbour graph construction (default \code{"PCA"}).
#'   This should be the same reduction used to build your UMAP, ensuring that
#'   velocity neighbours are consistent with the embedding. \code{n_pcs} is
#'   automatically capped to the number of components available. Set to
#'   \code{NULL} to let scVelo compute its own PCA internally (not recommended
#'   if you intend to project velocity onto an existing UMAP).
#' @param n_pcs Number of PCs to use for neighbour graph construction within
#'   scVelo (default \code{30L}).
#' @param n_neighbors Number of neighbours for the scVelo graph (default
#'   \code{30L}).
#' @param store_moments Logical; also store the moments-smoothed spliced/
#'   unspliced matrices as \code{assay(sce, "Ms")}/\code{assay(sce, "Mu")}
#'   (default \code{FALSE}). These are computed internally by scVelo
#'   regardless (needed to fit velocity itself); this only controls whether
#'   they're additionally returned and stored, since doing so roughly doubles
#'   the SCE's assay footprint. Enable for phase-portrait-style diagnostics
#'   (unspliced vs. spliced per gene).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#' @param ... Currently unused — kept for backward-compatible signature
#'   stability; nothing in the current implementation forwards these
#'   arguments anywhere.
#'
#' @return The input SCE with velocity results added (see Details).
#'   \code{metadata(sce)$velocity_run} is set to \code{TRUE} and
#'   \code{metadata(sce)$velocity_mode} records the mode used.
#'
#' @export
TRIPODES_run_velocity <- function(sce,
                                   spliced_assay   = "spliced",
                                   unspliced_assay = "unspliced",
                                   mode            = c("deterministic",
                                                       "stochastic",
                                                       "dynamical"),
                                   use_dimred      = "PCA",
                                   n_pcs           = 30L,
                                   n_neighbors     = 30L,
                                   store_moments   = FALSE,
                                   verbose         = TRUE,
                                   ...) {

  mode <- match.arg(mode)

  # ------------------------------------------------------------------
  # 1. Hard readiness check — errors if assays missing / wrong type
  # ------------------------------------------------------------------
  chk <- TRIPODES_check_velocity_ready(sce,
                                        spliced_assay   = spliced_assay,
                                        unspliced_assay = unspliced_assay,
                                        verbose         = FALSE)
  if (!chk$ready) {
    stop(
      "[TRIPODES] SCE did not pass velocity readiness check. ",
      "Run TRIPODES_check_velocity_ready(sce, verbose=TRUE) for details."
    )
  }

  # Soft warning if splicing assessment was not run
  if (!isTRUE(metadata(sce)$splicing_assessed)) {
    message(
      "[TRIPODES] [~] TRIPODES_assess_splicing() has not been run. ",
      "Consider running it first to review per-cell and per-gene splicing ",
      "metrics before estimating velocity."
    )
  }

  # ------------------------------------------------------------------
  # 3. Embedding check + n_pcs cap
  # ------------------------------------------------------------------
  dimred_names  <- reducedDimNames(sce)
  use_dimred_ok <- !is.null(use_dimred) && use_dimred %in% dimred_names

  if (!is.null(use_dimred) && !use_dimred_ok) {
    message(
      "[TRIPODES] [~] Dimensionality reduction '", use_dimred, "' not found in SCE. ",
      "Available: ", paste(dimred_names, collapse = ", "), ". ",
      "scVelo will compute its own PCA — velocity neighbours may not match your embedding."
    )
  }

  # Cap n_pcs to the number of components actually present in the PCA
  if (use_dimred_ok) {
    n_pcs_avail <- ncol(reducedDim(sce, use_dimred))
    if (n_pcs > n_pcs_avail) {
      message(
        "[TRIPODES] [~] n_pcs (", n_pcs, ") exceeds available components in '",
        use_dimred, "' (", n_pcs_avail, "). Capping to ", n_pcs_avail, "."
      )
      n_pcs <- n_pcs_avail
    }
  }

  # ------------------------------------------------------------------
  # 4. Prepare SCE for velocity
  #    (a) subset to velocity-reliable genes to avoid degenerate matrices
  #    (b) coerce spliced/unspliced to dgCMatrix float — sparse column slices
  #        from other formats return 2D (n,1) arrays in numpy, which breaks
  #        scVelo's stochastic least-squares solver
  # ------------------------------------------------------------------
  sce_vel <- sce

  if (isTRUE(metadata(sce)$splicing_assessed) &&
      "velocity_reliable" %in% colnames(rowData(sce))) {
    reliable <- which(rowData(sce)$velocity_reliable)
    if (length(reliable) >= 50L) {
      sce_vel <- sce[reliable, ]
      if (verbose) {
        .tripodes_msg(paste0(
          "  Subsetting to ", length(reliable), " / ", nrow(sce),
          " velocity-reliable genes."
        ))
      }
    } else {
      message(
        "[TRIPODES] [~] Too few velocity-reliable genes (", length(reliable),
        ") to subset — proceeding with all genes. ",
        "Consider reviewing TRIPODES_assess_splicing() output."
      )
    }
  }

  # Coerce to dgCMatrix float
  .to_dgc_float <- function(mat) {
    if (inherits(mat, "IterableMatrix")) return(mat)   # BPCells: leave as-is
    mat <- methods::as(mat, "CsparseMatrix")           # any sparse → column-compressed
    mat <- methods::as(mat, "dgCMatrix")               # ensure double (not integer)
    mat
  }
  assay(sce_vel, spliced_assay)   <- .to_dgc_float(assay(sce_vel, spliced_assay))
  assay(sce_vel, unspliced_assay) <- .to_dgc_float(assay(sce_vel, unspliced_assay))

  # ------------------------------------------------------------------
  # 5. Run scVelo via basilisk
  # ------------------------------------------------------------------
  if (verbose) {
    .tripodes_msg(paste0("Running RNA velocity (mode = '", mode, "')..."))
    if (mode == "dynamical")
      .tripodes_msg("  Dynamical mode fits full kinetic model per gene — this may take a while.")
  }

  script  <- system.file("python", "scvelo_run.py", package = "CAULDRON")

  # Extract PCA coordinates for neighbour construction (cells × n_pcs)
  pca_mat <- if (use_dimred_ok) {
    m <- reducedDim(sce_vel, use_dimred)
    m[, seq_len(n_pcs), drop = FALSE]
  } else NULL

  # Extract UMAP / tSNE for velocity projection (cheap — done after graph is built)
  dimred_names_vel <- reducedDimNames(sce_vel)
  umap_mat <- if ("UMAP" %in% dimred_names_vel) reducedDim(sce_vel, "UMAP") else NULL
  tsne_mat <- if ("tSNE" %in% dimred_names_vel) reducedDim(sce_vel, "tSNE") else NULL

  velo_result <- basilisk::basiliskRun(
    env           = .cauldron_scvelo_env,
    fun           = .tripodes_scvelo_run,
    script        = script,
    spliced_mat   = assay(sce_vel, spliced_assay),
    unspliced_mat = assay(sce_vel, unspliced_assay),
    genes         = rownames(sce_vel),
    cells         = colnames(sce_vel),
    mode          = mode,
    n_pcs         = n_pcs,
    n_neighbors   = n_neighbors,
    store_moments = store_moments,
    pca_mat       = pca_mat,
    umap_mat      = umap_mat,
    tsne_mat      = tsne_mat
  )

  # ------------------------------------------------------------------
  # 6. Store results in original SCE
  # ------------------------------------------------------------------

  # Velocity matrix (genes × cells) — written back into the full gene space;
  # genes excluded by velocity_reliable subsetting receive NA rows.
  if (!is.null(velo_result$velocity)) {
    vel_mat_sub <- velo_result$velocity
    if (nrow(sce_vel) < nrow(sce)) {
      vel_mat_full <- matrix(
        NA_real_, nrow = nrow(sce), ncol = ncol(sce),
        dimnames = list(rownames(sce), colnames(sce))
      )
      vel_mat_full[rownames(sce_vel), ] <- vel_mat_sub
      assay(sce, "velocity") <- vel_mat_full
    } else {
      assay(sce, "velocity") <- vel_mat_sub
    }
    if (verbose) .tripodes_msg("  Stored: assay 'velocity'")
  } else {
    warning("[TRIPODES] Velocity matrix not returned — check scVelo logs.")
  }

  # Moments-smoothed spliced/unspliced (Ms/Mu) — opt-in via store_moments,
  # since these are full dense gene × cell matrices. Same gene-space padding
  # as the velocity matrix above, since both come from the same
  # velocity_reliable gene subset.
  if (isTRUE(store_moments)) {
    .store_moment_assay <- function(mat_sub, assay_name) {
      if (is.null(mat_sub)) {
        warning("[TRIPODES] '", assay_name, "' matrix not returned — check scVelo logs.")
        return(invisible(NULL))
      }
      if (nrow(sce_vel) < nrow(sce)) {
        mat_full <- matrix(
          NA_real_, nrow = nrow(sce), ncol = ncol(sce),
          dimnames = list(rownames(sce), colnames(sce))
        )
        mat_full[rownames(sce_vel), ] <- mat_sub
        assay(sce, assay_name) <<- mat_full
      } else {
        assay(sce, assay_name) <<- mat_sub
      }
      if (verbose) .tripodes_msg(paste0("  Stored: assay '", assay_name, "'"))
    }

    .store_moment_assay(velo_result$Ms, "Ms")
    .store_moment_assay(velo_result$Mu, "Mu")
  }

  if (!is.null(velo_result$confidence)) {
    colData(sce)$velocity_confidence <- as.numeric(velo_result$confidence)
    if (verbose) .tripodes_msg("  Stored: colData 'velocity_confidence'")
  }

  if (!is.null(velo_result$length)) {
    colData(sce)$velocity_length <- as.numeric(velo_result$length)
    if (verbose) .tripodes_msg("  Stored: colData 'velocity_length'")
  }

  velo_dimred_name <- paste0("velocity_", use_dimred)
  if (!is.null(velo_result$velocity_pca)) {
    rownames(velo_result$velocity_pca) <- colnames(sce)
    reducedDim(sce, velo_dimred_name) <- velo_result$velocity_pca
    if (verbose) .tripodes_msg(paste0("  Stored: reducedDim '", velo_dimred_name, "'"))
  }

  if (!is.null(velo_result$velocity_umap)) {
    rownames(velo_result$velocity_umap) <- colnames(sce)
    reducedDim(sce, "velocity_UMAP") <- velo_result$velocity_umap
    if (verbose) .tripodes_msg("  Stored: reducedDim 'velocity_UMAP'")
  }

  if (!is.null(velo_result$velocity_tsne)) {
    rownames(velo_result$velocity_tsne) <- colnames(sce)
    reducedDim(sce, "velocity_tSNE") <- velo_result$velocity_tsne
    if (verbose) .tripodes_msg("  Stored: reducedDim 'velocity_tSNE'")
  }

  # ------------------------------------------------------------------
  # 7. Finalise
  # ------------------------------------------------------------------
  metadata(sce)$velocity_run  <- TRUE
  metadata(sce)$velocity_mode <- mode

  if (verbose) {
    .tripodes_msg(paste0("RNA velocity complete (mode = '", mode, "')."))
    if ("velocity_confidence" %in% names(colData(sce))) {
      conf <- colData(sce)$velocity_confidence
      .tripodes_msg(paste0(
        "  Median velocity confidence: ", round(stats::median(conf, na.rm = TRUE), 3)
      ))
    }
  }

  # Everything needed from these has already been written into `sce` above.
  # Both can be large (sce_vel holds a coerced-float copy of the assays;
  # velo_result holds the dense matrices pulled back from Python, doubled up
  # again with store_moments = TRUE) and R does not proactively garbage
  # collect, so without this they linger as reported (but unreferenced)
  # memory for the rest of the session.
  rm(sce_vel, velo_result)
  gc(full = TRUE)

  sce
}


# ── Internal helpers ───────────────────────────────────────────────────────────
# Named package-level function so that the CAULDRON namespace is in scope
# (ensures reticulate:: resolves correctly inside basiliskRun).

.tripodes_scvelo_run <- function(script, spliced_mat, unspliced_mat,
                                  genes, cells, mode, n_pcs, n_neighbors,
                                  store_moments, pca_mat, umap_mat, tsne_mat) {
  main <- reticulate::import("__main__")
  reticulate::py_set_attr(main, "r_spliced",     spliced_mat)
  reticulate::py_set_attr(main, "r_unspliced",   unspliced_mat)
  reticulate::py_set_attr(main, "r_genes",        genes)
  reticulate::py_set_attr(main, "r_cells",        cells)
  reticulate::py_set_attr(main, "r_mode",         mode)
  reticulate::py_set_attr(main, "r_n_pcs",        as.integer(n_pcs))
  reticulate::py_set_attr(main, "r_n_neighbors",  as.integer(n_neighbors))
  reticulate::py_set_attr(main, "r_store_moments", isTRUE(store_moments))
  reticulate::py_set_attr(main, "r_pca",          pca_mat)    # NULL → Python None
  reticulate::py_set_attr(main, "r_umap",         umap_mat)
  reticulate::py_set_attr(main, "r_tsne",         tsne_mat)

  reticulate::py_run_file(script)

  list(
    velocity      = reticulate::py_to_r(reticulate::py_get_attr(main, "result_velocity")),
    Ms            = reticulate::py_to_r(reticulate::py_get_attr(main, "result_Ms")),
    Mu            = reticulate::py_to_r(reticulate::py_get_attr(main, "result_Mu")),
    confidence    = reticulate::py_to_r(reticulate::py_get_attr(main, "result_confidence")),
    length        = reticulate::py_to_r(reticulate::py_get_attr(main, "result_length")),
    velocity_pca  = reticulate::py_to_r(reticulate::py_get_attr(main, "result_velocity_pca")),
    velocity_umap = reticulate::py_to_r(reticulate::py_get_attr(main, "result_velocity_umap")),
    velocity_tsne = reticulate::py_to_r(reticulate::py_get_attr(main, "result_velocity_tsne"))
  )
}
