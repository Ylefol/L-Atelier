###############################################################################
########################### Normalization #####################################
###############################################################################
# Omics-specific normalization methods. Currently mass spec normalization;
# this is where CPM/TPM/RPKM/quantile/TMM normalization would also live.

#' Normalize a Log2 Mass Spectrometry Intensity Matrix
#'
#' @description Normalizes a log2-scale DIA mass spectrometry intensity matrix
#' using either Variance Stabilizing Normalization (VSN, default) or per-sample
#' median centering.
#'
#' @param wide Numeric matrix, proteins (rows) × samples (cols). Values must be
#'   log2-transformed (as produced by \code{ELEUTHIA_load_massspec()}). Both
#'   methods expect log2 input; VSN back-transforms to raw intensities internally
#'   before fitting.
#' @param method Character. Normalization method: \code{"vsn"} (default) or
#'   \code{"median_centering"}.
#' @param verbose Logical. Print normalization summary. Default: TRUE.
#'
#' @return A numeric matrix of the same dimensions, normalized. NA values are
#'   preserved in their original positions.
#'
#' @details
#' \strong{VSN} (Huber et al. 2002, Bioinformatics): fits a generalized
#' log (glog/arsinh) calibration per sample via robust iterative regression on
#' the raw intensity scale. Simultaneously normalizes cross-sample shifts AND
#' stabilizes variance across the full intensity range (heteroscedasticity
#' correction) — proteins at low and high abundance are treated equally. Output
#' is on a log2-like scale (the "h" scale). Because the pipeline stores log2
#' data in \code{ms$wide}, VSN mode back-transforms (\code{2^wide}) before
#' fitting and re-applies the original NA mask after, since \code{vsn::justvsn()}
#' fitting is driven by observed values only.
#' Requires the \code{vsn} Bioconductor package:
#' \code{BiocManager::install("vsn")}.
#'
#' \strong{Median centering}: subtracts the per-sample median (non-NA values)
#' from each column. Removes run-to-run location shifts on the log2 scale but
#' does not correct heteroscedasticity. Retained as a lightweight fallback or
#' for comparison.
#'
#' For the analysis pipeline:
#' \preformatted{
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide, method = "vsn")
#' }
POSEIDON_normalize_massspec <- function(wide,
                                         method  = c("vsn", "median_centering"),
                                         verbose = TRUE) {

  if (!is.matrix(wide))
    stop("wide must be a numeric matrix (proteins x samples).")

  method <- match.arg(method)

  if (method == "vsn") {

    if (!requireNamespace("vsn", quietly = TRUE))
      stop("Package 'vsn' is required for VSN normalization.\n",
           "Install via: BiocManager::install('vsn')", call. = FALSE)

    na_mask <- is.na(wide)
    raw     <- 2^wide

    vsn_fit  <- vsn::justvsn(raw)
    norm_mat <- as.matrix(vsn_fit)
    dimnames(norm_mat) <- dimnames(wide)
    norm_mat[na_mask]  <- NA

    if (verbose) {
      cat("[POSEIDON] VSN normalization (mass spec):\n")
      cat("    Samples normalized:", ncol(norm_mat), "\n")
      cat("    Proteins:          ", nrow(norm_mat), "\n")
      cat("    NA values preserved:", sum(na_mask), "\n")
    }

  } else {

    col_medians <- apply(wide, 2L, stats::median, na.rm = TRUE)

    if (any(is.na(col_medians)))
      warning("Some samples have all-NA values; their medians are NA and they will not be centered.")

    norm_mat <- sweep(wide, 2L, col_medians, FUN = "-")

    if (verbose) {
      cat("[POSEIDON] Median centering (mass spec):\n")
      cat("    Samples normalized:", sum(!is.na(col_medians)), "\n")
      cat("    Median shift range: [",
          round(min(col_medians, na.rm = TRUE), 2), ", ",
          round(max(col_medians, na.rm = TRUE), 2), "]\n", sep = "")
    }

  }

  return(norm_mat)
}
