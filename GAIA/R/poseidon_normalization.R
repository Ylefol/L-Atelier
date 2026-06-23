###############################################################################
########################### Normalization #####################################
###############################################################################
# Omics-specific normalization methods. Currently mass spec median centering;
# this is where CPM/TPM/RPKM/quantile/TMM normalization would also live.

#' Median-Centre a Log2 Mass Spectrometry Intensity Matrix
#'
#' @description Applies per-sample median centering to a log2-scale intensity
#' matrix. Subtracts the per-sample median (computed on non-NA values) from each
#' column. This removes systematic run-to-run intensity shifts that persist after
#' DIA-NN MaxLFQ quantification, without distorting relative protein abundances.
#'
#' @param wide Numeric matrix, proteins (rows) × samples (cols). Values must be
#'   log2-transformed (as produced by \code{ELEUTHIA_load_massspec()}). Do NOT
#'   pass raw intensities.
#' @param verbose Logical. Print normalization summary. Default: TRUE.
#'
#' @return A numeric matrix of the same dimensions with per-sample median
#'   subtracted. NA values are unchanged.
#'
#' @details
#' Median centering is a standard post-processing step for DIA proteomics data.
#' After centering, each sample has a median of approximately zero (on the log2
#' scale), making cross-sample comparisons valid without distorting the relative
#' ordering of proteins within a sample.
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
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide)
#' }
POSEIDON_normalize_massspec <- function(wide, verbose = TRUE) {

  if (!is.matrix(wide))
    stop("wide must be a numeric matrix (proteins x samples).")

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

  return(norm_mat)
}
