#' Filter Mass Spectrometry Samples
#'
#' @description Filters a \code{massspec_data} object by sample type and/or
#' per-sample missingness. Updates \code{$wide} and \code{$sample_meta} to
#' reflect the retained samples.
#'
#' @param massspec_data A \code{massspec_data} object from
#'   \code{ELEUTHIA_load_massspec()}.
#' @param keep_types Character vector. SampleType values to retain. Default:
#'   \code{"SAMPLE"} (numeric-ID biological samples with metadata). Other
#'   valid values: \code{"SEP_SAMPLE"}, \code{"CONTROL"}, \code{"POOL"},
#'   \code{"OTHER"}.
#' @param max_sample_na_fraction Numeric (0–1). Samples whose log2 NA fraction
#'   across all proteins exceeds this threshold are removed. Default:
#'   \code{0.9} (remove only extreme outliers; DIA data normally has
#'   ~70% overall missingness).
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{massspec_data} object. \code{$protein_meta} is
#'   unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms_raw <- ELEUTHIA_load_massspec("report.pg_matrix.xlsx",
#'                                   metadata_file = "meta.csv")
#' ms <- HADES_filter_massspec(ms_raw, keep_types = "SAMPLE")
#' }
HADES_filter_massspec <- function(massspec_data,
                                   keep_types            = "SAMPLE",
                                   max_sample_na_fraction = 0.9,
                                   verbose               = TRUE) {

  if (!inherits(massspec_data, "massspec_data"))
    stop("massspec_data must be a massspec_data object from ELEUTHIA_load_massspec().")

  sid_col <- massspec_data$params$sample_col
  smeta   <- massspec_data$sample_meta
  n_start <- nrow(smeta)

  if (verbose) cat("[HADES] Filtering mass spec samples:\n")

  # ---------------------------------------------------------------------------
  # Filter by SampleType
  # ---------------------------------------------------------------------------

  if (!is.null(keep_types) && "SampleType" %in% colnames(smeta)) {
    n_before  <- nrow(smeta)
    smeta     <- smeta[smeta$SampleType %in% keep_types, ]
    n_removed <- n_before - nrow(smeta)
    if (verbose && n_removed > 0)
      cat("[HADES]   Removed", n_removed,
          "samples not in keep_types (", paste(keep_types, collapse = ", "), ")\n")
  }

  # ---------------------------------------------------------------------------
  # Filter by per-sample NA fraction
  # ---------------------------------------------------------------------------

  sample_ids <- as.character(smeta[[sid_col]])
  wide_sub   <- massspec_data$wide[, sample_ids, drop = FALSE]

  na_fracs   <- colMeans(is.na(wide_sub))
  pass_na    <- na_fracs <= max_sample_na_fraction

  if (any(!pass_na)) {
    removed_ids <- sample_ids[!pass_na]
    if (verbose)
      cat("[HADES]   Removed", sum(!pass_na),
          "samples with NA fraction >", max_sample_na_fraction, ":\n",
          "   ", paste(removed_ids, collapse = ", "), "\n")
    sample_ids <- sample_ids[pass_na]
    smeta      <- smeta[smeta[[sid_col]] %in% sample_ids, ]
  }

  n_retained <- nrow(smeta)
  if (verbose)
    cat("[HADES]   Retained:", n_retained, "of", n_start, "samples\n")

  wide_mat <- massspec_data$wide[, sample_ids, drop = FALSE]

  result <- list(
    wide         = wide_mat,
    sample_meta  = smeta,
    protein_meta = massspec_data$protein_meta,
    params       = massspec_data$params
  )
  class(result) <- c("massspec_data", "list")
  return(result)
}


#' Filter Mass Spectrometry Proteins by Missingness
#'
#' @description Removes proteins (rows) from a \code{massspec_data} object
#' whose missingness fraction across samples exceeds a threshold. DIA data
#' typically has high per-protein missingness, so a higher threshold than
#' used for Olink data is appropriate.
#'
#' @param massspec_data A \code{massspec_data} object, usually after
#'   \code{HADES_filter_massspec()}.
#' @param max_na_fraction Numeric (0–1). Maximum allowable proportion of NA
#'   values across samples. Proteins above this threshold are removed. Default:
#'   \code{0.5} (retain proteins detected in at least 50% of samples).
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{massspec_data} object. \code{$sample_meta} is
#'   unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_filter_massspec(ms_raw)
#' ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
#' }
HADES_filter_massspec_proteins <- function(massspec_data,
                                            max_na_fraction = 0.5,
                                            verbose         = TRUE) {

  if (!inherits(massspec_data, "massspec_data"))
    stop("massspec_data must be a massspec_data object from ELEUTHIA_load_massspec().")

  wide     <- massspec_data$wide
  n_before <- nrow(wide)

  na_fracs  <- rowMeans(is.na(wide))
  keep_mask <- na_fracs <= max_na_fraction

  n_removed  <- sum(!keep_mask)
  n_retained <- sum(keep_mask)

  if (verbose) {
    cat("[HADES] Filtering mass spec proteins:\n")
    cat("    Threshold  : NA fraction <=", max_na_fraction, "\n")
    cat("    Removed    :", n_removed, "proteins\n")
    cat("    Retained   :", n_retained, "of", n_before, "proteins\n")
  }

  wide_mat <- wide[keep_mask, , drop = FALSE]
  pmeta    <- massspec_data$protein_meta[keep_mask, , drop = FALSE]

  result <- list(
    wide         = wide_mat,
    sample_meta  = massspec_data$sample_meta,
    protein_meta = pmeta,
    params       = massspec_data$params
  )
  class(result) <- c("massspec_data", "list")
  return(result)
}


#' Detect Sample Outliers in Mass Spectrometry Data
#'
#' @description Flags multivariate outliers in a \code{massspec_data} object
#' using PCA-based Euclidean distance from the sample centroid. Mean imputation
#' and per-protein scaling are applied internally for the outlier computation
#' only — log2 values in \code{$wide} are never modified.
#'
#' @param massspec_data A \code{massspec_data} object after filtering and
#'   normalization.
#' @param n_pcs Integer. Number of principal components for the distance
#'   computation. Default: 5.
#' @param threshold Numeric. Multiplier applied to the spread measure above the
#'   centre. Default: 3.
#' @param method Character. \code{"sd"} (mean + threshold × sd) or
#'   \code{"mad"} (median + threshold × mad). Default: \code{"sd"}.
#' @param filter Logical. If TRUE, flagged samples are removed from
#'   \code{$wide} and \code{$sample_meta}. Default: FALSE (annotate only).
#' @param verbose Logical. Print detection summary. Default: TRUE.
#'
#' @return A \code{massspec_data} object. \code{$sample_meta} gains
#'   \code{outlier_distance} and \code{outlier} columns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_detect_outliers_massspec(ms, filter = FALSE)
#' ms$sample_meta[ms$sample_meta$outlier, c("SampleID", "outlier_distance")]
#' }
HADES_detect_outliers_massspec <- function(massspec_data,
                                            n_pcs     = 5L,
                                            threshold = 3,
                                            method    = c("sd", "mad"),
                                            filter    = FALSE,
                                            verbose   = TRUE) {

  if (!inherits(massspec_data, "massspec_data"))
    stop("massspec_data must be a massspec_data object from ELEUTHIA_load_massspec().")

  method  <- match.arg(method)
  sid_col <- massspec_data$params$sample_col
  wide    <- massspec_data$wide

  n_pcs <- min(as.integer(n_pcs), ncol(wide) - 1L, nrow(wide) - 1L)
  if (n_pcs < 1L)
    stop("Not enough samples/proteins for PCA-based outlier detection.")

  if (verbose) cat("[HADES] Detecting outliers (mass spec):\n")

  # Mean-impute NAs per protein and scale (internal copy only)
  mat        <- wide
  row_means  <- rowMeans(mat, na.rm = TRUE)
  na_idx     <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) mat[na_idx] <- row_means[na_idx[, 1L]]

  rv  <- apply(mat, 1L, var)
  mat <- mat[rv > 0, , drop = FALSE]

  mat_scaled <- t(scale(t(mat)))

  pca      <- stats::prcomp(t(mat_scaled), center = FALSE, scale. = FALSE)
  scores   <- pca$x[, seq_len(n_pcs), drop = FALSE]
  centroid <- colMeans(scores)
  distances <- sqrt(rowSums(sweep(scores, 2L, centroid)^2))

  if (method == "sd") {
    centre <- mean(distances)
    spread <- stats::sd(distances)
  } else {
    centre <- stats::median(distances)
    spread <- stats::mad(distances, constant = 1)
  }

  cutoff    <- centre + threshold * spread
  outliers  <- distances >= cutoff
  n_flagged <- sum(outliers)

  if (verbose) {
    var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    cat("[HADES]   PCs used:", n_pcs, " | variance explained:",
        paste0(var_pct[seq_len(n_pcs)], "%", collapse = ", "), "\n")
    cat("    Method  :", method, " | threshold:", threshold,
        " | cutoff:", round(cutoff, 3), "\n")
    cat("    Flagged :", n_flagged, "sample(s)\n")
    if (n_flagged > 0L)
      cat("   ", paste(names(distances)[outliers], collapse = ", "), "\n")
  }

  smeta   <- massspec_data$sample_meta
  dist_df <- data.frame(
    tmp_sid          = names(distances),
    outlier_distance = as.numeric(distances),
    outlier          = as.logical(outliers),
    stringsAsFactors = FALSE
  )
  colnames(dist_df)[1L] <- sid_col

  smeta <- smeta[, setdiff(colnames(smeta), c("outlier_distance", "outlier")),
                 drop = FALSE]
  smeta <- merge(smeta, dist_df, by = sid_col, all.x = TRUE, sort = FALSE)

  if (filter && n_flagged > 0L) {
    keep_ids <- smeta[[sid_col]][!smeta$outlier]
    smeta    <- smeta[!smeta$outlier, ]
    wide_mat <- massspec_data$wide[, colnames(massspec_data$wide) %in% keep_ids,
                                    drop = FALSE]
    if (verbose)
      cat("[HADES]   Removed", n_flagged, "outlier(s); retained", nrow(smeta), "\n")
  } else {
    wide_mat <- massspec_data$wide
    if (filter && n_flagged == 0L && verbose)
      cat("[HADES]   No outliers to remove.\n")
  }

  rownames(smeta) <- NULL

  result <- list(
    wide         = wide_mat,
    sample_meta  = smeta,
    protein_meta = massspec_data$protein_meta,
    params       = massspec_data$params
  )
  class(result) <- c("massspec_data", "list")
  return(result)
}
