#' Filter Olink Samples
#'
#' @description Applies QC-based filtering to an \code{olink_data} object:
#' optionally removes non-biological control samples and/or samples that
#' failed Olink's SampleQC. Updates \code{$data}, \code{$wide}, and
#' \code{$sample_meta} to reflect the filtered set.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}.
#' @param keep_controls Logical. Retain non-SAMPLE SampleTypes (e.g.
#'   NEGATIVE_CONTROL, PLATE_CONTROL) in the output. Default: FALSE (removes
#'   them, keeping only SampleType == "SAMPLE").
#' @param filter_qc Logical. Remove samples where SampleQC == "FAIL".
#'   Default: TRUE.
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{olink_data} object of the same class. \code{$data},
#'   \code{$wide}, and \code{$sample_meta} contain only the retained samples.
#'   \code{$assay_meta} is unchanged (warn_fraction is computed from the
#'   original load and reflects all samples).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol_raw <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
#' ol <- HADES_filter_olink(ol_raw, keep_controls = FALSE, filter_qc = TRUE)
#' print(ol)
#' }
HADES_filter_olink <- function(olink_data,
                                keep_controls = FALSE,
                                filter_qc     = TRUE,
                                verbose       = TRUE) {

  if (!inherits(olink_data, "olink_data")) {
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  }

  sid_col <- olink_data$params$sample_col
  smeta   <- olink_data$sample_meta
  n_start <- nrow(smeta)

  if (verbose) cat("HADES_filter_olink:\n")

  # ---------------------------------------------------------------------------
  # Filter by SampleType
  # ---------------------------------------------------------------------------

  if (!keep_controls) {
    if (!"SampleType" %in% colnames(smeta)) {
      warning("SampleType column not found in $sample_meta; ",
              "skipping control removal.")
    } else {
      n_before  <- nrow(smeta)
      smeta     <- smeta[smeta$SampleType == "SAMPLE", ]
      n_removed <- n_before - nrow(smeta)
      if (verbose && n_removed > 0) {
        cat("[HADES] Removed", n_removed, "non-SAMPLE samples (controls)\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Filter by SampleQC
  # ---------------------------------------------------------------------------

  if (filter_qc) {
    if (!"SampleQC" %in% colnames(smeta)) {
      warning("SampleQC column not found in $sample_meta; ",
              "skipping QC filter.")
    } else {
      failed    <- smeta[[sid_col]][smeta$SampleQC == "FAIL"]
      smeta     <- smeta[smeta$SampleQC != "FAIL", ]
      if (verbose && length(failed) > 0) {
        cat("[HADES] Removed", length(failed), "samples with SampleQC == 'FAIL':\n")
        cat("   ", paste(failed, collapse = ", "), "\n")
      }
    }
  }

  n_retained <- nrow(smeta)
  if (verbose) {
    cat("[HADES] Retained:", n_retained, "of", n_start, "samples\n")
  }

  # ---------------------------------------------------------------------------
  # Subset $data and $wide to retained sample IDs
  # ---------------------------------------------------------------------------

  keep_ids <- as.character(smeta[[sid_col]])

  data_df  <- olink_data$data[olink_data$data[[sid_col]] %in% keep_ids, ]
  rownames(data_df) <- NULL

  wide_mat <- olink_data$wide[, colnames(olink_data$wide) %in% keep_ids,
                               drop = FALSE]

  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = smeta,
    assay_meta  = olink_data$assay_meta,
    params      = olink_data$params
  )
  class(result) <- c("olink_data", "list")

  return(result)
}


#' Detect Sample Outliers in Olink Data
#'
#' @description Flags multivariate outliers in an \code{olink_data} object using
#' PCA-based Euclidean distance from the sample centroid. Outliers are annotated
#' in \code{$sample_meta} via two new columns (\code{outlier},
#' \code{outlier_distance}) and optionally removed.
#'
#' Mean imputation and per-protein scaling are applied internally for the outlier
#' computation only — the NPX values in \code{$data} and \code{$wide} are never
#' modified.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}
#'   or \code{HADES_filter_olink()}.
#' @param n_pcs Integer. Number of principal components to use for the distance
#'   computation. Default: 5.
#' @param threshold Numeric. Multiplier applied to the spread measure (SD or MAD)
#'   above the mean/median distance. Samples whose distance exceeds
#'   \code{centre + threshold * spread} are flagged. Default: 3.
#' @param method Character. Spread measure: \code{"sd"} uses
#'   \code{mean + threshold * sd}; \code{"mad"} uses
#'   \code{median + threshold * mad}. MAD is more robust when multiple outliers
#'   are present. Default: \code{"sd"}.
#' @param filter Logical. If \code{TRUE}, flagged samples are removed from
#'   \code{$data}, \code{$wide}, and \code{$sample_meta}. If \code{FALSE}
#'   (default), samples are annotated but retained.
#' @param verbose Logical. Print detection summary. Default: TRUE.
#'
#' @return An \code{olink_data} object. \code{$sample_meta} gains two columns:
#' \describe{
#'   \item{outlier_distance}{Euclidean distance from the centroid in PC space.}
#'   \item{outlier}{Logical. \code{TRUE} for samples exceeding the threshold.}
#' }
#' When \code{filter = TRUE}, flagged samples are also removed from all data
#' slots.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
#'
#' # Detect and annotate only (inspect before deciding)
#' ol <- HADES_detect_outliers_olink(ol)
#' ol$sample_meta[ol$sample_meta$outlier, c("SampleID", "outlier_distance")]
#'
#' # Remove detected outliers
#' ol <- HADES_detect_outliers_olink(ol, filter = TRUE)
#' }
HADES_detect_outliers_olink <- function(olink_data,
                                         n_pcs     = 5L,
                                         threshold = 3,
                                         method    = c("sd", "mad"),
                                         filter    = FALSE,
                                         verbose   = TRUE) {

  if (!inherits(olink_data, "olink_data")) {
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  }

  method  <- match.arg(method)
  sid_col <- olink_data$params$sample_col
  wide    <- olink_data$wide

  n_pcs <- min(as.integer(n_pcs), ncol(wide) - 1L, nrow(wide) - 1L)
  if (n_pcs < 1L) stop("Not enough samples/proteins for PCA-based outlier detection.")

  if (verbose) cat("[HADES] Detecting outliers:\n")

  # ---------------------------------------------------------------------------
  # Prepare matrix: mean-impute NAs per protein, then scale proteins
  # (internal only — $wide is never modified)
  # ---------------------------------------------------------------------------

  mat <- wide
  row_means <- rowMeans(mat, na.rm = TRUE)
  na_idx    <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    mat[na_idx] <- row_means[na_idx[, 1L]]
  }

  # Remove zero-variance proteins (would cause scale() to produce NaN)
  rv       <- apply(mat, 1L, var)
  mat      <- mat[rv > 0, , drop = FALSE]

  # Scale each protein to zero mean / unit variance, then transpose for PCA
  mat_scaled <- t(scale(t(mat)))

  # ---------------------------------------------------------------------------
  # PCA and distance computation
  # ---------------------------------------------------------------------------

  pca     <- stats::prcomp(t(mat_scaled), center = FALSE, scale. = FALSE)
  scores  <- pca$x[, seq_len(n_pcs), drop = FALSE]   # samples × n_pcs

  centroid  <- colMeans(scores)
  distances <- sqrt(rowSums(sweep(scores, 2L, centroid)^2))

  # ---------------------------------------------------------------------------
  # Flag outliers
  # ---------------------------------------------------------------------------

  if (method == "sd") {
    centre <- mean(distances)
    spread <- stats::sd(distances)
  } else {
    centre <- stats::median(distances)
    spread <- stats::mad(distances, constant = 1)
  }

  cutoff   <- centre + threshold * spread
  outliers <- distances >= cutoff

  n_flagged <- sum(outliers)

  if (verbose) {
    var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    cat("[HADES] PCs used:", n_pcs, " | variance explained:",
        paste0(var_pct[seq_len(n_pcs)], "%", collapse = ", "), "\n")
    cat("  Method  :", method, " | threshold:", threshold,
        " | cutoff:", round(cutoff, 3), "\n")
    cat("  Flagged :", n_flagged, "sample(s)\n")
    if (n_flagged > 0L) {
      flagged_ids <- names(distances)[outliers]
      cat("   ", paste(flagged_ids, collapse = ", "), "\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Annotate $sample_meta
  # ---------------------------------------------------------------------------

  smeta <- olink_data$sample_meta

  # Align distances to sample_meta row order
  dist_df <- data.frame(
    tmp_sid          = names(distances),
    outlier_distance = as.numeric(distances),
    outlier          = as.logical(outliers),
    stringsAsFactors = FALSE
  )
  colnames(dist_df)[1L] <- sid_col

  # Drop any pre-existing columns to avoid duplicates on re-run
  smeta <- smeta[, setdiff(colnames(smeta), c("outlier_distance", "outlier")),
                 drop = FALSE]

  smeta <- merge(smeta, dist_df, by = sid_col, all.x = TRUE, sort = FALSE)

  # ---------------------------------------------------------------------------
  # Optionally filter
  # ---------------------------------------------------------------------------

  if (filter && n_flagged > 0L) {
    keep_ids <- smeta[[sid_col]][!smeta$outlier]
    smeta    <- smeta[!smeta$outlier, ]

    data_df  <- olink_data$data[olink_data$data[[sid_col]] %in% keep_ids, ]
    rownames(data_df) <- NULL

    wide_mat <- olink_data$wide[, colnames(olink_data$wide) %in% keep_ids,
                                 drop = FALSE]

    if (verbose)
      cat("[HADES] Removed", n_flagged, "outlier sample(s); retained",
          nrow(smeta), "\n")
  } else {
    data_df  <- olink_data$data
    wide_mat <- olink_data$wide
    if (filter && n_flagged == 0L && verbose)
      cat("[HADES] No outliers to remove.\n")
  }

  rownames(smeta) <- NULL

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = smeta,
    assay_meta  = olink_data$assay_meta,
    params      = olink_data$params
  )
  class(result) <- c("olink_data", "list")

  return(result)
}


#' Filter Olink Proteins by Missingness
#'
#' @description Removes proteins (rows) from an \code{olink_data} object whose
#' NPX missingness fraction across samples exceeds a threshold. High-NA
#' proteins are dropped because missingness cannot be reliably attributed to
#' a true zero signal versus technical absence.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}
#'   or \code{HADES_filter_olink()}.
#' @param max_na_fraction Numeric (0-1). Maximum allowable proportion of NA
#'   NPX values across samples for a protein to be retained. Proteins with
#'   NA fraction above this threshold are removed. Default: 0.2 (20%).
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{olink_data} object. \code{$wide}, \code{$data},
#'   and \code{$assay_meta} contain only the retained proteins.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
#' }
HADES_filter_olink_proteins <- function(olink_data,
                                         max_na_fraction = 0.2,
                                         verbose         = TRUE) {

  if (!inherits(olink_data, "olink_data")) {
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  }

  wide      <- olink_data$wide
  n_samples <- ncol(wide)
  n_before  <- nrow(wide)

  na_fracs  <- rowSums(is.na(wide)) / n_samples
  keep_mask <- na_fracs <= max_na_fraction

  n_removed  <- sum(!keep_mask)
  n_retained <- sum(keep_mask)

  if (verbose) {
    cat("[HADES] Summary:\n")
    cat("  Threshold: NA fraction >", max_na_fraction, "\n")
    cat("  Removed :", n_removed,  "proteins\n")
    cat("  Retained:", n_retained, "of", n_before, "proteins\n")
  }

  keep_ids    <- rownames(wide)[keep_mask]
  sid_col     <- olink_data$params$sample_col

  wide_mat    <- wide[keep_mask, , drop = FALSE]

  data_df     <- olink_data$data[olink_data$data$OlinkID %in% keep_ids, ]
  rownames(data_df) <- NULL

  assay_meta  <- olink_data$assay_meta[
    olink_data$assay_meta$OlinkID %in% keep_ids, , drop = FALSE
  ]
  rownames(assay_meta) <- NULL

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = olink_data$sample_meta,
    assay_meta  = assay_meta,
    params      = olink_data$params
  )
  class(result) <- c("olink_data", "list")

  return(result)
}


#' Filter Olink Proteins by AssayQC Warn Fraction
#'
#' @description Removes proteins (rows) from an \code{olink_data} object whose
#' AssayQC warn fraction across samples exceeds a threshold. A high warn
#' fraction indicates systematic assay instability for that protein.
#'
#' \code{warn_fraction} is computed at load time by \code{ELEUTHIA_load_olink()}
#' and stored in \code{$assay_meta}. It reflects the proportion of samples
#' where that protein received \code{AssayQC == "WARN"}.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}.
#' @param max_warn_fraction Numeric (0–1). Maximum allowable proportion of samples
#'   with \code{AssayQC == "WARN"} for a protein to be retained. Proteins
#'   exceeding this threshold are removed. Default: \code{0.2} (20%). No
#'   established universal standard exists; verify this threshold against your
#'   panel and study design.
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{olink_data} object. \code{$wide}, \code{$data}, and
#'   \code{$assay_meta} contain only retained proteins. \code{$sample_meta} is
#'   unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
#' ol <- HADES_filter_olink_lod(ol, max_lod_fraction = 0.5)
#' ol <- HADES_filter_olink_assay_warn(ol, max_warn_fraction = 0.2)
#' }
HADES_filter_olink_assay_warn <- function(olink_data,
                                           max_warn_fraction = 0.2,
                                           verbose           = TRUE) {

  if (!inherits(olink_data, "olink_data"))
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")

  ameta    <- olink_data$assay_meta
  n_before <- nrow(ameta)

  if (!"warn_fraction" %in% colnames(ameta))
    stop("$assay_meta does not contain 'warn_fraction'. ",
         "Ensure data was loaded with ELEUTHIA_load_olink().")

  keep_mask  <- is.na(ameta$warn_fraction) | ameta$warn_fraction <= max_warn_fraction
  n_removed  <- sum(!keep_mask)
  n_retained <- sum(keep_mask)

  if (verbose) {
    cat("[HADES] AssayQC warn filter:\n")
    cat("  Threshold : warn fraction >", max_warn_fraction, "\n")
    cat("  Removed  :", n_removed,  "proteins\n")
    cat("  Retained :", n_retained, "of", n_before, "proteins\n")
  }

  keep_ids   <- ameta$OlinkID[keep_mask]
  ameta_filt <- ameta[keep_mask, , drop = FALSE]
  rownames(ameta_filt) <- NULL

  wide_mat  <- olink_data$wide[rownames(olink_data$wide) %in% keep_ids, , drop = FALSE]

  data_filt <- olink_data$data[olink_data$data$OlinkID %in% keep_ids, ]
  rownames(data_filt) <- NULL

  result <- list(
    data        = data_filt,
    wide        = wide_mat,
    sample_meta = olink_data$sample_meta,
    assay_meta  = ameta_filt,
    params      = olink_data$params
  )
  class(result) <- c("olink_data", "list")

  return(result)
}
