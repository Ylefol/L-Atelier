# ==============================================================================
# HADES - Olink QC and filtering
#
# Every filter below records what it did in the object's $qc_log slot: the
# arguments it actually ran with (defaults included), the sample / protein
# counts before and after, and the IDs it removed or flagged with a reason.
# ELEUTHIA_export_olink_qc() reads that log to write the QC summary, so no QC
# record has to be built by hand. Read the log with HADES_get_olink_qc_log().
# ==============================================================================


# --- Internal helpers: the QC log ---------------------------------------------

# One log entry. `action` is "loaded", "removed" or "flagged".
#   params           named list of the arguments the step ran with (resolved)
#   details          named list of values the step derived (e.g. a cutoff)
#   removed_samples  data.frame(sample_id, reason, value) or NULL
#   removed_proteins data.frame(OlinkID, Assay, reason, value) or NULL
#   flagged_samples  data.frame(sample_id, reason, value) or NULL (kept in data)
#' @noRd
.olink_log_entry <- function(step, action, params,
                             samples_before, samples_after,
                             proteins_before, proteins_after,
                             removed_samples  = NULL,
                             removed_proteins = NULL,
                             flagged_samples  = NULL,
                             details          = list()) {
  list(
    step             = step,
    action           = action,
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    params           = params,
    details          = details,
    samples_before   = as.integer(samples_before),
    samples_after    = as.integer(samples_after),
    proteins_before  = as.integer(proteins_before),
    proteins_after   = as.integer(proteins_after),
    removed_samples  = if (is.null(removed_samples))  .olink_empty_samples()  else removed_samples,
    removed_proteins = if (is.null(removed_proteins)) .olink_empty_proteins() else removed_proteins,
    flagged_samples  = if (is.null(flagged_samples))  .olink_empty_samples()  else flagged_samples
  )
}

#' @noRd
.olink_empty_samples <- function() {
  data.frame(sample_id = character(0), reason = character(0),
             value = numeric(0), stringsAsFactors = FALSE)
}

#' @noRd
.olink_empty_proteins <- function() {
  data.frame(OlinkID = character(0), Assay = character(0),
             reason = character(0), value = numeric(0),
             stringsAsFactors = FALSE)
}

# Returns the existing log with `entry` appended (starts a log if there is none).
#' @noRd
.olink_log_append <- function(olink_data, entry) {
  log <- olink_data$qc_log
  if (is.null(log)) log <- list()
  log[[length(log) + 1L]] <- entry
  log
}

# Removed-protein table (with Assay names where available).
#' @noRd
.olink_removed_proteins <- function(ids, assay_meta, reason, value = NA_real_) {
  if (length(ids) == 0L) return(.olink_empty_proteins())
  assay <- if (!is.null(assay_meta) && "Assay" %in% colnames(assay_meta)) {
    as.character(assay_meta$Assay[match(ids, assay_meta$OlinkID)])
  } else {
    rep(NA_character_, length(ids))
  }
  data.frame(OlinkID = as.character(ids), Assay = assay,
             reason = reason, value = as.numeric(value),
             stringsAsFactors = FALSE)
}

# "a=1, b=TRUE, c=NULL" for printing a parameter list.
#' @noRd
.olink_fmt_params <- function(params) {
  if (length(params) == 0L) return("")
  paste(vapply(names(params), function(nm) {
    v <- params[[nm]]
    paste0(nm, "=", if (is.null(v)) "NULL" else paste(format(v), collapse = "/"))
  }, character(1L)), collapse = ", ")
}


#' Get the QC Log of an Olink Object
#'
#' @description Returns the QC record that the Olink loading and filtering
#' functions keep inside an \code{olink_data} object (\code{$qc_log}) as tidy
#' data frames. Every step records the arguments it actually ran with
#' (defaults included), the sample and protein counts before and after, and the
#' IDs it removed or flagged, each with a reason.
#'
#' @param olink_data An \code{olink_data} object produced by
#'   \code{ELEUTHIA_load_olink()} and (optionally) the \code{HADES_*_olink*}
#'   filters.
#' @param what Character. Which table to return:
#'   \describe{
#'     \item{\code{"steps"}}{(default) One row per step: \code{step_id},
#'       \code{step}, \code{action} (\code{"loaded"}, \code{"removed"} or
#'       \code{"flagged"}), \code{parameters}, sample / protein counts before
#'       and after, number removed and flagged, \code{timestamp}.}
#'     \item{\code{"removed_samples"}}{Samples removed: \code{step_id},
#'       \code{step}, \code{sample_id}, \code{reason}, \code{value}.}
#'     \item{\code{"removed_proteins"}}{Proteins removed: \code{step_id},
#'       \code{step}, \code{OlinkID}, \code{Assay}, \code{reason},
#'       \code{value}.}
#'     \item{\code{"flagged_samples"}}{Samples flagged but kept in the data
#'       (outlier detection with \code{filter = FALSE}): same columns as
#'       \code{"removed_samples"}.}
#'   }
#'
#' @return A data.frame (zero rows if nothing was removed or flagged).
#'
#' @details \code{value} holds the number behind the decision where there is
#' one (the NA fraction, the warn fraction, the outlier distance) and is
#' \code{NA} otherwise. A step that is run again on already-filtered data
#' appears as a new row; the log is the history of calls, in order.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
#' HADES_get_olink_qc_log(ol)                     # the filter chain
#' HADES_get_olink_qc_log(ol, "removed_proteins") # which proteins and why
#' }
HADES_get_olink_qc_log <- function(olink_data,
                                    what = c("steps", "removed_samples",
                                             "removed_proteins", "flagged_samples")) {

  if (!inherits(olink_data, "olink_data"))
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")

  what <- match.arg(what)
  log  <- olink_data$qc_log
  if (is.null(log) || length(log) == 0L)
    stop("This olink_data object has no $qc_log. Create it with ",
         "ELEUTHIA_load_olink() and filter it with the HADES_*_olink* functions.")

  if (what == "steps") {
    rows <- lapply(seq_along(log), function(i) {
      e <- log[[i]]
      data.frame(
        step_id           = i,
        step              = e$step,
        action            = e$action,
        parameters        = .olink_fmt_params(e$params),
        samples_before    = e$samples_before,
        samples_after     = e$samples_after,
        proteins_before   = e$proteins_before,
        proteins_after    = e$proteins_after,
        n_samples_removed = nrow(e$removed_samples),
        n_proteins_removed = nrow(e$removed_proteins),
        n_samples_flagged = nrow(e$flagged_samples),
        timestamp         = e$timestamp,
        stringsAsFactors  = FALSE
      )
    })
    return(do.call(rbind, rows))
  }

  empty <- if (what == "removed_proteins") .olink_empty_proteins() else .olink_empty_samples()
  rows <- lapply(seq_along(log), function(i) {
    tbl <- log[[i]][[what]]
    if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
    cbind(step_id = i, step = log[[i]]$step, tbl, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    return(cbind(step_id = integer(0), step = character(0), empty,
                 stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Filter Olink Samples
#'
#' @description Applies QC-based filtering to an \code{olink_data} object:
#' optionally removes non-biological control samples and/or samples that
#' failed Olink's SampleQC. Updates \code{$data}, \code{$wide}, and
#' \code{$sample_meta} to reflect the filtered set. The step, its arguments and
#' the IDs it removed are appended to \code{$qc_log}.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}.
#' @param keep_controls Logical. Retain non-SAMPLE SampleTypes (e.g.
#'   NEGATIVE_CONTROL, PLATE_CONTROL) in the output. Default: FALSE (removes
#'   them, keeping only SampleType == "SAMPLE"). A missing (\code{NA})
#'   SampleType is treated as not being a SAMPLE and is removed.
#' @param filter_qc Logical. Remove samples where SampleQC == "FAIL".
#'   Default: TRUE. A missing (\code{NA}) SampleQC is not treated as a failure.
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return A filtered \code{olink_data} object of the same class. \code{$data},
#'   \code{$wide}, and \code{$sample_meta} contain only the retained samples.
#'   \code{$assay_meta} is unchanged (warn_fraction is computed from the
#'   original load and reflects all samples). \code{$qc_log} gains one entry.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol_raw <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
#' ol <- HADES_filter_olink(ol_raw, keep_controls = FALSE, filter_qc = TRUE)
#' print(ol)
#' HADES_get_olink_qc_log(ol, "removed_samples")
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
  n_prot  <- nrow(olink_data$wide)

  removed <- list()
  skipped <- character(0)

  if (verbose) cat("[HADES] HADES_filter_olink\n")

  # ---------------------------------------------------------------------------
  # Filter by SampleType
  # ---------------------------------------------------------------------------

  if (!keep_controls) {
    if (!"SampleType" %in% colnames(smeta)) {
      warning("SampleType column not found in $sample_meta; ",
              "skipping control removal.")
      skipped <- c(skipped, "SampleType column not found")
    } else {
      is_sample <- smeta$SampleType %in% "SAMPLE"
      if (any(!is_sample)) {
        type <- as.character(smeta$SampleType[!is_sample])
        type[is.na(type)] <- "NA"
        removed[[length(removed) + 1L]] <- data.frame(
          sample_id = as.character(smeta[[sid_col]][!is_sample]),
          reason    = paste0("SampleType != SAMPLE (", type, ")"),
          value     = NA_real_, stringsAsFactors = FALSE)
        if (verbose)
          cat("    Removed", sum(!is_sample), "non-SAMPLE samples (controls)\n")
      }
      smeta <- smeta[is_sample, , drop = FALSE]
    }
  }

  # ---------------------------------------------------------------------------
  # Filter by SampleQC
  # ---------------------------------------------------------------------------

  if (filter_qc) {
    if (!"SampleQC" %in% colnames(smeta)) {
      warning("SampleQC column not found in $sample_meta; ",
              "skipping QC filter.")
      skipped <- c(skipped, "SampleQC column not found")
    } else {
      failed_mask <- smeta$SampleQC %in% "FAIL"
      failed      <- as.character(smeta[[sid_col]][failed_mask])
      if (length(failed) > 0L) {
        removed[[length(removed) + 1L]] <- data.frame(
          sample_id = failed, reason = "SampleQC == FAIL",
          value = NA_real_, stringsAsFactors = FALSE)
        if (verbose) {
          cat("    Removed", length(failed), "samples with SampleQC == 'FAIL':\n")
          cat("   ", paste(failed, collapse = ", "), "\n")
        }
      }
      smeta <- smeta[!failed_mask, , drop = FALSE]
    }
  }

  n_retained <- nrow(smeta)
  if (verbose) {
    cat("    Retained:", n_retained, "of", n_start, "samples\n")
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
  # Log and assemble output
  # ---------------------------------------------------------------------------

  entry <- .olink_log_entry(
    step            = "HADES_filter_olink",
    action          = "removed",
    params          = list(keep_controls = keep_controls, filter_qc = filter_qc),
    samples_before  = n_start, samples_after  = n_retained,
    proteins_before = n_prot,  proteins_after = nrow(wide_mat),
    removed_samples = if (length(removed) > 0L) do.call(rbind, removed) else NULL,
    details         = if (length(skipped) > 0L) list(skipped = skipped) else list()
  )

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = smeta,
    assay_meta  = olink_data$assay_meta,
    params      = olink_data$params,
    qc_log      = .olink_log_append(olink_data, entry)
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
#' computation only - the NPX values in \code{$data} and \code{$wide} are never
#' modified.
#'
#' The outliers are recomputed from scratch on every call. A call with
#' \code{filter = TRUE} therefore removes the samples found with \emph{its own}
#' arguments, not those of an earlier detect-only call; pass the same
#' \code{n_pcs}, \code{threshold} and \code{method} to both if you want them to
#' agree. The arguments each call actually used are recorded in \code{$qc_log}.
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
#' slots. \code{$qc_log} gains one entry: action \code{"flagged"} when
#' \code{filter = FALSE} (the flagged samples stay in the data), or
#' \code{"removed"} when \code{filter = TRUE}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
#'
#' # Detect and annotate only (inspect before deciding)
#' ol <- HADES_detect_outliers_olink(ol, n_pcs = 10, threshold = 10)
#' ol$sample_meta[ol$sample_meta$outlier, c("SampleID", "outlier_distance")]
#'
#' # Remove detected outliers (same arguments, so the same samples)
#' ol <- HADES_detect_outliers_olink(ol, n_pcs = 10, threshold = 10, filter = TRUE)
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

  n_pcs_requested <- as.integer(n_pcs)
  n_pcs <- min(n_pcs_requested, ncol(wide) - 1L, nrow(wide) - 1L)
  if (n_pcs < 1L) stop("Not enough samples/proteins for PCA-based outlier detection.")

  if (verbose) cat("[HADES] Detecting outliers:\n")

  # ---------------------------------------------------------------------------
  # Prepare matrix: mean-impute NAs per protein, then scale proteins
  # (internal only - $wide is never modified)
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
  scores  <- pca$x[, seq_len(n_pcs), drop = FALSE]   # samples x n_pcs

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
    cat("    PCs used:", n_pcs, " | variance explained:",
        paste0(var_pct[seq_len(n_pcs)], "%", collapse = ", "), "\n")
    cat("    Method  :", method, " | threshold:", threshold,
        " | cutoff:", round(cutoff, 3), "\n")
    cat("    Flagged :", n_flagged, "sample(s)\n")
    if (n_flagged > 0L) {
      flagged_ids <- names(distances)[outliers]
      cat("   ", paste(flagged_ids, collapse = ", "), "\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Annotate $sample_meta
  # ---------------------------------------------------------------------------

  smeta <- olink_data$sample_meta
  n_samples_before <- nrow(smeta)

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

  flagged_tbl <- if (n_flagged > 0L) {
    data.frame(
      sample_id = names(distances)[outliers],
      reason    = sprintf("PCA distance %.3f >= cutoff %.3f (%s, threshold %s, %d PCs)",
                          distances[outliers], cutoff, method,
                          format(threshold), n_pcs),
      value     = as.numeric(distances[outliers]),
      stringsAsFactors = FALSE)
  } else {
    NULL
  }

  if (filter && n_flagged > 0L) {
    keep_ids <- smeta[[sid_col]][!smeta$outlier]
    smeta    <- smeta[!smeta$outlier, ]

    data_df  <- olink_data$data[olink_data$data[[sid_col]] %in% keep_ids, ]
    rownames(data_df) <- NULL

    wide_mat <- olink_data$wide[, colnames(olink_data$wide) %in% keep_ids,
                                 drop = FALSE]

    if (verbose)
      cat("    Removed", n_flagged, "outlier sample(s); retained",
          nrow(smeta), "\n")
  } else {
    data_df  <- olink_data$data
    wide_mat <- olink_data$wide
    if (filter && n_flagged == 0L && verbose)
      cat("    No outliers to remove.\n")
  }

  rownames(smeta) <- smeta[[sid_col]]

  # ---------------------------------------------------------------------------
  # Log and assemble output
  # ---------------------------------------------------------------------------

  params <- list(n_pcs = n_pcs, threshold = threshold, method = method,
                 filter = filter)
  if (n_pcs != n_pcs_requested) params$n_pcs_requested <- n_pcs_requested

  entry <- .olink_log_entry(
    step            = "HADES_detect_outliers_olink",
    action          = if (filter) "removed" else "flagged",
    params          = params,
    samples_before  = n_samples_before, samples_after  = nrow(smeta),
    proteins_before = nrow(olink_data$wide), proteins_after = nrow(wide_mat),
    removed_samples = if (filter) flagged_tbl else NULL,
    flagged_samples = if (!filter) flagged_tbl else NULL,
    details         = list(cutoff = cutoff, centre = centre, spread = spread)
  )

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = smeta,
    assay_meta  = olink_data$assay_meta,
    params      = olink_data$params,
    qc_log      = .olink_log_append(olink_data, entry)
  )
  class(result) <- c("olink_data", "list")

  return(result)
}


#' Filter Olink Proteins by Missingness
#'
#' @description Removes proteins (rows) from an \code{olink_data} object whose
#' NPX missingness fraction across samples exceeds a threshold. High-NA
#' proteins are dropped because missingness cannot be reliably attributed to
#' a true zero signal versus technical absence. The step and the removed protein
#' IDs (with their NA fractions) are appended to \code{$qc_log}.
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
    cat("[HADES] Olink protein missingness filter\n")
    cat("    Threshold: NA fraction >", max_na_fraction, "\n")
    cat("    Removed :", n_removed,  "proteins\n")
    cat("    Retained:", n_retained, "of", n_before, "proteins\n")
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

  entry <- .olink_log_entry(
    step            = "HADES_filter_olink_proteins",
    action          = "removed",
    params          = list(max_na_fraction = max_na_fraction),
    samples_before  = n_samples, samples_after  = n_samples,
    proteins_before = n_before,  proteins_after = n_retained,
    removed_proteins = .olink_removed_proteins(
      rownames(wide)[!keep_mask], olink_data$assay_meta,
      reason = paste0("NA fraction > ", max_na_fraction),
      value  = na_fracs[!keep_mask])
  )

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = olink_data$sample_meta,
    assay_meta  = assay_meta,
    params      = olink_data$params,
    qc_log      = .olink_log_append(olink_data, entry)
  )
  class(result) <- c("olink_data", "list")

  return(result)
}


#' Filter Olink Proteins by AssayQC Warn Fraction
#'
#' @description Removes proteins (rows) from an \code{olink_data} object whose
#' AssayQC warn fraction across samples exceeds a threshold. A high warn
#' fraction indicates systematic assay instability for that protein. The step
#' and the removed protein IDs (with their warn fractions) are appended to
#' \code{$qc_log}.
#'
#' \code{warn_fraction} is computed at load time by \code{ELEUTHIA_load_olink()}
#' and stored in \code{$assay_meta}. It reflects the proportion of samples
#' where that protein received \code{AssayQC == "WARN"}.
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}.
#' @param max_warn_fraction Numeric (0-1). Maximum allowable proportion of samples
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
    cat("[HADES] Olink AssayQC warn filter\n")
    cat("    Threshold: warn fraction >", max_warn_fraction, "\n")
    cat("    Removed  :", n_removed,  "proteins\n")
    cat("    Retained :", n_retained, "of", n_before, "proteins\n")
  }

  keep_ids   <- ameta$OlinkID[keep_mask]
  ameta_filt <- ameta[keep_mask, , drop = FALSE]
  rownames(ameta_filt) <- NULL

  wide_mat  <- olink_data$wide[rownames(olink_data$wide) %in% keep_ids, , drop = FALSE]

  data_filt <- olink_data$data[olink_data$data$OlinkID %in% keep_ids, ]
  rownames(data_filt) <- NULL

  entry <- .olink_log_entry(
    step            = "HADES_filter_olink_assay_warn",
    action          = "removed",
    params          = list(max_warn_fraction = max_warn_fraction),
    samples_before  = ncol(olink_data$wide), samples_after = ncol(wide_mat),
    proteins_before = nrow(olink_data$wide), proteins_after = nrow(wide_mat),
    removed_proteins = .olink_removed_proteins(
      ameta$OlinkID[!keep_mask], ameta,
      reason = paste0("AssayQC WARN fraction > ", max_warn_fraction),
      value  = ameta$warn_fraction[!keep_mask])
  )

  result <- list(
    data        = data_filt,
    wide        = wide_mat,
    sample_meta = olink_data$sample_meta,
    assay_meta  = ameta_filt,
    params      = olink_data$params,
    qc_log      = .olink_log_append(olink_data, entry)
  )
  class(result) <- c("olink_data", "list")

  return(result)
}
