#' Artemis - Utility Functions
#'
#' @description Miscellaneous utility functions for genomic analysis
#' that don't fit into more specific categories.


#' Calculate Distance to Nearest Feature
#'
#' @description For each query region, finds the nearest subject region and
#' reports the distance. Wraps GenomicRanges::distanceToNearest() with a
#' convenient data.frame interface.
#'
#' @param query_regions Data.frame with at minimum: chr, start, end.
#'   These are the regions you want to annotate (e.g., ChIP anchors).
#' @param subject_regions Data.frame with at minimum: chr, start, end.
#'   These are the reference features to find distances to (e.g., ATAC peaks).
#' @param query_id_col Character or NULL. Column name to use as query IDs.
#'   If NULL, auto-detects from common names (peak_id, name, region_id) or
#'   generates sequential IDs.
#' @param subject_id_col Character or NULL. Column name to use as subject IDs.
#'   Same auto-detection logic as query_id_col.
#' @param verbose Logical. Print summary statistics (default = TRUE).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{query_id}{ID of the query region}
#'   \item{query_chr, query_start, query_end}{Query region coordinates}
#'   \item{nearest_subject_id}{ID of the nearest subject region}
#'   \item{subject_chr, subject_start, subject_end}{Nearest subject coordinates}
#'   \item{distance}{Distance in bp (0 if overlapping)}
#'   \item{overlaps}{Logical. TRUE if distance == 0}
#' }
#'
#' @details
#' Distance is calculated as the minimum number of base pairs separating
#' the two regions. Overlapping regions have distance = 0.
#'
#' This function requires the GenomicRanges package from Bioconductor.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Find distance from ChIP anchors to nearest ATAC peak
#' atac_peaks <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' consensus_peaks <- ELEUTHIA_create_consensus_peaks(atac_peaks)
#'
#' distances <- ARTEMIS_distance_to_nearest(
#'   query_regions = chip_anchors,
#'   subject_regions = consensus_peaks
#' )
#'
#' # Check overlap rate
#' mean(distances$overlaps)  # Fraction of anchors at ATAC peaks
#'
#' # Examine distance distribution for non-overlapping
#' hist(distances$distance[!distances$overlaps])
#'
#' }
ARTEMIS_distance_to_nearest <- function(query_regions,
                                         subject_regions,
                                         query_id_col = NULL,
                                         subject_id_col = NULL,
                                         verbose = TRUE) {

  # Check for GenomicRanges

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required. Install from Bioconductor:\n",
         "  BiocManager::install('GenomicRanges')")
  }

  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required. Install from Bioconductor:\n",
         "  BiocManager::install('IRanges')")
  }

  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(query_regions)) {
    stop("query_regions must be a data.frame")
  }

  if (!is.data.frame(subject_regions)) {
    stop("subject_regions must be a data.frame")
  }

  required_cols <- c("chr", "start", "end")

  missing_query <- setdiff(required_cols, colnames(query_regions))
  if (length(missing_query) > 0) {
    stop("query_regions missing required columns: ",
         paste(missing_query, collapse = ", "))
  }

  missing_subject <- setdiff(required_cols, colnames(subject_regions))
  if (length(missing_subject) > 0) {
    stop("subject_regions missing required columns: ",
         paste(missing_subject, collapse = ", "))
  }

  if (nrow(query_regions) == 0) {
    stop("query_regions is empty")
  }

  if (nrow(subject_regions) == 0) {
    stop("subject_regions is empty")
  }

  # ---------------------------------------------------------------------------
  # Auto-detect or validate ID columns
  # ---------------------------------------------------------------------------
  common_id_cols <- c("peak_id", "region_id", "name", "id", "ID")

  # Query IDs
  if (is.null(query_id_col)) {
    found <- intersect(common_id_cols, colnames(query_regions))
    if (length(found) > 0) {
      query_id_col <- found[1]
      if (verbose) cat("Using '", query_id_col, "' as query ID column\n", sep = "")
    } else {
      query_regions$`.query_id` <- paste0("query_", seq_len(nrow(query_regions)))
      query_id_col <- ".query_id"
      if (verbose) cat("Generated sequential query IDs\n")
    }
  } else if (!query_id_col %in% colnames(query_regions)) {
    stop("query_id_col '", query_id_col, "' not found in query_regions")
  }

  # Subject IDs
  if (is.null(subject_id_col)) {
    found <- intersect(common_id_cols, colnames(subject_regions))
    if (length(found) > 0) {
      subject_id_col <- found[1]
      if (verbose) cat("Using '", subject_id_col, "' as subject ID column\n", sep = "")
    } else {
      subject_regions$`.subject_id` <- paste0("subject_", seq_len(nrow(subject_regions)))
      subject_id_col <- ".subject_id"
      if (verbose) cat("Generated sequential subject IDs\n")
    }
  } else if (!subject_id_col %in% colnames(subject_regions)) {
    stop("subject_id_col '", subject_id_col, "' not found in subject_regions")
  }

  # ---------------------------------------------------------------------------
  # Convert to GRanges
  # ---------------------------------------------------------------------------
  query_gr <- GenomicRanges::GRanges(
    seqnames = query_regions$chr,
    ranges = IRanges::IRanges(
      start = query_regions$start + 1,  # Convert from 0-based BED to 1-based GRanges
      end = query_regions$end
    )
  )

  subject_gr <- GenomicRanges::GRanges(
    seqnames = subject_regions$chr,
    ranges = IRanges::IRanges(
      start = subject_regions$start + 1,
      end = subject_regions$end
    )
  )

  if (verbose) {
    cat("\nFinding nearest features:\n")
    cat("  Query regions:", length(query_gr), "\n")
    cat("  Subject regions:", length(subject_gr), "\n")
  }

  # ---------------------------------------------------------------------------
  # Find nearest
  # ---------------------------------------------------------------------------
  nearest_hits <- GenomicRanges::distanceToNearest(query_gr, subject_gr)

  # Extract results
  query_idx <- S4Vectors::queryHits(nearest_hits)
  subject_idx <- S4Vectors::subjectHits(nearest_hits)
  distances <- S4Vectors::mcols(nearest_hits)$distance

  # Build result data.frame
  result <- data.frame(
    query_id = query_regions[[query_id_col]][query_idx],
    query_chr = query_regions$chr[query_idx],
    query_start = query_regions$start[query_idx],
    query_end = query_regions$end[query_idx],
    nearest_subject_id = subject_regions[[subject_id_col]][subject_idx],
    subject_chr = subject_regions$chr[subject_idx],
    subject_start = subject_regions$start[subject_idx],
    subject_end = subject_regions$end[subject_idx],
    distance = distances,
    overlaps = distances == 0,
    stringsAsFactors = FALSE
  )

  # Handle queries with no match on the same chromosome
  # distanceToNearest returns NA for these - we need to handle them
  if (length(query_idx) < nrow(query_regions)) {
    missing_idx <- setdiff(seq_len(nrow(query_regions)), query_idx)

    if (verbose && length(missing_idx) > 0) {
      cat("  Queries with no subject on same chromosome:", length(missing_idx), "\n")
    }

    missing_rows <- data.frame(
      query_id = query_regions[[query_id_col]][missing_idx],
      query_chr = query_regions$chr[missing_idx],
      query_start = query_regions$start[missing_idx],
      query_end = query_regions$end[missing_idx],
      nearest_subject_id = NA_character_,
      subject_chr = NA_character_,
      subject_start = NA_integer_,
      subject_end = NA_integer_,
      distance = NA_integer_,
      overlaps = FALSE,
      stringsAsFactors = FALSE
    )

    result <- rbind(result, missing_rows)
    # Restore original order
    result <- result[order(match(result$query_id, query_regions[[query_id_col]])), ]
    rownames(result) <- NULL
  }

  # ---------------------------------------------------------------------------
  # Summary statistics
  # ---------------------------------------------------------------------------
  if (verbose) {
    n_overlap <- sum(result$overlaps, na.rm = TRUE)
    n_valid <- sum(!is.na(result$distance))

    cat("\nResults:\n")
    cat("  Overlapping (distance = 0):", n_overlap,
        sprintf("(%.1f%%)\n", 100 * n_overlap / nrow(result)))
    cat("  Non-overlapping:", n_valid - n_overlap, "\n")

    if (n_valid > n_overlap) {
      non_overlap_dist <- result$distance[!result$overlaps & !is.na(result$distance)]
      cat("  Distance distribution (non-overlapping):\n")
      cat("    Min:", min(non_overlap_dist), "bp\n")
      cat("    Median:", median(non_overlap_dist), "bp\n")
      cat("    Mean:", round(mean(non_overlap_dist)), "bp\n")
      cat("    Max:", max(non_overlap_dist), "bp\n")

      # Binned summary
      bins <- c(0, 100, 500, 1000, 5000, 10000, Inf)
      bin_labels <- c("<100bp", "100-500bp", "500bp-1kb", "1-5kb", "5-10kb", ">10kb")
      binned <- cut(non_overlap_dist, breaks = bins, labels = bin_labels, right = FALSE)
      cat("    Binned:\n")
      for (i in seq_along(bin_labels)) {
        n_in_bin <- sum(binned == bin_labels[i], na.rm = TRUE)
        if (n_in_bin > 0) {
          cat(sprintf("      %s: %d (%.1f%%)\n",
                      bin_labels[i], n_in_bin,
                      100 * n_in_bin / length(non_overlap_dist)))
        }
      }
    }

    if (any(is.na(result$distance))) {
      cat("  No match on chromosome:", sum(is.na(result$distance)), "\n")
    }
  }

  return(result)
}
