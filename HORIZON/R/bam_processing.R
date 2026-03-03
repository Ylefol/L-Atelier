#' Sort and index a BAM file
#'
#' Coordinate-sorts the BAM produced by \code{\link{HORIZON_run_align}} and
#' creates a BAI index file. Sorting is required before feature counting.
#'
#' By default, the unsorted BAM is deleted after sorting to avoid doubling disk
#' usage (set \code{remove_unsorted = FALSE} to keep it).
#'
#' RAM note: total RAM used during sorting is approximately
#' \code{memory_per_thread * threads} MB. Reduce \code{threads} or
#' \code{memory_per_thread} on memory-constrained machines.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param threads Integer. Number of sort threads. Default 4.
#' @param memory_per_thread Integer. Memory per sort thread in MB. Total RAM
#'   used is approximately \code{memory_per_thread * threads} MB. Default 2000.
#' @param remove_unsorted Logical. Delete the unsorted BAM after sorting to
#'   free disk space. Default TRUE.
#'
#' @return Character. Path to the sorted, indexed BAM file, invisibly.
#' @export
HORIZON_sort_index_bam <- function(sample_sheet,
                                   sample_id,
                                   threads           = 4,
                                   memory_per_thread = 2000,
                                   remove_unsorted   = TRUE) {
  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")

  bam_in      <- file.path(out_dir, paste0(sample_id, ".bam"))
  dest_prefix <- file.path(out_dir, paste0(sample_id, "_sorted"))
  sorted_bam  <- paste0(dest_prefix, ".bam")

  if (!file.exists(bam_in)) stop("Unsorted BAM not found: ", bam_in)

  message("Sorting BAM for: ", sample_id)
  Rsamtools::sortBam(
    file        = bam_in,
    destination = dest_prefix,
    maxMemory   = memory_per_thread * threads
  )

  message("Indexing BAM for: ", sample_id)
  Rsamtools::indexBam(sorted_bam)

  if (remove_unsorted) {
    file.remove(bam_in)
    message("Unsorted BAM removed: ", bam_in)
  }

  message("BAM sort + index complete for: ", sample_id)
  invisible(sorted_bam)
}
