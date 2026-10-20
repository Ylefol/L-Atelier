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

  cat("Sorting BAM for: ", sample_id)
  Rsamtools::sortBam(
    file        = bam_in,
    destination = dest_prefix,
    maxMemory   = memory_per_thread * threads
  )

  cat("Indexing BAM for: ", sample_id)
  Rsamtools::indexBam(sorted_bam)

  if (remove_unsorted) {
    file.remove(bam_in)
    cat("Unsorted BAM removed: ", bam_in)
  }

  cat("BAM sort + index complete for: ", sample_id)
  invisible(sorted_bam)
}


#' Mark duplicate reads on a sorted BAM file
#'
#' Runs the standard samtools duplicate-marking pipeline (collate -> fixmate
#' -> sort -> markdup) on the coordinate-sorted BAM produced by
#' \code{\link{HORIZON_sort_index_bam}}. By default, duplicates are
#' FLAGGED (SAM FLAG 0x400) but not removed from the file.
#'
#' @details
#' Unlike \code{\link{HORIZON_process_bam}} (bam_atac_chip.R, the ATAC/ChIP
#' pipeline), which hardcodes \code{markdup -r} to remove duplicates
#' outright, this function defaults to marking only. For standard bulk
#' RNA-seq quantification, removing duplicate reads is not recommended:
#' unlike DNA-seq, where a duplicate read pair almost always reflects PCR
#' amplification bias, RNA-seq "duplicates" are frequently genuine
#' biological signal -- a highly expressed transcript naturally produces
#' many identical fragments, and removing them systematically under-counts
#' highly expressed genes. Flagging without removing is harmless for
#' standard counting (\code{\link{HORIZON_run_count}}/featureCounts counts
#' all reads by default regardless of the duplicate flag) while still
#' satisfying tools that need to exclude duplicates themselves downstream
#' -- e.g. \code{CYAN::CYAN_run_reditools(require_dup_marked = TRUE)} for
#' RNA editing detection, where PCR duplicates do need excluding from
#' editing-frequency estimates, but REDItools2 performs that exclusion
#' itself based on the flag, so removing reads beforehand is neither needed
#' nor wanted. Set \code{remove_duplicates = TRUE} only when a use case
#' specifically calls for it.
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param threads Integer. Threads for all samtools commands. Default 4.
#' @param memory_per_thread Integer. Memory in MB per thread for
#'   \code{samtools sort} (\code{-m}). Default 2000.
#' @param remove_duplicates Logical. If \code{TRUE}, duplicates are removed
#'   from the output BAM (\code{markdup -r}) instead of just flagged.
#'   Default \code{FALSE} -- standard RNA-seq quantification should not
#'   remove duplicate reads (see Details); flagging alone is sufficient for
#'   tools that need to exclude them (e.g. REDItools2 via CYAN).
#' @param remove_tmp Logical. Remove intermediate BAM files after
#'   processing. Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip processing when
#'   \code{_markdup.bam} already exists. Set \code{TRUE} to reprocess.
#'
#' @return Character. Path to the duplicate-marked BAM file, invisibly.
#' @export
HORIZON_mark_duplicates <- function(sample_sheet,
                                     sample_id,
                                     threads            = 4L,
                                     memory_per_thread  = 2000L,
                                     remove_duplicates  = FALSE,
                                     remove_tmp         = TRUE,
                                     force              = FALSE) {

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")

  input_bam    <- file.path(out_dir, paste0(sample_id, "_sorted.bam"))
  collate_bam  <- file.path(out_dir, paste0(sample_id, "_collated.bam"))
  fixmate_bam  <- file.path(out_dir, paste0(sample_id, "_fixmate.bam"))
  resorted_bam <- file.path(out_dir, paste0(sample_id, "_resorted.bam"))
  markdup_bam  <- file.path(out_dir, paste0(sample_id, "_markdup.bam"))
  markdup_stat <- file.path(out_dir, paste0(sample_id, "_markdup_stats.txt"))
  flagstat_f   <- file.path(out_dir, paste0(sample_id, "_markdup_flagstat.txt"))

  if (!file.exists(input_bam)) stop("Sorted BAM not found: ", input_bam)

  if (!isTRUE(force) && file.exists(markdup_bam)) {
    cat("Duplicate-marked BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to reprocess)\n")
    return(invisible(markdup_bam))
  }

  thr <- as.integer(threads)
  mem <- paste0(as.integer(memory_per_thread), "M")

  # Helper to run a samtools command via the Herper CLI env
  .run_samtools <- function(...) {
    exit <- .horizon_run_cli("samtools", c(...))
    if (exit != 0)
      stop("samtools exited with code ", exit, call. = FALSE)
    invisible(exit)
  }

  # Step 1: collate (name-sort) ------------------------------------------------
  cat("[", sample_id, "] samtools collate")
  .run_samtools("collate", "-u",
                "--threads", thr,
                input_bam,
                file.path(out_dir, paste0(sample_id, "_collate_tmp")),
                "-o", collate_bam)

  # Step 2: fixmate ------------------------------------------------------------
  cat("[", sample_id, "] samtools fixmate")
  .run_samtools("fixmate", "-m", "-u",
                "--threads", thr,
                collate_bam, fixmate_bam)
  if (isTRUE(remove_tmp)) file.remove(collate_bam)

  # Step 3: sort (coordinate) --------------------------------------------------
  cat("[", sample_id, "] samtools sort")
  .run_samtools("sort",
                "--threads", thr,
                "-m", mem,
                "-T", file.path(out_dir, paste0(sample_id, "_sort_tmp")),
                fixmate_bam, "-o", resorted_bam)
  if (isTRUE(remove_tmp)) file.remove(fixmate_bam)

  # Step 4: markdup + stats -----------------------------------------------------
  cat("[", sample_id, "] samtools markdup (",
          if (isTRUE(remove_duplicates)) "removing duplicates" else "flagging only, not removing",
          ")")
  .run_samtools("markdup",
                "--threads", thr,
                if (isTRUE(remove_duplicates)) "-r" else NULL,  # remove vs flag-only
                "-s",             # write stats to stderr (captured below)
                "-f", markdup_stat,
                resorted_bam, markdup_bam)
  if (isTRUE(remove_tmp)) file.remove(resorted_bam)

  # Step 5: index final BAM -----------------------------------------------------
  cat("[", sample_id, "] Indexing duplicate-marked BAM")
  .run_samtools("index", "--threads", thr, markdup_bam)

  # Step 6: flagstat -------------------------------------------------------------
  cat("[", sample_id, "] Running flagstat")
  .horizon_run_cli(
    "samtools",
    c("flagstat", "--threads", thr, markdup_bam),
    stdout = flagstat_f
  )

  cat("Duplicate marking complete for: ", sample_id,
          "\n  Output: ", markdup_bam)
  invisible(markdup_bam)
}
