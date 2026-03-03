#' Align reads to the genome using Rsubread
#'
#' Wraps \code{\link[Rsubread]{align}} for splice-aware RNA-seq alignment against
#' a pre-built Rsubread index. By default, uses the trimmed FASTQ files produced
#' by \code{\link{HORIZON_run_qc_trim}}. Raw FASTQs from the sample sheet can be
#' used instead by setting \code{use_trimmed = FALSE}.
#'
#' Output BAM is written to \code{<output_dir>/<sample_id>/aligned/}. The BAM
#' is unsorted at this stage; pass it to \code{\link{HORIZON_sort_index_bam}}
#' before counting.
#'
#' RAM note: Rsubread's RAM usage during alignment is largely determined by the
#' index size. Using a split index (\code{index_split = TRUE} in
#' \code{\link{HORIZON_build_index}}) substantially reduces peak RAM. Thread
#' count has a modest additional effect.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param index Character. Path to the Rsubread index basename as passed to
#'   \code{\link{HORIZON_build_index}}.
#' @param use_trimmed Logical. Use trimmed FASTQs from the \code{qc/} subfolder
#'   (output of \code{\link{HORIZON_run_qc_trim}}). Default TRUE. Set to FALSE
#'   to align the raw FASTQs listed in the sample sheet.
#' @param threads Integer. Number of alignment threads. Default 4.
#' @param max_mismatches Integer. Maximum mismatches allowed per read. Default 3.
#' @param ... Additional arguments passed to \code{\link[Rsubread]{align}}
#'   (e.g., \code{minFragLength}, \code{maxFragLength}, \code{unique},
#'   \code{nBestLocations}).
#'
#' @return Character. Path to the output BAM file, invisibly.
#' @export
HORIZON_run_align <- function(sample_sheet,
                              sample_id,
                              index,
                              use_trimmed    = TRUE,
                              threads        = 4,
                              max_mismatches = 3,
                              ...) {
  row    <- .get_sample_row(sample_sheet, sample_id)
  paired <- as.logical(row$paired_end)

  out_dir <- file.path(row$output_dir, sample_id, "aligned")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (use_trimmed) {
    qc_dir <- file.path(row$output_dir, sample_id, "qc")
    r1     <- file.path(qc_dir, paste0(sample_id, "_trimmed_R1.fastq.gz"))
    r2     <- file.path(qc_dir, paste0(sample_id, "_trimmed_R2.fastq.gz"))
  } else {
    r1 <- row$fastq_r1
    r2 <- row$fastq_r2
  }

  if (!file.exists(r1)) stop("R1 FASTQ not found: ", r1)
  if (paired && (is.na(r2) || !file.exists(r2))) stop("R2 FASTQ not found: ", r2)

  bam_out <- file.path(out_dir, paste0(sample_id, ".bam"))

  message("Aligning reads for: ", sample_id)
  message("  Index   : ", index)
  message("  R1      : ", r1)
  if (paired) message("  R2      : ", r2)
  message("  Output  : ", bam_out)

  Rsubread::align(
    index         = index,
    readfile1     = r1,
    readfile2     = if (paired) r2 else NULL,
    output_file   = bam_out,
    type          = "rna",
    nthreads      = threads,
    maxMismatches = max_mismatches,
    ...
  )

  message("Alignment complete for: ", sample_id)
  invisible(bam_out)
}
