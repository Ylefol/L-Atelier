#' Align FASTQ reads with Bowtie2
#'
#' Wraps \code{Rbowtie2::bowtie2()} to align paired-end FASTQ files against a
#' Bowtie2 index.  Applies \code{--no-mixed --no-discordant} by default
#' (appropriate for all paired-end chromatin assays) and filters low-MAPQ
#' alignments immediately after conversion to BAM.
#'
#' All chromatin assay processing assumes \strong{paired-end} sequencing.
#' Single-end data is not supported.
#'
#' \strong{Recommended \code{max_insert} values by assay:}
#' \itemize{
#'   \item ATAC-seq — \code{2000} (fragments span up to tri-nucleosome arrays)
#'   \item ChIP-seq / CUT&RUN / CUT&TAG — \code{1000} (shorter sonication/tagmentation fragments)
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param index Character. Bowtie2 index basename from
#'   \code{\link{HORIZON_build_bowtie2_index}}.
#' @param max_insert Integer. Maximum fragment insert size (\code{-X}).
#'   Default \code{2000L}.
#' @param use_trimmed Logical. If \code{TRUE} (default), use Rfastp-trimmed
#'   FASTQs from \code{<output_dir>/<sample_id>/qc/}.  If \code{FALSE}, use
#'   the raw FASTQs from the sample sheet (\code{fastq_r1} / \code{fastq_r2}).
#' @param threads Integer. Number of threads. Default 4.
#' @param min_mapq Integer. Minimum mapping quality for post-alignment
#'   filtering via \code{Rsamtools::filterBam}.  Reads below this threshold
#'   are discarded.  Default 30.
#' @param force Logical. If \code{FALSE} (default), skip alignment when the
#'   output \code{_mapq_filtered.bam} already exists.  Set \code{TRUE} to
#'   re-align and overwrite.
#' @param ... Additional flags passed to \code{Rbowtie2::bowtie2}.
#'
#' @return Character. Path to the MAPQ-filtered BAM file, invisibly.
#'   Pass to \code{\link{HORIZON_process_bam}}.
#' @export
HORIZON_run_bowtie2 <- function(sample_sheet,
                                 sample_id,
                                 index,
                                 max_insert  = 2000L,
                                 use_trimmed = TRUE,
                                 threads     = 4L,
                                 min_mapq    = 30L,
                                 force       = FALSE,
                                 ...) {

  if (!requireNamespace("Rbowtie2", quietly = TRUE))
    stop("Package 'Rbowtie2' is required. Install via BiocManager::install('Rbowtie2').",
         call. = FALSE)
  if (!requireNamespace("Rsamtools", quietly = TRUE))
    stop("Package 'Rsamtools' is required.", call. = FALSE)

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Resolve input FASTQs -------------------------------------------------------
  if (isTRUE(use_trimmed)) {
    qc_dir <- file.path(row$output_dir, sample_id, "qc")
    r1 <- file.path(qc_dir, paste0(sample_id, "_trimmed_R1.fastq.gz"))
    r2 <- file.path(qc_dir, paste0(sample_id, "_trimmed_R2.fastq.gz"))
    if (!file.exists(r1) || !file.exists(r2))
      stop("Trimmed FASTQs not found under ", qc_dir,
           ". Run HORIZON_run_qc_trim() first, or set use_trimmed=FALSE.",
           call. = FALSE)
  } else {
    r1 <- row$fastq_r1
    r2 <- row$fastq_r2
    if (!file.exists(r1) || !file.exists(r2))
      stop("Raw FASTQs not found: ", r1, " / ", r2, call. = FALSE)
  }

  # Skip if output already exists ----------------------------------------------
  filt_bam_check <- file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))
  if (!isTRUE(force) && file.exists(filt_bam_check)) {
    message("Aligned BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to re-align)")
    return(invisible(filt_bam_check))
  }

  # Check index ----------------------------------------------------------------
  if (!file.exists(paste0(index, ".1.bt2")))
    stop("Bowtie2 index not found: ", index,
         ". Run HORIZON_build_bowtie2_index() first.", call. = FALSE)

  # Paths ----------------------------------------------------------------------
  sam_path <- file.path(out_dir, paste0(sample_id, "_raw.sam"))
  raw_bam  <- file.path(out_dir, paste0(sample_id, "_raw.bam"))
  filt_bam <- file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))

  # Align: Rbowtie2 writes a SAM file -----------------------------------------
  message("Running Bowtie2 alignment for: ", sample_id,
          " [max_insert=", max_insert, "]")

  Rbowtie2::bowtie2(
    bt2Index  = index,
    samOutput = sam_path,
    seq1      = r1,
    seq2      = r2,
    paste0("--threads ", as.integer(threads)),
    "--no-mixed",        # discard reads where only one mate aligns
    "--no-discordant",   # discard pairs aligning in unexpected orientation
    "--no-unal",         # suppress unaligned reads from SAM output
    paste0("-X ", as.integer(max_insert)),
    ...
  )

  # SAM → BAM ------------------------------------------------------------------
  Rsamtools::asBam(sam_path, destination = tools::file_path_sans_ext(raw_bam),
                   overwrite = TRUE)
  file.remove(sam_path)

  # MAPQ filter ----------------------------------------------------------------
  message("Filtering MAPQ < ", min_mapq, " for: ", sample_id)
  filter_param <- Rsamtools::ScanBamParam(
    mapqFilter = as.integer(min_mapq)
  )
  Rsamtools::filterBam(raw_bam, destination = filt_bam,
                       param = filter_param)
  file.remove(raw_bam)
  # Remove accompanying index if created
  bai <- paste0(raw_bam, ".bai")
  if (file.exists(bai)) file.remove(bai)

  message("Bowtie2 alignment complete for: ", sample_id,
          "\n  Output: ", filt_bam)
  invisible(filt_bam)
}
