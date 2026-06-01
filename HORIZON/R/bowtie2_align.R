#' Align FASTQ reads with Bowtie2
#'
#' Aligns paired-end FASTQ files against a Bowtie2 index, piping output
#' directly through \code{samtools view} (MAPQ filtering) and
#' \code{samtools sort} (coordinate sort) without writing an intermediate
#' SAM file to disk. Output BAM is coordinate-sorted and ready for indexing.
#' Applies \code{--no-mixed --no-discordant} by default (appropriate for all
#' paired-end chromatin assays).
#'
#' Requires \code{bowtie2} and \code{samtools} to be present in the HORIZON
#' conda environment registered via \code{\link{HORIZON_set_conda_env}}.
#' Install via: \code{conda install -n horizon_cli -c bioconda bowtie2 samtools}.
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
#' @param threads Integer. Number of threads passed to bowtie2 (\code{--threads}).
#'   Default 4.
#' @param min_mapq Integer. Minimum mapping quality. Reads below this threshold
#'   are discarded by \code{samtools view -q}. Default 30.
#' @param force Logical. If \code{FALSE} (default), skip alignment when the
#'   output \code{_mapq_filtered.bam} already exists.  Set \code{TRUE} to
#'   re-align and overwrite.
#' @param ... Additional bowtie2 flags as character strings (e.g.
#'   \code{"--very-sensitive"}). Appended verbatim to the bowtie2 command.
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

  env_path <- getOption("horizon.conda_env", default = NULL)
  if (is.null(env_path))
    stop("HORIZON conda environment is not set.\n",
         "  Call HORIZON_set_conda_env('/path/to/conda/envs/horizon_cli') ",
         "before running pipeline functions.", call. = FALSE)

  bowtie2_bin  <- file.path(env_path, "bin", "bowtie2")
  samtools_bin <- file.path(env_path, "bin", "samtools")

  if (!file.exists(bowtie2_bin))
    stop("bowtie2 not found in conda env at: ", file.path(env_path, "bin"),
         "\n  Install via: conda install -n horizon_cli -c bioconda bowtie2",
         call. = FALSE)
  if (!file.exists(samtools_bin))
    stop("samtools not found in conda env at: ", file.path(env_path, "bin"),
         "\n  Install via: conda install -n horizon_cli -c bioconda samtools",
         call. = FALSE)

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
  filt_bam <- file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))
  if (!isTRUE(force) && file.exists(filt_bam)) {
    cat("Aligned BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to re-align)")
    return(invisible(filt_bam))
  }

  # Check index ----------------------------------------------------------------
  if (!file.exists(paste0(index, ".1.bt2")))
    stop("Bowtie2 index not found: ", index,
         ". Run HORIZON_build_bowtie2_index() first.", call. = FALSE)

  # Build and run pipeline: bowtie2 | samtools view (no SAM written to disk) --
  extra_flags <- paste(c(...), collapse = " ")

  cmd <- paste(
    shQuote(bowtie2_bin),
    "-x", shQuote(index),
    "-1", shQuote(r1),
    "-2", shQuote(r2),
    "--threads", as.integer(threads),
    "--no-mixed",
    "--no-discordant",
    "--no-unal",
    "-X", as.integer(max_insert),
    if (nzchar(extra_flags)) extra_flags else "",
    "|",
    shQuote(samtools_bin), "view -bS",
    "-q", as.integer(min_mapq),
    "|",
    shQuote(samtools_bin), "sort",
    "--threads", as.integer(threads),
    "-o", shQuote(filt_bam)
  )

  cat("Running Bowtie2 alignment for: ", sample_id,
          " [max_insert=", max_insert, ", min_mapq=", min_mapq, "]")

  ret <- system(cmd)
  if (ret != 0L)
    stop("Bowtie2/samtools pipeline failed for sample: ", sample_id,
         " (exit code ", ret, ")", call. = FALSE)

  cat("Bowtie2 alignment complete for: ", sample_id,
          "\n  Output: ", filt_bam)
  invisible(filt_bam)
}
