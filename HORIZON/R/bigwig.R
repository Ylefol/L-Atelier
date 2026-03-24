#' Convert a processed BAM to a BigWig using bamCoverage
#'
#' Wraps the deeptools \code{bamCoverage} command-line tool.  Supports two
#' normalization modes:
#' \itemize{
#'   \item \strong{RPKM} (default) — suitable for ATAC-seq and ChIP-seq
#'     without a spike-in control.
#'   \item \strong{Spike-in} — pass a numeric \code{scale_factor} from
#'     \code{\link{HORIZON_compute_spike_in_factors}}.  Uses
#'     \code{--normalizeUsing None --scaleFactor <value>}, which preserves
#'     absolute differences in chromatin occupancy across samples.
#' }
#'
#' The most downstream available BAM is used automatically via
#' \code{.latest_bam()} (priority: blacklist-filtered > processed > host >
#' mapq-filtered).
#'
#' \code{bamCoverage} is provided via the user-managed conda environment
#' registered with \code{\link{HORIZON_set_conda_env}}.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param effective_genome_size Integer, numeric, or character. The effective
#'   (mappable) genome size passed to \code{--effectiveGenomeSize}.  Pass a
#'   numeric value directly, or one of the following assembly names:
#'   \itemize{
#'     \item \code{"GRCh38"}   — 2913022398
#'     \item \code{"GRCm38"}   — 2652783500
#'     \item \code{"WBcel235"} — 100272607
#'     \item \code{"T2TCHM13"} — 3117292070
#'   }
#' @param scale_factor Numeric or \code{NULL}.  When \code{NULL} (default),
#'   RPKM normalisation is used.  When a numeric value is supplied (e.g. from
#'   \code{\link{HORIZON_compute_spike_in_factors}}), spike-in calibration
#'   normalisation is used instead (\code{--normalizeUsing None
#'   --scaleFactor <value>}).
#' @param threads Integer. Number of parallel threads (\code{-p}). Default 4.
#' @param bin_size Integer. Bin width in base pairs (\code{--binSize}).
#'   Default 10.
#' @param ignore_for_normalization Character vector. Chromosome names excluded
#'   from total read count during RPKM normalisation
#'   (\code{--ignoreForNormalization}).  Ignored when \code{scale_factor} is
#'   set.  Default \code{c("chrX", "chrM")}.
#' @param ignore_duplicates Logical. Skip duplicate reads
#'   (\code{--ignoreDuplicates}). Default \code{TRUE}.
#' @param extend_reads Logical. Extend reads to fragment size
#'   (\code{--extendReads}). Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip BigWig generation
#'   when the output file already exists.  Set \code{TRUE} to regenerate.
#'
#' @return Character. Path to the BigWig file, invisibly.
#' @export
HORIZON_bam_to_bigwig <- function(sample_sheet,
                                   sample_id,
                                   effective_genome_size,
                                   scale_factor             = NULL,
                                   threads                  = 4L,
                                   bin_size                 = 10L,
                                   ignore_for_normalization = c("chrX", "chrM"),
                                   ignore_duplicates        = TRUE,
                                   extend_reads             = TRUE,
                                   force                    = FALSE) {

  # Resolve effective genome size ----------------------------------------------
  .eff_genome_sizes <- c(
    GRCh38   = 2913022398,
    GRCm38   = 2652783500,
    WBcel235 = 100272607,
    T2TCHM13 = 3117292070
  )
  if (is.character(effective_genome_size)) {
    key <- effective_genome_size
    if (!key %in% names(.eff_genome_sizes))
      stop("Unknown genome assembly '", key, "'. ",
           "Supported names: ", paste(names(.eff_genome_sizes), collapse = ", "),
           ". Pass a numeric value directly for other genomes.",
           call. = FALSE)
    effective_genome_size <- .eff_genome_sizes[[key]]
  }

  if (!is.null(scale_factor) && (!is.numeric(scale_factor) || length(scale_factor) != 1L))
    stop("scale_factor must be a single numeric value or NULL.", call. = FALSE)
  if (!is.null(scale_factor) && is.na(scale_factor))
    stop("scale_factor is NA for sample '", sample_id,
         "'. Check HORIZON_compute_spike_in_factors() output.", call. = FALSE)

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")

  bam_path <- .latest_bam(sample_sheet, sample_id)

  norm_tag <- if (is.null(scale_factor)) "RPKM" else "spikein"
  bw_path  <- file.path(out_dir, paste0(sample_id, "_", norm_tag, ".bw"))

  if (!isTRUE(force) && file.exists(bw_path)) {
    message("BigWig already exists for: ", sample_id,
            " — skipping (use force=TRUE to regenerate)")
    return(invisible(bw_path))
  }

  if (!file.exists(paste0(bam_path, ".bai")))
    stop("BAM index not found for: ", bam_path,
         "\n  Run the relevant processing step to index the BAM.", call. = FALSE)

  # Build bamCoverage arguments ------------------------------------------------
  args <- c(
    "-p",                    as.integer(threads),
    "-b",                    bam_path,
    "-o",                    bw_path,
    "--effectiveGenomeSize", as.character(effective_genome_size),
    "--binSize",             as.integer(bin_size)
  )

  if (is.null(scale_factor)) {
    args <- c(args, "--normalizeUsing", "RPKM")
    if (length(ignore_for_normalization) > 0L)
      args <- c(args, "--ignoreForNormalization", ignore_for_normalization)
  } else {
    args <- c(args, "--normalizeUsing", "None",
              "--scaleFactor", as.character(scale_factor))
    message("[", sample_id, "] Using spike-in scale factor: ", round(scale_factor, 5L))
  }

  if (isTRUE(ignore_duplicates)) args <- c(args, "--ignoreDuplicates")
  if (isTRUE(extend_reads))      args <- c(args, "--extendReads")

  message("Running bamCoverage for: ", sample_id,
          " [normalization=", norm_tag, "]")
  exit_code <- .horizon_run_cli("bamCoverage", args)

  if (exit_code != 0)
    stop("bamCoverage failed for '", sample_id,
         "' (exit code ", exit_code, ").", call. = FALSE)

  message("BigWig complete for: ", sample_id, "\n  Output: ", bw_path)
  invisible(bw_path)
}
