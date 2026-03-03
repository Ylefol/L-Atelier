#' Convert a sorted BAM to an RPKM-normalised BigWig using bamCoverage
#'
#' Wraps the deeptools \code{bamCoverage} command-line tool to produce an
#' RPKM-normalised BigWig from a sorted, indexed BAM file.  The BigWig is
#' written alongside the BAM in the
#' \code{<output_dir>/<sample_id>/aligned/} directory.
#'
#' \code{bamCoverage} must be installed and available on \code{PATH}.
#' Install via conda (\code{conda install -c bioconda deeptools}) or
#' pip (\code{pip install deeptools}).
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param effective_genome_size Integer, numeric, or character. The effective
#'   (mappable) genome size passed to \code{--effectiveGenomeSize}.  Pass a
#'   numeric value directly, or one of the following genome assembly names for
#'   a built-in lookup:
#'   \itemize{
#'     \item \code{"GRCh38"}   — 2913022398
#'     \item \code{"GRCm38"}   — 2652783500
#'     \item \code{"WBcel235"} — 100272607
#'     \item \code{"T2TCHM13"} — 3117292070
#'   }
#' @param threads Integer. Number of parallel threads (\code{-p}). Default 4.
#' @param bin_size Integer. Bin width in base pairs (\code{--binSize}).
#'   Default 10, matching the bamCoverage default.
#' @param ignore_for_normalization Character vector. Chromosome names excluded
#'   when computing the total read count for RPKM normalisation
#'   (\code{--ignoreForNormalization}). Default \code{c("X", "M")} for
#'   Ensembl-style chromosome names; use \code{c("chrX", "chrM")} for
#'   UCSC-style genomes.
#' @param ignore_duplicates Logical. Skip duplicate reads
#'   (\code{--ignoreDuplicates}). Default TRUE.
#' @param extend_reads Logical. Extend reads to fragment size for paired-end
#'   data, or to the estimated fragment size for single-end
#'   (\code{--extendReads}). Default TRUE.
#'
#' @return Character. Path to the BigWig file, invisibly.
#' @export
HORIZON_bam_to_bigwig <- function(sample_sheet,
                                   sample_id,
                                   effective_genome_size,
                                   threads                  = 4L,
                                   bin_size                 = 10L,
                                   ignore_for_normalization = c("X", "M"),
                                   ignore_duplicates        = TRUE,
                                   extend_reads             = TRUE) {

  # Resolve effective genome size -----------------------------------------------
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

  if (!nzchar(Sys.which("bamCoverage")))
    stop("'bamCoverage' not found on PATH. ",
         "Install deeptools: conda install -c bioconda deeptools",
         call. = FALSE)

  row      <- .get_sample_row(sample_sheet, sample_id)
  out_dir  <- file.path(row$output_dir, sample_id, "aligned")
  bam_path <- file.path(out_dir, paste0(sample_id, "_sorted.bam"))
  bw_path  <- file.path(out_dir, paste0(sample_id, "_RPKM.bw"))

  if (!file.exists(bam_path))
    stop("Sorted BAM not found: ", bam_path, call. = FALSE)
  if (!file.exists(paste0(bam_path, ".bai")))
    stop("BAM index (.bai) not found. Run HORIZON_sort_index_bam() first.",
         call. = FALSE)

  args <- c(
    "-p",                   as.integer(threads),
    "-b",                   bam_path,
    "-o",                   bw_path,
    "--normalizeUsing",     "RPKM",
    "--effectiveGenomeSize", as.character(effective_genome_size),
    "--binSize",            as.integer(bin_size)
  )

  if (isTRUE(ignore_duplicates))
    args <- c(args, "--ignoreDuplicates")
  if (isTRUE(extend_reads))
    args <- c(args, "--extendReads")
  if (length(ignore_for_normalization) > 0)
    args <- c(args, "--ignoreForNormalization", ignore_for_normalization)

  message("Running bamCoverage for: ", sample_id)
  exit_code <- system2("bamCoverage", args = args, wait = TRUE)

  if (exit_code != 0)
    stop("bamCoverage failed for '", sample_id, "' (exit code ", exit_code, ").",
         call. = FALSE)

  message("BigWig complete for: ", sample_id)
  invisible(bw_path)
}
