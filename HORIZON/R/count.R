#' Count reads per feature using featureCounts
#'
#' Wraps \code{\link[Rsubread]{featureCounts}} to quantify reads mapping to
#' genomic features for a single sample. Strandedness is mapped to
#' featureCounts integer codes (0 = unstranded, 1 = forward, 2 = reverse).
#'
#' \strong{Two usage modes:}
#' \enumerate{
#'   \item \strong{Sample-sheet mode} (default): provide \code{sample_sheet} and
#'     \code{sample_id}. \code{paired_end} and \code{strandedness} are read from
#'     the sample sheet row, and the sorted BAM is located at the standard
#'     HORIZON path (\code{<output_dir>/<sample_id>/aligned/<sample_id>_sorted.bam}).
#'     Output is written to \code{<output_dir>/<sample_id>/counts/}.
#'   \item \strong{Direct BAM mode}: provide \code{bam_file} (path to a sorted
#'     BAM), plus \code{paired_end} and \code{strandedness} explicitly, since
#'     neither can be reliably inferred without a sample sheet row.
#'     \code{sample_id} is inferred from the filename if omitted. Output goes
#'     to \code{output_dir} (defaults to the directory containing the BAM).
#' }
#'
#' Two files are written to the resolved output directory:
#' \itemize{
#'   \item \code{<sample_id>_counts.txt} — tab-delimited table of gene ID and
#'     raw count (used by \code{\link{HORIZON_aggregate_counts}})
#'   \item \code{<sample_id>_featurecounts.rds} — full featureCounts result
#'     object including alignment statistics and annotation
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}. Required in sample-sheet
#'   mode; must be \code{NULL} when \code{bam_file} is supplied.
#' @param sample_id Character. Sample ID to process (sample-sheet mode).
#'   Optional in direct BAM mode — inferred from the filename if omitted.
#' @param bam_file Character or \code{NULL}. Path to a sorted BAM file. When
#'   supplied, \code{sample_sheet} must be \code{NULL}. Default \code{NULL}.
#' @param output_dir Character or \code{NULL}. Output directory for direct BAM
#'   mode. Ignored in sample-sheet mode. Defaults to the directory containing
#'   \code{bam_file}.
#' @param paired_end Logical or \code{NULL}. Whether the BAM is paired-end.
#'   Required in direct BAM mode; read from the sample sheet in sample-sheet
#'   mode.
#' @param strandedness Character or \code{NULL}. One of \code{"unstranded"},
#'   \code{"forward"}, \code{"reverse"}. Required in direct BAM mode; read
#'   from the sample sheet in sample-sheet mode.
#' @param annotation Character. Path to a GTF/GFF annotation file. Passed to
#'   \code{annot.ext} in \code{\link[Rsubread]{featureCounts}}.
#' @param feature_type Character. Feature type in the GTF to count reads
#'   against. Default \code{"exon"}.
#' @param attribute_type Character. GTF attribute used to group features into
#'   meta-features (genes). Default \code{"gene_id"}.
#' @param threads Integer. Number of threads. Default 4.
#' @param min_mapping_quality Integer. Minimum mapping quality score for a read
#'   to be counted. Default 10.
#' @param ... Additional arguments passed to \code{\link[Rsubread]{featureCounts}}
#'   (e.g., \code{allowMultiOverlap}, \code{countMultiMappingReads},
#'   \code{requireBothEndsMapped}, \code{minFragLength}, \code{maxFragLength}).
#'
#' @return Character. Path to the counts text file, invisibly.
#' @export
HORIZON_run_count <- function(sample_sheet        = NULL,
                              sample_id            = NULL,
                              bam_file             = NULL,
                              output_dir           = NULL,
                              paired_end           = NULL,
                              strandedness         = NULL,
                              annotation,
                              feature_type        = "exon",
                              attribute_type      = "gene_id",
                              threads             = 4,
                              min_mapping_quality = 10,
                              ...) {

  if (!is.null(bam_file) && !is.null(sample_sheet))
    stop("Provide either 'bam_file' or 'sample_sheet', not both.", call. = FALSE)
  if (is.null(bam_file) && is.null(sample_sheet))
    stop("One of 'bam_file' or 'sample_sheet' must be provided.", call. = FALSE)

  if (!is.null(bam_file)) {
    # Direct BAM mode
    if (!file.exists(bam_file))
      stop("BAM file not found: ", bam_file, call. = FALSE)
    if (is.null(sample_id))
      sample_id <- tools::file_path_sans_ext(
                     tools::file_path_sans_ext(basename(bam_file)))
    if (is.null(paired_end))
      stop("'paired_end' is required in direct BAM mode.", call. = FALSE)
    if (is.null(strandedness))
      stop("'strandedness' is required in direct BAM mode.", call. = FALSE)
    sorted_bam <- bam_file
    out_dir    <- if (!is.null(output_dir)) output_dir else dirname(bam_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  } else {
    # Sample-sheet mode
    if (is.null(sample_id))
      stop("'sample_id' is required when using sample-sheet mode.", call. = FALSE)
    row          <- .get_sample_row(sample_sheet, sample_id)
    paired_end   <- as.logical(row$paired_end)
    strandedness <- row$strandedness
    out_dir      <- file.path(row$output_dir, sample_id, "counts")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    sorted_bam   <- file.path(row$output_dir, sample_id, "aligned",
                              paste0(sample_id, "_sorted.bam"))
    if (!file.exists(sorted_bam)) {
      stop("Sorted BAM not found: ", sorted_bam,
           "\nRun HORIZON_sort_index_bam() first.")
    }
  }

  if (!file.exists(annotation)) stop("Annotation file not found: ", annotation)

  strand_int <- .strand_to_int(strandedness)

  counts_out <- file.path(out_dir, paste0(sample_id, "_counts.txt"))
  rds_out    <- file.path(out_dir, paste0(sample_id, "_featurecounts.rds"))

  cat("Counting reads for: ", sample_id)
  cat("  Strandedness: ", strandedness, " (", strand_int, ")")

  fc <- Rsubread::featureCounts(
    files               = sorted_bam,
    annot.ext           = annotation,
    isGTFAnnotationFile = TRUE,
    GTF.featureType     = feature_type,
    GTF.attrType        = attribute_type,
    isPairedEnd         = paired_end,
    strandSpecific      = strand_int,
    nthreads            = threads,
    minMQS              = min_mapping_quality,
    ...
  )

  # Write a clean tab-delimited counts file for downstream aggregation
  counts_df <- data.frame(
    gene_id = rownames(fc$counts),
    count   = as.integer(fc$counts[, 1]),
    stringsAsFactors = FALSE
  )
  write.table(counts_df, counts_out, sep = "\t", quote = FALSE, row.names = FALSE)

  # Save the full featureCounts result for diagnostics / alignment stats
  saveRDS(fc, rds_out)

  cat("Counts written to: ", counts_out)
  invisible(counts_out)
}
