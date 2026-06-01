#' Count reads per feature using featureCounts
#'
#' Wraps \code{\link[Rsubread]{featureCounts}} to quantify reads mapping to
#' genomic features for a single sample. Strandedness is read directly from the
#' sample sheet and mapped to featureCounts integer codes (0 = unstranded,
#' 1 = forward, 2 = reverse).
#'
#' Two files are written to \code{<output_dir>/<sample_id>/counts/}:
#' \itemize{
#'   \item \code{<sample_id>_counts.txt} — tab-delimited table of gene ID and
#'     raw count (used by \code{\link{HORIZON_aggregate_counts}})
#'   \item \code{<sample_id>_featurecounts.rds} — full featureCounts result
#'     object including alignment statistics and annotation
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
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
HORIZON_run_count <- function(sample_sheet,
                              sample_id,
                              annotation,
                              feature_type        = "exon",
                              attribute_type      = "gene_id",
                              threads             = 4,
                              min_mapping_quality = 10,
                              ...) {
  row    <- .get_sample_row(sample_sheet, sample_id)
  paired <- as.logical(row$paired_end)

  out_dir <- file.path(row$output_dir, sample_id, "counts")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sorted_bam <- file.path(row$output_dir, sample_id, "aligned",
                          paste0(sample_id, "_sorted.bam"))

  if (!file.exists(sorted_bam)) {
    stop("Sorted BAM not found: ", sorted_bam,
         "\nRun HORIZON_sort_index_bam() first.")
  }

  if (!file.exists(annotation)) stop("Annotation file not found: ", annotation)

  strandedness <- .strand_to_int(row$strandedness)

  counts_out <- file.path(out_dir, paste0(sample_id, "_counts.txt"))
  rds_out    <- file.path(out_dir, paste0(sample_id, "_featurecounts.rds"))

  cat("Counting reads for: ", sample_id)
  cat("  Strandedness: ", row$strandedness, " (", strandedness, ")")

  fc <- Rsubread::featureCounts(
    files               = sorted_bam,
    annot.ext           = annotation,
    isGTFAnnotationFile = TRUE,
    GTF.featureType     = feature_type,
    GTF.attrType        = attribute_type,
    isPairedEnd         = paired,
    strandSpecific      = strandedness,
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
