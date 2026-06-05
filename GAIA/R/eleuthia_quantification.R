#' Eleuthia - Quantification Functions
#'
#' @description Functions for quantifying reads against genomic features.
#' Includes both BAM-based (featureCounts) and BED-based quantification.


#' Quantify Reads in Peak Regions (BAM-based)
#'
#' @description Counts reads from BAM files in consensus peak regions using
#' Rsubread::featureCounts().
#'
#' @note This function requires BAM files and a sample sheet with a bam_loc column.
#' For most use cases, consider using ELEUTHIA_quantify_bed() instead, which
#' quantifies directly from BED fragment files and is more lightweight.
#'
#' @param sample_sheet A validated sample sheet data.frame with bam_loc column.
#'   Note: The standard sample sheet uses bed_loc; you may need to add bam_loc
#'   manually if using this function.
#' @param consensus_peaks A consensus peak data.frame from
#'   ELEUTHIA_create_consensus_peaks().
#' @param omics Character string. Omics type to quantify ("ATACseq" or "CHIPseq").
#' @param paired_end Logical. Are the reads paired-end? (default = TRUE for ATAC-seq).
#' @param count_fragments Logical. If paired-end, count fragments instead of reads
#'   (default = TRUE).
#' @param min_mapq Integer. Minimum mapping quality to count a read (default = 10).
#' @param nthreads Integer. Number of threads for parallel processing (default = 4).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{counts}{Matrix of read counts (peaks x samples)}
#'   \item{annotation}{Data.frame with peak annotations}
#'   \item{targets}{Data.frame with sample information}
#'   \item{stat}{Data.frame with counting statistics}
#' }
#'
#' @details
#' This function wraps Rsubread::featureCounts() for ease of use with the
#' GAIA pipeline. It:
#' \enumerate{
#'   \item Extracts BAM file paths from the sample sheet
#'   \item Converts consensus peaks to SAF format
#'   \item Runs featureCounts with appropriate parameters for ATAC/ChIP-seq
#'   \item Returns a structured list with counts and metadata
#' }
#'
#' For ATAC-seq data with paired-end reads, the function counts fragments
#' (read pairs) by default rather than individual reads.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load and create consensus peaks
#' peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' consensus <- ELEUTHIA_create_consensus_peaks(peak_list, min_overlap = 2)
#'
#' # Quantify reads
#' counts <- ELEUTHIA_quantify_peaks(sample_sheet, consensus, "ATACseq")
#'
#' # Access count matrix
#' count_matrix <- counts$counts
#'
#' }
ELEUTHIA_quantify_peaks <- function(sample_sheet,
                                     consensus_peaks,
                                     omics,
                                     paired_end = TRUE,
                                     count_fragments = TRUE,
                                     min_mapq = 10,
                                     nthreads = 4,
                                     verbose = TRUE) {

  # Check if Rsubread is available
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop("Package 'Rsubread' is required for quantification. ",
         "Install with: BiocManager::install('Rsubread')")
  }

  # Check for bam_loc column
  if (!"bam_loc" %in% colnames(sample_sheet)) {
    stop("sample_sheet must have a 'bam_loc' column for BAM-based quantification.\n",
         "Note: The standard sample sheet uses 'bed_loc' for BED fragment files.\n",
         "Consider using ELEUTHIA_quantify_bed() instead, which works with BED files.\n",
         "If you need BAM-based quantification, add a 'bam_loc' column to your sample sheet.")
  }

  # Get subset for this omics type
  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type: ", omics)
  }

  # Get BAM file paths
  bam_files <- subset_df$bam_loc
  names(bam_files) <- subset_df$sample_id

  # Check BAM files exist
  missing_bams <- bam_files[!file.exists(bam_files)]
  if (length(missing_bams) > 0) {
    stop("BAM files not found:\n  ",
         paste(names(missing_bams), ":", missing_bams, collapse = "\n  "))
  }

  if (verbose) {
    cat("[ELEUTHIA] Quantifying", length(bam_files), "BAM files against",
        nrow(consensus_peaks), "consensus peaks...\n")
  }

  # Convert peaks to SAF format
  saf <- ELEUTHIA_peaks_to_saf(consensus_peaks)

  # Run featureCounts
  if (verbose) cat("[ELEUTHIA] Running featureCounts...\n")

  fc_result <- Rsubread::featureCounts(
    files = bam_files,
    annot.ext = saf,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = paired_end,
    countReadPairs = count_fragments,
    minMQS = min_mapq,
    nthreads = nthreads,
    verbose = verbose
  )

  # Extract and clean up count matrix
  counts <- fc_result$counts

  # Clean column names (remove path, keep sample_id)
  colnames(counts) <- subset_df$sample_id

  # Create targets data.frame — include replicate/batch columns only if present
  meta_cols <- intersect(c("sample_id", "group", "bio_rep", "tech_rep", "batch"),
                         colnames(subset_df))
  targets <- subset_df[, meta_cols, drop = FALSE]
  rownames(targets) <- targets$sample_id

  if (verbose) {
    cat("[ELEUTHIA] Quantification complete.\n")
    cat("    Count matrix dimensions:", nrow(counts), "peaks x", ncol(counts), "samples\n")
    cat("    Total counts:", sum(counts), "\n")
    cat("    Mean counts per peak:", round(mean(rowSums(counts)), 1), "\n")
  }

  return(list(
    counts = counts,
    annotation = consensus_peaks,
    targets = targets,
    stat = fc_result$stat
  ))
}


#' Summarize Quantification Statistics
#'
#' @description Prints a summary of quantification results. Works with both
#' BAM-based (featureCounts) and BED-based quantification outputs.
#'
#' @param quant_result Result from ELEUTHIA_quantify_peaks() or ELEUTHIA_quantify_bed().
#'
#' @return Invisibly returns a summary data.frame with per-sample statistics.
#'
#' @export
#'
ELEUTHIA_summarize_quantification <- function(quant_result) {

  cat("================================================================================\n")
  cat("[ELEUTHIA] QUANTIFICATION SUMMARY\n")
  cat("================================================================================\n\n")

  counts <- quant_result$counts
  n_regions <- nrow(counts)
  n_samples <- ncol(counts)
  sample_names <- colnames(counts)

  # Check if this is featureCounts output (has stat element)
  if (!is.null(quant_result$stat)) {
    stat <- quant_result$stat
    assigned_row <- which(stat$Status == "Assigned")

    if (length(assigned_row) > 0) {
      assigned <- as.numeric(stat[assigned_row, -1])
      total <- colSums(stat[, -1, drop = FALSE])
      pct_assigned <- round(100 * assigned / total, 1)

      cat("[ELEUTHIA] Assignment rates (featureCounts):\n")
      fc_sample_names <- colnames(stat)[-1]
      for (i in seq_along(fc_sample_names)) {
        cat(sprintf("    %-30s: %s assigned (%.1f%%)\n",
                    fc_sample_names[i],
                    format(assigned[i], big.mark = ","),
                    pct_assigned[i]))
      }

      cat("    Overall:\n")
      cat("    Mean assignment rate:", round(mean(pct_assigned), 1), "%\n")
      cat("    Total assigned reads:", format(sum(assigned), big.mark = ","), "\n")
    }
  } else {
    # BED-based quantification - summarize from counts matrix
    cat("    Regions:", format(n_regions, big.mark = ","), "\n")
    cat("    Samples:", n_samples, "\n\n")

    cat("    Per-sample fragment counts:\n")
    sample_totals <- colSums(counts)
    sample_means <- colMeans(counts)
    sample_nonzero <- apply(counts, 2, function(x) sum(x > 0))

    for (i in seq_along(sample_names)) {
      cat(sprintf("    %-30s: %10s total | %6.1f mean | %s regions with signal\n",
                  sample_names[i],
                  format(sample_totals[i], big.mark = ","),
                  sample_means[i],
                  format(sample_nonzero[i], big.mark = ",")))
    }

    cat("[ELEUTHIA] Overall:\n")
    cat("    Total fragments counted:", format(sum(counts), big.mark = ","), "\n")
    cat("    Mean fragments per region:", round(mean(rowSums(counts)), 1), "\n")
    cat("    Median fragments per region:", round(median(rowSums(counts)), 1), "\n")
    cat("    Regions with zero counts:", format(sum(rowSums(counts) == 0), big.mark = ","),
        sprintf("(%.1f%%)\n", 100 * sum(rowSums(counts) == 0) / n_regions))
  }

  cat("================================================================================\n")

  # Create summary data.frame
  summary_df <- data.frame(
    sample = sample_names,
    total_counts = colSums(counts),
    mean_per_region = colMeans(counts),
    regions_with_signal = apply(counts, 2, function(x) sum(x > 0)),
    row.names = NULL
  )

  invisible(summary_df)
}
