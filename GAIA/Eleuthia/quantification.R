#' Eleuthia - Quantification Functions
#'
#' @description Functions for quantifying reads from BAM files against
#' genomic features (peaks, genes, etc.) using Rsubread.


#' Quantify Reads in Peak Regions
#'
#' @description Counts reads from BAM files in consensus peak regions using
#' Rsubread::featureCounts().
#'
#' @param sample_sheet A validated sample sheet data.frame with bam_loc column.
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
    cat("Quantifying", length(bam_files), "BAM files against",
        nrow(consensus_peaks), "consensus peaks...\n")
  }

  # Convert peaks to SAF format
  saf <- ELEUTHIA_peaks_to_saf(consensus_peaks)

  # Run featureCounts
  if (verbose) cat("Running featureCounts...\n")

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

  # Create targets data.frame (sample metadata)
  targets <- subset_df[, c("sample_id", "group", "bio_rep", "tech_rep", "batch")]
  rownames(targets) <- targets$sample_id

  if (verbose) {
    cat("\nQuantification complete.\n")
    cat("  Count matrix dimensions:", nrow(counts), "peaks x", ncol(counts), "samples\n")
    cat("  Total counts:", sum(counts), "\n")
    cat("  Mean counts per peak:", round(mean(rowSums(counts)), 1), "\n")
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
#' @description Prints a summary of featureCounts statistics including
#' assignment rates per sample.
#'
#' @param quant_result Result from ELEUTHIA_quantify_peaks().
#'
#' @return Invisibly returns a summary data.frame.
#'
#' @export
#'
ELEUTHIA_summarize_quantification <- function(quant_result) {

  stat <- quant_result$stat

  cat("================================================================================\n")
  cat("QUANTIFICATION SUMMARY\n")
  cat("================================================================================\n\n")

  # Get assigned and total reads per sample
  assigned_row <- which(stat$Status == "Assigned")

  if (length(assigned_row) > 0) {
    assigned <- as.numeric(stat[assigned_row, -1])
    total <- colSums(stat[, -1, drop = FALSE])
    pct_assigned <- round(100 * assigned / total, 1)

    cat("Assignment rates:\n")
    sample_names <- colnames(stat)[-1]
    for (i in seq_along(sample_names)) {
      cat(sprintf("  %-30s: %s assigned (%.1f%%)\n",
                  sample_names[i],
                  format(assigned[i], big.mark = ","),
                  pct_assigned[i]))
    }

    cat("\nOverall:\n")
    cat("  Mean assignment rate:", round(mean(pct_assigned), 1), "%\n")
    cat("  Total assigned reads:", format(sum(assigned), big.mark = ","), "\n")
  }

  cat("\n================================================================================\n")

  invisible(stat)
}
