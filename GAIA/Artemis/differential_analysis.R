#' Artemis - Differential Analysis Functions
#'
#' @description General-purpose differential analysis functions using
#' established statistical methods (DESeq2, etc.) for count-based data.

library(DESeq2)


#' Differential Analysis for Count Data
#'
#' @description Performs differential analysis on count data using DESeq2.
#' Works with any count-based data: RNA-seq, ATAC-seq, CUT&TAG, ChIP-seq, etc.
#'
#' @param quant_result A quantification result containing:
#'   \itemize{
#'     \item counts: Matrix of raw integer counts (features x samples)
#'     \item targets: Data.frame with sample metadata including group column
#'     \item annotation: (Optional) Data.frame with feature annotations
#'   }
#' @param reference Character. The reference/baseline group label (e.g., "WT").
#'   This is the denominator in fold change calculations.
#' @param experiment Character. The experimental/comparison group label
#'   (e.g., "KO"). This is the numerator in fold change calculations.
#' @param group_col Character. Column name in targets containing group labels
#'   (default = "group").
#' @param batch_col Character or NULL. Optional column name for batch variable
#'   to include in design formula (default = NULL).
#' @param alpha Numeric. FDR threshold for summary statistics (default = 0.05).
#' @param verbose Logical. Print progress and summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{Data.frame with per-feature DESeq2 results: feature_id,
#'     (annotation columns if available), baseMean, log2FoldChange, lfcSE,
#'     stat, pvalue, padj}
#'   \item{dds}{The DESeqDataSet object for further analysis}
#'   \item{summary}{List with summary statistics: n_features, n_sig_up,
#'     n_sig_down, n_tested, alpha}
#'   \item{comparison}{Character describing the comparison (experiment vs reference)}
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Validates that exactly two groups are present
#'   \item Validates that both groups have at least 2 replicates
#'   \item Creates DESeq2 dataset with design ~ group (or ~ batch + group)
#'   \item Runs DESeq2 with default settings
#'   \item Returns per-feature statistics with multiple testing correction
#' }
#'
#' The reference level is set so that positive log2FC means higher in experiment
#' compared to reference.
#'
#' @section Data Types:
#' This function is designed to work with any count-based quantification:
#' \itemize{
#'   \item RNA-seq gene/transcript counts
#'   \item ATAC-seq counts at peaks/regions
#'   \item CUT&TAG/ChIP-seq counts at binding sites
#'   \item Any other integer count matrix
#' }
#'
#' @export
#'
#' @examples
#' # RNA-seq differential expression
#' de_result <- ARTEMIS_differential_counts(
#'   rna_counts,
#'   reference = "WT",
#'   experiment = "KO"
#' )
#'
#' # ATAC-seq differential accessibility at regions
#' da_result <- ARTEMIS_differential_counts(
#'   atac_at_peaks,
#'   reference = "Control",
#'   experiment = "Treatment"
#' )
#'
#' # With batch correction
#' result <- ARTEMIS_differential_counts(
#'   counts_data,
#'   reference = "WT",
#'   experiment = "KO",
#'   batch_col = "batch"
#' )
#'
#' # View significant features
#' sig_features <- subset(result$results, padj < 0.05)
#'
ARTEMIS_differential_counts <- function(quant_result,
                                         reference,
                                         experiment,
                                         group_col = "group",
                                         batch_col = NULL,
                                         alpha = 0.05,
                                         verbose = TRUE) {

  # ===== Input Validation =====

  if (!is.list(quant_result) || !all(c("counts", "targets") %in% names(quant_result))) {
    stop("quant_result must be a list with 'counts' and 'targets' elements")
  }

  counts <- quant_result$counts
  targets <- quant_result$targets
  annotation <- quant_result$annotation

  # Check counts are integers

  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  if (!all(counts == floor(counts))) {
    stop("Counts must be integers. DESeq2 requires raw counts, not normalized or ",
         "batch-corrected values. If your data has been processed, use the original ",
         "raw counts instead.")
  }
  # Ensure integer storage mode
  storage.mode(counts) <- "integer"

  # Check group column exists
  if (!group_col %in% colnames(targets)) {
    stop("group_col '", group_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  groups <- targets[[group_col]]
  unique_groups <- unique(groups)

  # ===== Two-Group Validation =====

  # Check reference and experiment exist
  if (!reference %in% unique_groups) {
    stop("reference '", reference, "' not found in data. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))
  }
  if (!experiment %in% unique_groups) {
    stop("experiment '", experiment, "' not found in data. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))
  }

  # Check exactly two groups (only reference and experiment)
  other_groups <- setdiff(unique_groups, c(reference, experiment))
  if (length(other_groups) > 0) {
    stop("Data contains more than two groups. This function requires exactly two groups.\n",
         "  Specified: '", reference, "' (reference) and '", experiment, "' (experiment)\n",
         "  Additional groups found: ", paste(other_groups, collapse = ", "), "\n",
         "  Please subset your data to include only the two groups of interest.")
  }

  # Check replicates
  n_reference <- sum(groups == reference)
  n_experiment <- sum(groups == experiment)

  if (n_reference < 2) {
    stop("Reference group '", reference, "' has only ", n_reference, " sample(s). ",
         "DESeq2 requires at least 2 replicates per group.")
  }
  if (n_experiment < 2) {
    stop("Experiment group '", experiment, "' has only ", n_experiment, " sample(s). ",
         "DESeq2 requires at least 2 replicates per group.")
  }

  # Check batch column if provided
  if (!is.null(batch_col)) {
    if (!batch_col %in% colnames(targets)) {
      stop("batch_col '", batch_col, "' not found in targets. ",
           "Available columns: ", paste(colnames(targets), collapse = ", "))
    }
  }

  # ===== Setup DESeq2 =====

  if (verbose) {
    cat("Differential analysis (DESeq2):\n")
    cat("  Reference: ", reference, " (n = ", n_reference, ")\n", sep = "")
    cat("  Experiment: ", experiment, " (n = ", n_experiment, ")\n", sep = "")
    cat("  Features: ", nrow(counts), "\n", sep = "")
    if (!is.null(batch_col)) {
      cat("  Batch correction: Yes (", batch_col, ")\n", sep = "")
    }
    cat("\n")
  }

  # Create colData for DESeq2
  col_data <- data.frame(
    group = factor(groups, levels = c(reference, experiment)),
    row.names = colnames(counts)
  )

  # Add batch if specified
  if (!is.null(batch_col)) {
    col_data$batch <- factor(targets[[batch_col]])
    design_formula <- ~ batch + group
  } else {
    design_formula <- ~ group
  }

  # Create DESeqDataSet
  if (verbose) cat("Creating DESeq2 dataset...\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design = design_formula
  )

  # ===== Run DESeq2 =====

  if (verbose) cat("Running DESeq2...\n")

  # Suppress messages from DESeq2 unless verbose
  if (verbose) {
    dds <- DESeq(dds)
  } else {
    dds <- suppressMessages(DESeq(dds))
  }

  # Extract results
  # Contrast: experiment vs reference (so positive = higher in experiment)
  res <- results(dds, contrast = c("group", experiment, reference), alpha = alpha)

  if (verbose) cat("Extracting results...\n\n")

  # ===== Format Output =====

  # Create results data.frame
  results_df <- as.data.frame(res)

  # Add feature identifiers
  if (!is.null(rownames(counts))) {
    results_df$feature_id <- rownames(counts)
  } else if (!is.null(annotation) && "peak_id" %in% colnames(annotation)) {
    results_df$feature_id <- annotation$peak_id
  } else if (!is.null(annotation) && "gene_id" %in% colnames(annotation)) {
    results_df$feature_id <- annotation$gene_id
  } else {
    results_df$feature_id <- paste0("feature_", seq_len(nrow(results_df)))
  }

  # Add genomic coordinates if available (for peak/region data)
  if (!is.null(annotation)) {
    if ("chr" %in% colnames(annotation)) results_df$chr <- annotation$chr
    if ("start" %in% colnames(annotation)) results_df$start <- annotation$start
    if ("end" %in% colnames(annotation)) results_df$end <- annotation$end
    # Add gene symbol if available (for RNA-seq)
    if ("gene_name" %in% colnames(annotation)) results_df$gene_name <- annotation$gene_name
  }

  # Reorder columns for clarity
  id_cols <- intersect(c("feature_id", "gene_name", "chr", "start", "end"), colnames(results_df))
  stat_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  results_df <- results_df[, c(id_cols, stat_cols)]

  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$padj), ]
  rownames(results_df) <- NULL

  # ===== Summary Statistics =====

  n_tested <- sum(!is.na(results_df$padj))
  n_sig <- sum(results_df$padj < alpha, na.rm = TRUE)
  n_sig_up <- sum(results_df$padj < alpha & results_df$log2FoldChange > 0, na.rm = TRUE)
  n_sig_down <- sum(results_df$padj < alpha & results_df$log2FoldChange < 0, na.rm = TRUE)

  summary_stats <- list(
    n_features = nrow(results_df),
    n_tested = n_tested,
    n_significant = n_sig,
    n_sig_up = n_sig_up,
    n_sig_down = n_sig_down,
    alpha = alpha
  )

  if (verbose) {
    cat("Results summary (FDR < ", alpha, "):\n", sep = "")
    cat("  Features tested: ", n_tested, "\n", sep = "")
    cat("  Significant: ", n_sig, " (", round(100 * n_sig / n_tested, 1), "%)\n", sep = "")
    cat("    - Increased in ", experiment, ": ", n_sig_up, "\n", sep = "")
    cat("    - Decreased in ", experiment, ": ", n_sig_down, "\n", sep = "")
    cat("\n")

    # Show top hits
    if (n_sig > 0) {
      cat("Top significant features:\n")
      top_n <- min(5, n_sig)
      top_hits <- head(results_df[!is.na(results_df$padj) & results_df$padj < alpha, ], top_n)
      for (i in seq_len(nrow(top_hits))) {
        direction <- if (top_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
        cat("  ", top_hits$feature_id[i], ": log2FC = ",
            round(top_hits$log2FoldChange[i], 2), " (", direction, "), ",
            "padj = ", format.pval(top_hits$padj[i], digits = 2), "\n", sep = "")
      }
    }
  }

  return(list(
    results = results_df,
    dds = dds,
    summary = summary_stats,
    comparison = paste0(experiment, " vs ", reference)
  ))
}
