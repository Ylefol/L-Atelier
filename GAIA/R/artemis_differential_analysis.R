#' Artemis - Differential Analysis Functions
#'
#' @description General-purpose differential analysis functions using
#' established statistical methods (DESeq2, etc.) for count-based data.


# ==============================================================================
# NORMALIZATION (for multi-comparison workflows)
# ==============================================================================

#' Normalize Count Data with DESeq2
#'
#' Creates a DESeq2 object with size factors estimated on the FULL dataset.
#' Use this when performing multiple pairwise comparisons from the same
#' experiment - normalizing once ensures consistent size factors across all
#' subsequent comparisons.
#'
#' @param counts Integer matrix of raw counts (features x samples). Rownames required.
#' @param targets Data.frame with sample metadata. Rownames must match colnames
#'   of counts. Must include a group column.
#' @param group_col Character. Column in targets for group labels. Default: "group".
#' @param batch_col Character or NULL. Optional batch column to include in design.
#'   Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_norm"} containing:
#'   \describe{
#'     \item{dds}{DESeqDataSet with size factors estimated on full dataset}
#'     \item{norm_counts}{Normalized count matrix}
#'     \item{targets}{The targets data.frame}
#'     \item{size_factors}{Named vector of size factors per sample}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Creates a DESeqDataSet from the full count matrix
#'   \item Estimates size factors using all samples together
#'   \item Stores the normalized counts and size factors
#' }
#'
#' The resulting object can be passed to \code{ARTEMIS_differential_counts()}
#' which will subset samples for each comparison while preserving the global
#' size factors.
#'
#' @section Why normalize first:
#' When performing multiple DEA comparisons from the same experiment (e.g.,
#' A vs B, A vs C, B vs C), normalizing each comparison independently can
#' introduce inconsistencies. Estimating size factors once on all samples
#' ensures that the normalized expression values are comparable across all
#' comparisons.
#'
#' @examples
#' \dontrun{
#' # Step 1: Normalize the full dataset
#' norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet)
#'
#' # Step 2: Run multiple DEA comparisons (all use same size factors)
#' dea_A_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "A")
#' dea_C_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "C")
#'
#' # Step 3: Combine results for downstream analysis
#' all_genes <- ARTEMIS_select_de_genes(list(dea_A_vs_B, dea_C_vs_B))
#'
#' }
#' @export
ARTEMIS_normalize_counts <- function(counts,
                                      targets,
                                      group_col = "group",
                                      batch_col = NULL,
                                      verbose = TRUE) {

  # --- Validate inputs ---
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  if (is.null(rownames(counts))) {
    stop("counts matrix must have rownames (feature IDs)")
  }

  # Check counts are integers
  if (!all(counts == floor(counts))) {
    stop("Counts must be integers. DESeq2 requires raw counts.")
  }
  storage.mode(counts) <- "integer"

  # Check group column
  if (!group_col %in% colnames(targets)) {
    stop("group_col '", group_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  # Align samples
  common_samples <- intersect(colnames(counts), rownames(targets))
  if (length(common_samples) == 0) {
    stop("No matching samples between counts colnames and targets rownames.\n",
         "Counts columns: ", paste(head(colnames(counts), 3), collapse = ", "), "...\n",
         "Targets rownames: ", paste(head(rownames(targets), 3), collapse = ", "), "...")
  }

  counts <- counts[, common_samples, drop = FALSE]
  targets <- targets[common_samples, , drop = FALSE]

  if (verbose) {
    cat("[ARTEMIS] Count Normalization \n")
    cat("    Samples:", ncol(counts), "\n")
    cat("    Features:", nrow(counts), "\n")
    cat("    Groups:", paste(unique(targets[[group_col]]), collapse = ", "), "\n")
    if (!is.null(batch_col)) {
      cat("    Batch correction: Yes (", batch_col, ")\n", sep = "")
    }
    cat("\n")
  }

  # --- Create DESeq2 object ---
  col_data <- data.frame(
    condition = factor(targets[[group_col]]),
    row.names = colnames(counts)
  )

  # Add batch if specified
  if (!is.null(batch_col)) {
    if (!batch_col %in% colnames(targets)) {
      stop("batch_col '", batch_col, "' not found in targets")
    }
    col_data$batch <- factor(targets[[batch_col]])
    design_formula <- ~ batch + condition
  } else {
    design_formula <- ~ condition
  }

  if (verbose) cat("[ARTEMIS] Creating DESeq2 dataset...\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design = design_formula
  )

  # --- Estimate size factors on FULL dataset ---
  if (verbose) cat("    Estimating size factors on full dataset...\n")
  dds <- estimateSizeFactors(dds)

  size_factors <- sizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)

  if (verbose) {
    cat("    Size factor range:", round(min(size_factors), 3), "-",
        round(max(size_factors), 3), "\n\n")
  }

  # --- Return artemis_norm object ---
  result <- list(
    dds = dds,
    norm_counts = norm_counts,
    targets = targets,
    size_factors = size_factors,
    parameters = list(
      group_col = group_col,
      batch_col = batch_col
    )
  )

  class(result) <- c("artemis_norm", "list")

  if (verbose) {
    cat("    Normalization complete.\n")
    cat("    Pass this object to ARTEMIS_differential_counts() for DEA.\n")
  }

  return(result)
}


#' @method print artemis_norm
#' @export
print.artemis_norm <- function(x, ...) {
  cat("Normalized Count Data\n")
  cat("------------------------------\n")
  cat("Samples:", ncol(x$norm_counts), "\n")
  cat("Features:", nrow(x$norm_counts), "\n")
  cat("Groups:", paste(unique(x$targets[[x$parameters$group_col]]), collapse = ", "), "\n")
  cat("Size factors: ", round(min(x$size_factors), 3), " - ",
      round(max(x$size_factors), 3), "\n", sep = "")
  if (!is.null(x$parameters$batch_col)) {
    cat("Batch variable:", x$parameters$batch_col, "\n")
  }
  invisible(x)
}


# ==============================================================================
# DIFFERENTIAL ANALYSIS
# ==============================================================================

#' Differential Analysis for Count Data
#'
#' @description Performs differential analysis on count data using DESeq2.
#' Works with any count-based data: RNA-seq, ATAC-seq, CUT&TAG, ChIP-seq, etc.
#'
#' @param quant_result Either:
#'   \itemize{
#'     \item A quantification result list with \code{counts} (raw integer matrix)
#'       and \code{targets} (sample metadata), optionally \code{annotation}
#'     \item An \code{artemis_norm} object from \code{ARTEMIS_normalize_counts()}.
#'       When using pre-normalized data, the function subsets to the two groups
#'       while preserving the size factors estimated on the full dataset.
#'   }
#' @param reference Character. The reference/baseline group label (e.g., "WT").
#'   This is the denominator in fold change calculations.
#' @param experiment Character. The experimental/comparison group label
#'   (e.g., "KO"). This is the numerator in fold change calculations.
#' @param group_col Character. Column name in targets containing group labels
#'   (default = "group"). Ignored if quant_result is artemis_norm (uses stored value).
#' @param batch_col Character or NULL. Optional column name for batch variable
#'   to include in design formula (default = NULL). Ignored if quant_result is
#'   artemis_norm (uses stored value).
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
#'   \item{norm_counts}{Normalized count matrix (useful for downstream analysis)}
#' }
#'
#' @details
#' This function supports two workflows:
#'
#' \strong{Single comparison:} Pass raw counts directly. Size factors are
#' estimated on the two groups being compared.
#'
#' \strong{Multiple comparisons:} First call \code{ARTEMIS_normalize_counts()}
#' on the full dataset, then pass the result to this function for each pairwise
#' comparison. This ensures consistent size factors across all comparisons.
#'
#' The function:
#' \enumerate{
#'   \item Validates that exactly two groups are present (subsets if using artemis_norm)
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
#' \dontrun{
#' # Single comparison (size factors estimated on these two groups)
#' de_result <- ARTEMIS_differential_counts(
#'   quant_result,
#'   reference = "WT",
#'   experiment = "KO"
#' )
#'
#' # Multiple comparisons from same experiment (recommended)
#' # Step 1: Normalize full dataset once
#' norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet)
#'
#' # Step 2: Run each comparison (uses pre-computed size factors)
#' dea_A_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "A")
#' dea_C_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "C")
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
#' }
ARTEMIS_differential_counts <- function(quant_result,
                                         reference,
                                         experiment,
                                         group_col = "group",
                                         batch_col = NULL,
                                         alpha = 0.05,
                                         verbose = TRUE) {

  #Check for pre-normalized data

  use_prenorm <- inherits(quant_result, "artemis_norm")

  if (use_prenorm) {
    # Using pre-normalized data - route to specialized handler
    return(.run_deseq_from_norm(
      norm_data = quant_result,
      reference = reference,
      experiment = experiment,
      alpha = alpha,
      verbose = verbose
    ))
  }

  #Standard pathway: raw counts

  if (!is.list(quant_result) || !all(c("counts", "targets") %in% names(quant_result))) {
    stop("quant_result must be either:\n",
         "  - A list with 'counts' and 'targets' elements, or\n",
         "  - An artemis_norm object from ARTEMIS_normalize_counts()")
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

  # Two-Group Validation

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
    stop("Data contains more than two groups. When running multiple comparisons\n",
         "  from the same dataset, size factors must be estimated once on all samples\n",
         "  to ensure consistency across comparisons.\n",
         "  Use ARTEMIS_normalize_counts() on the full dataset first, then pass\n",
         "  the result to ARTEMIS_differential_counts() for each pairwise comparison.\n",
         "  Specified: '", reference, "' (reference) and '", experiment, "' (experiment)\n",
         "  Additional groups found: ", paste(other_groups, collapse = ", "))
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

  # Setup DESeq2

  if (verbose) {
    cat("[ARTEMIS] Differential analysis (DESeq2):\n")
    cat("    Reference: ", reference, " (n = ", n_reference, ")\n", sep = "")
    cat("    Experiment: ", experiment, " (n = ", n_experiment, ")\n", sep = "")
    cat("    Features: ", nrow(counts), "\n", sep = "")
    if (!is.null(batch_col)) {
      cat("    Batch correction: Yes (", batch_col, ")\n", sep = "")
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
  if (verbose) cat("[ARTEMIS] Creating DESeq2 dataset...\n")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design = design_formula
  )

  # Run DESeq2 

  if (verbose) cat("[ARTEMIS] Running DESeq2...\n")

  # Suppress messages from DESeq2 unless verbose
  if (verbose) {
    dds <- DESeq(dds)
  } else {
    dds <- suppressMessages(DESeq(dds))
  }

  # Extract results
  # Contrast: experiment vs reference (so positive = higher in experiment)
  res <- results(dds, contrast = c("group", experiment, reference), alpha = alpha)

  if (verbose) cat("[ARTEMIS] Extracting results...\n\n")

  # Format Output 

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

  # Summary Statistics

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
    cat("[ARTEMIS] Results summary (FDR < ", alpha, "):\n", sep = "")
    cat("    Features tested: ", n_tested, "\n", sep = "")
    cat("    Significant: ", n_sig, " (", round(100 * n_sig / n_tested, 1), "%)\n", sep = "")
    cat("        - Increased in ", experiment, ": ", n_sig_up, "\n", sep = "")
    cat("        - Decreased in ", experiment, ": ", n_sig_down, "\n", sep = "")
    cat("\n")

    # Show top hits
    if (n_sig > 0) {
      cat("    Top significant features:\n")
      top_n <- min(5, n_sig)
      top_hits <- head(results_df[!is.na(results_df$padj) & results_df$padj < alpha, ], top_n)
      for (i in seq_len(nrow(top_hits))) {
        direction <- if (top_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
        cat("        ", top_hits$feature_id[i], ": log2FC = ",
            round(top_hits$log2FoldChange[i], 2), " (", direction, "), ",
            "padj = ", format.pval(top_hits$padj[i], digits = 2), "\n", sep = "")
      }
    }
  }

  return(list(
    results = results_df,
    dds = dds,
    summary = summary_stats,
    comparison = paste0(experiment, " vs ", reference),
    norm_counts = counts(dds, normalized = TRUE)
  ))
}


# ------------------------------------------------------------------------------
# Internal helper for pre-normalized data
# ------------------------------------------------------------------------------

.run_deseq_from_norm <- function(norm_data, reference, experiment, alpha, verbose) {

  targets <- norm_data$targets
  group_col <- norm_data$parameters$group_col
  batch_col <- norm_data$parameters$batch_col
  dds_full <- norm_data$dds

  groups <- targets[[group_col]]
  unique_groups <- unique(groups)

  # Validate groups exist

  if (!reference %in% unique_groups) {
    stop("reference '", reference, "' not found in data. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))
  }
  if (!experiment %in% unique_groups) {
    stop("experiment '", experiment, "' not found in data. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))
  }

  # Subset to the two groups of interest
  subset_mask <- groups %in% c(reference, experiment)
  subset_samples <- rownames(targets)[subset_mask]
  subset_targets <- targets[subset_samples, , drop = FALSE]
  subset_groups <- subset_targets[[group_col]]

  # Check replicates
  n_reference <- sum(subset_groups == reference)
  n_experiment <- sum(subset_groups == experiment)

  if (n_reference < 2) {
    stop("Reference group '", reference, "' has only ", n_reference, " sample(s). ",
         "DESeq2 requires at least 2 replicates per group.")
  }
  if (n_experiment < 2) {
    stop("Experiment group '", experiment, "' has only ", n_experiment, " sample(s). ",
         "DESeq2 requires at least 2 replicates per group.")
  }

  if (verbose) {
    cat("[ARTEMIS] Differential analysis (DESeq2 with pre-normalized data):\n")
    cat("    Reference: ", reference, " (n = ", n_reference, ")\n", sep = "")
    cat("    Experiment: ", experiment, " (n = ", n_experiment, ")\n", sep = "")
    cat("    Features: ", nrow(dds_full), "\n", sep = "")
    cat("    Using pre-computed size factors from full dataset\n")
    if (!is.null(batch_col)) {
      cat("    Batch correction: Yes (", batch_col, ")\n", sep = "")
    }
    cat("\n")
  }

  # Subset the DESeqDataSet
  if (verbose) cat("[ARTEMIS] Subsetting to comparison groups...\n")
  dds_subset <- dds_full[, subset_samples]

  # Update the condition factor to only have the two levels
  colData(dds_subset)$condition <- factor(
    subset_groups,
    levels = c(reference, experiment)
  )

  # Update design to use condition (simple two-group comparison)
  if (!is.null(batch_col)) {
    design(dds_subset) <- ~ batch + condition
  } else {
    design(dds_subset) <- ~ condition
  }

  # Run DESeq on subset - size factors are already set from the full dataset
  if (verbose) cat("[ARTEMIS] Running DESeq2 (size factors preserved from full dataset)...\n")

  # Estimate dispersions and run Wald test (size factors already set)
  if (verbose) {
    dds_subset <- estimateDispersions(dds_subset)
    dds_subset <- nbinomWaldTest(dds_subset)
  } else {
    dds_subset <- suppressMessages(estimateDispersions(dds_subset))
    dds_subset <- suppressMessages(nbinomWaldTest(dds_subset))
  }

  # Extract results
  res <- results(dds_subset, contrast = c("condition", experiment, reference), alpha = alpha)

  if (verbose) cat("[ARTEMIS] Extracting results...\n\n")

  # Format output
  results_df <- as.data.frame(res)

  # Add feature identifiers
  results_df$feature_id <- rownames(results_df)

  # Reorder columns
  id_cols <- "feature_id"
  stat_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  results_df <- results_df[, c(id_cols, stat_cols)]

  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$padj), ]
  rownames(results_df) <- NULL

  # Summary statistics
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
    cat("[ARTEMIS] Results summary (FDR < ", alpha, "):\n", sep = "")
    cat("    Features tested: ", n_tested, "\n", sep = "")
    cat("    Significant: ", n_sig, " (", round(100 * n_sig / n_tested, 1), "%)\n", sep = "")
    cat("        - Increased in ", experiment, ": ", n_sig_up, "\n", sep = "")
    cat("        - Decreased in ", experiment, ": ", n_sig_down, "\n", sep = "")
    cat("\n")

    # Show top hits
    if (n_sig > 0) {
      cat("    Top significant features:\n")
      top_n <- min(5, n_sig)
      top_hits <- head(results_df[!is.na(results_df$padj) & results_df$padj < alpha, ], top_n)
      for (i in seq_len(nrow(top_hits))) {
        direction <- if (top_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
        cat("        ", top_hits$feature_id[i], ": log2FC = ",
            round(top_hits$log2FoldChange[i], 2), " (", direction, "), ",
            "padj = ", format.pval(top_hits$padj[i], digits = 2), "\n", sep = "")
      }
    }
  }

  return(list(
    results = results_df,
    dds = dds_subset,
    summary = summary_stats,
    comparison = paste0(experiment, " vs ", reference),
    norm_counts = counts(dds_subset, normalized = TRUE)
  ))
}


# ==============================================================================
# GENE SELECTION
# ==============================================================================

#' Select DE Genes Across Differential Expression Results
#'
#' Collects the union (or intersection) of significant genes across one or more
#' differential expression comparisons. Works with time series DEA results,
#' standard DEA results from \code{ARTEMIS_differential_counts()}, or raw
#' data.frames containing DE statistics.
#'
#' @param de_results DE results in any of the following formats:
#'   \itemize{
#'     \item \code{artemis_ts_de} object (from time series analysis)
#'     \item DEA result from \code{ARTEMIS_differential_counts()}
#'     \item A data.frame with gene IDs, log2FoldChange, and p-values
#'     \item A list of any combination of the above
#'   }
#' @param l2fc_thresh Numeric. Absolute log2 fold-change threshold. Only genes
#'   with |log2FC| >= this value are selected. Default: 1.0.
#' @param p_col Character. P-value column to use: "padj" or "pvalue".
#'   Default: "padj".
#' @param p_thresh Numeric. P-value threshold. Default: 0.05.
#' @param gene_col Character or NULL. Column name containing gene IDs. If NULL,
#'   auto-detects from common names (gene_id, feature_id, gene, etc.). Default: NULL.
#' @param union Logical. If TRUE, take union of genes across comparisons.
#'   If FALSE, take intersection. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A character vector of selected gene IDs with attributes:
#'   \describe{
#'     \item{gene_summary}{Data.frame showing which comparisons each gene
#'       was significant in}
#'     \item{per_comparison}{Named list of gene sets per comparison}
#'     \item{parameters}{Selection parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' # From time series DEA
#' genes <- ARTEMIS_select_de_genes(list(cond_de, temp_de), l2fc_thresh = 2)
#'
#' # From standard DEA (ARTEMIS_differential_counts)
#' genes <- ARTEMIS_select_de_genes(dea_result, l2fc_thresh = 1)
#'
#' # From multiple DEA results
#' genes <- ARTEMIS_select_de_genes(
#'   list(dea_result1, dea_result2, dea_result3),
#'   l2fc_thresh = 1,
#'   union = TRUE
#' )
#'
#' }
#' @export
ARTEMIS_select_de_genes <- function(de_results,
                                     l2fc_thresh = 1.0,
                                     p_col = "padj",
                                     p_thresh = 0.05,
                                     gene_col = NULL,
                                     union = TRUE,
                                     verbose = TRUE) {

  # ==========================================================================
  # Normalize inputs to a named list of data.frames
  # Supports: artemis_ts_de, DEA result (from differential_counts), data.frame,
  #           or lists of any of these
  # ==========================================================================

  .extract_results_df <- function(x, name = NULL) {
    # Returns a named list of data.frames

    if (inherits(x, "artemis_ts_de")) {
      # Time series DEA: nested structure $results[[comparison]]$results
      out <- list()
      for (exp_name in names(x$results)) {
        out[[exp_name]] <- x$results[[exp_name]]$results
      }
      return(out)

    } else if (is.data.frame(x)) {
      # Raw data.frame - use provided name or generate one
      nm <- if (!is.null(name)) name else "comparison_1"
      return(setNames(list(x), nm))

    } else if (is.list(x) && "results" %in% names(x) && is.data.frame(x$results)) {
      # Single DEA result from ARTEMIS_differential_counts
      # Structure: $results (data.frame), $dds, $summary
      nm <- if (!is.null(x$comparison)) x$comparison else
            if (!is.null(name)) name else "comparison_1"
      return(setNames(list(x$results), nm))

    } else {
      return(NULL)
    }
  }

  .detect_gene_col <- function(df) {
    # Auto-detect gene ID column
    candidates <- c("gene_id", "feature_id", "gene", "ensembl_id", "symbol")
    for (col in candidates) {
      if (col %in% colnames(df)) return(col)
    }
    # Fall back to first column if it looks like gene IDs
    if (ncol(df) > 0 && is.character(df[[1]])) return(colnames(df)[1])
    return(NULL)
  }

  # Normalize input
  all_dfs <- list()

  if (is.list(de_results) && !is.data.frame(de_results)) {
    # Could be: list of objects, artemis_ts_de, or single DEA result

    # Check if it's a single DEA result (has $results data.frame)
    extracted <- .extract_results_df(de_results)
    if (!is.null(extracted)) {
      all_dfs <- c(all_dfs, extracted)
    } else {
      # It's a list of objects - process each

      for (i in seq_along(de_results)) {
        item <- de_results[[i]]
        nm <- names(de_results)[i]
        if (is.null(nm) || nm == "") nm <- paste0("comparison_", i)
        extracted <- .extract_results_df(item, nm)
        if (!is.null(extracted)) {
          all_dfs <- c(all_dfs, extracted)
        } else {
          warning("Skipping unrecognized element at position ", i)
        }
      }
    }
  } else if (is.data.frame(de_results)) {
    all_dfs <- .extract_results_df(de_results)
  } else {
    stop("'de_results' must be a data.frame, DEA result, artemis_ts_de object, or list of these")
  }

  if (length(all_dfs) == 0) {
    stop("No valid DE results found in input")
  }

  if (verbose) {
    cat("[ARTEMIS] Selecting DE Genes \n")
    cat("    Comparisons found:", length(all_dfs), "\n")
    cat("    Thresholds: |log2FC| >=", l2fc_thresh, ",", p_col, "<=", p_thresh, "\n")
  }

  # ==========================================================================
  # Process each comparison
  # ==========================================================================

  per_comparison <- list()

  for (exp_name in names(all_dfs)) {
    res_df <- all_dfs[[exp_name]]

    # Detect gene column
    g_col <- gene_col
    if (is.null(g_col)) {
      g_col <- .detect_gene_col(res_df)
    }
    if (is.null(g_col) || !g_col %in% colnames(res_df)) {
      warning("Could not find gene ID column in results for ", exp_name, ". Skipping.")
      next
    }

    if (!p_col %in% colnames(res_df)) {
      warning("Column '", p_col, "' not found in results for ", exp_name, ". Skipping.")
      next
    }

    if (!"log2FoldChange" %in% colnames(res_df)) {
      warning("Column 'log2FoldChange' not found in results for ", exp_name, ". Skipping.")
      next
    }

    # Filter by p-value and L2FC
    sig_mask <- !is.na(res_df[[p_col]]) &
      res_df[[p_col]] <= p_thresh &
      !is.na(res_df$log2FoldChange) &
      abs(res_df$log2FoldChange) >= l2fc_thresh

    sig_genes <- res_df[[g_col]][sig_mask]
    per_comparison[[exp_name]] <- sig_genes

    if (verbose) {
      cat("    ", exp_name, ":", length(sig_genes), "genes\n")
    }
  }

  if (length(per_comparison) == 0) {
    warning("No significant genes found in any comparison")
    result <- character(0)
    attr(result, "per_comparison") <- per_comparison
    attr(result, "parameters") <- list(
      l2fc_thresh = l2fc_thresh, p_col = p_col,
      p_thresh = p_thresh, union = union
    )
    return(result)
  }

  # ==========================================================================
  # Combine: union or intersection
  # ==========================================================================

  if (union) {
    selected <- unique(unlist(per_comparison))
  } else {
    selected <- Reduce(intersect, per_comparison)
  }

  if (verbose) {
    cat("\n", if (union) "Union" else "Intersection", ":", length(selected), "genes\n")
  }

  # Build gene summary: which comparisons each gene appeared in
  all_genes <- unique(unlist(per_comparison))
  gene_presence <- sapply(per_comparison, function(g) all_genes %in% g)
  if (!is.matrix(gene_presence)) {
    gene_presence <- matrix(gene_presence, ncol = 1,
                             dimnames = list(all_genes, names(per_comparison)))
  }
  gene_summary <- data.frame(
    gene = all_genes,
    n_comparisons = rowSums(gene_presence),
    stringsAsFactors = FALSE
  )
  gene_summary <- gene_summary[order(-gene_summary$n_comparisons), ]

  # Attach metadata as attributes
  attr(selected, "gene_summary") <- gene_summary
  attr(selected, "per_comparison") <- per_comparison
  attr(selected, "parameters") <- list(
    l2fc_thresh = l2fc_thresh, p_col = p_col,
    p_thresh = p_thresh, union = union
  )

  return(selected)
}


# ==============================================================================
# MATRIX PREPARATION FOR CLUSTERING
# ==============================================================================

#' Prepare Expression Matrix for PART Clustering
#'
#' Extracts normalized counts from DEA results and subsets to selected genes,
#' optionally z-score scaling each gene. This bridges DEA results to PART input.
#'
#' @param de_result DEA result in one of the following formats:
#'   \itemize{
#'     \item DEA result from \code{ARTEMIS_differential_counts()} (extracts
#'       normalized counts from the DESeqDataSet)
#'     \item \code{artemis_ts_de} object (extracts norm_counts)
#'     \item A numeric matrix of normalized expression (genes x samples)
#'   }
#' @param genes Character vector of gene IDs to include (e.g., from
#'   \code{ARTEMIS_select_de_genes()}).
#' @param samples Character vector of sample IDs (column names) to include,
#'   in the desired order. NULL = all samples. Default: NULL.
#' @param scale Logical. Z-score scale each gene (row) across samples.
#'   Recommended for PART when using euclidean distance. Default: TRUE.
#' @param log_transform Logical. Apply log2(x + 1) transformation before scaling.
#'   Default: FALSE (matches original TiSA behavior). Set TRUE if you want
#'   log-scale distances for normalized counts.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Numeric matrix ready for \code{ARTEMIS_part()}.
#'
#' @examples
#' \dontrun{
#' # From standard DEA result (default: no log transform, matches TiSA)
#' dea <- ARTEMIS_differential_counts(quant, "WT", "KO")
#' genes <- ARTEMIS_select_de_genes(dea, l2fc_thresh = 1)
#' part_mat <- ARTEMIS_prepare_part_matrix(dea, genes)
#' part_result <- ARTEMIS_part(part_mat, seed = 42)
#'
#' # From time series DEA
#' genes <- ARTEMIS_select_de_genes(cond_de, l2fc_thresh = 2)
#' part_mat <- ARTEMIS_prepare_part_matrix(cond_de, genes)
#'
#' # With log transformation (optional, if you prefer log-scale distances)
#' part_mat <- ARTEMIS_prepare_part_matrix(dea, genes, log_transform = TRUE)
#'
#' }
#' @export
ARTEMIS_prepare_part_matrix <- function(de_result,
                                         genes,
                                         samples = NULL,
                                         scale = TRUE,
                                         log_transform = FALSE,
                                         verbose = TRUE) {

  # ==========================================================================
  # Extract normalized counts based on input type
  # ==========================================================================

  counts <- NULL
  source_type <- NULL

  # Normalized data from ARTEMIS_normalize_counts
  if (inherits(de_result, "artemis_norm")) {
    counts <- de_result$norm_counts
    source_type <- "artemis_norm"

  # DEA result from ARTEMIS_differential_counts (has $dds or $norm_counts)
  } else if (is.list(de_result) && "dds" %in% names(de_result) &&
      inherits(de_result$dds, "DESeqDataSet")) {
    # Prefer $norm_counts if available, otherwise extract from dds
    if (!is.null(de_result$norm_counts)) {
      counts <- de_result$norm_counts
    } else {
      counts <- DESeq2::counts(de_result$dds, normalized = TRUE)
    }
    source_type <- "DEA result (DESeqDataSet)"

  # Time series DEA result (has $norm_counts)
  } else if (inherits(de_result, "artemis_ts_de")) {
    if (is.null(de_result$norm_counts)) {
      stop("artemis_ts_de object does not contain norm_counts. ",
           "This may be from an older version. Please re-run normalization.")
    }
    counts <- de_result$norm_counts
    source_type <- "artemis_ts_de"

  # Raw matrix or data.frame
  } else if (is.matrix(de_result) || is.data.frame(de_result)) {
    counts <- as.matrix(de_result)
    source_type <- "matrix"
    # Assume raw matrix is already appropriately transformed
    if (log_transform && verbose) {
      cat("[ARTEMIS] Note: log_transform=TRUE but input is raw matrix. ",
          "Set log_transform=FALSE if already transformed.\n")
    }

  } else {
    stop("'de_result' must be a DEA result (from ARTEMIS_differential_counts), ",
         "an artemis_norm object, an artemis_ts_de object, or a numeric matrix")
  }

  if (verbose) {
    cat("[ARTEMIS] Extracting counts from:", source_type, "\n")
  }

  # ==========================================================================
  # Subset to available genes
  # ==========================================================================

  available <- genes[genes %in% rownames(counts)]
  missing <- length(genes) - length(available)

  if (length(available) == 0) {
    stop("None of the specified genes found in count matrix rownames.\n",
         "First few gene IDs requested: ", paste(head(genes, 3), collapse = ", "), "\n",
         "First few rownames in counts: ", paste(head(rownames(counts), 3), collapse = ", "))
  }

  mat <- counts[available, , drop = FALSE]

  # ==========================================================================
  # Subset samples
  # ==========================================================================

  if (!is.null(samples)) {
    avail_samps <- samples[samples %in% colnames(mat)]
    if (length(avail_samps) == 0) {
      stop("None of the specified samples found in count matrix colnames")
    }
    mat <- mat[, avail_samps, drop = FALSE]
  }

  # ==========================================================================
  # Log transform (if requested and source is normalized counts)
  # ==========================================================================

  if (log_transform && source_type != "matrix") {
    mat <- log2(mat + 1)
    if (verbose) cat("[ARTEMIS] Applied log2(x + 1) transformation\n")
  }

  # ==========================================================================
  # Z-score scale
  # ==========================================================================

  if (scale) {
    mat <- t(scale(t(mat), center = TRUE, scale = TRUE))
    # Remove constant rows (NaN from zero variance)
    nan_rows <- rowSums(is.nan(mat)) > 0
    if (any(nan_rows)) {
      mat <- mat[!nan_rows, , drop = FALSE]
      if (verbose) cat("[ARTEMIS] Removed", sum(nan_rows), "constant genes after scaling\n")
    }
  }

  # ==========================================================================
  # Summary
  # ==========================================================================

  if (verbose) {
    cat("[ARTEMIS] Prepared matrix:", nrow(mat), "genes x", ncol(mat), "samples")
    if (missing > 0) cat(" (", missing, " genes not found)", sep = "")
    if (scale) cat(" [z-scored]")
    cat("\n")
  }

  return(mat)
}
