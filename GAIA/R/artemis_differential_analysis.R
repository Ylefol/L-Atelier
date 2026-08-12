#' Artemis - Differential Analysis Functions
#'
#' @description General-purpose differential analysis functions using
#' established statistical methods (DESeq2, etc.) for count-based data.


# ==============================================================================
# NORMALIZATION (for multi-comparison workflows)
# ==============================================================================

# Resolve a user-supplied `size_factors` argument (NULL is handled by the
# caller before this is invoked) into a full, sample-ordered numeric vector.
# A single scalar recycles across every sample; a named vector must cover
# every sample_id exactly.
.artemis_resolve_size_factors <- function(size_factors, sample_ids) {
  if (!is.numeric(size_factors)) {
    stop("size_factors must be NULL, a named numeric vector, or a single ",
         "numeric scalar.")
  }
  if (any(size_factors <= 0)) {
    stop("size_factors must be strictly positive.")
  }

  if (length(size_factors) == 1L) {
    return(setNames(rep(size_factors, length(sample_ids)), sample_ids))
  }

  if (is.null(names(size_factors))) {
    stop("size_factors must be named (names matching sample IDs) when ",
         "supplying more than one value. To apply the same factor to every ",
         "sample, pass a single scalar instead.")
  }

  missing_samples <- setdiff(sample_ids, names(size_factors))
  if (length(missing_samples) > 0) {
    stop("size_factors is missing values for: ",
         paste(missing_samples, collapse = ", "))
  }

  size_factors[sample_ids]
}


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
#' @param batch_col Character or NULL. Optional column name for a single blocking
#'   variable to include in the design formula as \code{~ batch_col + condition}.
#'   Default: NULL. \strong{Limitation:} only one blocking variable is supported;
#'   multi-variable designs are not yet implemented.
#' @param size_factors NULL, a named numeric vector, or a single numeric scalar.
#'   Controls how per-sample size factors are set. Default: NULL.
#'   \describe{
#'     \item{NULL}{DESeq2 estimates size factors itself via
#'       \code{estimateSizeFactors()} (median-of-ratios). This is correct when
#'       the input counts are raw and untouched by any prior normalization.}
#'     \item{Named numeric vector}{Names must match \code{rownames(targets)}
#'       (equivalently \code{colnames(counts)}); every sample must have a
#'       value. Assigned directly via \code{sizeFactors(dds) <-}, skipping
#'       DESeq2's own estimation entirely. Use this to hand DESeq2 an
#'       externally computed scale, e.g. spike-in-derived factors from
#'       \code{HORIZON_compute_spike_in_factors()}'s \code{size_factor_deseq2}
#'       column.}
#'     \item{Single scalar}{Recycled across every sample (e.g. \code{1}). Use
#'       this to disable DESeq2-level normalization entirely -- appropriate
#'       when the input counts have already been put on a comparable scale
#'       upstream (e.g. by physical read downsampling) and no further
#'       per-sample correction should be applied.}
#'   }
#'   All values must be strictly positive (DESeq2's own requirement).
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_norm"} containing:
#'   \describe{
#'     \item{dds}{DESeqDataSet with size factors set on full dataset}
#'     \item{norm_counts}{Normalized count matrix}
#'     \item{targets}{The targets data.frame}
#'     \item{size_factors}{Named vector of size factors per sample}
#'     \item{parameters}{List of parameters used, including
#'       \code{size_factors_source} ("deseq2_estimated" or "external")}
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Creates a DESeqDataSet from the full count matrix
#'   \item Sets size factors using all samples together -- either estimated by
#'     DESeq2 or supplied via \code{size_factors}
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
#' @section Externally supplied size factors:
#' DESeq2's default (median-of-ratios) normalization assumes most features are
#' non-differential between samples and infers relative library scale from the
#' count matrix itself. When an experiment includes a spike-in control (ATAC,
#' CHIP/CUT&RUN, R-loop), that assumption can be exactly what you don't want --
#' a real, global shift in signal violates it. \code{size_factors} lets you
#' substitute a ground-truth scale (from the spike-in) in place of DESeq2's own
#' estimate. Do not combine externally supplied size factors with counts that
#' have \emph{already} been normalized upstream (e.g. via physical BAM
#' downsampling to a spike-in ratio) -- that applies the same correction twice.
#' Pick one: raw counts + \code{size_factors} (spike-in derived), or
#' already-normalized counts + \code{size_factors = 1} (disables further
#' correction).
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
#' # Spike-in normalized ATAC/CHIP: supply factors instead of letting DESeq2
#' # estimate its own (raw, non-downsampled counts required for this to be a
#' # single, correct normalization rather than a double one).
#' spikein_factors <- HORIZON_compute_spike_in_factors(spikein_bams)
#' sf <- setNames(spikein_factors$size_factor_deseq2, spikein_factors$sample_id)
#' norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, size_factors = sf)
#'
#' # Counts already normalized upstream -- disable DESeq2-level normalization
#' norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, size_factors = 1)
#' }
#' @export
ARTEMIS_normalize_counts <- function(counts,
                                      targets,
                                      group_col = "group",
                                      batch_col = NULL,
                                      size_factors = NULL,
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

  # --- Set size factors on FULL dataset ---
  if (is.null(size_factors)) {
    if (verbose) cat("    Estimating size factors on full dataset (DESeq2 median-of-ratios)...\n")
    dds <- estimateSizeFactors(dds)
    size_factors_source <- "deseq2_estimated"
  } else {
    resolved_sf <- .artemis_resolve_size_factors(size_factors, colnames(counts))
    if (verbose) cat("    Using externally supplied size factors (skipping DESeq2 estimation)...\n")
    sizeFactors(dds) <- resolved_sf
    size_factors_source <- "external"
  }

  size_factors_final <- sizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)

  if (verbose) {
    cat("    Size factor range:", round(min(size_factors_final), 3), "-",
        round(max(size_factors_final), 3), "\n\n")
  }

  # --- Return artemis_norm object ---
  result <- list(
    dds = dds,
    norm_counts = norm_counts,
    targets = targets,
    size_factors = size_factors_final,
    parameters = list(
      group_col = group_col,
      batch_col = batch_col,
      size_factors_source = size_factors_source
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
  sf_source <- x$parameters$size_factors_source
  cat("Size factor source: ",
      if (is.null(sf_source)) "unknown (older artemis_norm object)"
      else switch(sf_source,
                  deseq2_estimated = "DESeq2 (median-of-ratios)",
                  external          = "externally supplied",
                  "unknown"),
      "\n", sep = "")
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
#' Requires an \code{artemis_norm} object from \code{ARTEMIS_normalize_counts()}
#' as input, ensuring size factors and batch correction are set consistently
#' across all comparisons before any DEA is run.
#'
#' @param norm_data An \code{artemis_norm} object from
#'   \code{ARTEMIS_normalize_counts()}. Size factors and batch correction
#'   settings are inherited from this object.
#' @param reference Character. The reference/baseline group label (e.g., "WT").
#'   This is the denominator in fold change calculations.
#' @param experiment Character. The experimental/comparison group label
#'   (e.g., "KO"). This is the numerator in fold change calculations.
#' @param alpha Numeric. FDR threshold for summary statistics (default = 0.05).
#' @param verbose Logical. Print progress and summary (default = TRUE).
#'
#' @return An \code{artemis_dea} object (list) containing:
#' \describe{
#'   \item{results}{Data.frame with per-feature DESeq2 results: feature_id,
#'     baseMean, log2FoldChange, lfcSE, stat, pvalue, padj}
#'   \item{dds}{The DESeqDataSet object for further analysis}
#'   \item{summary}{List with summary statistics: n_features, n_sig_up,
#'     n_sig_down, n_tested, alpha}
#'   \item{comparison}{Character describing the comparison (experiment vs reference)}
#'   \item{norm_counts}{Normalized count matrix (useful for downstream analysis)}
#' }
#' Calling \code{print()} on the returned object reprints the comparison and
#' significance summary (and top 5 hits) without rerunning DESeq2 -- see
#' \code{\link{print.artemis_dea}}.
#'
#' @details
#' Always call \code{ARTEMIS_normalize_counts()} on the full dataset first, then
#' pass the result here for each pairwise comparison. This ensures:
#' \itemize{
#'   \item Size factors are estimated on all samples, not just the two being compared
#'   \item Batch correction is configured in one place and applied consistently
#'   \item Normalized counts are available for visualization before running DEA
#' }
#'
#' The function subsets to the two groups of interest, preserves the pre-computed
#' size factors, estimates dispersions, and runs the Wald test. Positive log2FC
#' means higher in experiment relative to reference.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Step 1: Normalize full dataset (with optional batch correction)
#' norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, batch_col = "batch")
#'
#' # Step 2: Run each pairwise comparison
#' dea_KO_vs_WT  <- ARTEMIS_differential_counts(norm_data, "WT", "KO")
#' dea_OE_vs_WT  <- ARTEMIS_differential_counts(norm_data, "WT", "OE")
#'
#' # View significant features
#' sig_features <- subset(dea_KO_vs_WT$results, padj < 0.05)
#' }
ARTEMIS_differential_counts <- function(norm_data,
                                         reference,
                                         experiment,
                                         alpha = 0.05,
                                         verbose = TRUE) {

  if (!inherits(norm_data, "artemis_norm")) {
    stop("norm_data must be an artemis_norm object from ARTEMIS_normalize_counts().\n",
         "  Call ARTEMIS_normalize_counts() on your full dataset first, then pass\n",
         "  the result here for each pairwise comparison.")
  }

  return(.run_deseq_from_norm(
    norm_data = norm_data,
    reference = reference,
    experiment = experiment,
    alpha = alpha,
    verbose = verbose
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

  result <- list(
    results = results_df,
    dds = dds_subset,
    summary = summary_stats,
    comparison = paste0(experiment, " vs ", reference),
    norm_counts = counts(dds_subset, normalized = TRUE)
  )
  class(result) <- c("artemis_dea", "list")
  result
}


#' Print method for ARTEMIS differential analysis results
#'
#' Reprints the comparison and significance summary (and top 5 hits) from an
#' \code{artemis_dea} object without rerunning DESeq2 -- the same summary
#' \code{\link{ARTEMIS_differential_counts}} prints live when \code{verbose =
#' TRUE}, but callable on the already-computed result.
#'
#' @param x An \code{artemis_dea} object from
#'   \code{\link{ARTEMIS_differential_counts}}.
#' @param ... Additional arguments (ignored)
#'
#' @method print artemis_dea
#' @export
print.artemis_dea <- function(x, ...) {
  cat("Differential Analysis Results\n")
  cat("------------------------------\n")
  cat("Comparison: ", x$comparison, "\n", sep = "")
  cat("Features tested: ", x$summary$n_tested, "\n", sep = "")
  cat("Significant (FDR < ", x$summary$alpha, "): ", x$summary$n_significant,
      " (", round(100 * x$summary$n_significant / x$summary$n_tested, 1), "%)\n", sep = "")
  cat("    - Increased in ", sub(" vs .*", "", x$comparison), ": ", x$summary$n_sig_up, "\n", sep = "")
  cat("    - Decreased in ", sub(" vs .*", "", x$comparison), ": ", x$summary$n_sig_down, "\n", sep = "")

  if (x$summary$n_significant > 0) {
    cat("\nTop significant features:\n")
    top_n <- min(5, x$summary$n_significant)
    top_hits <- head(x$results[!is.na(x$results$padj) & x$results$padj < x$summary$alpha, ], top_n)
    for (i in seq_len(nrow(top_hits))) {
      direction <- if (top_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
      cat("    ", top_hits$feature_id[i], ": log2FC = ",
          round(top_hits$log2FoldChange[i], 2), " (", direction, "), ",
          "padj = ", format.pval(top_hits$padj[i], digits = 2), "\n", sep = "")
    }
  }
  invisible(x)
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
