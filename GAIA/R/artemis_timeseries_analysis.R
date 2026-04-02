# GAIA/Artemis/timeseries_analysis.R
# Time Series Differential Expression Analysis
#
# Proper time series DEA requires normalizing ALL samples together ONCE,
# then subsetting for individual comparisons. This ensures consistent
# size factors across all comparisons.
#
# Workflow:
# 1. ARTEMIS_normalize_timeseries() - normalize full dataset (REQUIRED FIRST)
# 2. ARTEMIS_timeseries_conditional() - compare groups at each timepoint
# 3. ARTEMIS_timeseries_temporal()    - compare timepoints within groups
# 4. ARTEMIS_select_de_genes()        - collect significant genes across comparisons
# 5. ARTEMIS_prepare_part_matrix()    - prepare matrix for PART clustering


# ==============================================================================
# NORMALIZATION (MUST BE RUN FIRST)
# ==============================================================================

#' Normalize Time Series Data with DESeq2
#'
#' Creates a DESeq2 object with size factors estimated on the FULL dataset.
#' This is a REQUIRED first step before running time series DEA. Normalizing
#' all samples together ensures consistent size factors across all subsequent
#' comparisons.
#'
#' @param counts Integer matrix of raw counts (genes x samples). Rownames required.
#' @param targets Data.frame with sample metadata. Rownames must match colnames
#'   of counts. Must include columns for group and timepoint.
#' @param group_col Character. Column in targets for group labels. Default: "group".
#' @param time_col Character. Column in targets for timepoint labels. Default: "timepoint".
#' @param batch_col Character or NULL. Optional batch column to include in design.
#'   Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_norm"} containing:
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
#' The resulting object is then passed to \code{ARTEMIS_timeseries_conditional()}
#' or \code{ARTEMIS_timeseries_temporal()} which will subset and run DESeq on
#' subsets while preserving the global size factors.
#'
#' @examples
#' \dontrun{
#' # Step 1: Normalize the full dataset
#' ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)
#'
#' # Step 2: Run conditional DEA (uses pre-computed size factors)
#' cond_de <- ARTEMIS_timeseries_conditional(
#'   ts_norm, reference = "Control", experiment = "Treatment"
#' )
#'
#' }
#' @export
ARTEMIS_normalize_timeseries <- function(counts,
                                          targets,
                                          group_col = "group",
                                          time_col = "timepoint",
                                          batch_col = NULL,
                                          verbose = TRUE) {


  # --- Validate inputs ---
  .validate_ts_inputs(counts, targets, group_col, time_col)

  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  # Check counts are integers

if (!all(counts == floor(counts))) {
    stop("Counts must be integers. DESeq2 requires raw counts.")
  }
  storage.mode(counts) <- "integer"

  # Align samples
  common_samples <- intersect(colnames(counts), rownames(targets))
  counts <- counts[, common_samples, drop = FALSE]
  targets <- targets[common_samples, , drop = FALSE]

  if (verbose) {
    cat("[ARTEMIS] Time Series Normalization \n")
    cat("    Samples:", ncol(counts), "\n")
    cat("    Genes:", nrow(counts), "\n")
    cat("    Groups:", paste(unique(targets[[group_col]]), collapse = ", "), "\n")
    cat("    Timepoints:", paste(sort(unique(targets[[time_col]])), collapse = ", "), "\n")
    if (!is.null(batch_col)) {
      cat("    Batch correction: Yes (", batch_col, ")\n", sep = "")
    }
    cat("\n")
  }

  # --- Create DESeq2 object ---
  # Design uses group as the condition (will be overridden for each comparison)
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

  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)

  if (verbose) {
    cat("    Size factors range:", round(min(size_factors), 3), "-",
        round(max(size_factors), 3), "\n")
    cat("    Normalization complete.\n\n")
  }

  result <- list(
    dds = dds,
    norm_counts = norm_counts,
    targets = targets,
    size_factors = size_factors,
    parameters = list(
      group_col = group_col,
      time_col = time_col,
      batch_col = batch_col,
      n_samples = ncol(counts),
      n_genes = nrow(counts)
    )
  )

  class(result) <- c("artemis_ts_norm", "list")
  return(result)
}


#' Print method for artemis_ts_norm
#' @param x An artemis_ts_norm object
#' @param ... Additional arguments (ignored)
#' @method print artemis_ts_norm
#' @export
print.artemis_ts_norm <- function(x, ...) {
  cat("Time series normalized data:\n")
  cat("------------------------------\n")
  cat("Samples:", x$parameters$n_samples, "\n")
  cat("Genes:", x$parameters$n_genes, "\n")
  cat("Size factors:", round(min(x$size_factors), 3), "-",
      round(max(x$size_factors), 3), "\n")
  invisible(x)
}


# ==============================================================================
# CONDITIONAL DEA
# ==============================================================================

#' Time Series Conditional Differential Expression
#'
#' Compares two groups at each timepoint independently. For example, Treatment
#' vs Control at TP1, TP2, TP3. Uses size factors from the pre-normalized
#' dataset to ensure consistent normalization.
#'
#' @param ts_norm An \code{artemis_ts_norm} object from
#'   \code{ARTEMIS_normalize_timeseries()}. This ensures all comparisons use
#'   consistent size factors estimated on the full dataset.
#' @param reference Character. Reference/baseline group (denominator in FC).
#' @param experiment Character. Experimental group (numerator in FC).
#' @param alpha Numeric. FDR threshold. Default: 0.05.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_de"} containing:
#'   \describe{
#'     \item{results}{Named list of per-timepoint DESeq2 result data.frames}
#'     \item{summary}{Data.frame: timepoint, n_tested, n_sig_up, n_sig_down}
#'     \item{type}{"conditional"}
#'     \item{comparison}{Character description}
#'     \item{norm_counts}{Normalized count matrix from ts_norm}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' # First normalize
#' ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)
#'
#' # Then run conditional DEA
#' cond_de <- ARTEMIS_timeseries_conditional(
#'   ts_norm, reference = "Control", experiment = "Treatment"
#' )
#'
#' }
#' @export
ARTEMIS_timeseries_conditional <- function(ts_norm,
                                            reference,
                                            experiment,
                                            alpha = 0.05,
                                            verbose = TRUE) {

  # --- Validate ---
  if (!inherits(ts_norm, "artemis_ts_norm")) {
    stop("'ts_norm' must be an artemis_ts_norm object from ARTEMIS_normalize_timeseries().\n",
         "Time series DEA requires normalizing the full dataset first to ensure ",
         "consistent size factors across all comparisons.")
  }

  targets <- ts_norm$targets
  group_col <- ts_norm$parameters$group_col
  time_col <- ts_norm$parameters$time_col

  if (!reference %in% targets[[group_col]]) {
    stop("reference '", reference, "' not found in '", group_col, "' column")
  }
  if (!experiment %in% targets[[group_col]]) {
    stop("experiment '", experiment, "' not found in '", group_col, "' column")
  }

  # Get timepoints (sorted)
  timepoints <- sort(unique(targets[[time_col]]))

  if (verbose) {
    cat("[ARTEMIS] Time Series Conditional DEA \n")
    cat("    Comparison:", experiment, "vs", reference, "\n")
    cat("    Timepoints:", paste(timepoints, collapse = ", "), "\n")
    cat("    Using pre-computed size factors from full dataset\n\n")
  }

  # --- Run DEA per timepoint ---
  results <- list()
  summary_rows <- list()

  for (tp in timepoints) {
    tp_label <- as.character(tp)

    # Subset to this timepoint and the two groups
    tp_mask <- targets[[time_col]] == tp &
      targets[[group_col]] %in% c(reference, experiment)
    tp_samples <- rownames(targets)[tp_mask]

    # Check we have enough samples
    tp_targets <- targets[tp_samples, , drop = FALSE]
    n_ref <- sum(tp_targets[[group_col]] == reference)
    n_exp <- sum(tp_targets[[group_col]] == experiment)

    if (n_ref < 2 || n_exp < 2) {
      if (verbose) {
        cat("        TP", tp_label, ": Skipping (", reference, "=", n_ref,
            ", ", experiment, "=", n_exp, " samples)\n")
      }
      next
    }

    if (verbose) cat("        TP", tp_label, ": ")

    # Run DESeq on subset with pre-computed size factors
    de_result <- tryCatch({
      .run_deseq_subset(
        ts_norm = ts_norm,
        samples = tp_samples,
        group_labels = tp_targets[[group_col]],
        reference = reference,
        experiment = experiment,
        alpha = alpha
      )
    }, error = function(e) {
      if (verbose) cat("    Error -", e$message, "\n")
      NULL
    })

    if (is.null(de_result)) next

    # Store result
    exp_name <- paste0(experiment, "_vs_", reference, "_TP", tp_label)
    de_result$timepoint <- tp
    de_result$experiment_name <- exp_name
    results[[exp_name]] <- de_result

    # Summary stats
    summary_rows[[exp_name]] <- data.frame(
      timepoint = tp,
      experiment_name = exp_name,
      n_tested = de_result$summary$n_tested,
      n_sig_up = de_result$summary$n_sig_up,
      n_sig_down = de_result$summary$n_sig_down,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat("    ", de_result$summary$n_sig_up, "up,",
          de_result$summary$n_sig_down, "down\n")
    }
  }

  if (length(results) == 0) {
    stop("No successful comparisons. Check that each timepoint has >= 2 ",
         "replicates per group.")
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  result <- list(
    results = results,
    summary = summary_df,
    type = "conditional",
    comparison = paste(experiment, "vs", reference),
    norm_counts = ts_norm$norm_counts,
    parameters = list(
      reference = reference,
      experiment = experiment,
      group_col = group_col,
      time_col = time_col,
      batch_col = ts_norm$parameters$batch_col,
      alpha = alpha
    )
  )

  class(result) <- c("artemis_ts_de", "list")
  return(result)
}


# ==============================================================================
# TEMPORAL DEA
# ==============================================================================

#' Time Series Temporal Differential Expression
#'
#' Compares timepoints across ALL samples, pooling groups together. This extracts
#' the pure temporal effect - genes that change over time regardless of condition.
#' For example, TP3 vs TP1 compares all TP3 samples (IgM + LPS) against all TP1
#' samples (IgM + LPS). The group/condition variation becomes biological noise
#' that averages out, leaving only the time-dependent changes.
#'
#' Uses size factors from the pre-normalized dataset to ensure consistent
#' normalization across all comparisons.
#'
#' @param ts_norm An \code{artemis_ts_norm} object from
#'   \code{ARTEMIS_normalize_timeseries()}.
#' @param comparisons Character. "consecutive" (TP2 vs TP1, TP3 vs TP2) or
#'   "all" (all pairwise timepoint combinations). Default: "consecutive".
#' @param alpha Numeric. FDR threshold for counting significant genes.
#'   Default: 0.05.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_de"} containing:
#'   \describe{
#'     \item{results}{Named list of per-comparison DESeq2 results}
#'     \item{summary}{Data.frame with comparison details and significance counts}
#'     \item{type}{"temporal"}
#'     \item{norm_counts}{Normalized count matrix from ts_norm}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @details
#' Unlike conditional DEA (which compares groups at each timepoint), temporal
#' DEA pools all groups together and uses timepoint as the "condition" for
#' DESeq2. This design captures genes with consistent temporal behavior across
#' all experimental conditions.
#'
#' @examples
#' \dontrun{
#' # First normalize
#' ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)
#'
#' # Consecutive timepoint comparisons (TP2 vs TP1, TP3 vs TP2, etc.)
#' temp_de <- ARTEMIS_timeseries_temporal(ts_norm)
#'
#' # All pairwise timepoint comparisons
#' temp_de <- ARTEMIS_timeseries_temporal(ts_norm, comparisons = "all")
#'
#' }
#' @export
ARTEMIS_timeseries_temporal <- function(ts_norm,
                                         comparisons = "consecutive",
                                         alpha = 0.05,
                                         verbose = TRUE) {

  # --- Validate ---
  if (!inherits(ts_norm, "artemis_ts_norm")) {
    stop("'ts_norm' must be an artemis_ts_norm object from ARTEMIS_normalize_timeseries().\n",
         "Time series DEA requires normalizing the full dataset first to ensure ",
         "consistent size factors across all comparisons.")
  }

  comparisons <- match.arg(comparisons, c("consecutive", "all"))

  targets <- ts_norm$targets
  time_col <- ts_norm$parameters$time_col

  timepoints <- sort(unique(targets[[time_col]]))

  if (length(timepoints) < 2) {
    stop("Need at least 2 timepoints for temporal analysis. Found: ",
         length(timepoints))
  }

  # Generate comparison pairs
  if (comparisons == "consecutive") {
    pairs <- data.frame(
      reference = timepoints[-length(timepoints)],
      experiment = timepoints[-1],
      stringsAsFactors = FALSE
    )
  } else {
    pairs <- expand.grid(
      reference = timepoints,
      experiment = timepoints,
      stringsAsFactors = FALSE
    )
    pairs <- pairs[pairs$reference < pairs$experiment, ]
  }

  if (verbose) {
    cat("[ARTEMIS] Time Series Temporal DEA \n")
    cat("    Timepoints:", paste(timepoints, collapse = ", "), "\n")
    cat("    Comparisons:", nrow(pairs), "(", comparisons, ")\n")
    cat("    Note: Pooling all groups together to extract pure temporal effect\n")
    cat("    Using pre-computed size factors from full dataset\n\n")
  }

  # --- Run DEA per timepoint comparison (pooling ALL groups) ---
  results <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(pairs))) {
    tp_ref <- pairs$reference[i]
    tp_exp <- pairs$experiment[i]

    ref_label <- paste0("TP_", tp_ref)
    exp_label <- paste0("TP_", tp_exp)

    # Get ALL samples from both timepoints (regardless of group)
    tp_mask <- targets[[time_col]] %in% c(tp_ref, tp_exp)
    tp_samples <- rownames(targets)[tp_mask]
    tp_targets <- targets[tp_samples, , drop = FALSE]

    # Create condition labels based on timepoint (not group!)
    tp_labels <- paste0("TP_", tp_targets[[time_col]])

    n_ref <- sum(tp_labels == ref_label)
    n_exp <- sum(tp_labels == exp_label)

    if (n_ref < 2 || n_exp < 2) {
      if (verbose) {
        cat("    ", exp_label, " vs ", ref_label,
            ": Skipping (n=", n_ref, ",", n_exp, ")\n", sep = "")
      }
      next
    }

    if (verbose) cat("    ", exp_label, " vs ", ref_label, ": ", sep = "")

    de_result <- tryCatch({
      .run_deseq_subset(
        ts_norm = ts_norm,
        samples = tp_samples,
        group_labels = tp_labels,
        reference = ref_label,
        experiment = exp_label,
        alpha = alpha
      )
    }, error = function(e) {
      if (verbose) cat("    Error -", e$message, "\n")
      NULL
    })

    if (is.null(de_result)) next

    exp_name <- paste0(exp_label, "_vs_", ref_label)
    de_result$timepoint_ref <- tp_ref
    de_result$timepoint_exp <- tp_exp
    de_result$experiment_name <- exp_name
    results[[exp_name]] <- de_result

    summary_rows[[exp_name]] <- data.frame(
      tp_reference = tp_ref,
      tp_experiment = tp_exp,
      experiment_name = exp_name,
      n_samples_ref = n_ref,
      n_samples_exp = n_exp,
      n_tested = de_result$summary$n_tested,
      n_sig_up = de_result$summary$n_sig_up,
      n_sig_down = de_result$summary$n_sig_down,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat("    ", de_result$summary$n_sig_up, "up,",
          de_result$summary$n_sig_down, "down",
          "(n=", n_ref, "+", n_exp, " samples)\n", sep = "")
    }
  }

  if (length(results) == 0) {
    stop("No successful comparisons. Check that each timepoint has >= 2 samples.")
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  result <- list(
    results = results,
    summary = summary_df,
    type = "temporal",
    comparison = paste("Temporal:", comparisons),
    norm_counts = ts_norm$norm_counts,
    parameters = list(
      time_col = time_col,
      batch_col = ts_norm$parameters$batch_col,
      comparisons = comparisons,
      alpha = alpha
    )
  )

  class(result) <- c("artemis_ts_de", "list")
  return(result)
}


# ==============================================================================
# INTERNAL: Run DESeq on subset with pre-computed size factors
# ==============================================================================

#' Run DESeq2 on a subset of samples using pre-computed size factors
#'
#' This is the core function that ensures consistent normalization across
#' all time series comparisons by using size factors estimated on the full
#' dataset.
#'
#' @param ts_norm The artemis_ts_norm object with pre-computed size factors
#' @param samples Character vector of sample names to include
#' @param group_labels Character vector of group labels for these samples
#' @param reference Reference group label
#' @param experiment Experiment group label
#' @param alpha FDR threshold
#'
#' @return List with results data.frame and summary stats
#' @noRd
.run_deseq_subset <- function(ts_norm, samples, group_labels,
                               reference, experiment, alpha) {

  # Subset the DESeq2 object
  dds_subset <- ts_norm$dds[, samples]

  # Update the condition factor for this comparison
  dds_subset$condition <- factor(group_labels, levels = c(reference, experiment))

  # The size factors are already set from the full dataset normalization
  # DESeq will use these existing size factors

  # Run DESeq (dispersion estimation + statistical testing)
  dds_subset <- DESeq(dds_subset, quiet = TRUE)

  # Extract results
  # Note: We do NOT pass alpha to results() to match TiSA behavior.
  # The alpha parameter in results() controls independent filtering optimization
  # (which genes get NA padj), NOT the significance threshold for counting.
  # TiSA uses the default (alpha=0.1), then filters at the user's threshold.
  # Passing alpha=0.05 would cause more aggressive filtering and fewer testable genes.
  res <- results(dds_subset,
                  contrast = c("condition", experiment, reference))

  # Format output
  results_df <- as.data.frame(res)
  results_df$feature_id <- rownames(results_df)

  # Reorder columns
  results_df <- results_df[, c("feature_id", "baseMean", "log2FoldChange",
                                "lfcSE", "stat", "pvalue", "padj")]

  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$padj), ]
  rownames(results_df) <- NULL

  # Summary statistics
  n_tested <- sum(!is.na(results_df$padj))
  n_sig_up <- sum(results_df$padj < alpha & results_df$log2FoldChange > 0, na.rm = TRUE)
  n_sig_down <- sum(results_df$padj < alpha & results_df$log2FoldChange < 0, na.rm = TRUE)

  list(
    results = results_df,
    dds = dds_subset,
    summary = list(
      n_tested = n_tested,
      n_sig_up = n_sig_up,
      n_sig_down = n_sig_down
    )
  )
}


# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print method for artemis_ts_de
#' @param x An artemis_ts_de object
#' @param ... Additional arguments (ignored)
#' @method print artemis_ts_de
#' @export
print.artemis_ts_de <- function(x, ...) {
  cat("Time series DE (", x$type, "):", length(x$results), "comparisons\n")
  cat("------------------------------\n")
  total_up <- sum(x$summary$n_sig_up)
  total_down <- sum(x$summary$n_sig_down)
  cat("Total significant: ", total_up, " up, ", total_down, " down\n", sep = "")
  invisible(x)
}


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Validate time series inputs
#' @noRd
.validate_ts_inputs <- function(counts, targets, group_col, time_col) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame")
  }
  if (!is.data.frame(targets)) {
    stop("'targets' must be a data.frame")
  }
  if (!group_col %in% colnames(targets)) {
    stop("'", group_col, "' column not found in targets. Available: ",
         paste(colnames(targets), collapse = ", "))
  }
  if (!time_col %in% colnames(targets)) {
    stop("'", time_col, "' column not found in targets. Available: ",
         paste(colnames(targets), collapse = ", "))
  }

  # Ensure targets has rownames matching count columns
  if (is.null(rownames(targets))) {
    stop("'targets' must have rownames matching colnames of 'counts'")
  }
  common <- intersect(colnames(counts), rownames(targets))
  if (length(common) == 0) {
    stop("No matching sample names between counts colnames and targets rownames")
  }
}
