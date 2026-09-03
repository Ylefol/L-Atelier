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
#'
#' Mirrors the per-comparison format of \code{\link{print.artemis_dea}}
#' (features tested, significance breakdown, top hits), looped once per
#' comparison bundled in this object -- since \code{ARTEMIS_differential_counts}
#' and the time series DEA functions in this file (DESeq2- and limma-backed
#' alike) all answer the same question (differential expression between two
#' conditions), their console summaries read the same way, one comparison at
#' a time, regardless of which engine produced them.
#'
#' @param x An artemis_ts_de object
#' @param ... Additional arguments (ignored)
#' @method print artemis_ts_de
#' @export
print.artemis_ts_de <- function(x, ...) {
  cat("Time Series Differential Analysis Results (", x$type, ")\n", sep = "")
  cat(x$comparison, " -- ", length(x$results), " comparison(s)\n", sep = "")

  alpha <- x$parameters$alpha

  for (nm in names(x$results)) {
    de <- x$results[[nm]]
    s  <- de$summary

    cat("\n------------------------------\n")
    cat("Comparison: ", gsub("_", " ", nm), "\n", sep = "")
    cat("Features tested: ", s$n_tested, "\n", sep = "")

    n_sig <- s$n_sig_up + s$n_sig_down
    pct   <- if (s$n_tested > 0) round(100 * n_sig / s$n_tested, 1) else 0
    cat("Significant (FDR < ", alpha, "): ", n_sig, " (", pct, "%)\n", sep = "")
    cat("    - Up: ", s$n_sig_up, "\n", sep = "")
    cat("    - Down: ", s$n_sig_down, "\n", sep = "")

    if (n_sig > 0) {
      df       <- de$results
      sig_mask <- !is.na(df$padj) & df$padj < alpha
      top_hits <- head(df[sig_mask, ], min(5, sum(sig_mask)))

      cat("\n  Top significant features:\n")
      for (i in seq_len(nrow(top_hits))) {
        direction <- if (top_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
        cat("    ", top_hits$feature_id[i], ": log2FC = ",
            round(top_hits$log2FoldChange[i], 2), " (", direction, "), ",
            "padj = ", format.pval(top_hits$padj[i], digits = 2), "\n", sep = "")
      }
    }
  }
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


# ==============================================================================
# LIMMA TIME SERIES WORKFLOW
# For log2-scale data: Olink NPX, log2 mass spec, microarray
# Do NOT use for raw count data (RNA-seq, ATAC-seq) — use the DESeq2 workflow above
# ==============================================================================

#' Normalize Time Series Data for limma
#'
#' Validates and stores a log2-scale matrix with sample metadata as the first
#' step in the limma-based time series workflow. Unlike the DESeq2 workflow
#' (\code{ARTEMIS_normalize_timeseries}), no transformation is applied here:
#' the matrix is assumed to be already on a log2 scale (Olink NPX, log2
#' mass spec intensities, microarray log2 values). This step aligns samples,
#' validates required columns, and packages the inputs into an
#' \code{artemis_ts_norm_limma} object consumed by
#' \code{ARTEMIS_timeseries_conditional_limma()} and
#' \code{ARTEMIS_timeseries_temporal_limma()}.
#'
#' @param matrix Numeric matrix, features x samples (log2-scale). Rownames
#'   required. Columns must match \code{rownames(targets)}.
#' @param targets Data.frame with sample metadata. Rownames must match
#'   \code{colnames(matrix)}. Must include columns for group and timepoint.
#' @param group_col Character. Column in targets for group labels. Default:
#'   "group".
#' @param time_col Character. Column in targets for timepoint labels. Default:
#'   "timepoint".
#' @param batch_col Character or NULL. Batch column in targets. When provided,
#'   it is automatically included as a covariate in all downstream comparisons
#'   (prepended to any \code{covariates} passed to conditional/temporal
#'   functions). Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_norm_limma"} containing:
#'   \describe{
#'     \item{matrix}{The aligned log2-scale matrix (features x samples)}
#'     \item{targets}{The aligned targets data.frame}
#'     \item{parameters}{List: group_col, time_col, batch_col, n_samples,
#'       n_features}
#'   }
#'
#' @details
#' \strong{Olink NPX}: values are already log2-normalized by the Olink
#' platform. Pass \code{ol$wide} directly — no prior transformation needed.
#'
#' \strong{Mass spectrometry}: \code{ELEUTHIA_load_massspec()} stores
#' log2-transformed intensities in \code{massspec_data$wide}. Verify that
#' sample-level normalization (e.g., median centering) has been applied
#' upstream if required by your experiment design before passing the matrix
#' here — unlike Olink, DIA-NN output does not guarantee cross-sample
#' comparability without an explicit normalization step.
#'
#' @examples
#' \dontrun{
#' # Olink longitudinal
#' ts_norm <- ARTEMIS_normalize_timeseries_limma(
#'   ol$wide, ol$sample_meta,
#'   group_col = "Group", time_col = "Timepoint"
#' )
#'
#' # Mass spec with batch
#' ts_norm <- ARTEMIS_normalize_timeseries_limma(
#'   ms$wide, ms$sample_meta,
#'   group_col = "Group", time_col = "Visit", batch_col = "Plate"
#' )
#' }
#' @export
ARTEMIS_normalize_timeseries_limma <- function(matrix,
                                                targets,
                                                group_col = "group",
                                                time_col  = "timepoint",
                                                batch_col = NULL,
                                                verbose   = TRUE) {

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Install with: BiocManager::install('limma')")
  }

  .validate_ts_inputs(matrix, targets, group_col, time_col)

  if (!is.matrix(matrix)) matrix <- as.matrix(matrix)

  # Warn if matrix looks like raw counts
  non_na <- matrix[!is.na(matrix)]
  if (length(non_na) > 0 &&
      all(non_na == floor(non_na)) &&
      all(non_na >= 0)) {
    warning("Matrix contains only non-negative integers — it may be raw counts. ",
            "ARTEMIS_normalize_timeseries_limma() expects a log2-scale matrix. ",
            "For count data use the DESeq2 workflow (ARTEMIS_normalize_timeseries).")
  }

  if (!is.null(batch_col) && !batch_col %in% colnames(targets)) {
    stop("batch_col '", batch_col, "' not found in targets")
  }

  # Align samples
  common_samples <- intersect(colnames(matrix), rownames(targets))
  matrix  <- matrix[, common_samples, drop = FALSE]
  targets <- targets[common_samples, , drop = FALSE]

  if (verbose) {
    cat("[ARTEMIS] Time Series Normalization (limma)\n")
    cat("    Samples  :", ncol(matrix), "\n")
    cat("    Features :", nrow(matrix), "\n")
    cat("    Groups   :", paste(unique(targets[[group_col]]), collapse = ", "), "\n")
    cat("    Timepoints:", paste(sort(unique(targets[[time_col]])), collapse = ", "), "\n")
    if (!is.null(batch_col)) {
      cat("    Batch col :", batch_col, "(auto-included as covariate in comparisons)\n")
    }
    cat("    Note: matrix stored as-is (log2-scale assumed).\n\n")
  }

  result <- list(
    matrix     = matrix,
    targets    = targets,
    parameters = list(
      group_col  = group_col,
      time_col   = time_col,
      batch_col  = batch_col,
      n_samples  = ncol(matrix),
      n_features = nrow(matrix)
    )
  )
  class(result) <- c("artemis_ts_norm_limma", "list")
  return(result)
}


#' Print method for artemis_ts_norm_limma
#' @param x An artemis_ts_norm_limma object
#' @param ... Additional arguments (ignored)
#' @method print artemis_ts_norm_limma
#' @export
print.artemis_ts_norm_limma <- function(x, ...) {
  cat("Time series limma data (log2-scale):\n")
  cat("--------------------------------------\n")
  cat("Samples :", x$parameters$n_samples, "\n")
  cat("Features:", x$parameters$n_features, "\n")
  invisible(x)
}


# ==============================================================================
# LIMMA CONDITIONAL DEA
# ==============================================================================

#' Time Series Conditional Differential Expression (limma)
#'
#' Compares two groups at each timepoint independently using limma linear
#' models. For example, Treatment vs Control at TP1, TP2, TP3. Operates on
#' the log2-scale matrix stored in an \code{artemis_ts_norm_limma} object —
#' no transformation is applied during the comparison.
#'
#' @param ts_norm An \code{artemis_ts_norm_limma} object from
#'   \code{ARTEMIS_normalize_timeseries_limma()}.
#' @param reference Character. Reference/baseline group (denominator in FC).
#' @param experiment Character. Experimental group (numerator in FC).
#' @param block_col Character or NULL. Column in targets for the blocking
#'   variable (repeated measures within each timepoint comparison, e.g. a
#'   crossover design). When provided, \code{limma::duplicateCorrelation()}
#'   estimates the within-block correlation per timepoint. Default: NULL.
#' @param covariates Character vector or NULL. Additional covariate columns
#'   from targets to include in the design. If \code{batch_col} was set in
#'   \code{ARTEMIS_normalize_timeseries_limma()}, it is automatically
#'   prepended. Default: NULL.
#' @param alpha Numeric. FDR threshold for counting significant features.
#'   Default: 0.05.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_de"} (same class as the
#'   DESeq2 workflow) containing:
#'   \describe{
#'     \item{results}{Named list of per-timepoint result lists, each with
#'       \code{$results} (data.frame: feature_id, log2FoldChange, AveExpr,
#'       t, pvalue, padj, B, sig), \code{$fit}, \code{$summary}}
#'     \item{summary}{Data.frame: timepoint, experiment_name, n_tested,
#'       n_sig_up, n_sig_down}
#'     \item{type}{"conditional"}
#'     \item{comparison}{Character "experiment vs reference"}
#'     \item{norm_counts}{The log2 matrix (named for compatibility with
#'       \code{ARTEMIS_prepare_part_matrix()} and \code{ARTEMIS_select_de_genes()})}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' ts_norm <- ARTEMIS_normalize_timeseries_limma(ol$wide, ol$sample_meta,
#'   group_col = "Group", time_col = "Timepoint")
#'
#' cond_de <- ARTEMIS_timeseries_conditional_limma(
#'   ts_norm, reference = "Control", experiment = "Sepsis"
#' )
#' }
#' @export
ARTEMIS_timeseries_conditional_limma <- function(ts_norm,
                                                  reference,
                                                  experiment,
                                                  block_col  = NULL,
                                                  covariates = NULL,
                                                  alpha      = 0.05,
                                                  verbose    = TRUE) {

  if (!inherits(ts_norm, "artemis_ts_norm_limma")) {
    stop("'ts_norm' must be an artemis_ts_norm_limma object from ",
         "ARTEMIS_normalize_timeseries_limma().")
  }

  targets   <- ts_norm$targets
  group_col <- ts_norm$parameters$group_col
  time_col  <- ts_norm$parameters$time_col

  if (!reference %in% targets[[group_col]]) {
    stop("reference '", reference, "' not found in '", group_col, "' column")
  }
  if (!experiment %in% targets[[group_col]]) {
    stop("experiment '", experiment, "' not found in '", group_col, "' column")
  }
  if (!is.null(block_col) && !block_col %in% colnames(targets)) {
    stop("block_col '", block_col, "' not found in targets")
  }

  # Auto-prepend batch_col to covariates
  effective_covariates <- unique(c(ts_norm$parameters$batch_col, covariates))
  if (length(effective_covariates) == 0) effective_covariates <- NULL

  if (!is.null(effective_covariates)) {
    missing_cov <- setdiff(effective_covariates, colnames(targets))
    if (length(missing_cov) > 0) {
      stop("Covariate columns not found in targets: ",
           paste(missing_cov, collapse = ", "))
    }
  }

  timepoints <- sort(unique(targets[[time_col]]))

  if (verbose) {
    cat("[ARTEMIS] Time Series Conditional DEA (limma)\n")
    cat("    Comparison:", experiment, "vs", reference, "\n")
    cat("    Timepoints:", paste(timepoints, collapse = ", "), "\n")
    if (!is.null(block_col)) {
      cat("    Block     :", block_col,
          "(repeated measures via duplicateCorrelation)\n")
    }
    if (!is.null(effective_covariates)) {
      cat("    Covariates:", paste(effective_covariates, collapse = ", "), "\n")
    }
    cat("\n")
  }

  results      <- list()
  summary_rows <- list()

  for (tp in timepoints) {
    tp_label <- as.character(tp)

    tp_mask    <- targets[[time_col]] == tp &
                  targets[[group_col]] %in% c(reference, experiment)
    tp_samples <- rownames(targets)[tp_mask]
    tp_targets <- targets[tp_samples, , drop = FALSE]

    n_ref <- sum(tp_targets[[group_col]] == reference)
    n_exp <- sum(tp_targets[[group_col]] == experiment)

    if (n_ref < 2 || n_exp < 2) {
      if (verbose) {
        cat("        TP", tp_label, ": Skipping (",
            reference, "=", n_ref, ", ", experiment, "=", n_exp, " samples)\n")
      }
      next
    }

    if (verbose) cat("        TP", tp_label, ": ")

    mat_subset   <- ts_norm$matrix[, tp_samples, drop = FALSE]
    block_vec    <- if (!is.null(block_col)) tp_targets[[block_col]] else NULL
    covariate_df <- if (!is.null(effective_covariates))
                      tp_targets[, effective_covariates, drop = FALSE] else NULL

    de_result <- tryCatch({
      .run_limma_ts_subset(
        mat_subset   = mat_subset,
        group_labels = tp_targets[[group_col]],
        reference    = reference,
        experiment   = experiment,
        alpha        = alpha,
        block_vec    = block_vec,
        covariate_df = covariate_df
      )
    }, error = function(e) {
      if (verbose) cat("Error -", e$message, "\n")
      NULL
    })

    if (is.null(de_result)) next

    exp_name <- paste0(experiment, "_vs_", reference, "_TP", tp_label)
    de_result$timepoint       <- tp
    de_result$experiment_name <- exp_name
    results[[exp_name]] <- de_result

    summary_rows[[exp_name]] <- data.frame(
      timepoint       = tp,
      experiment_name = exp_name,
      n_tested        = de_result$summary$n_tested,
      n_sig_up        = de_result$summary$n_sig_up,
      n_sig_down      = de_result$summary$n_sig_down,
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
    results     = results,
    summary     = summary_df,
    type        = "conditional",
    comparison  = paste(experiment, "vs", reference),
    norm_counts = ts_norm$matrix,
    parameters  = list(
      reference   = reference,
      experiment  = experiment,
      group_col   = group_col,
      time_col    = time_col,
      block_col   = block_col,
      covariates  = effective_covariates,
      alpha       = alpha
    )
  )
  class(result) <- c("artemis_ts_de", "list")
  return(result)
}


# ==============================================================================
# LIMMA TEMPORAL DEA
# ==============================================================================

#' Time Series Temporal Differential Expression (limma)
#'
#' Compares timepoints across ALL samples, pooling groups together, to extract
#' the pure temporal effect — features that change over time regardless of
#' condition. For example, TP3 vs TP1 pools all TP3 samples against all TP1
#' samples; the condition variation averages out as noise. When subjects are
#' measured at multiple timepoints, use \code{block_col} (e.g. SubjectID) to
#' account for repeated measures via \code{limma::duplicateCorrelation()}.
#'
#' @param ts_norm An \code{artemis_ts_norm_limma} object from
#'   \code{ARTEMIS_normalize_timeseries_limma()}.
#' @param comparisons Character. "consecutive" (TP2 vs TP1, TP3 vs TP2) or
#'   "all" (all pairwise timepoint combinations). Default: "consecutive".
#' @param block_col Character or NULL. Column in targets for the blocking
#'   variable for repeated measures (e.g. SubjectID when the same individuals
#'   are measured at multiple timepoints). When provided,
#'   \code{limma::duplicateCorrelation()} estimates the within-subject
#'   correlation per comparison. Default: NULL.
#' @param covariates Character vector or NULL. Additional covariate columns
#'   from targets. If \code{batch_col} was set in
#'   \code{ARTEMIS_normalize_timeseries_limma()}, it is automatically
#'   prepended. Default: NULL.
#' @param alpha Numeric. FDR threshold for counting significant features.
#'   Default: 0.05.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_ts_de"} containing:
#'   \describe{
#'     \item{results}{Named list of per-comparison result lists, each with
#'       \code{$results} (data.frame: feature_id, log2FoldChange, AveExpr,
#'       t, pvalue, padj, B, sig), \code{$fit}, \code{$summary}}
#'     \item{summary}{Data.frame with comparison details and significance
#'       counts}
#'     \item{type}{"temporal"}
#'     \item{comparison}{Character description}
#'     \item{norm_counts}{The log2 matrix}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' ts_norm <- ARTEMIS_normalize_timeseries_limma(ol$wide, ol$sample_meta,
#'   group_col = "Group", time_col = "Timepoint")
#'
#' # Consecutive comparisons, blocking on subject (paired design)
#' temp_de <- ARTEMIS_timeseries_temporal_limma(
#'   ts_norm, block_col = "SubjectID"
#' )
#'
#' # All pairwise timepoint comparisons
#' temp_de <- ARTEMIS_timeseries_temporal_limma(
#'   ts_norm, comparisons = "all", block_col = "SubjectID"
#' )
#' }
#' @export
ARTEMIS_timeseries_temporal_limma <- function(ts_norm,
                                               comparisons = "consecutive",
                                               block_col   = NULL,
                                               covariates  = NULL,
                                               alpha       = 0.05,
                                               verbose     = TRUE) {

  if (!inherits(ts_norm, "artemis_ts_norm_limma")) {
    stop("'ts_norm' must be an artemis_ts_norm_limma object from ",
         "ARTEMIS_normalize_timeseries_limma().")
  }

  comparisons <- match.arg(comparisons, c("consecutive", "all"))

  targets  <- ts_norm$targets
  time_col <- ts_norm$parameters$time_col

  if (!is.null(block_col) && !block_col %in% colnames(targets)) {
    stop("block_col '", block_col, "' not found in targets")
  }

  # Auto-prepend batch_col to covariates
  effective_covariates <- unique(c(ts_norm$parameters$batch_col, covariates))
  if (length(effective_covariates) == 0) effective_covariates <- NULL

  if (!is.null(effective_covariates)) {
    missing_cov <- setdiff(effective_covariates, colnames(targets))
    if (length(missing_cov) > 0) {
      stop("Covariate columns not found in targets: ",
           paste(missing_cov, collapse = ", "))
    }
  }

  timepoints <- sort(unique(targets[[time_col]]))

  if (length(timepoints) < 2) {
    stop("Need at least 2 timepoints for temporal analysis. Found: ",
         length(timepoints))
  }

  if (comparisons == "consecutive") {
    pairs <- data.frame(
      reference  = timepoints[-length(timepoints)],
      experiment = timepoints[-1],
      stringsAsFactors = FALSE
    )
  } else {
    pairs <- expand.grid(
      reference  = timepoints,
      experiment = timepoints,
      stringsAsFactors = FALSE
    )
    pairs <- pairs[pairs$reference < pairs$experiment, ]
  }

  if (verbose) {
    cat("[ARTEMIS] Time Series Temporal DEA (limma)\n")
    cat("    Timepoints  :", paste(timepoints, collapse = ", "), "\n")
    cat("    Comparisons :", nrow(pairs), "(", comparisons, ")\n")
    cat("    Note: pooling all groups to extract pure temporal effect\n")
    if (!is.null(block_col)) {
      cat("    Block       :", block_col,
          "(repeated measures via duplicateCorrelation)\n")
    }
    if (!is.null(effective_covariates)) {
      cat("    Covariates  :", paste(effective_covariates, collapse = ", "), "\n")
    }
    cat("\n")
  }

  results      <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(pairs))) {
    tp_ref <- pairs$reference[i]
    tp_exp <- pairs$experiment[i]

    ref_label <- paste0("TP_", tp_ref)
    exp_label <- paste0("TP_", tp_exp)

    # All samples from both timepoints (all groups pooled)
    tp_mask    <- targets[[time_col]] %in% c(tp_ref, tp_exp)
    tp_samples <- rownames(targets)[tp_mask]
    tp_targets <- targets[tp_samples, , drop = FALSE]

    # Condition label is timepoint, not group
    tp_labels <- paste0("TP_", tp_targets[[time_col]])

    n_ref <- sum(tp_labels == ref_label)
    n_exp <- sum(tp_labels == exp_label)

    if (n_ref < 2 || n_exp < 2) {
      if (verbose) {
        cat("    ", exp_label, " vs ", ref_label,
            ": Skipping (n=", n_ref, ", ", n_exp, ")\n", sep = "")
      }
      next
    }

    if (verbose) cat("    ", exp_label, " vs ", ref_label, ": ", sep = "")

    mat_subset   <- ts_norm$matrix[, tp_samples, drop = FALSE]
    block_vec    <- if (!is.null(block_col)) tp_targets[[block_col]] else NULL
    covariate_df <- if (!is.null(effective_covariates))
                      tp_targets[, effective_covariates, drop = FALSE] else NULL

    de_result <- tryCatch({
      .run_limma_ts_subset(
        mat_subset   = mat_subset,
        group_labels = tp_labels,
        reference    = ref_label,
        experiment   = exp_label,
        alpha        = alpha,
        block_vec    = block_vec,
        covariate_df = covariate_df
      )
    }, error = function(e) {
      if (verbose) cat("Error -", e$message, "\n")
      NULL
    })

    if (is.null(de_result)) next

    exp_name <- paste0(exp_label, "_vs_", ref_label)
    de_result$timepoint_ref   <- tp_ref
    de_result$timepoint_exp   <- tp_exp
    de_result$experiment_name <- exp_name
    results[[exp_name]] <- de_result

    summary_rows[[exp_name]] <- data.frame(
      tp_reference    = tp_ref,
      tp_experiment   = tp_exp,
      experiment_name = exp_name,
      n_samples_ref   = n_ref,
      n_samples_exp   = n_exp,
      n_tested        = de_result$summary$n_tested,
      n_sig_up        = de_result$summary$n_sig_up,
      n_sig_down      = de_result$summary$n_sig_down,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(de_result$summary$n_sig_up, "up,",
          de_result$summary$n_sig_down, "down",
          "(n=", n_ref, "+", n_exp, " samples)\n")
    }
  }

  if (length(results) == 0) {
    stop("No successful comparisons. Check that each timepoint has >= 2 samples.")
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  result <- list(
    results     = results,
    summary     = summary_df,
    type        = "temporal",
    comparison  = paste("Temporal:", comparisons),
    norm_counts = ts_norm$matrix,
    parameters  = list(
      time_col   = time_col,
      block_col  = block_col,
      covariates = effective_covariates,
      comparisons = comparisons,
      alpha      = alpha
    )
  )
  class(result) <- c("artemis_ts_de", "list")
  return(result)
}


# ==============================================================================
# INTERNAL: Run limma on a subset of samples
# ==============================================================================

#' Core limma fitting for one time series comparison
#'
#' @param mat_subset Matrix (features x samples) already subsetted to the
#'   relevant samples.
#' @param group_labels Character vector of group/timepoint labels for each
#'   column of mat_subset.
#' @param reference Reference group label.
#' @param experiment Experiment group label.
#' @param alpha FDR threshold for the sig flag.
#' @param block_vec Character/factor vector of blocking variable values (one
#'   per sample), or NULL.
#' @param covariate_df Data.frame of covariate values (one row per sample),
#'   or NULL.
#' @return List: $results (data.frame), $fit (MArrayLM), $summary
#' @noRd
.run_limma_ts_subset <- function(mat_subset, group_labels, reference, experiment,
                                  alpha, block_vec = NULL, covariate_df = NULL) {

  safe_ref <- make.names(reference)
  safe_exp <- make.names(experiment)

  grp_factor  <- factor(group_labels, levels = c(reference, experiment))
  design_data <- data.frame(grp_factor = grp_factor, stringsAsFactors = FALSE)

  if (!is.null(covariate_df)) {
    for (cov in colnames(covariate_df)) {
      val <- covariate_df[[cov]]
      if (is.character(val)) val <- factor(val)
      design_data[[cov]] <- val
    }
    formula_str <- paste("~ 0 + grp_factor +",
                          paste(colnames(covariate_df), collapse = " + "))
  } else {
    formula_str <- "~ 0 + grp_factor"
  }

  design   <- stats::model.matrix(stats::as.formula(formula_str),
                                   data = design_data)
  grp_cols <- grep("^grp_factor", colnames(design))
  colnames(design)[grp_cols] <- c(safe_ref, safe_exp)

  contrast_str <- paste0(safe_exp, " - ", safe_ref)
  contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  if (!is.null(block_vec)) {
    corfit <- limma::duplicateCorrelation(mat_subset, design, block = block_vec)
    fit    <- limma::lmFit(mat_subset, design,
                            block       = block_vec,
                            correlation = corfit$consensus)
  } else {
    fit <- limma::lmFit(mat_subset, design)
  }

  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")

  results_df <- data.frame(
    feature_id     = rownames(tt),
    log2FoldChange = tt$logFC,
    AveExpr        = tt$AveExpr,
    t              = tt$t,
    pvalue         = tt$P.Value,
    padj           = tt$adj.P.Val,
    B              = tt$B,
    stringsAsFactors = FALSE
  )
  results_df$sig <- !is.na(results_df$padj) & results_df$padj < alpha
  results_df <- results_df[order(results_df$pvalue), ]
  rownames(results_df) <- NULL

  n_tested   <- sum(!is.na(results_df$padj))
  n_sig_up   <- sum(results_df$sig & results_df$log2FoldChange > 0, na.rm = TRUE)
  n_sig_down <- sum(results_df$sig & results_df$log2FoldChange < 0, na.rm = TRUE)

  list(
    results = results_df,
    fit     = fit2,
    summary = list(
      n_tested   = n_tested,
      n_sig_up   = n_sig_up,
      n_sig_down = n_sig_down
    )
  )
}
