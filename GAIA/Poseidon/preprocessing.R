###############################################################################
########### Preprocessing & Normalization for Multi-Omics Data ###########
###############################################################################

#' Scale and Center Data Matrix
#'
#' @description Centers and scales columns of a data matrix.
#' Wrapper around base::scale with consistent handling.
#'
#' @param X Data matrix (samples x features)
#' @param center Logical, should columns be centered? (default = TRUE)
#' @param scale Logical, should columns be scaled to unit variance? (default = TRUE)
#'
#' @return Scaled matrix with centering and scaling attributes
#'
#' @export
POSEIDON_scale_data <- function(X, center = TRUE, scale = TRUE) {
    return(scale(X, center = center, scale = scale))
}

#' Filter Features by Variance
#'
#' @description Removes features with zero or near-zero variance.
#' Essential preprocessing step before CCA to avoid numerical issues.
#'
#' @param X Data matrix or list of data matrices (samples x features)
#' @param variance_threshold Minimum variance to retain feature (default = 0)
#' @param verbose Print filtering summary (default = TRUE)
#'
#' @return Filtered data matrix or list of matrices
#'
#' @export
POSEIDON_filter_zero_variance <- function(X, variance_threshold = 0, verbose = TRUE) {

    # Handle single matrix
    if (!is.list(X)) {
        vars <- apply(X, 2, var, na.rm = TRUE)
        keep_features <- vars > variance_threshold

        if (verbose) {
            cat(sprintf("Removing %d features with variance <= %.2e\n",
                       sum(!keep_features), variance_threshold))
            cat(sprintf("Keeping %d features with variance > %.2e\n",
                       sum(keep_features), variance_threshold))
        }

        return(X[, keep_features, drop = FALSE])
    }

    # Handle list of matrices
    if (verbose) {
        cat("Filtering features across multiple datasets...\n")
    }

    n_datasets <- length(X)
    vars_list <- lapply(X, function(x) apply(x, 2, var, na.rm = TRUE))

    filtered_X <- list()
    for (i in 1:n_datasets) {
        keep_features <- vars_list[[i]] > variance_threshold

        if (verbose) {
            cat(sprintf("  Dataset %d: Removing %d features, keeping %d\n",
                       i, sum(!keep_features), sum(keep_features)))
        }

        filtered_X[[i]] <- X[[i]][, keep_features, drop = FALSE]
    }

    return(filtered_X)
}

#' Filter Features by Shared Variance Across Datasets
#'
#' @description For multi-dataset integration, removes genes that have
#' zero variance in ANY dataset. Ensures all datasets have valid features.
#'
#' @param X_list List of data matrices (all must have same features in rows or columns)
#' @param by_row Logical, are features in rows? (default = TRUE for gene expression)
#' @param variance_threshold Minimum variance (default = 0)
#' @param verbose Print filtering summary (default = TRUE)
#'
#' @return List of filtered matrices with only features having variance > 0 in ALL datasets
#'
#' @export
POSEIDON_filter_shared_variance <- function(X_list, by_row = TRUE,
                                   variance_threshold = 0, verbose = TRUE) {

    n_datasets <- length(X_list)

    # Compute variance for each dataset
    if (by_row) {
        vars_list <- lapply(X_list, function(x) apply(x, 1, var, na.rm = TRUE))
    } else {
        vars_list <- lapply(X_list, function(x) apply(x, 2, var, na.rm = TRUE))
    }

    # Find features with variance > threshold in ALL datasets
    keep_features <- Reduce("&", lapply(vars_list, function(v) v > variance_threshold))

    if (verbose) {
        cat(sprintf("Filtering features with variance > %.2e across %d datasets:\n",
                   variance_threshold, n_datasets))
        cat(sprintf("  Removing %d features with zero variance in ANY dataset\n",
                   sum(!keep_features)))
        cat(sprintf("  Keeping %d features with variance > %.2e in ALL datasets\n",
                   sum(keep_features), variance_threshold))
    }

    # Filter all datasets
    if (by_row) {
        filtered_list <- lapply(X_list, function(x) x[keep_features, , drop = FALSE])
    } else {
        filtered_list <- lapply(X_list, function(x) x[, keep_features, drop = FALSE])
    }

    return(filtered_list)
}


#' Filter Top Variable Features
#'
#' @description Selects the top N most variable features from each dataset.
#' Essential for reducing dimensionality before methods like sCCA that don't
#' scale well with very high feature counts (>10k features).
#'
#' Works with any numeric data: counts, normalized values, log-transformed, etc.
#'
#' @param X Data matrix (samples x features) or a named list of data matrices.
#'   For sCCA input, this should be the transposed format (samples in rows).
#' @param n_top Integer. Number of top variable features to keep per dataset.
#'   If a single value, applied to all datasets. Can also be a named list
#'   matching dataset names for dataset-specific values (default = 5000).
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered data in same format as input (matrix or list of matrices).
#'   Column names are preserved to track which features were retained.
#'
#' @details
#' Variance is calculated per feature (column) using \code{var()}.
#' Features are ranked by variance and the top N are retained.
#'
#' For multi-omics sCCA, typical values:
#' \itemize{
#'   \item ATAC-seq: 5,000 - 10,000 peaks
#'   \item ChIP-seq: 1,000 - 5,000 regions (often fewer to start)
#'   \item RNA-seq: 5,000 - 10,000 genes
#' }
#'
#' If n_top exceeds the number of features in a dataset, all features are kept
#' with a warning.
#'
#' @export
#'
#' @examples
#' # Single matrix
#' X_reduced <- POSEIDON_filter_top_variable(X, n_top = 5000)
#'
#' # List of matrices (for sCCA)
#' scca_data <- list(ATAC = atac_matrix, ChIP = chip_matrix, RNA = rna_matrix)
#' scca_reduced <- POSEIDON_filter_top_variable(scca_data, n_top = 5000)
#'
#' # Different n_top per dataset
#' scca_reduced <- POSEIDON_filter_top_variable(
#'   scca_data,
#'   n_top = list(ATAC = 10000, ChIP = 2000, RNA = 5000)
#' )
#'
POSEIDON_filter_top_variable <- function(X, n_top = 5000, verbose = TRUE) {

  # Helper function for single matrix
  filter_matrix <- function(mat, n, name = NULL) {
    n_features <- ncol(mat)

    # Check if n_top exceeds available features
    if (n >= n_features) {
      if (verbose) {
        prefix <- if (!is.null(name)) paste0(name, ": ") else ""
        warning(prefix, "n_top (", n, ") >= number of features (", n_features,
                "). Keeping all features.")
      }
      return(mat)
    }

    # Calculate variance per feature
    vars <- apply(mat, 2, var, na.rm = TRUE)

    # Handle any NA variances (shouldn't happen but be safe)
    vars[is.na(vars)] <- 0

    # Get indices of top N by variance
    top_idx <- order(vars, decreasing = TRUE)[1:n]

    # Sort indices to preserve original order
    top_idx <- sort(top_idx)

    if (verbose) {
      prefix <- if (!is.null(name)) paste0("  ", name, ": ") else ""
      cat(sprintf("%s%d -> %d features (top %.1f%% by variance)\n",
                  prefix, n_features, n, 100 * n / n_features))
      cat(sprintf("%s  Variance range kept: %.2e to %.2e\n",
                  prefix, min(vars[top_idx]), max(vars[top_idx])))
    }

    return(mat[, top_idx, drop = FALSE])
  }

  # Handle single matrix
  if (!is.list(X)) {
    if (verbose) {
      cat("Filtering to top", n_top, "variable features...\n")
    }
    return(filter_matrix(X, n_top))
  }

  # Handle list of matrices
  if (verbose) {
    cat("Filtering to top variable features per dataset...\n")
  }

  dataset_names <- names(X)
  n_datasets <- length(X)

  # Handle n_top as list or single value
  if (is.list(n_top)) {
    # Ensure n_top has entries for all datasets
    if (!all(dataset_names %in% names(n_top))) {
      missing <- setdiff(dataset_names, names(n_top))
      stop("n_top list missing entries for: ", paste(missing, collapse = ", "))
    }
    n_top_vec <- unlist(n_top[dataset_names])
  } else {
    # Single value applied to all
    n_top_vec <- rep(n_top, n_datasets)
    names(n_top_vec) <- dataset_names
  }

  # Filter each dataset
  filtered_X <- list()
  for (i in seq_along(X)) {
    name <- dataset_names[i]
    filtered_X[[name]] <- filter_matrix(X[[i]], n_top_vec[i], name)
  }

  if (verbose) {
    cat("Done.\n")
  }

  return(filtered_X)
}


#' Remove Low Variance Features
#'
#' @description Removes the bottom percentage of features by variance.
#' A conservative noise-removal approach that doesn't bias toward high-variance
#' features - it simply removes features that are clearly uninformative.
#'
#' Preferred over \code{POSEIDON_filter_top_variable()} for methods like sCCA
#' where you want to preserve cross-dataset correlations, not just variance.
#'
#' @param X Data matrix (samples x features) or a named list of data matrices.
#' @param bottom_pct Numeric. Percentage of lowest-variance features to remove
#'   (default = 0.5, removes bottom 50%). Value between 0 and 1.
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered data in same format as input (matrix or list of matrices).
#'   Column names are preserved to track which features were retained.
#'
#' @details
#' Unlike \code{POSEIDON_filter_top_variable()} which selects features based on
#' high variance (potentially biasing results), this function only removes
#' features that are clearly uninformative (very low variance).
#'
#' This is more appropriate for integration methods like sCCA where the goal
#' is to find cross-dataset correlations, not necessarily high-variance features.
#' A feature with moderate variance but strong cross-dataset correlation is
#' valuable and should be retained.
#'
#' Suggested values:
#' \itemize{
#'   \item 0.25 - Conservative, removes only bottom 25%
#'   \item 0.50 - Moderate, removes bottom half (default)
#'   \item 0.75 - Aggressive, keeps only top 25%
#' }
#'
#' @export
#'
#' @examples
#' # Remove bottom 50% of features by variance
#' X_filtered <- POSEIDON_filter_low_variance(X, bottom_pct = 0.5)
#'
#' # For sCCA with multiple datasets
#' scca_data <- list(ATAC = atac_matrix, ChIP = chip_matrix, RNA = rna_matrix)
#' scca_filtered <- POSEIDON_filter_low_variance(scca_data, bottom_pct = 0.5)
#'
POSEIDON_filter_low_variance <- function(X, bottom_pct = 0.5, verbose = TRUE) {

  # Validate bottom_pct

  if (bottom_pct <= 0 || bottom_pct >= 1) {
    stop("bottom_pct must be between 0 and 1 (exclusive)")
  }

  # Helper function for single matrix
  filter_matrix <- function(mat, pct, name = NULL) {
    n_features <- ncol(mat)
    n_remove <- floor(n_features * pct)
    n_keep <- n_features - n_remove

    if (n_keep < 2) {
      warning("bottom_pct too high - would remove all but ", n_keep, " features. ",
              "Keeping at least 2 features.")
      n_keep <- max(2, n_keep)
      n_remove <- n_features - n_keep
    }

    # Calculate variance per feature
    vars <- apply(mat, 2, var, na.rm = TRUE)

    # Handle any NA variances
    vars[is.na(vars)] <- 0

    # Get indices to KEEP (everything except bottom n_remove)
    keep_idx <- order(vars, decreasing = TRUE)[1:n_keep]

    # Sort indices to preserve original order
    keep_idx <- sort(keep_idx)

    # Get variance threshold (the cutoff point)
    var_threshold <- vars[order(vars, decreasing = TRUE)[n_keep]]

    if (verbose) {
      prefix <- if (!is.null(name)) paste0("  ", name, ": ") else ""
      cat(sprintf("%s%d -> %d features (removed bottom %.0f%%)\n",
                  prefix, n_features, n_keep, pct * 100))
      cat(sprintf("%s  Variance threshold: %.2e (features below this removed)\n",
                  prefix, var_threshold))
    }

    return(mat[, keep_idx, drop = FALSE])
  }

  # Handle single matrix
  if (!is.list(X)) {
    if (verbose) {
      cat(sprintf("Removing bottom %.0f%% of features by variance...\n", bottom_pct * 100))
    }
    return(filter_matrix(X, bottom_pct))
  }

  # Handle list of matrices
  if (verbose) {
    cat(sprintf("Removing bottom %.0f%% of features by variance per dataset...\n",
                bottom_pct * 100))
  }

  dataset_names <- names(X)

  # Filter each dataset
  filtered_X <- list()
  for (i in seq_along(X)) {
    name <- dataset_names[i]
    filtered_X[[name]] <- filter_matrix(X[[i]], bottom_pct, name)
  }

  if (verbose) {
    cat("Done.\n")
  }

  return(filtered_X)
}


#' Validate Multi-Dataset Structure
#'
#' @description Checks that all datasets in a list have:
#' - Same number of samples (rows)
#' - No missing values
#' - Valid dimensions
#'
#' @param X_list List of data matrices
#' @param sample_names Optional vector of expected sample names
#'
#' @return Invisible TRUE if valid, stops with error otherwise
#'
#' @export
POSEIDON_validate_multi_dataset <- function(X_list, sample_names = NULL) {

    if (!is.list(X_list)) {
        stop("X must be a list of data matrices")
    }

    n_datasets <- length(X_list)
    n_samples <- nrow(X_list[[1]])

    for (i in 1:n_datasets) {
        # Check for missing values
        if (any(is.na(X_list[[i]]))) {
            stop(sprintf("Dataset %d contains missing values. Please remove/impute before analysis.", i))
        }

        # Check sample size consistency
        if (nrow(X_list[[i]]) != n_samples) {
            stop(sprintf("Dataset %d has %d samples, expected %d. All datasets must have the same samples.",
                        i, nrow(X_list[[i]]), n_samples))
        }

        # Check for valid dimensions
        if (ncol(X_list[[i]]) == 0) {
            stop(sprintf("Dataset %d has zero features after filtering.", i))
        }
    }

    # Check sample names if provided
    if (!is.null(sample_names)) {
        for (i in 1:n_datasets) {
            if (!all(rownames(X_list[[i]]) == sample_names)) {
                warning(sprintf("Dataset %d has different sample names/order than expected.", i))
            }
        }
    }

    invisible(TRUE)
}


#' Filter Low-Count Features (Peaks/Genes)
#'
#' @description Filters features (peaks, genes, etc.) with insufficient counts
#' across samples. Commonly used to remove unreliable features before
#' differential analysis or integration.
#'
#' @param quant_result A list containing at minimum:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item annotation: Data.frame with feature annotations (optional)
#'   }
#' @param min_count Integer. Minimum count threshold (default = 10).
#' @param min_samples Integer. Minimum number of samples that must meet
#'   min_count threshold (default = 2).
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return The input list with filtered counts and annotation.
#'
#' @details
#' A feature is retained if at least \code{min_samples} samples have
#' counts >= \code{min_count}. This filtering removes features that are:
#' \itemize{
#'   \item Not detected in most samples
#'   \item Too low-count for reliable statistical analysis
#'   \item Likely to introduce noise in downstream analyses
#' }
#'
#' Common thresholds:
#' \itemize{
#'   \item ATAC-seq peaks: min_count = 10, min_samples = 3
#'   \item RNA-seq genes: min_count = 10, min_samples = 2-3
#' }
#'
#' @export
#'
#' @examples
#' # Filter peaks with at least 10 counts in at least 3 samples
#' counts_filtered <- POSEIDON_filter_low_counts(counts, min_count = 10, min_samples = 3)
#'
POSEIDON_filter_low_counts <- function(quant_result,
                                        min_count = 10,
                                        min_samples = 2,
                                        verbose = TRUE) {

  counts <- quant_result$counts

  # Count how many samples meet threshold for each feature
  n_above_threshold <- rowSums(counts >= min_count)

  # Filter

  keep <- n_above_threshold >= min_samples

  n_before <- nrow(counts)
  n_after <- sum(keep)

  if (verbose) {
    cat("Filtering low-count features:\n")
    cat("  Before:", n_before, "features\n")
    cat("  After:", n_after, "features\n")
    cat("  Removed:", n_before - n_after, "features\n")
    cat("  (min_count =", min_count, ", min_samples =", min_samples, ")\n")
  }

  # Update result
  quant_result$counts <- counts[keep, , drop = FALSE]

  # Update annotation if present
  if (!is.null(quant_result$annotation)) {
    quant_result$annotation <- quant_result$annotation[keep, , drop = FALSE]
  }

  return(quant_result)
}


#' Batch Correction for Count Data
#'
#' @description Removes batch effects from count data using ComBat or limma.
#' Preserves biological variation while removing technical batch effects.
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata including batch information
#'   }
#' @param batch_col Character string. Column name in targets containing batch
#'   information (default = "batch").
#' @param group_col Character string. Optional column name for biological groups
#'   to preserve during correction (default = "group").
#' @param method Character string. Batch correction method:
#'   "combat" (default) or "limma".
#' @param log_transform Logical. Should data be log-transformed before ComBat?
#'   ComBat assumes approximately normal data. If TRUE, applies log2(counts + 1)
#'   before correction and back-transforms after (default = TRUE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return The input list with batch-corrected counts matrix.
#'
#' @details
#' Two methods are available:
#'
#' \strong{ComBat} (sva package):
#' \itemize{
#'   \item Empirical Bayes framework for batch correction
#'   \item Recommended for most applications
#'   \item Can preserve biological groups during correction
#'   \item Assumes approximately normal data (use log_transform = TRUE for counts)
#' }
#'
#' \strong{limma::removeBatchEffect}:
#' \itemize{
#'   \item Linear model-based batch removal
#'   \item Faster than ComBat
#'   \item Also preserves biological groups
#'   \item Works well with log-transformed data
#' }
#'
#' For count data (ATAC-seq, RNA-seq), log transformation is recommended
#' before batch correction since these methods assume approximately normal
#' distributions.
#'
#' @export
#'
#' @examples
#' # Basic batch correction
#' counts_corrected <- POSEIDON_correct_batch(atac_counts)
#'
#' # Preserve group differences while correcting batch
#' counts_corrected <- POSEIDON_correct_batch(atac_counts,
#'                                             batch_col = "batch",
#'                                             group_col = "group")
#'
#' # Use limma instead of ComBat
#' counts_corrected <- POSEIDON_correct_batch(atac_counts, method = "limma")
#'
POSEIDON_correct_batch <- function(quant_result,
                                    batch_col = "batch",
                                    group_col = "group",
                                    method = "combat",
                                    log_transform = TRUE,
                                    verbose = TRUE) {

  # Validate method

  method <- tolower(method)
  if (!method %in% c("combat", "limma")) {
    stop("method must be either 'combat' or 'limma'")
  }

  # Check required packages
  if (method == "combat" && !requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBat batch correction. ",
         "Install with: BiocManager::install('sva')")
  }

  if (method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for limma batch correction. ",
         "Install with: BiocManager::install('limma')")
  }

  # Extract data
  counts <- quant_result$counts
  targets <- quant_result$targets

  # Check batch column exists
  if (!batch_col %in% colnames(targets)) {
    stop("Batch column '", batch_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  batch <- targets[[batch_col]]

  # Check we have multiple batches
  n_batches <- length(unique(batch))
  if (n_batches < 2) {
    warning("Only one batch found. No batch correction needed.")
    return(quant_result)
  }

  if (verbose) {
    cat("Batch correction using", method, "\n")
    cat("  Batches:", paste(unique(batch), collapse = ", "), "\n")
    cat("  Samples per batch:\n")
    batch_counts <- table(batch)
    for (b in names(batch_counts)) {
      cat("    Batch", b, ":", batch_counts[b], "samples\n")
    }
  }

  # Check for group to preserve
  mod <- NULL
  if (!is.null(group_col) && group_col %in% colnames(targets)) {
    group <- factor(targets[[group_col]])
    mod <- model.matrix(~ group)
    if (verbose) {
      cat("  Preserving group:", group_col, "\n")
      cat("    Groups:", paste(levels(group), collapse = ", "), "\n")
    }
  }

  # Log transform if requested
  if (log_transform) {
    if (verbose) cat("  Log-transforming data (log2(x + 1))...\n")
    counts_input <- log2(counts + 1)
  } else {
    counts_input <- counts
  }

  # Apply batch correction
  if (verbose) cat("  Running batch correction...\n")

  if (method == "combat") {
    # ComBat batch correction
    counts_corrected <- sva::ComBat(
      dat = counts_input,
      batch = batch,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
  } else {
    # limma batch correction
    if (!is.null(mod)) {
      # Remove batch effect while preserving group
      design <- mod
      counts_corrected <- limma::removeBatchEffect(
        counts_input,
        batch = batch,
        design = design
      )
    } else {
      counts_corrected <- limma::removeBatchEffect(
        counts_input,
        batch = batch
      )
    }
  }

  # Back-transform if we log-transformed
  if (log_transform) {
    if (verbose) cat("  Back-transforming (2^x - 1)...\n")
    counts_corrected <- 2^counts_corrected - 1
    # Ensure no negative values (can happen with correction)
    counts_corrected[counts_corrected < 0] <- 0
  }

  if (verbose) {
    cat("Batch correction complete.\n")
    cat("  Output dimensions:", nrow(counts_corrected), "x",
        ncol(counts_corrected), "\n")
  }

  # Update result
  quant_result$counts <- counts_corrected
  quant_result$batch_corrected <- TRUE
  quant_result$batch_method <- method

  return(quant_result)
}


#' Visualize Batch Effects (Before/After Correction)
#'
#' @description Creates PCA plots to visualize batch effects before and
#' optionally after batch correction.
#'
#' @param quant_result A list containing counts and targets.
#' @param batch_col Character string. Column name for batch (default = "batch").
#' @param group_col Character string. Column name for group (default = "group").
#' @param log_transform Logical. Log-transform counts for PCA (default = TRUE).
#' @param title Character string. Plot title (default = "PCA - Batch Effect").
#'
#' @return A ggplot object showing PCA colored by batch and shaped by group.
#'
#' @export
#'
POSEIDON_plot_batch_pca <- function(quant_result,
                                     batch_col = "batch",
                                     group_col = "group",
                                     log_transform = TRUE,
                                     title = "PCA - Batch Effect") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  counts <- quant_result$counts
  targets <- quant_result$targets

  # Log transform if requested
  if (log_transform) {
    counts <- log2(counts + 1)
  }

  # Run PCA on transposed data (samples as rows)
  pca_result <- prcomp(t(counts), scale. = TRUE, center = TRUE)

  # Create plot data
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample_id = rownames(pca_result$x)
  )

  # Add batch and group info
  pca_data$batch <- factor(targets[[batch_col]])

  if (!is.null(group_col) && group_col %in% colnames(targets)) {
    pca_data$group <- factor(targets[[group_col]])
  } else {
    pca_data$group <- factor("all")
  }

  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

  # Create plot
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2,
                                               color = batch,
                                               shape = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      color = "Batch",
      shape = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}


#' Filter Samples by Group
#'
#' @description Filters a count dataset to include only samples belonging
#' to specified group(s). Works with the standard quant_result structure.
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata
#'   }
#' @param group_col Character string. Column name in targets containing
#'   group information (default = "group").
#' @param groups Character vector. Group value(s) to keep.
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered quant_result with only samples from specified group(s).
#'
#' @export
#'
#' @examples
#' # Filter for WT samples only
#' wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")
#'
#' # Filter for multiple groups
#' subset_data <- POSEIDON_filter_by_group(data, groups = c("WT", "Control"))
#'
POSEIDON_filter_by_group <- function(quant_result,
                                      group_col = "group",
                                      groups,
                                      verbose = TRUE) {

  targets <- quant_result$targets
  counts <- quant_result$counts

  # Check group column exists
 if (!group_col %in% colnames(targets)) {
    stop("Group column '", group_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  # Find samples in specified groups
  keep_samples <- targets[[group_col]] %in% groups
  n_before <- nrow(targets)
  n_after <- sum(keep_samples)

  if (n_after == 0) {
    stop("No samples found in group(s): ", paste(groups, collapse = ", "))
  }

  if (verbose) {
    cat("Filtering by group:\n")
    cat("  Keeping groups:", paste(groups, collapse = ", "), "\n")
    cat("  Samples before:", n_before, "\n")
    cat("  Samples after:", n_after, "\n")
  }

  # Filter counts and targets
  quant_result$counts <- counts[, keep_samples, drop = FALSE]
  quant_result$targets <- targets[keep_samples, , drop = FALSE]

  # Update annotation if present (no change needed, features stay the same)

  return(quant_result)
}


#' Filter Samples by Metadata Value
#'
#' @description General function to filter samples by any metadata column.
#' More flexible than filter_by_group.
#'
#' @param quant_result A list containing counts and targets.
#' @param column Character string. Column name in targets to filter by.
#' @param values Vector. Value(s) to keep.
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered quant_result.
#'
#' @export
#'
#' @examples
#' # Filter by batch
#' batch1_data <- POSEIDON_filter_by_metadata(data, column = "batch", values = 1)
#'
#' # Filter by biological replicate
#' bio1_data <- POSEIDON_filter_by_metadata(data, column = "bio_rep", values = c(1, 2))
#'
POSEIDON_filter_by_metadata <- function(quant_result,
                                         column,
                                         values,
                                         verbose = TRUE) {

  targets <- quant_result$targets
  counts <- quant_result$counts

  if (!column %in% colnames(targets)) {
    stop("Column '", column, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  keep_samples <- targets[[column]] %in% values
  n_before <- nrow(targets)
  n_after <- sum(keep_samples)

  if (n_after == 0) {
    stop("No samples found with ", column, " in: ", paste(values, collapse = ", "))
  }

  if (verbose) {
    cat("Filtering by", column, ":\n")
    cat("  Keeping values:", paste(values, collapse = ", "), "\n")
    cat("  Samples before:", n_before, "\n")
    cat("  Samples after:", n_after, "\n")
  }

  quant_result$counts <- counts[, keep_samples, drop = FALSE]
  quant_result$targets <- targets[keep_samples, , drop = FALSE]

  return(quant_result)
}


#' Match Samples Across Datasets
#'
#' @description Identifies and filters datasets to include only samples
#' that are present in all datasets, based on a matching column (e.g., bio_rep).
#'
#' @param data_list Named list of quant_result objects (each with counts and targets).
#' @param match_col Character string. Column in targets to match on
#'   (default = "bio_rep").
#' @param verbose Logical. Print matching summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{datasets}{Named list of filtered quant_results with matched samples}
#'   \item{matched_values}{Vector of matched values (e.g., bio_rep IDs)}
#'   \item{sample_map}{Data.frame showing sample correspondence across datasets}
#' }
#'
#' @details
#' This function finds the intersection of values in match_col across all
#' datasets and filters each dataset to include only samples with those values.
#'
#' Note: This does NOT reorder samples. Use POSEIDON_align_samples() after
#' matching to ensure consistent sample ordering across datasets.
#'
#' @export
#'
#' @examples
#' data_list <- list(
#'   ATAC = atac_counts,
#'   ChIP = chip_counts,
#'   RNA = rna_data
#' )
#' matched <- POSEIDON_match_samples(data_list, match_col = "bio_rep")
#'
POSEIDON_match_samples <- function(data_list,
                                    match_col = "bio_rep",
                                    verbose = TRUE) {

  if (!is.list(data_list) || length(data_list) < 2) {
    stop("data_list must be a list with at least 2 datasets")
  }

  dataset_names <- names(data_list)
  if (is.null(dataset_names)) {
    dataset_names <- paste0("Dataset_", seq_along(data_list))
    names(data_list) <- dataset_names
  }

  # Get values of match_col for each dataset
  values_per_dataset <- lapply(data_list, function(d) {
    if (!match_col %in% colnames(d$targets)) {
      stop("Column '", match_col, "' not found in targets")
    }
    unique(d$targets[[match_col]])
  })

  if (verbose) {
    cat("Sample matching by:", match_col, "\n\n")
    cat("Values per dataset:\n")
    for (nm in dataset_names) {
      cat("  ", nm, ":", paste(sort(values_per_dataset[[nm]]), collapse = ", "), "\n")
    }
  }

  # Find intersection
  matched_values <- Reduce(intersect, values_per_dataset)

  if (length(matched_values) == 0) {
    stop("No matching values found across all datasets")
  }

  if (verbose) {
    cat("\nMatched values:", paste(sort(matched_values), collapse = ", "), "\n")
    cat("Number of matched", match_col, ":", length(matched_values), "\n\n")
  }

  # Filter each dataset to matched values
  filtered_list <- list()

  for (nm in dataset_names) {
    filtered_list[[nm]] <- POSEIDON_filter_by_metadata(
      data_list[[nm]],
      column = match_col,
      values = matched_values,
      verbose = verbose
    )
    if (verbose) cat("\n")
  }

  # Create sample map showing correspondence
  sample_map <- data.frame(match_value = matched_values)
  colnames(sample_map) <- match_col

  for (nm in dataset_names) {
    targets <- filtered_list[[nm]]$targets
    # Get sample counts per match value
    samples_per_value <- sapply(matched_values, function(v) {
      sum(targets[[match_col]] == v)
    })
    sample_map[[paste0(nm, "_n")]] <- samples_per_value
  }

  if (verbose) {
    cat("Sample map:\n")
    print(sample_map)
  }

  return(list(
    datasets = filtered_list,
    matched_values = matched_values,
    sample_map = sample_map
  ))
}


#' Average Technical Replicates
#'
#' @description Aggregates technical replicates by averaging counts within
#' each biological replicate. Produces one sample per unique combination of
#' grouping variables (excluding tech_rep).
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata
#'   }
#' @param bio_rep_col Character string. Column identifying biological replicates
#'   (default = "bio_rep").
#' @param tech_rep_col Character string. Column identifying technical replicates
#'   to average over (default = "tech_rep").
#' @param group_by Character vector. Additional columns to group by when
#'   averaging. Samples are averaged within each unique combination of
#'   bio_rep_col + group_by columns. Default includes "group" and "batch"
#'   if present in targets.
#' @param method Character string. Aggregation method: "mean" (default) or
#'   "median".
#' @param verbose Logical. Print summary (default = TRUE).
#'
#' @return A quant_result with:
#'   \itemize{
#'     \item counts: Averaged count matrix (features x biological samples)
#'     \item targets: Updated metadata with one row per biological sample
#'     \item n_tech_reps: Number of technical replicates averaged per sample
#'   }
#'
#' @details
#' This function averages technical replicates to produce one observation per
#' biological sample. This is appropriate when:
#' \itemize{
#'   \item Preparing data for integration methods requiring matched samples
#'   \item Reducing technical noise while preserving biological variation
#'   \item Creating a cleaner dataset for correlation-based analyses
#' }
#'
#' The function automatically detects which columns to use for grouping:
#' \itemize{
#'   \item Always groups by bio_rep_col
#'   \item Includes "group" if present (preserves experimental groups)
#'   \item Includes "batch" if present (keeps batches separate)
#'   \item Additional columns can be specified via group_by parameter
#' }
#'
#' Can be used before or after filtering by group - the function respects
#' whatever samples are present in the input.
#'
#' @export
#'
#' @examples
#' # Basic averaging
#' averaged <- POSEIDON_average_tech_reps(atac_counts)
#'
#' # Average after filtering for WT
#' wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")
#' wt_averaged <- POSEIDON_average_tech_reps(wt_data)
#'
#' # Average before filtering (groups will be preserved)
#' averaged <- POSEIDON_average_tech_reps(atac_counts)
#' wt_averaged <- POSEIDON_filter_by_group(averaged, groups = "WT")
#'
#' # Custom grouping
#' averaged <- POSEIDON_average_tech_reps(data, group_by = c("group", "batch", "treatment"))
#'
POSEIDON_average_tech_reps <- function(quant_result,
                                        bio_rep_col = "bio_rep",
                                        tech_rep_col = "tech_rep",
                                        group_by = NULL,
                                        method = "mean",
                                        verbose = TRUE) {

  counts <- quant_result$counts
  targets <- quant_result$targets

  # Validate columns exist
  if (!bio_rep_col %in% colnames(targets)) {
    stop("bio_rep_col '", bio_rep_col, "' not found in targets")
  }

  if (!tech_rep_col %in% colnames(targets)) {
    stop("tech_rep_col '", tech_rep_col, "' not found in targets")
  }

  # Determine grouping columns
  # Always include bio_rep_col, plus any standard columns that exist
  standard_group_cols <- c("group", "batch")
  auto_group_by <- standard_group_cols[standard_group_cols %in% colnames(targets)]

  if (is.null(group_by)) {
    group_cols <- c(bio_rep_col, auto_group_by)
  } else {
    group_cols <- unique(c(bio_rep_col, group_by))
  }

  # Remove tech_rep_col from grouping if accidentally included
  group_cols <- setdiff(group_cols, tech_rep_col)

  if (verbose) {
    cat("Averaging technical replicates:\n")
    cat("  Grouping by:", paste(group_cols, collapse = ", "), "\n")
    cat("  Averaging over:", tech_rep_col, "\n")
    cat("  Method:", method, "\n")
  }

  # Create grouping key for each sample
  group_key <- apply(targets[, group_cols, drop = FALSE], 1, paste, collapse = "_")

  # Get unique groups
  unique_groups <- unique(group_key)
  n_groups <- length(unique_groups)

  if (verbose) {
    cat("  Samples before:", ncol(counts), "\n")
    cat("  Samples after:", n_groups, "\n")
  }

  # Aggregate counts
  agg_func <- if (method == "mean") rowMeans else function(x) apply(x, 1, median)

  agg_counts <- matrix(
    nrow = nrow(counts),
    ncol = n_groups,
    dimnames = list(rownames(counts), NULL)
  )

  agg_targets <- data.frame(matrix(
    nrow = n_groups,
    ncol = length(group_cols),
    dimnames = list(NULL, group_cols)
  ), stringsAsFactors = FALSE)

  n_tech_reps <- numeric(n_groups)

  for (i in seq_along(unique_groups)) {
    grp <- unique_groups[i]
    sample_idx <- which(group_key == grp)

    # Average counts
    if (length(sample_idx) == 1) {
      agg_counts[, i] <- counts[, sample_idx]
    } else {
      agg_counts[, i] <- agg_func(counts[, sample_idx, drop = FALSE])
    }

    # Take first row's metadata for grouping columns
    agg_targets[i, ] <- targets[sample_idx[1], group_cols]

    # Record number of tech reps
    n_tech_reps[i] <- length(sample_idx)
  }

  # Create sample IDs for averaged data
  new_sample_ids <- apply(agg_targets, 1, paste, collapse = "_")
  colnames(agg_counts) <- new_sample_ids
  rownames(agg_targets) <- new_sample_ids
  agg_targets$sample_id <- new_sample_ids
  agg_targets$n_tech_reps <- n_tech_reps

  if (verbose) {
    cat("\nTechnical replicates per biological sample:\n")
    tech_rep_summary <- table(n_tech_reps)
    for (n in names(tech_rep_summary)) {
      cat("  ", tech_rep_summary[n], "sample(s) with", n, "tech rep(s)\n")
    }
  }

  # Build result
  result <- list(
    counts = agg_counts,
    targets = agg_targets
  )

  # Preserve annotation if present
  if (!is.null(quant_result$annotation)) {
    result$annotation <- quant_result$annotation
  }

  # Add metadata
  result$averaged <- TRUE
  result$averaging_method <- method
  result$group_cols <- group_cols

  return(result)
}


#' Filter Features with Pairwise Identical Values (Debug Filter)
#'
#' @description Aggressively filters features where ANY pair of samples has
#' identical values. Such features will produce NaN when scaled in CV folds
#' containing that pair (due to zero variance → division by zero).
#'
#' @param X A matrix (samples x features) or a list of matrices.
#' @param tolerance Numeric. Values within this tolerance are considered
#'   identical (default = 1e-10).
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered matrix or list of matrices with problematic features removed.
#'
#' @details
#' This is an AGGRESSIVE filter intended for debugging small-sample scenarios.
#' It removes features where any two samples have identical (or near-identical)
#' values, which would cause scaling to fail in CV when those samples end up
#' in the same training fold.
#'
#' For a feature to pass this filter, ALL pairwise comparisons between samples
#' must show some variance. This can remove a substantial number of features,
#' especially with:
#' \itemize{
#'   \item Count data with many zeros
#'   \item Averaged technical replicates that converge to similar values
#'   \item Small sample sizes
#' }
#'
#' WARNING: This filter is NOT recommended for production analyses as it
#' may remove biologically meaningful features. Use only for debugging or
#' when small sample sizes make CV scaling problematic.
#'
#' @export
#'
#' @examples
#' # Filter a single matrix
#' X_filtered <- POSEIDON_filter_pairwise_identical(X)
#'
#' # Filter a list of matrices (e.g., for multi-omics)
#' data_filtered <- POSEIDON_filter_pairwise_identical(scca_data)
#'
POSEIDON_filter_pairwise_identical <- function(X,
                                                tolerance = 1e-10,
                                                verbose = TRUE) {

  # Handle list of matrices
 if (is.list(X) && !is.data.frame(X)) {
    if (verbose) cat("Filtering pairwise identical features across", length(X), "datasets:\n")

    result <- lapply(names(X), function(nm) {
      if (verbose) cat("\n  Dataset:", nm, "\n")
      POSEIDON_filter_pairwise_identical(X[[nm]], tolerance = tolerance, verbose = verbose)
    })
    names(result) <- names(X)
    return(result)
  }

  # Single matrix processing
  n_samples <- nrow(X)
  n_features <- ncol(X)

  if (n_samples < 2) {
    warning("Need at least 2 samples to check pairwise values")
    return(X)
  }

  # For each feature, check if any pair of samples has identical values
  # A feature is problematic if min(pairwise differences) == 0 for any pair

  keep_feature <- apply(X, 2, function(col) {
    # Check all pairwise differences
    for (i in 1:(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        if (abs(col[i] - col[j]) <= tolerance) {
          return(FALSE)  # Found identical pair, remove this feature
        }
      }
    }
    return(TRUE)  # All pairs are different, keep this feature
  })

  n_removed <- sum(!keep_feature)
  n_kept <- sum(keep_feature)

  if (verbose) {
    cat(sprintf("    Features before: %d\n", n_features))
    cat(sprintf("    Features removed (pairwise identical): %d (%.1f%%)\n",
                n_removed, 100 * n_removed / n_features))
    cat(sprintf("    Features kept: %d\n", n_kept))
  }

  if (n_kept == 0) {
    warning("All features removed! Every feature has at least one pair of identical values.")
  }

  return(X[, keep_feature, drop = FALSE])
}
