###############################################################################
########### Feature Filtering & Scaling for Multi-Dataset Integration #######
###############################################################################
# Generic matrix / list-of-matrices operations (no quant_result structure
# assumed) used to prepare features for integration methods like sCCA.

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
            cat(sprintf("[POSEIDON] Removing %d features with variance <= %.2e\n",
                       sum(!keep_features), variance_threshold))
            cat(sprintf("[POSEIDON] Keeping %d features with variance > %.2e\n",
                       sum(keep_features), variance_threshold))
        }

        return(X[, keep_features, drop = FALSE])
    }

    # Handle list of matrices
    if (verbose) {
        cat("[POSEIDON] Filtering features across multiple datasets...\n")
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
        cat(sprintf("[POSEIDON] Filtering features with variance > %.2e across %d datasets:\n",
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
#' \dontrun{
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
#' }
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
      cat(sprintf("[POSEIDON] %s%d -> %d features (top %.1f%% by variance)\n",
                  prefix, n_features, n, 100 * n / n_features))
      cat(sprintf("[POSEIDON] %s  Variance range kept: %.2e to %.2e\n",
                  prefix, min(vars[top_idx]), max(vars[top_idx])))
    }

    return(mat[, top_idx, drop = FALSE])
  }

  # Handle single matrix
  if (!is.list(X)) {
    if (verbose) {
      cat("[POSEIDON] Filtering to top", n_top, "variable features...\n")
    }
    return(filter_matrix(X, n_top))
  }

  # Handle list of matrices
  if (verbose) {
    cat("[POSEIDON] Filtering to top variable features per dataset...\n")
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
    cat("[POSEIDON] Done.\n")
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
#' \dontrun{
#' # Remove bottom 50% of features by variance
#' X_filtered <- POSEIDON_filter_low_variance(X, bottom_pct = 0.5)
#'
#' # For sCCA with multiple datasets
#' scca_data <- list(ATAC = atac_matrix, ChIP = chip_matrix, RNA = rna_matrix)
#' scca_filtered <- POSEIDON_filter_low_variance(scca_data, bottom_pct = 0.5)
#'
#' }
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
      cat(sprintf("[POSEIDON] %s%d -> %d features (removed bottom %.0f%%)\n",
                  prefix, n_features, n_keep, pct * 100))
      cat(sprintf("[POSEIDON] %s  Variance threshold: %.2e (features below this removed)\n",
                  prefix, var_threshold))
    }

    return(mat[, keep_idx, drop = FALSE])
  }

  # Handle single matrix
  if (!is.list(X)) {
    if (verbose) {
      cat(sprintf("[POSEIDON] Removing bottom %.0f%% of features by variance...\n", bottom_pct * 100))
    }
    return(filter_matrix(X, bottom_pct))
  }

  # Handle list of matrices
  if (verbose) {
    cat(sprintf("[POSEIDON] Removing bottom %.0f%% of features by variance per dataset...\n",
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
    cat("[POSEIDON] Done.\n")
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
#' \dontrun{
#' # Filter a single matrix
#' X_filtered <- POSEIDON_filter_pairwise_identical(X)
#'
#' # Filter a list of matrices (e.g., for multi-omics)
#' data_filtered <- POSEIDON_filter_pairwise_identical(scca_data)
#'
#' }
POSEIDON_filter_pairwise_identical <- function(X,
                                                tolerance = 1e-10,
                                                verbose = TRUE) {

  # Handle list of matrices
 if (is.list(X) && !is.data.frame(X)) {
    if (verbose) cat("[POSEIDON] Filtering pairwise identical features across", length(X), "datasets:\n")

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
    warning("[POSEIDON] All features removed! Every feature has at least one pair of identical values.")
  }

  return(X[, keep_feature, drop = FALSE])
}
