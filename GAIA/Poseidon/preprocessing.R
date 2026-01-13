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
