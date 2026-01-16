###############################################################################
########### Cross-Validation Framework for Sparse CCA ###########
###############################################################################

# Import dependencies from other GAIA modules using box
source('~/A_Projects/ZERO_DAWN/GAIA/Hephaestus/ConvCCA.R')
source('~/A_Projects/ZERO_DAWN/GAIA/Hephaestus/RelPMDCCA.R')

#' Create K-Fold Cross-Validation Indices
#'
#' @description Creates k-fold CV splits, ensuring same samples are held out
#' across all datasets (maintains sample pairing)
#'
#' @param n_samples Number of samples
#' @param k Number of folds (default = 5)
#' @param seed Random seed for reproducibility
#'
#' @return List of k vectors, each containing test indices for that fold
#'
#' @export
MINERVA_create_folds <- function(n_samples, k = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Shuffle sample indices
  indices <- sample(1:n_samples)

  # Split into k roughly equal folds
  folds <- split(indices, cut(seq_along(indices), breaks = k, labels = FALSE))

  return(folds)
}


#' Fit Sparse CCA Model (Unified Interface)
#'
#' @description Wrapper function that fits any sCCA method with consistent interface
#'
#' @param X List of datasets (for multi-dataset) or use X_1, X_2 for two datasets
#' @param X_1 First dataset (alternative to X for two-dataset methods)
#' @param X_2 Second dataset (alternative to X for two-dataset methods)
#' @param method Method to use: "ConvCCA", "RelPMDCCA", or "PMDCCA"
#' @param tau List of tau values (one per dataset) or tauW_1, tauW_2 for two datasets
#' @param lambda Lambda parameter (RelPMDCCA only, default = 10)
#' @param nIter Maximum iterations
#' @param penalty Penalty function: "LASSO" (default) or "SCAD"
#' @param ... Additional arguments passed to the method
#'
#' @return List with canonical vectors W (list of vectors, one per dataset)
#'
#' @export
MINERVA_scca_fit <- function(X = NULL, X_1 = NULL, X_2 = NULL,
                     method = c("ConvCCA", "RelPMDCCA", "PMDCCA"),
                     tau = NULL, lambda = 10, nIter = 100, penalty = "LASSO", ...) {

  method <- match.arg(method)

  # Determine if multi-dataset or two-dataset
  is_multi <- !is.null(X)
  n_datasets <- if (is_multi) length(X) else 2

  # Fit the appropriate method
  if (method == "ConvCCA") {
    if (is_multi) {
      result <- HEPHAESTUS_multi_convCCA(X = X, tau = tau, nIter = nIter,
                              penalty = penalty, ...)
      return(list(W = result$W))
    } else {
      result <- HEPHAESTUS_convCCA(X_1 = X_1, X_2 = X_2,
                        tauW_1 = tau[[1]], tauW_2 = tau[[2]],
                        nIter = nIter, penalty = penalty, ...)
      return(list(W = list(result$W_1, result$W_2)))
    }

  } else if (method == "RelPMDCCA") {
    if (is_multi) {
      result <- HEPHAESTUS_multi_relPMDCCA(X = X, lambda = lambda, tau = tau,
                                nIter = nIter, penalty = penalty, ...)
      return(list(W = result$W))
    } else {
      result <- HEPHAESTUS_relPMDCCA(X_1 = X_1, X_2 = X_2, lambda = lambda,
                          tauW_1 = tau[[1]], tauW_2 = tau[[2]],
                          nIter = nIter, penalty = penalty, ...)
      return(list(W = list(result$W_1, result$W_2)))
    }

  } else if (method == "PMDCCA") {
    # PMDCCA from PMA package (two datasets only)
    if (is_multi) {
      stop("PMDCCA only supports two datasets. Use ConvCCA or RelPMDCCA for multi-dataset.")
    }
    if (!requireNamespace("PMA", quietly = TRUE)) {
      stop("Package 'PMA' is required for PMDCCA. Install with: install.packages('PMA')")
    }
    result <- PMA::CCA(x = X_1, z = X_2,
                       typex = "standard", typez = "standard",
                       penaltyx = tau[[1]], penaltyz = tau[[2]],
                       K = 1, ...)
    return(list(W = list(result$u, result$v)))

  } else {
    stop("Unknown method: ", method)
  }
}


#' Evaluate Sparse CCA Model
#'
#' @description Computes evaluation metrics for sCCA results
#'
#' @param X_list List of datasets (test data)
#' @param W_list List of canonical vectors
#' @param metric Metric to compute: "correlation" (default) or "both"
#'   (both computes correlation and sparsity metrics)
#'
#' @return Named list of metrics
#'
#' @export
MINERVA_evaluate_scca <- function(X_list, W_list, metric = "correlation") {

  # Validate metric parameter
  if (!metric %in% c("correlation", "both")) {
    stop("metric must be either 'correlation' or 'both'")
  }

  n_datasets <- length(X_list)

  # Compute canonical scores for each dataset
  scores <- lapply(1:n_datasets, function(i) {
    X_list[[i]] %*% W_list[[i]]
  })

  results <- list()

  # Canonical correlation (average pairwise correlation)
  if (metric %in% c("correlation", "both")) {
    cors <- c()
    for (i in 1:(n_datasets - 1)) {
      for (j in (i + 1):n_datasets) {
        cors <- c(cors, cor(scores[[i]], scores[[j]]))
      }
    }
    results$correlation <- mean(cors)
    results$correlation_sd <- sd(cors)
    results$correlations_all <- cors
  }

  # Sparsity metrics (only when metric="both")
  if (metric == "both") {
    sparsity <- sapply(W_list, function(w) {
      sum(abs(w) > 1e-6) / length(w)  # Proportion of non-zero features
    })
    results$sparsity <- sparsity
    results$sparsity_mean <- mean(sparsity)
  }

  return(results)
}


#' K-Fold Cross-Validation for Sparse CCA
#'
#' @description Performs k-fold cross-validation to select optimal tau parameters
#'
#' @param X List of datasets (for multi-dataset) or use X_1, X_2 for two datasets
#' @param X_1 First dataset (alternative for two-dataset methods)
#' @param X_2 Second dataset (alternative for two-dataset methods)
#' @param method Method to use: "ConvCCA", "RelPMDCCA", or "PMDCCA"
#' @param tau_grid Data frame or list of tau combinations to test.
#'   For multi-dataset: data.frame(tau1 = ..., tau2 = ..., tau3 = ...)
#'   For two-dataset: data.frame(tau1 = ..., tau2 = ...)
#' @param lambda Lambda parameter (RelPMDCCA only, default = 10)
#' @param k Number of CV folds (default = 5)
#' @param nIter Maximum iterations for each fit
#' @param penalty Penalty function: "LASSO" (default) or "SCAD"
#' @param metric Metric to use: "correlation" (default, optimizes for correlation)
#'   or "both" (computes both correlation and sparsity, optimizes for correlation)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress (default = TRUE)
#'
#' @return List containing:
#'   \item{best_tau}{Optimal tau values (list)}
#'   \item{best_score}{Best CV score}
#'   \item{cv_results}{Data frame with all CV results}
#'   \item{best_model}{Model fit on full data with best tau}
#'
#' @export
MINERVA_cv_scca <- function(X = NULL, X_1 = NULL, X_2 = NULL,
                    method = c("ConvCCA", "RelPMDCCA", "PMDCCA"),
                    tau_grid,
                    lambda = 10,
                    k = 5,
                    nIter = 100,
                    penalty = "LASSO",
                    metric = "correlation",
                    seed = 123,
                    verbose = TRUE) {

  method <- match.arg(method)

  # Validate metric parameter
  if (!metric %in% c("correlation", "both")) {
    stop("metric must be either 'correlation' or 'both'")
  }

  # Determine dataset structure
  is_multi <- !is.null(X)
  n_datasets <- if (is_multi) length(X) else 2
  n_samples <- if (is_multi) nrow(X[[1]]) else nrow(X_1)

  # Convert tau_grid to data frame if needed
  if (is.list(tau_grid) && !is.data.frame(tau_grid)) {
    tau_grid <- expand.grid(tau_grid)
  }

  # Create CV folds
  folds <- MINERVA_create_folds(n_samples, k = k, seed = seed)

  # Calculate minimum training set size and warn if problematic
  min_train_size <- n_samples - max(sapply(folds, length))

  if (min_train_size <= 2) {
    warning(
      "\n",
      "=======================================================================\n",
      "SMALL SAMPLE WARNING\n",
      "=======================================================================\n",
      sprintf("With %d samples and %d-fold CV, training sets will have only %d samples.\n",
              n_samples, k, min_train_size),
      "\n",
      "This can cause issues:\
",
      "  1. Insufficient statistical power for reliable correlation estimates\n",
      "  2. Features with identical values across training samples will produce\n",
      "     NaN after scaling, causing CCA to fail\n",
      "  3. Results may be highly unstable and not generalizable\n",
      "\n",
      "Consider:\n",
      "  - Using fewer folds (if not already LOOCV)\n",
      "  - Acquiring more samples\n",
      "  - Pre-filtering features with low variance across all samples\n",
      "=======================================================================\n",
      immediate. = TRUE
    )
  }

  if (verbose) {
    cat(sprintf("Running %d-fold CV for %s with %d tau combinations\n",
                k, method, nrow(tau_grid)))
    cat(sprintf("Training samples per fold: %d\n", min_train_size))
    cat(sprintf("Total fits: %d\n\n", k * nrow(tau_grid)))
  }

  # Storage for results
  cv_scores <- matrix(NA, nrow = nrow(tau_grid), ncol = k)
  
  # Grid search over tau combinations
  for (i in 1:nrow(tau_grid)) {
    # Extract tau combination
    if (is.data.frame(tau_grid)){
      # Convert row to list, stripping names and ensuring numeric
      tau_combo <- as.list(as.numeric(unname(tau_grid[i,])))
    } else {
      tau_combo <- as.list(as.numeric(tau_grid[i,]))
    }

    if (verbose && i %% 10 == 0) {
      cat(sprintf("Testing tau combination %d/%d...\n", i, nrow(tau_grid)))
    }

    # K-fold CV for this tau combination
    for (fold_idx in 1:k) {
      test_indices <- folds[[fold_idx]]
      train_indices <- setdiff(1:n_samples, test_indices)

      # Split data
      if (is_multi) {
        X_train <- lapply(X, function(x) x[train_indices, , drop = FALSE])
        X_test <- lapply(X, function(x) x[test_indices, , drop = FALSE])
      } else {
        X_1_train <- X_1[train_indices, , drop = FALSE]
        X_1_test <- X_1[test_indices, , drop = FALSE]
        X_2_train <- X_2[train_indices, , drop = FALSE]
        X_2_test <- X_2[test_indices, , drop = FALSE]
      }

      # Fit model on training data
      tryCatch({
        if (is_multi) {
          fit <- MINERVA_scca_fit(X = X_train, method = method, tau = tau_combo,
                         lambda = lambda, nIter = nIter, penalty = penalty)
        } else {
          fit <- MINERVA_scca_fit(X_1 = X_1_train, X_2 = X_2_train,
                         method = method, tau = tau_combo,
                         lambda = lambda, nIter = nIter, penalty = penalty)
        }

        # Evaluate on test data
        if (is_multi) {
          eval_result <- MINERVA_evaluate_scca(X_test, fit$W, metric = metric)
        } else {
          eval_result <- MINERVA_evaluate_scca(list(X_1_test, X_2_test), fit$W, metric = metric)
        }

        # Extract score for CV optimization
        # When metric="both", still optimize based on correlation
        cv_metric <- if (metric == "both") "correlation" else metric
        cv_scores[i, fold_idx] <- eval_result[[cv_metric]]

      }, error = function(e) {
        if (verbose) {
          cat(sprintf("  Error in fold %d: %s\n", fold_idx, e$message))
        }
        cv_scores[i, fold_idx] <- NA
      })
    }
  }

  # Compute mean CV score for each tau combination
  mean_scores <- rowMeans(cv_scores, na.rm = TRUE)
  sd_scores <- apply(cv_scores, 1, sd, na.rm = TRUE)

  # Check if we got any valid scores
  if (all(is.na(mean_scores)) || all(is.nan(mean_scores))){
    stop("All CV scores are NA or NaN. Check if model fitting is failing for all folds.")
  }

  # Find best tau combination
  best_idx <- which.max(mean_scores)

  if (length(best_idx) == 0){
    stop("which.max returned empty index. mean_scores:", paste(mean_scores, collapse=", "))
  }

  # Extract best tau - handle single row properly
  if (is.data.frame(tau_grid)){
    best_tau_row <- tau_grid[best_idx, , drop = FALSE]
    best_tau <- as.list(as.numeric(best_tau_row[1, ]))
  } else {
    best_tau <- as.list(as.numeric(tau_grid[best_idx, ]))
  }

  best_score <- mean_scores[best_idx]

  if (verbose) {
    cat("\n=== CV Results ===\n")
    cat(sprintf("Best %s: %.4f (± %.4f)\n", metric, best_score, sd_scores[best_idx]))
    cat("Best tau values:\n")
    for (j in 1:length(best_tau)) {
      cat(sprintf("  Dataset %d: %.3f\n", j, best_tau[[j]]))
    }
  }

  # Fit final model on full data with best tau
  if (verbose) cat("\nFitting final model on full data...\n")

  if (is_multi) {
    final_model <- MINERVA_scca_fit(X = X, method = method, tau = best_tau,
                           lambda = lambda, nIter = nIter, penalty = penalty)
  } else {
    final_model <- MINERVA_scca_fit(X_1 = X_1, X_2 = X_2, method = method, tau = best_tau,
                           lambda = lambda, nIter = nIter, penalty = penalty)
  }

  # Compile results
  cv_results_df <- cbind(tau_grid,
                         mean_score = mean_scores,
                         sd_score = sd_scores)
  cv_results_df <- cv_results_df[order(-cv_results_df$mean_score), ]

  return(list(
    best_tau = best_tau,
    best_score = best_score,
    cv_results = cv_results_df,
    best_model = final_model,
    method = method,
    metric = metric
  ))
}


#' Compute Sparsity Quality Score
#'
#' @description Scores sparsity based on how well it falls within ideal range.
#' Penalizes both extreme sparsity (too few features) and lack of sparsity (too many features).
#'
#' @param sparsity Numeric. Proportion of non-zero features (0 to 1)
#' @param ideal_min Numeric. Minimum ideal sparsity (default = 0.05, or 5%)
#' @param ideal_max Numeric. Maximum ideal sparsity (default = 0.30, or 30%)
#' @param penalty_steepness Numeric. How quickly to penalize outside ideal range (default = 2)
#'
#' @return Numeric quality score between 0 and 1, where 1 = ideal sparsity
#'
#' @details
#' The quality score is:
#' \itemize{
#'   \item 1.0 if sparsity is within [ideal_min, ideal_max]
#'   \item Decreases as sparsity moves outside this range
#'   \item Uses exponential decay based on penalty_steepness
#' }
#'
#' @export
MINERVA_sparsity_quality <- function(sparsity,
                                      ideal_min = 0.05,
                                      ideal_max = 0.30,
                                      penalty_steepness = 2) {

  if (sparsity >= ideal_min && sparsity <= ideal_max) {
    # Within ideal range
    return(1.0)
  } else if (sparsity < ideal_min) {
    # Too sparse (too few features)
    distance <- ideal_min - sparsity
    return(exp(-penalty_steepness * distance / ideal_min))
  } else {
    # Not sparse enough (too many features)
    distance <- sparsity - ideal_max
    max_distance <- 1 - ideal_max
    return(exp(-penalty_steepness * distance / max_distance))
  }
}


#' Pilot Comparison of sCCA Methods
#'
#' @description Quick comparison of all sCCA methods with default tau values
#' to determine which method works best for the data
#'
#' @param X List of datasets (for multi-dataset) or use X_1, X_2
#' @param X_1 First dataset (for two-dataset)
#' @param X_2 Second dataset (for two-dataset)
#' @param test_proportion Proportion of data to hold out for testing (default = 0.2)
#' @param default_tau_convCCA Default tau for ConvCCA (default = 0.3 for all datasets)
#' @param default_tau_relPMDCCA Default tau for RelPMDCCA (default = 0.7 for all datasets)
#' @param lambda Lambda for RelPMDCCA (default = 10)
#' @param nIter Iterations for each method
#' @param seed Random seed
#' @param use_composite_score Use composite score combining correlation and sparsity quality (default = TRUE)
#' @param ideal_sparsity_min Minimum ideal sparsity (default = 0.05, or 5%)
#' @param ideal_sparsity_max Maximum ideal sparsity (default = 0.30, or 30%)
#' @param correlation_weight Weight for correlation in composite score (default = 0.7)
#' @param verbose Print results (default = TRUE)
#'
#' @return List containing:
#'   \item{best_method}{Name of best-performing method}
#'   \item{results}{Data frame with results for all methods}
#'   \item{models}{List of fitted models}
#'
#' @export
MINERVA_pilot_compare_methods <- function(X = NULL, X_1 = NULL, X_2 = NULL,
                                  test_proportion = 0.2,
                                  default_tau_convCCA = NULL,
                                  default_tau_relPMDCCA = NULL,
                                  lambda = 10,
                                  nIter = 100,
                                  seed = 123,
                                  use_composite_score = TRUE,
                                  ideal_sparsity_min = 0.05,
                                  ideal_sparsity_max = 0.30,
                                  correlation_weight = 0.7,
                                  verbose = TRUE) {

  set.seed(seed)

  # Determine structure
  is_multi <- !is.null(X)
  n_datasets <- if (is_multi) length(X) else 2
  n_samples <- if (is_multi) nrow(X[[1]]) else nrow(X_1)

  # Set default taus if not provided
  if (is.null(default_tau_convCCA)) {
    default_tau_convCCA <- as.list(rep(0.1, n_datasets))
  }
  if (is.null(default_tau_relPMDCCA)) {
    default_tau_relPMDCCA <- as.list(rep(0.9, n_datasets))
  }

  # Create train/test split
  n_test <- floor(n_samples * test_proportion)
  test_indices <- sample(1:n_samples, n_test)
  train_indices <- setdiff(1:n_samples, test_indices)

  if (is_multi) {
    X_train <- lapply(X, function(x) x[train_indices, , drop = FALSE])
    X_test <- lapply(X, function(x) x[test_indices, , drop = FALSE])
  } else {
    X_1_train <- X_1[train_indices, , drop = FALSE]
    X_1_test <- X_1[test_indices, , drop = FALSE]
    X_2_train <- X_2[train_indices, , drop = FALSE]
    X_2_test <- X_2[test_indices, , drop = FALSE]
  }

  if (verbose) {
    cat("=== Pilot Comparison of sCCA Methods ===\n")
    cat(sprintf("Training on %d samples, testing on %d samples\n\n",
                length(train_indices), length(test_indices)))
  }

  # Test each method
  methods <- c("ConvCCA", "RelPMDCCA")
  results <- data.frame(
    method = character(),
    correlation = numeric(),
    sparsity_mean = numeric(),
    sparsity_quality = numeric(),
    composite_score = numeric(),
    time_sec = numeric(),
    stringsAsFactors = FALSE
  )

  models <- list()

  for (method in methods) {
    if (verbose) cat(sprintf("Testing %s...\n", method))

    # Choose appropriate tau
    tau <- if (method == "ConvCCA") default_tau_convCCA else default_tau_relPMDCCA

    start_time <- Sys.time()

    tryCatch({
      # Fit model
      if (is_multi) {
        fit <- MINERVA_scca_fit(X = X_train, method = method, tau = tau,
                       lambda = lambda, nIter = nIter, penalty = "LASSO")
        eval_result <- MINERVA_evaluate_scca(X_test, fit$W, metric = "both")
      } else {
        fit <- MINERVA_scca_fit(X_1 = X_1_train, X_2 = X_2_train, method = method,
                       tau = tau, lambda = lambda, nIter = nIter, penalty = "LASSO")
        eval_result <- MINERVA_evaluate_scca(list(X_1_test, X_2_test), fit$W, metric = "both")
      }

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      # Compute sparsity quality score
      spars_quality <- MINERVA_sparsity_quality(
        eval_result$sparsity_mean,
        ideal_min = ideal_sparsity_min,
        ideal_max = ideal_sparsity_max
      )

      # Compute composite score
      sparsity_weight <- 1 - correlation_weight
      comp_score <- correlation_weight * eval_result$correlation +
                    sparsity_weight * spars_quality

      # Store results
      results <- rbind(results, data.frame(
        method = method,
        correlation = eval_result$correlation,
        sparsity_mean = eval_result$sparsity_mean,
        sparsity_quality = spars_quality,
        composite_score = comp_score,
        time_sec = elapsed
      ))

      models[[method]] <- fit

      if (verbose) {
        cat(sprintf("  Correlation: %.4f\n", eval_result$correlation))
        cat(sprintf("  Sparsity: %.2f%%\n", eval_result$sparsity_mean * 100))
        cat(sprintf("  Sparsity Quality: %.3f\n", spars_quality))
        if (use_composite_score) {
          cat(sprintf("  Composite Score: %.4f\n", comp_score))
        }
        cat(sprintf("  Time: %.2f sec\n\n", elapsed))
      }

    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  ERROR: %s\n\n", e$message))
      }
    })
  }

  # Determine best method
  if (use_composite_score) {
    best_idx <- which.max(results$composite_score)
  } else {
    best_idx <- which.max(results$correlation)
  }
  best_method <- results$method[best_idx]

  if (verbose) {
    cat("=== Recommendation ===\n")
    if (use_composite_score) {
      cat(sprintf("Best method: %s (composite score = %.4f, correlation = %.4f, sparsity quality = %.3f)\n",
                  best_method, results$composite_score[best_idx],
                  results$correlation[best_idx], results$sparsity_quality[best_idx]))
    } else {
      cat(sprintf("Best method: %s (correlation = %.4f)\n",
                  best_method, results$correlation[best_idx]))
    }
  }

  return(list(
    best_method = best_method,
    results = results,
    models = models
  ))
}


#' Convenience wrapper for ConvCCA cross-validation
#'
#' @param ... Arguments passed to MINERVA_cv_scca
#' @export
MINERVA_cv_convCCA <- function(...) {
  MINERVA_cv_scca(method = "ConvCCA", ...)
}


#' Convenience wrapper for RelPMDCCA cross-validation
#'
#' @param ... Arguments passed to MINERVA_cv_scca
#' @export
MINERVA_cv_relPMDCCA <- function(...) {
  MINERVA_cv_scca(method = "RelPMDCCA", ...)
}
