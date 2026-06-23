###############################################################################
############################# Imputation #######################################
###############################################################################

#' Multiple Imputation by Chained Equations (MICE)
#'
#' @description Imputes missing values in a numeric data frame using MICE with
#' predictive mean matching. Supports log1p-transformation of right-skewed
#' variables before imputation, with automatic back-transformation in the output.
#' Returns all completed datasets plus a mean-imputed dataset for use with
#' \code{ARTEMIS_consensus_cluster()}.
#'
#' @param data Numeric data frame. Should contain only the variables to be
#'   imputed/used for clustering — do not pass the full study dataset.
#' @param n_imputations Integer. Number of imputed datasets to generate.
#'   Oslo (Wick et al.) used 100; fewer imputations increase variability. Default = 100.
#' @param method Character. Imputation method passed to \code{mice::mice()}.
#'   "pmm" (predictive mean matching) is appropriate for continuous variables. Default = "pmm".
#' @param max_iter Integer. Maximum MICE iterations per imputation. Default = 20.
#' @param log_vars Character vector. Columns to log1p-transform before imputation
#'   (back-transformed in all output datasets). Intended for right-skewed lab values
#'   (e.g., CRP, lactate, bilirubin). Default = NULL.
#' @param seed Integer. Random seed for reproducibility. Default = 42.
#' @param verbose Logical. Print progress. Default = TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{datasets}{List of \code{n_imputations} completed data frames on original scale}
#'   \item{mean_dataset}{Data frame with column means across all imputations (used for centroid computation)}
#'   \item{log_vars}{Character vector of variables that were log1p-transformed}
#'   \item{n_imputations}{Number of imputations generated}
#'   \item{n_obs}{Number of observations}
#'   \item{n_vars}{Number of variables}
#' }
#'
#' @details
#' Requires the \code{mice} package (install with \code{install.packages("mice")}).
#'
#' @export
POSEIDON_mice_impute <- function(data, n_imputations = 100, method = "pmm",
                                  max_iter = 20, log_vars = NULL,
                                  seed = 42, verbose = TRUE) {

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required. Install with: install.packages('mice')")
  }

  if (!is.data.frame(data)) data <- as.data.frame(data)

  n_obs  <- nrow(data)
  n_var  <- ncol(data)
  n_miss <- colSums(is.na(data))

  if (verbose) {
    cat("[POSEIDON] MICE Imputation\n")
    cat("    Variables:", n_var, " | Samples:", n_obs,
        " | Imputations:", n_imputations, "\n")
    if (any(n_miss > 0)) {
      cat("    Missing values:\n")
      for (v in names(n_miss[n_miss > 0])) {
        cat("        ", v, ": ", n_miss[v],
            " (", round(100 * n_miss[v] / n_obs, 1), "%)\n", sep = "")
      }
    } else {
      cat("    No missing values detected.\n")
    }
  }

  # Log1p-transform skewed variables before imputation
  data_trans      <- data
  log_vars_present <- intersect(log_vars, colnames(data))
  if (length(log_vars_present) > 0) {
    if (verbose)
      cat("    log1p-transforming:", paste(log_vars_present, collapse = ", "), "\n")
    for (v in log_vars_present) data_trans[[v]] <- log1p(data_trans[[v]])
  }

  if (verbose) cat("    Running MICE (this may take a few minutes)...\n")
  set.seed(seed)
  mids <- mice::mice(data_trans, m = n_imputations, method = method,
                     maxit = max_iter, printFlag = FALSE, seed = seed)

  # Extract and back-transform completed datasets
  if (verbose) cat("    Extracting completed datasets...\n")
  datasets <- vector("list", n_imputations)
  for (i in seq_len(n_imputations)) {
    completed <- mice::complete(mids, action = i)
    for (v in log_vars_present) completed[[v]] <- expm1(completed[[v]])
    datasets[[i]] <- completed
  }

  # Mean dataset: column-wise mean across all imputations (original scale)
  mean_dataset <- data
  imp_cols <- names(n_miss[n_miss > 0])
  for (v in imp_cols) {
    imp_mat <- do.call(cbind, lapply(datasets, function(d) d[[v]]))
    mean_dataset[[v]] <- rowMeans(imp_mat)
  }

  if (verbose) cat("[POSEIDON] Imputation complete.\n\n")

  return(list(
    datasets     = datasets,
    mean_dataset = mean_dataset,
    log_vars     = log_vars_present,
    n_imputations = n_imputations,
    n_obs        = n_obs,
    n_vars       = n_var
  ))
}
