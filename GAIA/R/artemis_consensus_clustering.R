###############################################################################
########### Consensus Clustering for Patient Stratification ###########
###############################################################################

#' Consensus Clustering Across Multiple Imputed Datasets
#'
#' @description Runs ConsensusClusterPlus on each imputed dataset and averages
#' the pairwise consensus matrices across all imputations. This extends the
#' standard single-dataset approach to account for imputation uncertainty,
#' following the method of Wick et al. (Oslo sepsis clustering, 2024) who
#' averaged matrices from 100 imputed datasets.
#'
#' @param imputed_obj Either the output of \code{POSEIDON_mice_impute()} or a
#'   plain named list of data frames (one per imputation).
#' @param vars Character vector. Variables to use for clustering. If NULL, all
#'   columns are used. Default = NULL.
#' @param log_vars Character vector. Variables to log1p-transform before scaling.
#'   Should match the variables used in \code{POSEIDON_mice_impute()}. Default = NULL.
#' @param scale_data Logical. Scale variables to mean=0, SD=1 before clustering.
#'   Default = TRUE.
#' @param k_max Integer. Maximum number of clusters to test. Default = 10.
#' @param n_resample Integer. Number of subsampling resamples per dataset in
#'   ConsensusClusterPlus. Default = 100.
#' @param resample_prop Numeric. Fraction of samples drawn per resample. Default = 0.8.
#' @param clusterAlg Character. Internal clustering algorithm: "km" (k-means) or
#'   "hc" (hierarchical). Default = "km".
#' @param distance Character. Distance metric. Default = "euclidean".
#' @param seed Integer. Base random seed; each imputation uses seed+i. Default = 42.
#' @param verbose Logical. Print progress. Default = TRUE.
#'
#' @return A list with class "artemis_cc" containing:
#' \describe{
#'   \item{consensus_matrices}{List (indexed by k) of averaged n×n consensus matrices}
#'   \item{cluster_assignments}{List (indexed by k) of named integer cluster assignment vectors}
#'   \item{wcss}{Named vector of within-cluster sum-of-squares per k}
#'   \item{cdf_area}{Named vector of mean consensus value (area proxy) per k}
#'   \item{delta_area}{Named vector of relative CDF area change between consecutive k}
#'   \item{cdf_values}{List (indexed by k) of sorted upper-triangle values for CDF plotting}
#'   \item{k_max}{Maximum k tested}
#'   \item{k_range}{Integer vector 2:k_max}
#'   \item{n_imputations}{Number of imputed datasets used}
#'   \item{n_samples}{Number of samples}
#'   \item{sample_names}{Sample identifiers}
#'   \item{vars_used}{Variables used for clustering}
#'   \item{scaled_mean_data}{Mean preprocessed data matrix (for centroid calculation)}
#' }
#'
#' @details
#' Requires package \code{ConsensusClusterPlus} (BiocManager::install("ConsensusClusterPlus")).
#'
#' Cluster assignments are derived from the averaged consensus matrix via
#' hierarchical clustering (UPGMA) on (1 - consensus_matrix), then cutree at k.
#' WCSS is computed on the mean log1p-transformed, scaled dataset.
#'
#' @export
ARTEMIS_consensus_cluster <- function(imputed_obj, vars = NULL, log_vars = NULL,
                                       scale_data = TRUE, k_max = 10,
                                       n_resample = 100, resample_prop = 0.8,
                                       clusterAlg = "km", distance = "euclidean",
                                       seed = 42, verbose = TRUE) {

  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("Package 'ConsensusClusterPlus' is required. ",
         "Install with: BiocManager::install('ConsensusClusterPlus')")
  }

  # Accept POSEIDON_mice_impute output or plain list of data frames
  if (is.list(imputed_obj) && !is.null(imputed_obj$datasets)) {
    datasets <- imputed_obj$datasets
  } else if (is.list(imputed_obj) && is.data.frame(imputed_obj[[1]])) {
    datasets <- imputed_obj
  } else {
    stop("imputed_obj must be a list of data frames or output of POSEIDON_mice_impute()")
  }

  n_imp <- length(datasets)

  if (!is.null(vars)) {
    miss <- setdiff(vars, colnames(datasets[[1]]))
    if (length(miss) > 0)
      stop("Variables not found in data: ", paste(miss, collapse = ", "))
    datasets <- lapply(datasets, function(d) d[, vars, drop = FALSE])
  }

  n_samples    <- nrow(datasets[[1]])
  n_vars       <- ncol(datasets[[1]])
  sample_names <- rownames(datasets[[1]])
  vars_used    <- colnames(datasets[[1]])

  if (verbose) {
    cat("[ARTEMIS] Consensus Clustering\n")
    cat("    Imputed datasets:", n_imp, "\n")
    cat("    Samples:", n_samples, "| Variables:", n_vars, "\n")
    cat("    k range: 2 to", k_max, "\n")
    cat("    Resamples per dataset:", n_resample,
        "| Subsample fraction:", resample_prop, "\n\n")
  }

  # Internal preprocessing: log1p then scale
  .preprocess <- function(df) {
    mat <- as.matrix(df)
    if (!is.null(log_vars)) {
      lv <- intersect(log_vars, colnames(mat))
      for (v in lv) mat[, v] <- log1p(mat[, v])
    }
    if (scale_data) mat <- scale(mat)
    return(mat)
  }

  proc_datasets <- lapply(datasets, .preprocess)

  # Consensus matrix accumulators (list indexed 1:k_max, positions 1 unused)
  consensus_accum <- lapply(seq_len(k_max), function(i) {
    matrix(0, nrow = n_samples, ncol = n_samples)
  })

  tmp_dir <- tempfile(pattern = "ccp_")
  dir.create(tmp_dir, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  for (i in seq_len(n_imp)) {
    if (verbose && (i == 1 || i %% 10 == 0))
      cat("    Processing imputation", i, "/", n_imp, "\n")

    mat_i <- t(proc_datasets[[i]])  # CCP expects features x samples
    set.seed(seed + i)

    # Redirect graphics device: CCP draws to the active device internally even
    # when plot=NULL, flooding the IDE plot pane with 100s of renders per call.
    grDevices::pdf(file = nullfile())
    cc_i <- tryCatch(
      suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(
          d          = mat_i,
          maxK       = k_max,
          reps       = n_resample,
          pItem      = resample_prop,
          clusterAlg = clusterAlg,
          distance   = distance,
          title      = file.path(tmp_dir, paste0("imp", i)),
          plot       = NULL,
          writeTable = FALSE,
          verbose    = FALSE
        )
      ),
      finally = grDevices::dev.off()
    )

    for (k in 2:k_max) {
      if (!is.null(cc_i[[k]]$consensusMatrix))
        consensus_accum[[k]] <- consensus_accum[[k]] + cc_i[[k]]$consensusMatrix
    }
  }

  if (verbose) cat("\n    Averaging matrices and computing summaries...\n")

  mean_proc <- Reduce("+", proc_datasets) / n_imp

  consensus_matrices  <- vector("list", k_max)
  cluster_assignments <- vector("list", k_max)
  wcss_vals  <- rep(NA_real_, k_max)
  cdf_area   <- rep(NA_real_, k_max)
  cdf_values <- vector("list", k_max)

  for (k in 2:k_max) {
    cm <- consensus_accum[[k]] / n_imp
    if (!is.null(sample_names)) rownames(cm) <- colnames(cm) <- sample_names
    consensus_matrices[[k]] <- cm

    # HC on (1 - consensus) → cutree at k  (same derivation as CCP)
    hc   <- hclust(as.dist(1 - cm), method = "average")
    asgn <- cutree(hc, k = k)
    if (!is.null(sample_names)) names(asgn) <- sample_names
    cluster_assignments[[k]] <- asgn

    # WCSS on mean preprocessed data
    wcss <- 0
    for (cl in seq_len(k)) {
      idx <- which(asgn == cl)
      if (length(idx) > 1) {
        cl_data  <- mean_proc[idx, , drop = FALSE]
        centroid <- colMeans(cl_data)
        wcss     <- wcss + sum(sweep(cl_data, 2, centroid)^2)
      }
    }
    wcss_vals[k] <- wcss

    upper          <- cm[upper.tri(cm)]
    cdf_values[[k]] <- sort(upper)
    cdf_area[k]    <- mean(upper)
  }

  # Relative delta area between consecutive k
  delta_area    <- rep(NA_real_, k_max)
  delta_area[2] <- cdf_area[2]
  for (k in 3:k_max) {
    if (!is.na(cdf_area[k]) && !is.na(cdf_area[k - 1]) && cdf_area[k - 1] > 0)
      delta_area[k] <- (cdf_area[k] - cdf_area[k - 1]) / cdf_area[k - 1]
  }

  k_range <- 2:k_max
  result <- structure(
    list(
      consensus_matrices  = consensus_matrices,
      cluster_assignments = cluster_assignments,
      wcss                = setNames(wcss_vals[k_range],  as.character(k_range)),
      cdf_area            = setNames(cdf_area[k_range],   as.character(k_range)),
      delta_area          = setNames(delta_area[k_range], as.character(k_range)),
      cdf_values          = cdf_values,
      k_max               = k_max,
      k_range             = k_range,
      n_imputations       = n_imp,
      n_samples           = n_samples,
      sample_names        = sample_names,
      vars_used           = vars_used,
      scaled_mean_data    = mean_proc
    ),
    class = "artemis_cc"
  )

  if (verbose) cat("[ARTEMIS] Consensus Clustering complete.\n\n")
  return(result)
}


#' Cluster Stability Analysis for Consensus Clustering
#'
#' @description Assesses the stability of a chosen k by running ConsensusClusterPlus
#' \code{n_runs} times on the same imputed dataset with different random seeds. For
#' each patient, counts how often they are assigned to the same (modal) cluster
#' across all runs.
#'
#' @param imputed_obj Output of \code{POSEIDON_mice_impute()} or a list of data frames.
#' @param k Integer. Number of clusters to assess.
#' @param n_runs Integer. Number of independent CC runs. Default = 10.
#' @param vars Character vector. Variables to use. If NULL, all columns used. Default = NULL.
#' @param log_vars Character vector. Variables to log1p-transform before scaling. Default = NULL.
#' @param scale_data Logical. Scale data before clustering. Default = TRUE.
#' @param n_resample Integer. Resamples per CC run. Default = 100.
#' @param resample_prop Numeric. Subsample fraction. Default = 0.8.
#' @param clusterAlg Character. "km" or "hc". Default = "km".
#' @param distance Character. Distance metric. Default = "euclidean".
#' @param imputation_idx Integer. Which imputed dataset to use. Default = 1.
#' @param seed Integer. Base seed; run i uses seed+i. Default = 42.
#' @param verbose Logical. Print progress. Default = TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{assignment_matrix}{Integer matrix (n_samples × n_runs) of cluster assignments per run}
#'   \item{stability_counts}{Integer vector — for each patient, number of runs in their modal cluster}
#'   \item{pct_stable}{Numeric — percentage of patients in the same cluster across ALL runs}
#'   \item{k}{Number of clusters assessed}
#'   \item{n_runs}{Number of runs performed}
#'   \item{sample_names}{Sample identifiers}
#' }
#'
#' @details
#' Requires package \code{ConsensusClusterPlus}.
#'
#' @export
ARTEMIS_consensus_stability <- function(imputed_obj, k, n_runs = 10,
                                         vars = NULL, log_vars = NULL,
                                         scale_data = TRUE,
                                         n_resample = 100, resample_prop = 0.8,
                                         clusterAlg = "km", distance = "euclidean",
                                         imputation_idx = 1, seed = 42,
                                         verbose = TRUE) {

  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("Package 'ConsensusClusterPlus' is required. ",
         "Install with: BiocManager::install('ConsensusClusterPlus')")
  }

  if (is.list(imputed_obj) && !is.null(imputed_obj$datasets)) {
    base_data <- imputed_obj$datasets[[imputation_idx]]
  } else if (is.list(imputed_obj) && is.data.frame(imputed_obj[[1]])) {
    base_data <- imputed_obj[[imputation_idx]]
  } else {
    stop("imputed_obj must be a list of data frames or output of POSEIDON_mice_impute()")
  }

  if (!is.null(vars)) {
    base_data <- base_data[, vars, drop = FALSE]
  }

  sample_names <- rownames(base_data)
  n_samples    <- nrow(base_data)

  if (verbose) {
    cat("[ARTEMIS] Stability Analysis (k=", k, ", n_runs=", n_runs, ")\n", sep = "")
    cat("    Using imputation index:", imputation_idx, "\n\n")
  }

  # Preprocess once
  mat <- as.matrix(base_data)
  if (!is.null(log_vars)) {
    lv <- intersect(log_vars, colnames(mat))
    for (v in lv) mat[, v] <- log1p(mat[, v])
  }
  if (scale_data) mat <- scale(mat)
  mat_ccp <- t(mat)  # features x samples

  assignment_matrix <- matrix(NA_integer_, nrow = n_samples, ncol = n_runs)
  if (!is.null(sample_names)) rownames(assignment_matrix) <- sample_names
  colnames(assignment_matrix) <- paste0("run", seq_len(n_runs))

  tmp_dir <- tempfile(pattern = "ccp_stab_")
  dir.create(tmp_dir, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  for (r in seq_len(n_runs)) {
    if (verbose) cat("    Run", r, "/", n_runs, "\n")

    set.seed(seed + r)
    grDevices::pdf(file = nullfile())
    cc_r <- tryCatch(
      suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(
          d          = mat_ccp,
          # CCP's clusterTrackingPlot is called unconditionally; when maxK=2
          # colorM has 1 row, column subsetting drops the dimension to a vector,
          # and nrow(vector)=NULL causes 1:NULL to error. Using max(k,3) avoids
          # this while still yielding a valid cc_r[[k]] result for k=2.
          maxK       = max(k, 3L),
          reps       = n_resample,
          pItem      = resample_prop,
          clusterAlg = clusterAlg,
          distance   = distance,
          title      = file.path(tmp_dir, paste0("run", r)),
          plot       = NULL,
          writeTable = FALSE,
          verbose    = FALSE
        )
      ),
      finally = grDevices::dev.off()
    )

    # Derive assignments from consensus matrix at target k
    cm <- cc_r[[k]]$consensusMatrix
    hc <- hclust(as.dist(1 - cm), method = "average")
    assignment_matrix[, r] <- cutree(hc, k = k)
  }

  # For each patient, count how often they land in their modal cluster
  stability_counts <- apply(assignment_matrix, 1, function(row) {
    max(table(row))
  })

  pct_stable <- round(100 * mean(stability_counts == n_runs), 1)

  if (verbose) {
    cat("\n[ARTEMIS] Stability summary (k=", k, "):\n", sep = "")
    cat("    Patients stable across all", n_runs, "runs:",
        sum(stability_counts == n_runs), "(", pct_stable, "%)\n")
  }

  return(list(
    assignment_matrix = assignment_matrix,
    stability_counts  = stability_counts,
    pct_stable        = pct_stable,
    k                 = k,
    n_runs            = n_runs,
    sample_names      = sample_names
  ))
}
