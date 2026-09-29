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


#' Item Consensus per Sample, per Cluster
#'
#' @description For a chosen k, computes each sample's mean pairwise consensus
#' (from the averaged consensus matrix) against the members of every cluster,
#' not just its own assigned cluster. This is the classic Monti/Wilkerson
#' consensus-clustering "item consensus" concept (the same quantity
#' \code{ConsensusClusterPlus::calcICL()} computes for a single CCP run),
#' applied here to \code{ARTEMIS_consensus_cluster()}'s multi-imputation-
#' averaged consensus matrix.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param k Integer. Cluster solution to evaluate.
#' @param verbose Logical. Print progress/summary. Default = TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{profile}{n_samples x k matrix — mean pairwise consensus of each
#'     sample against each cluster's members (self excluded from its own
#'     cluster's average).}
#'   \item{summary}{Data frame with one row per sample: \code{cluster}
#'     (assigned cluster), \code{assigned_consensus} (consensus with its own
#'     cluster), \code{best_other_consensus} (highest consensus with any
#'     other cluster), and \code{margin} (assigned minus best_other). A low
#'     or negative margin flags a sample nearly as consensus-linked to a
#'     different cluster as to its own — a fringe/unstable call, as opposed
#'     to a sample solidly tied to its own cluster.}
#'   \item{k}{The k evaluated.}
#' }
#'
#' @details
#' A cluster with only one member has no other same-cluster sample to
#' average against, so that sample's own-cluster consensus (and therefore
#' its margin) is returned as \code{NA}.
#'
#' @export
ARTEMIS_item_consensus <- function(cc_result, k, verbose = TRUE) {

  if (!inherits(cc_result, "artemis_cc"))
    stop("cc_result must be the output of ARTEMIS_consensus_cluster()")
  if (!k %in% cc_result$k_range)
    stop("k=", k, " is not in cc_result$k_range (",
         paste(range(cc_result$k_range), collapse = "-"), ")")

  cm   <- cc_result$consensus_matrices[[k]]
  asgn <- cc_result$cluster_assignments[[k]]
  if (is.null(cm) || is.null(asgn))
    stop("cc_result has no consensus matrix / cluster assignment stored for k=", k)

  samples      <- names(asgn)
  clust_levels <- seq_len(k)

  if (verbose) {
    cat("[ARTEMIS] Item Consensus (k=", k, ")\n", sep = "")
    cat("    Samples:", length(samples), "\n\n")
  }

  profile <- matrix(NA_real_, nrow = length(samples), ncol = k,
                     dimnames = list(samples, paste0("C", clust_levels)))

  for (cl in clust_levels) {
    members     <- samples[asgn == cl]
    non_members <- setdiff(samples, members)
    if (length(members) == 0) next

    if (length(non_members) > 0)
      profile[non_members, cl] <- rowMeans(cm[non_members, members, drop = FALSE])

    if (length(members) > 1) {
      for (s in members)
        profile[s, cl] <- mean(cm[s, setdiff(members, s)])
    }
    # length(members) == 1: leave as NA -- no other same-cluster sample
  }

  assigned_consensus <- vapply(samples, function(s) profile[s, asgn[s]], numeric(1))
  best_other <- vapply(samples, function(s) {
    other_cols <- setdiff(clust_levels, asgn[s])
    if (length(other_cols) == 0) return(NA_real_)
    max(profile[s, other_cols], na.rm = TRUE)
  }, numeric(1))

  summary_df <- data.frame(
    sample               = samples,
    cluster              = unname(asgn[samples]),
    assigned_consensus   = unname(assigned_consensus),
    best_other_consensus = unname(best_other),
    margin               = unname(assigned_consensus - best_other),
    row.names            = NULL,
    stringsAsFactors     = FALSE
  )

  if (verbose) {
    cat("    Mean assigned-cluster consensus:",
        round(mean(summary_df$assigned_consensus, na.rm = TRUE), 3), "\n")
    cat("    Mean margin:", round(mean(summary_df$margin, na.rm = TRUE), 3), "\n")
    cat("    Samples with negative margin (closer to a different cluster):",
        sum(summary_df$margin < 0, na.rm = TRUE), "/", nrow(summary_df), "\n\n")
  }

  list(profile = profile, summary = summary_df, k = k)
}


#' Neighbor-Based Confidence for Projected (Out-of-Sample) Samples
#'
#' @description For samples projected into the Train factor space (e.g. via
#' \code{HEPHAESTUS_project_mofa_samples()} plus the same standardization
#' Train used), looks up each sample's nearest Train neighbors in that real
#' K-dimensional factor space and computes two DIFFERENT, complementary
#' quantities from that neighbor set:
#' \itemize{
#'   \item Aggregated neighbor \code{ARTEMIS_item_consensus()} margin/consensus
#'     — how solidly, INDIVIDUALLY, those specific neighbors were themselves
#'     tied to their own cluster. A projected sample has no resampling
#'     history of its own, so it cannot have a true item consensus — this is
#'     a neighbor-borrowed stand-in for that.
#'   \item Neighbor cluster purity — a SEPARATE question the above cannot
#'     answer: do those same neighbors even AGREE on a cluster? A sample can
#'     have high aggregated neighbor margin (its neighbors are each
#'     individually solid members of their own clusters) while its neighbor
#'     SET is still split across two different clusters — i.e. the sample
#'     itself sits ambiguously between them. Purity is a plain majority-vote
#'     count over the neighbors' cluster labels, with no resampling involved
#'     at all (unlike margin) — a different mechanism, not just a
#'     differently-named copy of the same idea.
#' }
#' Empirically (checked against a real projected cohort), these two are only
#' weakly correlated (Spearman rho around 0.3) — genuinely complementary, not
#' redundant.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}. Supplies the
#'   Train reference factor space (\code{cc_result$scaled_mean_data}).
#' @param item_consensus Output of \code{ARTEMIS_item_consensus()}, computed
#'   for the SAME k as the clustering being evaluated against.
#' @param new_scores Matrix or data frame of projected samples' factor
#'   scores, ALREADY standardized with Train's own center/SD (the same
#'   scaling Train's clustering used — not re-derived here). Rows = samples;
#'   columns must include every column of \code{cc_result$scaled_mean_data}.
#' @param n_neighbors Integer. Number of nearest Train neighbors per new
#'   sample. Default = 10.
#' @param aggregate Character. How to summarize the neighbors' margin/
#'   consensus: \code{"mean"} (average confidence of the neighborhood) or
#'   \code{"min"} (worst-case neighbor, most conservative). Default = "mean".
#' @param verbose Logical. Print progress/summary. Default = TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{per_sample}{Named list (by new-sample ID) of neighbor IDs,
#'     neighbor distances, the aggregated neighbor margin/consensus, and the
#'     neighbor-set cluster purity/majority cluster.}
#'   \item{summary}{Data frame, one row per new sample: \code{neighbor_margin},
#'     \code{neighbor_assigned_consensus}, \code{neighbor_cluster_purity}
#'     (fraction of the \code{n_neighbors} nearest Train neighbors sharing the
#'     majority cluster — NOT affected by \code{aggregate}, since it is a
#'     property of the whole neighbor set, not an aggregate of per-neighbor
#'     values), \code{neighbor_majority_cluster} (which cluster that is).}
#'   \item{n_neighbors}{Neighbors used.}
#'   \item{aggregate}{Aggregation method used (applies to margin/consensus only).}
#' }
#'
#' @details
#' Distance is Euclidean, computed in \code{cc_result$scaled_mean_data}'s own
#' column space — the same space/metric \code{ARTEMIS_consensus_cluster()}
#' and nearest-centroid assignment already use, not any 2D visualization
#' space (a UMAP or PCA embedding can distort neighbor relationships and
#' should not be used to compute this).
#'
#' @export
ARTEMIS_project_confidence <- function(cc_result, item_consensus, new_scores,
                                        n_neighbors = 10,
                                        aggregate = c("mean", "min"),
                                        verbose = TRUE) {

  aggregate <- match.arg(aggregate)

  if (!inherits(cc_result, "artemis_cc"))
    stop("cc_result must be the output of ARTEMIS_consensus_cluster()")
  if (is.null(item_consensus$summary))
    stop("item_consensus must be the output of ARTEMIS_item_consensus()")

  train_mat  <- cc_result$scaled_mean_data
  train_vars <- colnames(train_mat)

  new_scores <- as.matrix(new_scores)
  miss_vars  <- setdiff(train_vars, colnames(new_scores))
  if (length(miss_vars) > 0)
    stop("new_scores is missing column(s) present in cc_result$scaled_mean_data: ",
         paste(miss_vars, collapse = ", "))
  new_scores <- new_scores[, train_vars, drop = FALSE]

  ic_summary <- item_consensus$summary
  rownames(ic_summary) <- ic_summary$sample
  ref_pool <- intersect(rownames(train_mat), ic_summary$sample)
  if (length(ref_pool) < nrow(train_mat))
    warning(nrow(train_mat) - length(ref_pool),
            " Train sample(s) in cc_result have no entry in item_consensus$summary -- excluded from the neighbor pool")
  if (n_neighbors > length(ref_pool))
    stop("n_neighbors (", n_neighbors, ") exceeds the usable Train reference pool (",
         length(ref_pool), ")")

  train_mat   <- train_mat[ref_pool, , drop = FALSE]
  ic_margin   <- stats::setNames(ic_summary[ref_pool, "margin"], ref_pool)
  ic_assigned <- stats::setNames(ic_summary[ref_pool, "assigned_consensus"], ref_pool)
  ic_cluster  <- stats::setNames(ic_summary[ref_pool, "cluster"], ref_pool)

  new_ids <- rownames(new_scores)
  if (is.null(new_ids)) new_ids <- paste0("new", seq_len(nrow(new_scores)))

  if (verbose) {
    cat("[ARTEMIS] Projection Confidence\n")
    cat("    New samples:", nrow(new_scores), "| Train reference pool:", length(ref_pool), "\n")
    cat("    Neighbors per sample:", n_neighbors, "| Aggregate:", aggregate, "\n\n")
  }

  agg_fun <- if (aggregate == "mean") {
    function(x) mean(x, na.rm = TRUE)
  } else {
    function(x) min(x, na.rm = TRUE)
  }

  per_sample <- lapply(seq_len(nrow(new_scores)), function(i) {
    d  <- sqrt(rowSums(sweep(train_mat, 2, new_scores[i, ], "-")^2))
    nn <- names(sort(d))[seq_len(n_neighbors)]
    nn_tab <- table(ic_cluster[nn])
    list(
      neighbors                   = nn,
      neighbor_distances          = d[nn],
      neighbor_margin             = agg_fun(ic_margin[nn]),
      neighbor_assigned_consensus = agg_fun(ic_assigned[nn]),
      neighbor_cluster_purity     = as.numeric(max(nn_tab)) / length(nn),
      neighbor_majority_cluster   = names(nn_tab)[which.max(nn_tab)]
    )
  })
  names(per_sample) <- new_ids

  summary_df <- data.frame(
    sample                       = new_ids,
    neighbor_margin               = vapply(per_sample, `[[`, numeric(1), "neighbor_margin"),
    neighbor_assigned_consensus   = vapply(per_sample, `[[`, numeric(1), "neighbor_assigned_consensus"),
    neighbor_cluster_purity       = vapply(per_sample, `[[`, numeric(1), "neighbor_cluster_purity"),
    neighbor_majority_cluster     = vapply(per_sample, `[[`, character(1), "neighbor_majority_cluster"),
    row.names                     = NULL,
    stringsAsFactors              = FALSE
  )

  if (verbose) {
    cat("    Neighbor-margin across new samples -- median:",
        round(stats::median(summary_df$neighbor_margin), 3),
        "| range:", round(min(summary_df$neighbor_margin), 3), "to",
        round(max(summary_df$neighbor_margin), 3), "\n")
    cat("    Neighbor-cluster-purity across new samples -- median:",
        round(stats::median(summary_df$neighbor_cluster_purity), 3),
        "| samples below 0.6 (neighbor set meaningfully mixed):",
        sum(summary_df$neighbor_cluster_purity < 0.6), "/", nrow(summary_df), "\n\n")
  }

  list(per_sample = per_sample, summary = summary_df,
       n_neighbors = n_neighbors, aggregate = aggregate)
}
