#' Cluster Mixed Data Using Gower Distance
#'
#' @description Clustering for datasets containing both quantitative and
#' qualitative variables. Uses Gower distance which handles mixed types
#' natively, then applies PAM (partitioning around medoids) or hierarchical
#' clustering.
#'
#' @param data Data frame with mixed variable types. Should be pre-cleaned
#'   (no high-NA columns, appropriate factor conversions done).
#' @param k Integer or "auto". Number of clusters. If "auto", determines
#'   optimal k using silhouette scores over k_range. Default = "auto".
#' @param k_range Integer vector. Range of k values to test when k = "auto".
#'   Default = 2:10.
#' @param method Character. Clustering method: "pam" (partitioning around
#'   medoids) or "hierarchical". Default = "pam".
#' @param hclust_method Character. Linkage method for hierarchical clustering.
#'   One of "ward.D2", "complete", "average", "single". Default = "ward.D2".
#' @param stand Logical. Standardize variables before computing Gower distance?
#'   Default = FALSE (Gower handles scale differences internally).
#' @param verbose Logical. Print progress messages. Default = TRUE.
#'
#' @return A list with class "artemis_cluster" containing:
#' \describe{
#'   \item{clusters}{Integer vector of cluster assignments}
#'   \item{k}{Number of clusters used}
#'   \item{k_selection}{If k="auto", data frame with silhouette scores per k}
#'   \item{k_evaluation}{If k="auto", full evaluation data (class "artemis_k_evaluation")
#'     containing assignments at each k and transition data for Sankey visualization.
#'     Pass to AETHER_plot_cluster_sankey() or AETHER_plot_k_selection(). NULL if k
#'     was specified directly.}
#'   \item{silhouette}{Silhouette information for final clustering}
#'   \item{silhouette_avg}{Average silhouette width (cluster quality measure)}
#'   \item{medoids}{For PAM: indices of medoid samples}
#'   \item{distance}{The Gower distance matrix}
#'   \item{method}{Clustering method used}
#'   \item{data_used}{The data frame used for clustering}
#'   \item{call}{The function call}
#' }
#'
#' @details
#' Gower distance handles mixed data by computing appropriate distances for
#' each variable type:
#' \itemize{
#'   \item Quantitative: Manhattan distance, normalized by range
#'   \item Qualitative: Simple matching (0 if same, 1 if different)
#'   \item Ordinal: Ranked, then treated as quantitative
#' }
#'
#' PAM is generally preferred over k-means for mixed data because:
#' \itemize{
#'   \item Works with any distance metric (not just Euclidean)
#'   \item More robust to outliers (uses medoids, not centroids)
#'   \item Medoids are actual data points, aiding interpretation
#' }
#'
#' @export
#'
#' @examples
#' # Auto-select k
#' result <- ARTEMIS_cluster_mixed(my_data)
#'
#' # Visualize k selection and cluster transitions (when k="auto")
#' AETHER_plot_k_selection(result$k_evaluation)
#' AETHER_plot_cluster_sankey(result$k_evaluation)
#'
#' # Specify k directly (no k_evaluation generated)
#' result <- ARTEMIS_cluster_mixed(my_data, k = 3)
#'
#' # Use hierarchical clustering
#' result <- ARTEMIS_cluster_mixed(my_data, method = "hierarchical")
#'
ARTEMIS_cluster_mixed <- function(data,
                                   k = "auto",
                                   k_range = 2:10,
                                   method = "pam",
                                   hclust_method = "ward.D2",
                                   stand = FALSE,
                                   verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (nrow(data) < 3) {
    stop("Need at least 3 observations for clustering")
  }

  if (!method %in% c("pam", "hierarchical")) {
    stop("method must be 'pam' or 'hierarchical'")
  }

  if (!hclust_method %in% c("ward.D", "ward.D2", "single", "complete",
                            "average", "mcquitty", "median", "centroid")) {
    stop("Invalid hclust_method. See ?hclust for options.")
  }

  # Check for cluster package

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Install with: install.packages('cluster')")
  }

  # Check for NA values
  if (any(is.na(data))) {
    n_na <- sum(is.na(data))
    stop("Data contains ", n_na, " NA values. Please handle missing data before clustering.\n",
         "Use HADES_drop_high_na_variables() or ARTEMIS_famd() with na_action='impute'.")
  }

  # ---------------------------------------------------------------------------
  # Convert character columns to factors (required by daisy)
  # ---------------------------------------------------------------------------
  char_cols <- sapply(data, is.character)
  if (any(char_cols)) {
    data[char_cols] <- lapply(data[char_cols], as.factor)
  }

  if (verbose) {
    cat("ARTEMIS Mixed Data Clustering\n")
    cat("==============================\n")
    cat("Input:", nrow(data), "samples,", ncol(data), "variables\n")

    # Report variable types
    n_quanti <- sum(sapply(data, is.numeric))
    n_quali <- sum(sapply(data, is.factor))
    n_char_converted <- sum(char_cols)
    cat("  Quantitative:", n_quanti, "\n")
    cat("  Qualitative:", n_quali)
    if (n_char_converted > 0) {
      cat(" (", n_char_converted, " converted from character)", sep = "")
    }
    cat("\n\n")
  }

  # ---------------------------------------------------------------------------
  # Step 1: Compute Gower distance
  # ---------------------------------------------------------------------------
  if (verbose) cat("Computing Gower distance matrix...\n")

  gower_dist <- cluster::daisy(data, metric = "gower", stand = stand)

  # Check for issues
  if (any(is.na(gower_dist))) {
    warning("Gower distance matrix contains NA values. This may cause issues.")
  }

  # ---------------------------------------------------------------------------
  # Step 2: Determine k if auto
  # ---------------------------------------------------------------------------
  k_selection <- NULL
  k_evaluation <- NULL  # Full evaluation data for Sankey visualization

  if (identical(k, "auto")) {
    if (verbose) cat("Determining optimal k using silhouette scores...\n")

    # Limit k_range to reasonable values
    max_k <- min(max(k_range), floor(nrow(data) / 2), 15)
    k_range <- k_range[k_range <= max_k & k_range >= 2]

    if (length(k_range) == 0) {
      stop("No valid k values in k_range for this dataset size")
    }

    # Get sample names
    sample_names <- rownames(data)
    if (is.null(sample_names)) {
      sample_names <- paste0("S", seq_len(nrow(data)))
    }

    # Storage for full evaluation
    metrics_list <- list()
    assignments <- data.frame(sample = sample_names, stringsAsFactors = FALSE)

    for (i in seq_along(k_range)) {
      k_test <- k_range[i]

      if (method == "pam") {
        clust_test <- cluster::pam(gower_dist, k = k_test, diss = TRUE)
        clusters_test <- clust_test$clustering
        sil_info_test <- cluster::silhouette(clust_test)
      } else {
        hc <- hclust(gower_dist, method = hclust_method)
        clusters_test <- cutree(hc, k = k_test)
        sil_info_test <- cluster::silhouette(clusters_test, gower_dist)
      }

      sil_score <- mean(sil_info_test[, "sil_width"])

      # Store metrics
      cluster_sizes <- table(clusters_test)
      size_str <- paste(as.numeric(cluster_sizes), collapse = "/")

      metrics_list[[as.character(k_test)]] <- data.frame(
        k = k_test,
        silhouette_avg = sil_score,
        cluster_sizes = size_str,
        min_cluster_size = min(cluster_sizes),
        max_cluster_size = max(cluster_sizes),
        stringsAsFactors = FALSE
      )

      # Store assignments
      col_name <- paste0("k", k_test)
      assignments[[col_name]] <- clusters_test

      if (verbose) {
        cat("  k =", k_test, ": avg silhouette =", round(sil_score, 3),
            " sizes:", size_str, "\n")
      }
    }

    # Build metrics data frame
    k_selection <- do.call(rbind, metrics_list)
    rownames(k_selection) <- NULL

    # Build transitions data for Sankey diagram
    transitions_list <- list()
    for (i in seq_len(length(k_range) - 1)) {
      k_from <- k_range[i]
      k_to <- k_range[i + 1]

      col_from <- paste0("k", k_from)
      col_to <- paste0("k", k_to)

      trans_table <- table(
        from = paste0("k", k_from, "_C", assignments[[col_from]]),
        to = paste0("k", k_to, "_C", assignments[[col_to]])
      )

      trans_df <- as.data.frame(trans_table, stringsAsFactors = FALSE)
      names(trans_df) <- c("source", "target", "value")
      trans_df <- trans_df[trans_df$value > 0, ]
      trans_df$k_from <- k_from
      trans_df$k_to <- k_to

      transitions_list[[i]] <- trans_df
    }

    transitions <- do.call(rbind, transitions_list)
    rownames(transitions) <- NULL

    # Select k with highest silhouette
    best_idx <- which.max(k_selection$silhouette_avg)
    k <- k_range[best_idx]

    # Build k_evaluation object (compatible with AETHER_plot_cluster_sankey)
    k_evaluation <- list(
      metrics = k_selection,
      assignments = assignments,
      transitions = transitions,
      optimal_k = k,
      k_range = k_range,
      method = method
    )
    class(k_evaluation) <- c("artemis_k_evaluation", "list")

    if (verbose) {
      cat("\nOptimal k =", k, "(silhouette =", round(k_selection$silhouette_avg[best_idx], 3), ")\n\n")
    }
  } else {
    # Validate provided k
    k <- as.integer(k)
    if (k < 2) stop("k must be >= 2")
    if (k >= nrow(data)) stop("k must be less than number of samples")
  }

  # ---------------------------------------------------------------------------
  # Step 3: Perform clustering
  # ---------------------------------------------------------------------------
  if (verbose) cat("Clustering with k =", k, "using", method, "method...\n")

  if (method == "pam") {
    clust_result <- cluster::pam(gower_dist, k = k, diss = TRUE)
    clusters <- clust_result$clustering
    medoids <- clust_result$medoids
    sil_info <- cluster::silhouette(clust_result)
  } else {
    hc <- hclust(gower_dist, method = hclust_method)
    clusters <- cutree(hc, k = k)
    medoids <- NULL
    sil_info <- cluster::silhouette(clusters, gower_dist)
  }

  sil_avg <- mean(sil_info[, "sil_width"])

  # ---------------------------------------------------------------------------
  # Step 4: Summarize results
  # ---------------------------------------------------------------------------
  cluster_sizes <- table(clusters)

  if (verbose) {
    cat("\nClustering complete.\n")
    cat("  Average silhouette width:", round(sil_avg, 3), "\n")
    cat("  Cluster sizes:\n")
    for (i in seq_along(cluster_sizes)) {
      cat("    Cluster", names(cluster_sizes)[i], ":", cluster_sizes[i], "samples\n")
    }

    # Interpret silhouette
    cat("\n  Silhouette interpretation:\n")
    if (sil_avg > 0.7) {
      cat("    Strong structure found\n")
    } else if (sil_avg > 0.5) {
      cat("    Reasonable structure found\n")
    } else if (sil_avg > 0.25) {
      cat("    Weak structure - clusters may overlap\n")
    } else {
      cat("    No substantial structure - consider different k or variables\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Build result
  # ---------------------------------------------------------------------------
  result <- list(
    clusters = clusters,
    k = k,
    k_selection = k_selection,
    k_evaluation = k_evaluation,  # Full evaluation data for Sankey (NULL if k was specified)
    silhouette = sil_info,
    silhouette_avg = sil_avg,
    medoids = medoids,
    distance = gower_dist,
    method = method,
    hclust_method = if (method == "hierarchical") hclust_method else NULL,
    cluster_sizes = as.integer(cluster_sizes),
    data_used = data,
    call = match.call()
  )

  class(result) <- c("artemis_cluster", "list")

  return(result)
}


#' Print method for ARTEMIS clustering results
#'
#' @param x An artemis_cluster object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.artemis_cluster <- function(x, ...) {
  cat("ARTEMIS Clustering Result\n")
  cat("=========================\n")
  cat("Method:", x$method)
  if (!is.null(x$hclust_method)) cat(" (", x$hclust_method, ")", sep = "")
  cat("\n")
  cat("Samples:", nrow(x$data_used), "\n")
  cat("Variables:", ncol(x$data_used), "\n")
  cat("Clusters (k):", x$k, "\n")
  cat("\nCluster sizes:\n")
  for (i in seq_along(x$cluster_sizes)) {
    cat("  Cluster", i, ":", x$cluster_sizes[i], "\n")
  }
  cat("\nAverage silhouette width:", round(x$silhouette_avg, 3), "\n")

  if (!is.null(x$k_selection)) {
    cat("\nk selection (silhouette scores):\n")
    best_k <- x$k_selection$k[which.max(x$k_selection$silhouette_avg)]
    for (i in seq_len(nrow(x$k_selection))) {
      marker <- if (x$k_selection$k[i] == best_k) " <-- best" else ""
      cat("  k =", x$k_selection$k[i], ":",
          round(x$k_selection$silhouette_avg[i], 3), marker, "\n")
    }
  }

  invisible(x)
}


#' Assess Cluster Stability via Bootstrap
#'
#' @description Evaluates how stable cluster assignments are under bootstrap
#' resampling. High stability suggests the clustering captures real structure
#' rather than noise.
#'
#' @param cluster_result An artemis_cluster object from ARTEMIS_cluster_mixed().
#' @param n_boot Integer. Number of bootstrap iterations. Default = 100.
#' @param verbose Logical. Print progress. Default = TRUE.
#'
#' @return A list with class "artemis_stability" containing:
#' \describe{
#'   \item{jaccard_mean}{Mean Jaccard similarity across bootstrap samples}
#'   \item{jaccard_per_cluster}{Jaccard similarity for each cluster}
#'   \item{sample_stability}{Per-sample stability score (proportion of times
#'     assigned to same cluster as in original)}
#'   \item{n_boot}{Number of bootstrap iterations performed}
#'   \item{interpretation}{Text interpretation of stability}
#' }
#'
#' @details
#' For each bootstrap iteration:
#' 1. Resample data with replacement
#' 2. Re-cluster using same parameters
#' 3. Compare cluster assignments for samples present in both
#'
#' Jaccard similarity measures overlap between original and bootstrap clusters.
#' Values > 0.75 indicate stable clusters; < 0.5 indicates instability.
#'
#' @export
#'
ARTEMIS_cluster_stability <- function(cluster_result,
                                       n_boot = 100,
                                       verbose = TRUE) {

  if (!inherits(cluster_result, "artemis_cluster")) {
    stop("cluster_result must be an artemis_cluster object")
  }

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required.")
  }

  data <- cluster_result$data_used
  k <- cluster_result$k
  method <- cluster_result$method
  hclust_method <- cluster_result$hclust_method
  original_clusters <- cluster_result$clusters

  n_samples <- nrow(data)
  sample_stability <- numeric(n_samples)
  jaccard_all <- numeric(n_boot)

  if (verbose) {
    cat("Assessing cluster stability with", n_boot, "bootstrap iterations...\n")
    pb_interval <- max(1, floor(n_boot / 10))
  }

  for (b in seq_len(n_boot)) {
    # Bootstrap sample (with replacement)
    boot_idx <- sample(n_samples, replace = TRUE)
    boot_data <- data[boot_idx, , drop = FALSE]

    # Compute Gower distance
    gower_dist <- cluster::daisy(boot_data, metric = "gower")

    # Cluster
    if (method == "pam") {
      clust_boot <- cluster::pam(gower_dist, k = k, diss = TRUE)
      boot_clusters <- clust_boot$clustering
    } else {
      hc <- hclust(gower_dist, method = hclust_method)
      boot_clusters <- cutree(hc, k = k)
    }

    # Map bootstrap clusters to original clusters (best match by overlap)
    # This handles the label switching problem
    mapping <- .match_cluster_labels(original_clusters[boot_idx], boot_clusters)
    boot_clusters_mapped <- mapping[boot_clusters]

    # Calculate Jaccard for this bootstrap
    # For each unique original sample index, check if cluster matches
    unique_idx <- unique(boot_idx)
    matches <- 0
    total <- 0

    for (orig_idx in unique_idx) {
      orig_cluster <- original_clusters[orig_idx]
      boot_positions <- which(boot_idx == orig_idx)
      boot_cluster <- boot_clusters_mapped[boot_positions[1]]

      if (orig_cluster == boot_cluster) {
        matches <- matches + 1
        sample_stability[orig_idx] <- sample_stability[orig_idx] + 1
      }
      total <- total + 1
    }

    jaccard_all[b] <- matches / total

    if (verbose && b %% pb_interval == 0) {
      cat("  Iteration", b, "/", n_boot, "\n")
    }
  }

  # Normalize sample stability
  sample_stability <- sample_stability / n_boot

  # Calculate per-cluster Jaccard
  jaccard_per_cluster <- tapply(sample_stability, original_clusters, mean)

  # Overall stats
  jaccard_mean <- mean(jaccard_all)

  # Interpretation
  if (jaccard_mean > 0.85) {
    interpretation <- "Highly stable clusters"
  } else if (jaccard_mean > 0.75) {
    interpretation <- "Stable clusters"
  } else if (jaccard_mean > 0.6) {
    interpretation <- "Moderately stable clusters"
  } else if (jaccard_mean > 0.5) {
    interpretation <- "Somewhat unstable clusters"
  } else {
    interpretation <- "Unstable clusters - structure may not be reliable"
  }

  if (verbose) {
    cat("\nStability Results:\n")
    cat("  Mean Jaccard similarity:", round(jaccard_mean, 3), "\n")
    cat("  Interpretation:", interpretation, "\n")
    cat("\n  Per-cluster stability:\n")
    for (i in seq_along(jaccard_per_cluster)) {
      cat("    Cluster", names(jaccard_per_cluster)[i], ":",
          round(jaccard_per_cluster[i], 3), "\n")
    }
  }

  result <- list(
    jaccard_mean = jaccard_mean,
    jaccard_all = jaccard_all,
    jaccard_per_cluster = jaccard_per_cluster,
    sample_stability = sample_stability,
    n_boot = n_boot,
    interpretation = interpretation
  )

  class(result) <- c("artemis_stability", "list")

  return(result)
}


#' Match cluster labels between two clusterings
#'
#' @description Internal helper to handle label switching problem.
#' Finds the best mapping from cluster labels in clust2 to clust1.
#'
#' @param clust1 Original cluster assignments
#' @param clust2 New cluster assignments to map
#'
#' @return Named integer vector: mapping from clust2 labels to clust1 labels
#'
#' @keywords internal
.match_cluster_labels <- function(clust1, clust2) {
  k <- length(unique(clust1))
  labels1 <- sort(unique(clust1))
  labels2 <- sort(unique(clust2))

  # Build confusion matrix
  conf_mat <- table(clust2, clust1)

  # Greedy matching (good enough for this purpose)
  mapping <- integer(max(labels2))
  used <- logical(length(labels1))

  for (i in seq_len(nrow(conf_mat))) {
    # Find best unused match for this cluster
    overlaps <- conf_mat[i, ]
    overlaps[used] <- -1
    best_match <- which.max(overlaps)
    mapping[labels2[i]] <- labels1[best_match]
    used[best_match] <- TRUE
  }

  return(mapping)
}


#' Characterize Clusters by Variable Discriminative Power
#'
#' @description Determines which variables best distinguish between clusters
#' by testing each variable's association with cluster membership and
#' calculating effect sizes.
#'
#' @param cluster_result An artemis_cluster object from ARTEMIS_cluster_mixed().
#' @param data Optional data frame. If NULL (default), uses data stored in
#'   cluster_result. Provide if you want to test additional variables not
#'   used in clustering.
#' @param p_adjust Character. Method for p-value adjustment. One of
#'   "BH" (Benjamini-Hochberg), "bonferroni", "holm", "none". Default = "BH".
#' @param verbose Logical. Print progress and summary. Default = TRUE.
#'
#' @return A list with class "artemis_characterization" containing:
#' \describe{
#'   \item{results}{Data frame with per-variable statistics, sorted by effect size}
#'   \item{top_variables}{Names of variables with large effect sizes}
#'   \item{cluster_profiles}{Summary statistics per cluster for top variables}
#'   \item{n_clusters}{Number of clusters}
#'   \item{n_variables}{Number of variables tested}
#' }
#'
#' @details
#' For quantitative variables:
#' \itemize{
#'   \item Test: Kruskal-Wallis rank sum test (non-parametric ANOVA)
#'   \item Effect size: Eta-squared (η²) = H / (n-1), where H is the test statistic
#'   \item Interpretation: 0.01 = small, 0.06 = medium, 0.14 = large
#' }
#'
#' For qualitative variables:
#' \itemize{
#'   \item Test: Chi-squared test of independence
#'   \item Effect size: Cramér's V = sqrt(χ² / (n * min(r-1, c-1)))
#'   \item Interpretation: 0.1 = small, 0.3 = medium, 0.5 = large
#' }
#'
#' Variables are ranked by effect size. Those with high effect sizes are
#' the most discriminative - they differ substantially between clusters.
#'
#' @export
#'
#' @examples
#' # Characterize clusters
#' char <- ARTEMIS_characterize_clusters(cluster_result)
#'
#' # View top discriminating variables
#' head(char$results, 10)
#'
#' # Get variable names with large effects
#' char$top_variables
#'
ARTEMIS_characterize_clusters <- function(cluster_result,
                                           data = NULL,
                                           p_adjust = "BH",
                                           verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(cluster_result, "artemis_cluster")) {
    stop("cluster_result must be an artemis_cluster object")
  }

  if (is.null(data)) {
    data <- cluster_result$data_used
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  clusters <- cluster_result$clusters
  k <- cluster_result$k

  if (length(clusters) != nrow(data)) {
    stop("Number of cluster assignments doesn't match number of rows in data")
  }

  if (!p_adjust %in% c("BH", "bonferroni", "holm", "none")) {
    stop("p_adjust must be one of: 'BH', 'bonferroni', 'holm', 'none'")
  }

  if (verbose) {
    cat("ARTEMIS Cluster Characterization\n")
    cat("=================================\n")
    cat("Clusters:", k, "\n")
    cat("Variables:", ncol(data), "\n")
    cat("Samples:", nrow(data), "\n\n")
  }

  # ---------------------------------------------------------------------------
  # Test each variable
  # ---------------------------------------------------------------------------
  if (verbose) cat("Testing variable associations with clusters...\n")

  results_list <- list()

  for (i in seq_len(ncol(data))) {
    var_name <- names(data)[i]
    var_data <- data[[i]]

    # Determine variable type
    if (is.numeric(var_data)) {
      var_type <- "quantitative"
      test_result <- .test_quanti_cluster(var_data, clusters)
    } else {
      var_type <- "qualitative"
      # Ensure factor
      if (!is.factor(var_data)) {
        var_data <- as.factor(var_data)
      }
      test_result <- .test_quali_cluster(var_data, clusters)
    }

    results_list[[i]] <- data.frame(
      variable = var_name,
      type = var_type,
      test = test_result$test,
      statistic = test_result$statistic,
      p_value = test_result$p_value,
      effect_size = test_result$effect_size,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, results_list)

  # ---------------------------------------------------------------------------
  # Adjust p-values
  # ---------------------------------------------------------------------------
  if (p_adjust != "none") {
    results$p_adjusted <- p.adjust(results$p_value, method = p_adjust)
  } else {
    results$p_adjusted <- results$p_value
  }

  # ---------------------------------------------------------------------------
  # Add effect size interpretation
  # ---------------------------------------------------------------------------
  results$effect_interpretation <- mapply(function(es, type) {
    if (is.na(es)) return("NA")

    if (type == "quantitative") {
      # Eta-squared interpretation
      if (es >= 0.14) return("Large")
      else if (es >= 0.06) return("Medium")
      else if (es >= 0.01) return("Small")
      else return("Negligible")
    } else {
      # Cramér's V interpretation
      if (es >= 0.5) return("Large")
      else if (es >= 0.3) return("Medium")
      else if (es >= 0.1) return("Small")
      else return("Negligible")
    }
  }, results$effect_size, results$type)

  # ---------------------------------------------------------------------------
  # Sort by effect size (descending)
  # ---------------------------------------------------------------------------
  results <- results[order(-results$effect_size), ]
  rownames(results) <- NULL

  # ---------------------------------------------------------------------------
  # Identify top variables
  # ---------------------------------------------------------------------------
  # Variables with at least medium effect size (handle NA safely)
  top_quanti <- results$variable[results$type == "quantitative" &
                                   !is.na(results$effect_size) &
                                   results$effect_size >= 0.06]
  top_quali <- results$variable[results$type == "qualitative" &
                                  !is.na(results$effect_size) &
                                  results$effect_size >= 0.3]
  top_variables <- c(top_quanti, top_quali)

  # Also include any with "Large" interpretation (handle NA safely)
  large_effect <- results$variable[!is.na(results$effect_interpretation) &
                                     results$effect_interpretation == "Large"]
  top_variables <- unique(c(large_effect, top_variables))

  # Remove any NA values that might have slipped through
  top_variables <- top_variables[!is.na(top_variables)]

  # ---------------------------------------------------------------------------
  # Generate cluster profiles for top variables
  # ---------------------------------------------------------------------------
  cluster_profiles <- list()

  for (var_name in head(top_variables, 10)) {  # Top 10 max
    var_data <- data[[var_name]]

    if (is.numeric(var_data)) {
      # Quantitative: mean and SD per cluster
      profile <- tapply(var_data, clusters, function(x) {
        c(mean = mean(x, na.rm = TRUE),
          sd = sd(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE))
      })
      cluster_profiles[[var_name]] <- do.call(rbind, profile)
    } else {
      # Qualitative: proportion table
      cluster_profiles[[var_name]] <- prop.table(table(clusters, var_data), margin = 1)
    }
  }

  # ---------------------------------------------------------------------------
  # Summary output
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\nTop discriminating variables:\n")
    cat("─────────────────────────────────────────────────────────────────\n")

    n_show <- min(15, nrow(results))
    for (i in seq_len(n_show)) {
      row <- results[i, ]
      sig_marker <- if (!is.na(row$p_adjusted) && row$p_adjusted < 0.05) "*" else " "
      cat(sprintf("  %2d. %-25s %s  effect=%.3f (%s)%s\n",
                  i, row$variable, row$type,
                  row$effect_size, row$effect_interpretation, sig_marker))
    }
    if (nrow(results) > 15) {
      cat("  ... and", nrow(results) - 15, "more variables\n")
    }

    cat("\n* = significant after adjustment (p < 0.05)\n")
    cat("\nVariables with large effect:", length(large_effect), "\n")
    cat("Variables with medium+ effect:", length(top_variables), "\n")
  }

  # ---------------------------------------------------------------------------
  # Build result
  # ---------------------------------------------------------------------------
  result <- list(
    results = results,
    top_variables = top_variables,
    cluster_profiles = cluster_profiles,
    n_clusters = k,
    n_variables = ncol(data),
    clusters = clusters
  )

  class(result) <- c("artemis_characterization", "list")

  return(result)
}


#' Test quantitative variable against clusters
#' @keywords internal
.test_quanti_cluster <- function(x, clusters) {
  # Kruskal-Wallis test (non-parametric)
  tryCatch({
    kw <- kruskal.test(x ~ clusters)

    # Eta-squared approximation: H / (n - 1)
    # This is epsilon-squared, a common effect size for Kruskal-Wallis
    n <- length(x)
    eta_sq <- kw$statistic / (n - 1)

    list(
      test = "kruskal",
      statistic = kw$statistic,
      p_value = kw$p.value,
      effect_size = as.numeric(eta_sq)
    )
  }, error = function(e) {
    list(
      test = "kruskal",
      statistic = NA,
      p_value = NA,
      effect_size = NA
    )
  })
}


#' Test qualitative variable against clusters
#' @keywords internal
.test_quali_cluster <- function(x, clusters) {
  # Chi-squared test
  tryCatch({
    tbl <- table(clusters, x)

    # Check for minimum expected counts
    chi <- chisq.test(tbl)

    # Cramér's V
    n <- sum(tbl)
    min_dim <- min(nrow(tbl) - 1, ncol(tbl) - 1)
    if (min_dim > 0) {
      cramers_v <- sqrt(chi$statistic / (n * min_dim))
    } else {
      cramers_v <- NA
    }

    list(
      test = "chi-sq",
      statistic = chi$statistic,
      p_value = chi$p.value,
      effect_size = as.numeric(cramers_v)
    )
  }, error = function(e) {
    list(
      test = "chi-sq",
      statistic = NA,
      p_value = NA,
      effect_size = NA
    )
  })
}


#' Print method for cluster characterization
#' @export
print.artemis_characterization <- function(x, n = 20, ...) {
  cat("ARTEMIS Cluster Characterization\n")
  cat("=================================\n")
  cat("Clusters:", x$n_clusters, "\n")
  cat("Variables tested:", x$n_variables, "\n")
  cat("Variables with large effect:", sum(x$results$effect_interpretation == "Large"), "\n")
  cat("Variables with medium effect:", sum(x$results$effect_interpretation == "Medium"), "\n")

  cat("\nTop", min(n, nrow(x$results)), "variables by effect size:\n\n")

  # Format as table
  show_df <- head(x$results, n)
  show_df$effect_size <- round(show_df$effect_size, 3)
  show_df$p_adjusted <- format.pval(show_df$p_adjusted, digits = 2)

  print(show_df[, c("variable", "type", "effect_size", "effect_interpretation", "p_adjusted")],
        row.names = FALSE)

  invisible(x)
}


#' Compare Variable Importance from Clustering and FAMD
#'
#' @description Combines variable importance metrics from cluster characterization
#' (effect sizes) and FAMD (contributions) to identify variables that may be
#' candidates for removal during variable refinement.
#'
#' Variables that score low on BOTH metrics are likely not contributing
#' meaningfully to either cluster discrimination or overall data structure.
#'
#' @param cluster_char An artemis_characterization object from
#'   ARTEMIS_characterize_clusters().
#' @param famd_result An artemis_famd object from ARTEMIS_famd().
#' @param n_dims Integer. Number of top FAMD dimensions to consider for
#'   contribution averaging. Default = 3.
#' @param low_threshold Numeric (0-1). Percentile below which a variable is
#'   considered "low" on a metric. Default = 0.25 (bottom 25%).
#' @param verbose Logical. Print summary. Default = TRUE.
#'
#' @return A list with class "artemis_variable_comparison" containing:
#' \describe{
#'   \item{comparison}{Data frame with all variables, their metrics, ranks, and flags}
#'   \item{low_signal}{Character vector of variables low on BOTH metrics (removal candidates)}
#'   \item{high_signal}{Character vector of variables high on at least one metric (keep)}
#'   \item{summary}{Summary statistics}
#' }
#'
#' @details
#' The function computes percentile ranks for each variable on both metrics:
#' \itemize{
#'   \item Effect size from clustering (eta² for quantitative, Cramér's V for qualitative)
#'   \item Average contribution to top FAMD dimensions
#' }
#'
#' For qualitative variables, FAMD contributions are aggregated from category-level
#' to variable-level by summing contributions across all categories of each variable.
#'
#' Variables in the bottom percentile (controlled by low_threshold) on BOTH metrics
#' are flagged as candidates for removal. This is conservative - a variable only
#' needs to be important on ONE metric to be retained.
#'
#' @export
#'
#' @examples
#' # Run clustering and characterization
#' clust <- ARTEMIS_cluster_mixed(data)
#' char <- ARTEMIS_characterize_clusters(clust)
#' famd <- ARTEMIS_famd(data)
#'
#' # Compare importance metrics
#' comparison <- ARTEMIS_compare_variable_importance(char, famd)
#'
#' # See which variables to consider removing
#' comparison$low_signal
#'
#' # Iteratively refine
#' data_refined <- data[, !names(data) %in% comparison$low_signal]
#'
ARTEMIS_compare_variable_importance <- function(cluster_char,
                                                  famd_result,
                                                  n_dims = 3,
                                                  low_threshold = 0.25,
                                                  verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(cluster_char, "artemis_characterization")) {
    stop("cluster_char must be an artemis_characterization object from ARTEMIS_characterize_clusters()")
  }

  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  if (low_threshold < 0 || low_threshold > 1) {
    stop("low_threshold must be between 0 and 1")
  }

  # Limit n_dims to available dimensions
  max_dims <- ncol(famd_result$var$quanti$contrib)
  n_dims <- min(n_dims, max_dims)

  if (verbose) {
    cat("ARTEMIS Variable Importance Comparison\n")
    cat("=======================================\n")
  }

  # ---------------------------------------------------------------------------
  # Extract effect sizes from cluster characterization
  # ---------------------------------------------------------------------------
  effect_df <- cluster_char$results[, c("variable", "type", "effect_size")]
  names(effect_df)[3] <- "cluster_effect"

  # ---------------------------------------------------------------------------
  # Extract and aggregate FAMD contributions
  # ---------------------------------------------------------------------------
  # Quantitative variables: straightforward
  quanti_contrib <- famd_result$var$quanti$contrib[, 1:n_dims, drop = FALSE]
  quanti_avg <- rowMeans(quanti_contrib)

  quanti_df <- data.frame(
    variable = names(quanti_avg),
    famd_contrib = as.numeric(quanti_avg),
    stringsAsFactors = FALSE
  )

  # Qualitative variables: aggregate from category-level to variable-level
  quali_contrib <- famd_result$var$quali$contrib[, 1:n_dims, drop = FALSE]
  category_names <- rownames(quali_contrib)

  # Get original variable names
  quali_vars <- names(famd_result$data_used)[sapply(famd_result$data_used, is.factor)]

  # Map categories to parent variables and sum contributions
  quali_agg <- list()
  for (v in quali_vars) {
    # Match categories belonging to this variable (handles prefixed levels)
    pattern <- paste0("^", v, "_")
    idx <- grep(pattern, category_names)

    if (length(idx) > 0) {
      # Sum contributions across categories and dimensions, then average by dims
      total_contrib <- sum(rowMeans(quali_contrib[idx, , drop = FALSE]))
      quali_agg[[v]] <- total_contrib
    }
  }

  if (length(quali_agg) > 0) {
    quali_df <- data.frame(
      variable = names(quali_agg),
      famd_contrib = unlist(quali_agg),
      stringsAsFactors = FALSE
    )
  } else {
    quali_df <- data.frame(
      variable = character(0),
      famd_contrib = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Combine quantitative and qualitative
  famd_df <- rbind(quanti_df, quali_df)

  # ---------------------------------------------------------------------------
  # Merge effect sizes and FAMD contributions
  # ---------------------------------------------------------------------------
  comparison <- merge(effect_df, famd_df, by = "variable", all = TRUE)

  # Handle any NAs (variables in one but not other - shouldn't happen normally)
  comparison$cluster_effect[is.na(comparison$cluster_effect)] <- 0
  comparison$famd_contrib[is.na(comparison$famd_contrib)] <- 0

  # ---------------------------------------------------------------------------
  # Compute percentile ranks
  # ---------------------------------------------------------------------------
  comparison$cluster_rank <- rank(comparison$cluster_effect) / nrow(comparison)
  comparison$famd_rank <- rank(comparison$famd_contrib) / nrow(comparison)

  # Average rank (simple combined score)
  comparison$avg_rank <- (comparison$cluster_rank + comparison$famd_rank) / 2

  # ---------------------------------------------------------------------------
  # Flag low-signal variables
  # ---------------------------------------------------------------------------
  comparison$low_cluster <- comparison$cluster_rank <= low_threshold
  comparison$low_famd <- comparison$famd_rank <= low_threshold
  comparison$low_both <- comparison$low_cluster & comparison$low_famd

  # Sort by average rank (lowest first = removal candidates at top)
  comparison <- comparison[order(comparison$avg_rank), ]
  rownames(comparison) <- NULL

  # ---------------------------------------------------------------------------
  # Extract variable lists
  # ---------------------------------------------------------------------------
  low_signal <- comparison$variable[comparison$low_both]
  high_signal <- comparison$variable[!comparison$low_both]

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------
  summary_stats <- list(
    n_variables = nrow(comparison),
    n_low_cluster = sum(comparison$low_cluster),
    n_low_famd = sum(comparison$low_famd),
    n_low_both = sum(comparison$low_both),
    low_threshold = low_threshold,
    n_dims_used = n_dims
  )

  if (verbose) {
    cat("Variables analyzed:", summary_stats$n_variables, "\n")
    cat("Dimensions used for FAMD:", n_dims, "\n")
    cat("Low threshold (percentile):", low_threshold * 100, "%\n\n")

    cat("Variables below threshold:\n")
    cat("  Low cluster effect only:", sum(comparison$low_cluster & !comparison$low_famd), "\n")
    cat("  Low FAMD contrib only:", sum(comparison$low_famd & !comparison$low_cluster), "\n")
    cat("  Low on BOTH (removal candidates):", summary_stats$n_low_both, "\n\n")

    if (length(low_signal) > 0) {
      cat("Candidates for removal:\n")
      n_show <- min(20, length(low_signal))
      for (i in seq_len(n_show)) {
        v <- low_signal[i]
        row <- comparison[comparison$variable == v, ]
        cat(sprintf("  %2d. %-25s cluster=%.3f  famd=%.2f\n",
                    i, v, row$cluster_effect, row$famd_contrib))
      }
      if (length(low_signal) > 20) {
        cat("  ... and", length(low_signal) - 20, "more\n")
      }
    } else {
      cat("No variables flagged as low on both metrics.\n")
      cat("Consider lowering low_threshold or reviewing manually.\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Build result
  # ---------------------------------------------------------------------------
  result <- list(
    comparison = comparison,
    low_signal = low_signal,
    high_signal = high_signal,
    summary = summary_stats
  )

  class(result) <- c("artemis_variable_comparison", "list")

  return(result)
}


#' Print method for variable comparison
#' @export
print.artemis_variable_comparison <- function(x, n = 30, ...) {
  cat("ARTEMIS Variable Importance Comparison\n")
  cat("=======================================\n")
  cat("Variables:", x$summary$n_variables, "\n")
  cat("Low threshold:", x$summary$low_threshold * 100, "percentile\n")
  cat("Removal candidates (low on both):", x$summary$n_low_both, "\n\n")

  cat("Variable rankings (sorted by combined score):\n\n")

  show_df <- head(x$comparison, n)
  show_df$cluster_effect <- round(show_df$cluster_effect, 3)
  show_df$famd_contrib <- round(show_df$famd_contrib, 2)
  show_df$avg_rank <- round(show_df$avg_rank, 2)

  # Add flag column
  show_df$flag <- ifelse(show_df$low_both, "REMOVE?", "")

  print(show_df[, c("variable", "type", "cluster_effect", "famd_contrib", "avg_rank", "flag")],
        row.names = FALSE)

  if (nrow(x$comparison) > n) {
    cat("\n... and", nrow(x$comparison) - n, "more variables\n")
  }

  invisible(x)
}


#' Evaluate Clustering Across a Range of k Values
#'
#' @description Systematically evaluates clustering quality across multiple k
#' values and tracks how samples move between clusters. This helps identify
#' the optimal number of clusters and visualize cluster stability.
#'
#' @param data Data frame with mixed variable types, OR a pre-computed distance
#'   matrix (class "dist").
#' @param k_range Integer vector. Range of k values to evaluate. Default = 2:8.
#' @param method Character. Clustering method: "pam" or "hierarchical".
#'   Default = "pam".
#' @param hclust_method Character. Linkage method for hierarchical clustering.
#'   Default = "ward.D2".
#' @param stand Logical. Standardize variables before computing Gower distance?
#'   Default = FALSE.
#' @param verbose Logical. Print progress messages. Default = TRUE.
#'
#' @return A list with class "artemis_k_evaluation" containing:
#' \describe{
#'   \item{metrics}{Data frame with k, silhouette_avg, and cluster size info}
#'   \item{assignments}{Data frame with sample assignments at each k (wide format)}
#'   \item{transitions}{Data frame tracking sample flows between consecutive k values}
#'   \item{optimal_k}{Suggested optimal k based on silhouette}
#'   \item{distance}{The distance matrix used}
#'   \item{k_range}{The k values evaluated}
#' }
#'
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Determining optimal number of clusters (via silhouette scores)
#'   \item Understanding cluster stability (do samples stay together as k increases?)
#'   \item Visualizing sample flows with AETHER_plot_cluster_sankey()
#' }
#'
#' The transitions output tracks how clusters at k split into clusters at k+1,
#' which can reveal natural substructure or artificial splits.
#'
#' @export
#'
#' @examples
#' # Evaluate k from 2 to 6
#' k_eval <- ARTEMIS_evaluate_k_range(my_data, k_range = 2:6)
#'
#' # View metrics
#' print(k_eval)
#'
#' # Visualize transitions
#' AETHER_plot_cluster_sankey(k_eval)
#'
ARTEMIS_evaluate_k_range <- function(data,
                                      k_range = 2:8,
                                      method = "pam",
                                      hclust_method = "ward.D2",
                                      stand = FALSE,
                                      verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Install with: install.packages('cluster')")
  }

  # Check if data is already a distance matrix
  if (inherits(data, "dist")) {
    gower_dist <- data
    n_samples <- attr(data, "Size")
    sample_names <- attr(data, "Labels")
    if (is.null(sample_names)) {
      sample_names <- paste0("S", seq_len(n_samples))
    }
  } else {
    # Convert to data frame if needed
    if (!is.data.frame(data)) {
      data <- as.data.frame(data)
    }

    n_samples <- nrow(data)
    sample_names <- rownames(data)
    if (is.null(sample_names)) {
      sample_names <- paste0("S", seq_len(n_samples))
    }

    # Check for NA values
    if (any(is.na(data))) {
      stop("Data contains NA values. Please handle missing data before clustering.")
    }

    # Convert character columns to factors
    char_cols <- sapply(data, is.character)
    if (any(char_cols)) {
      data[char_cols] <- lapply(data[char_cols], as.factor)
    }

    if (verbose) {
      cat("ARTEMIS k-Range Evaluation\n")
      cat("==========================\n")
      cat("Input:", n_samples, "samples,", ncol(data), "variables\n\n")
    }

    # Compute Gower distance
    if (verbose) cat("Computing Gower distance matrix...\n")
    gower_dist <- cluster::daisy(data, metric = "gower", stand = stand)
  }

  # Validate k_range
  max_k <- min(max(k_range), floor(n_samples / 2), 15)
  k_range <- sort(unique(k_range[k_range <= max_k & k_range >= 2]))

  if (length(k_range) == 0) {
    stop("No valid k values in k_range for this dataset size")
  }

  if (verbose) {
    cat("Evaluating k =", paste(k_range, collapse = ", "), "\n\n")
  }

  # ---------------------------------------------------------------------------
  # Evaluate each k
  # ---------------------------------------------------------------------------
  metrics_list <- list()
  assignments <- data.frame(sample = sample_names, stringsAsFactors = FALSE)

  for (k in k_range) {
    if (verbose) cat("  k =", k, "... ")

    # Cluster
    if (method == "pam") {
      clust_result <- cluster::pam(gower_dist, k = k, diss = TRUE)
      clusters <- clust_result$clustering
    } else {
      hc <- hclust(gower_dist, method = hclust_method)
      clusters <- cutree(hc, k = k)
    }

    # Compute silhouette
    sil_info <- cluster::silhouette(clusters, gower_dist)
    sil_avg <- mean(sil_info[, "sil_width"])

    # Store metrics
    cluster_sizes <- table(clusters)
    size_str <- paste(as.numeric(cluster_sizes), collapse = "/")

    metrics_list[[as.character(k)]] <- data.frame(
      k = k,
      silhouette_avg = sil_avg,
      cluster_sizes = size_str,
      min_cluster_size = min(cluster_sizes),
      max_cluster_size = max(cluster_sizes),
      stringsAsFactors = FALSE
    )

    # Store assignments
    col_name <- paste0("k", k)
    assignments[[col_name]] <- clusters

    if (verbose) {
      cat("silhouette =", round(sil_avg, 3),
          " sizes:", size_str, "\n")
    }
  }

  metrics <- do.call(rbind, metrics_list)
  rownames(metrics) <- NULL

  # ---------------------------------------------------------------------------
  # Build transitions data for Sankey diagram
  # ---------------------------------------------------------------------------
  transitions_list <- list()

  for (i in seq_len(length(k_range) - 1)) {
    k_from <- k_range[i]
    k_to <- k_range[i + 1]

    col_from <- paste0("k", k_from)
    col_to <- paste0("k", k_to)

    # Count transitions
    trans_table <- table(
      from = paste0("k", k_from, "_C", assignments[[col_from]]),
      to = paste0("k", k_to, "_C", assignments[[col_to]])
    )

    trans_df <- as.data.frame(trans_table, stringsAsFactors = FALSE)
    names(trans_df) <- c("source", "target", "value")
    trans_df <- trans_df[trans_df$value > 0, ]

    trans_df$k_from <- k_from
    trans_df$k_to <- k_to

    transitions_list[[i]] <- trans_df
  }

  transitions <- do.call(rbind, transitions_list)
  rownames(transitions) <- NULL

  # ---------------------------------------------------------------------------
  # Determine optimal k
  # ---------------------------------------------------------------------------
  optimal_k <- metrics$k[which.max(metrics$silhouette_avg)]

  if (verbose) {
    cat("\n")
    cat("Optimal k =", optimal_k,
        "(silhouette =", round(max(metrics$silhouette_avg), 3), ")\n")
  }

  # ---------------------------------------------------------------------------
  # Build result
  # ---------------------------------------------------------------------------
  result <- list(
    metrics = metrics,
    assignments = assignments,
    transitions = transitions,
    optimal_k = optimal_k,
    distance = gower_dist,
    k_range = k_range,
    method = method
  )

  class(result) <- c("artemis_k_evaluation", "list")

  return(result)
}


#' Print method for k evaluation
#' @export
print.artemis_k_evaluation <- function(x, ...) {
  cat("ARTEMIS k-Range Evaluation\n")
  cat("==========================\n")
  cat("k values tested:", paste(x$k_range, collapse = ", "), "\n")
  cat("Optimal k:", x$optimal_k,
      "(silhouette =", round(max(x$metrics$silhouette_avg), 3), ")\n\n")

  cat("Metrics by k:\n")
  print_df <- x$metrics
  print_df$silhouette_avg <- round(print_df$silhouette_avg, 3)
  print_df$optimal <- ifelse(print_df$k == x$optimal_k, " <--", "")
  print(print_df, row.names = FALSE)

  invisible(x)
}
