# GAIA/Artemis/part_clustering.R
# PART Clustering - Partitioning Algorithm based on Recursive Thresholding
#
# Standalone implementation of the PART algorithm for detecting multi-scale
# cluster structure in genomic data. Uses the Gap statistic to determine
# optimal k at each recursive step and a dendrogram height threshold for
# the splitting criterion.
#
# Reference: Nilsen et al. (2013) "Identifying clusters in genomics data
#            by recursive partitioning"
# Original: clusterGenomics package (Gro Nilsen)
# Ported from: to_process/TiSA/R/clusterGenomics.R + PART.R


# ==============================================================================
# PUBLIC API
# ==============================================================================

#' PART Clustering
#'
#' Performs clustering using the PART algorithm (Partitioning Algorithm based
#' on Recursive Thresholding). The method recursively evaluates whether to
#' split clusters using the Gap statistic and a dendrogram height threshold.
#' It detects multi-scale cluster structure while guarding against spurious
#' splits.
#'
#' @param mat Numeric matrix with features (genes) as rows and samples as
#'   columns. Rownames are used as feature identifiers.
#' @param q Numeric (0-1). Tuning parameter controlling splitting aggressiveness.
#'   Threshold = (1-q) quantile of dendrogram heights. Lower q = less splitting.
#'   Default: 0.25.
#' @param min_size Integer. Minimum features per cluster. Clusters smaller than
#'   this are marked as outliers (cluster C0). Default: 8.
#' @param B Integer. Number of bootstrap reference datasets for the Gap
#'   statistic. Higher = more stable but slower. Default: 100.
#' @param Kmax Integer. Maximum clusters to consider per recursive step.
#'   Default: 10.
#' @param dist_method Character. Distance metric: "euclidean", "sq.euclidean"
#'   (squared Euclidean), "correlation" (1 - Pearson), "manhattan".
#'   Default: "euclidean".
#' @param linkage Character. Hierarchical clustering linkage method:
#'   "average", "complete", "ward.D2", "single", etc. Default: "average".
#' @param scale Logical. Z-score scale each row (feature) before clustering.
#'   Default: FALSE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_part"} containing:
#'   \describe{
#'     \item{clusters}{Named character vector of cluster assignments (C1, C2, ...; C0 for outliers)}
#'     \item{cluster_map}{Data.frame with gene, cluster, cluster_color columns}
#'     \item{n_clusters}{Number of clusters (excluding outliers)}
#'     \item{n_outliers}{Number of outlier features}
#'     \item{cluster_sizes}{Named integer vector of cluster sizes}
#'     \item{cluster_colors}{Named character vector of cluster hex colors}
#'     \item{dendrogram}{hclust object for the full dataset}
#'     \item{data}{The (possibly scaled) input matrix, rows ordered by cluster}
#'     \item{parameters}{List of all parameters used}
#'     \item{computation_time}{Elapsed time in seconds}
#'   }
#'
#' @details
#' The PART algorithm works as follows:
#' \enumerate{
#'   \item Compute a stopping threshold from the (1-q) quantile of dendrogram heights
#'   \item Apply the Gap statistic to determine optimal k
#'   \item If k > 1: accept the split, recursively apply PART to each cluster
#'   \item If k = 1 but dendrogram height exceeds threshold: try splitting into 2, recurse
#'   \item If k = 1 and height < threshold: stop (terminal cluster)
#'   \item Clusters smaller than min_size are marked as outliers (C0)
#' }
#'
#' After clustering, clusters are reordered hierarchically so that similar
#' clusters are adjacent, and labeled C1, C2, ... in order.
#'
#' @examples
#' \dontrun{
#' # Cluster gene expression matrix
#' part_result <- ARTEMIS_part(expr_matrix, scale = TRUE, seed = 42)
#'
#' # More aggressive splitting with larger minimum cluster size
#' part_result <- ARTEMIS_part(expr_matrix, q = 0.5, min_size = 50)
#'
#' # Use correlation distance
#' part_result <- ARTEMIS_part(expr_matrix, dist_method = "correlation")
#'
#' }
#' @export
ARTEMIS_part <- function(mat,
                          q = 0.25,
                          min_size = 8,
                          B = 100,
                          Kmax = 10,
                          dist_method = "euclidean",
                          linkage = "average",
                          scale = FALSE,
                          seed = NULL,
                          verbose = TRUE) {

  # --- Validate inputs ---
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("'mat' must be a matrix or data.frame")
  }
  if (is.data.frame(mat)) mat <- as.matrix(mat)

  if (is.null(rownames(mat))) {
    stop("'mat' must have rownames (feature identifiers)")
  }
  if (nrow(mat) < 2 * min_size) {
    stop("Too few features (", nrow(mat), ") for min_size = ", min_size,
         ". Need at least ", 2 * min_size, " features.")
  }

  if (q <= 0 || q >= 1) stop("'q' must be between 0 and 1 (exclusive)")
  if (min_size < 2) stop("'min_size' must be at least 2")

  if (verbose) {
    cat("[ARTEMIS] PART Clustering \n")
    cat("    Features:", nrow(mat), "| Samples:", ncol(mat), "\n")
    cat("    Parameters: q =", q, "| min_size =", min_size,
        "| B =", B, "| Kmax =", Kmax, "\n")
    cat("    Distance:", dist_method, "| Linkage:", linkage, "\n")
  }

  # --- Set seed ---
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) cat("    Seed:", seed, "\n")
  }

  # --- Scale if requested ---
  if (scale) {
    mat <- t(scale(t(mat), center = TRUE, scale = TRUE))
    # Handle constant rows (NaN from zero variance)
    nan_rows <- rowSums(is.nan(mat)) > 0
    if (any(nan_rows)) {
      if (verbose) cat("    Removing", sum(nan_rows),
                        "constant features (zero variance after scaling)\n")
      mat <- mat[!nan_rows, , drop = FALSE]
    }
  }

  # --- Build parameter list for internal functions ---
  params <- list(
    q = q,
    min_size = min_size,
    B = B,
    Kmax = Kmax,
    Kmax_rec = 5,
    dist_method = dist_method,
    cl_method = "hclust",
    linkage = linkage,
    cor_method = "pearson",
    nstart = 10,
    ref_gen = "PC",
    min_dist = NULL
  )

  # --- Compute threshold ---
  if (verbose) cat("[ARTEMIS] Computing threshold...\n")
  params$min_dist <- .part_get_threshold(mat, q, params)

  # --- Run recursive PART ---
  if (verbose) cat("[ARTEMIS] Running PART recursion...\n")
  start_time <- proc.time()

  cluster_matrix <- .part_recursive(
    X = mat,
    Kmax = Kmax,
    ind = rep(1, nrow(mat)),
    cl_lab = NULL,
    params = params,
    verbose = verbose,
    depth = 0
  )

  # --- Assign labels ---
  labels <- .part_get_labels(cluster_matrix, min_size)

  # Count outliers
  outlier_idx <- which(labels == 0)
  n_outliers <- length(outlier_idx)
  if (n_outliers > 0) {
    n_clusters <- length(unique(labels[-outlier_idx]))
  } else {
    n_clusters <- length(unique(labels))
  }

  if (verbose) {
    cat("[ARTEMIS] Raw PART found", n_clusters, "clusters")
    if (n_outliers > 0) cat(" +", n_outliers, "outliers")
    cat("\n")
  }

  # --- Reorder clusters hierarchically ---
  if (verbose) cat("[ARTEMIS] Reordering clusters...\n")
  reordered <- .part_reorder_clusters(labels, mat, dist_method, linkage)

  elapsed <- (proc.time() - start_time)["elapsed"]

  if (verbose) {
    cat("    Done in", round(elapsed, 1), "seconds\n")
    cat("    Final:", reordered$n_clusters, "clusters")
    if (reordered$n_outliers > 0) {
      cat(",", reordered$n_outliers, "outliers (C0)")
    }
    cat("\n")
  }

  # --- Build result object ---
  # Order matrix rows by cluster
  ordered_genes <- rownames(reordered$cluster_map)
  ordered_mat <- mat[ordered_genes, , drop = FALSE]

  # Full dendrogram
  dX <- .part_get_dist(mat, dist_method, params$cor_method)
  full_dendro <- stats::hclust(dX, method = linkage)

  result <- list(
    clusters = reordered$clusters,
    cluster_map = reordered$cluster_map,
    n_clusters = reordered$n_clusters,
    n_outliers = reordered$n_outliers,
    cluster_sizes = reordered$cluster_sizes,
    cluster_colors = reordered$cluster_colors,
    dendrogram = full_dendro,
    data = ordered_mat,
    parameters = list(
      q = q,
      min_size = min_size,
      B = B,
      Kmax = Kmax,
      dist_method = dist_method,
      linkage = linkage,
      scale = scale,
      seed = seed,
      threshold = params$min_dist
    ),
    computation_time = elapsed
  )

  class(result) <- c("artemis_part", "list")
  return(result)
}


#' Print PART result summary
#'
#' @param part_result An \code{artemis_part} object.
#' @param verbose Logical. Print to console. Default: TRUE.
#'
#' @return Invisible NULL.
#' @export
ARTEMIS_part_summary <- function(part_result, verbose = TRUE) {

  if (!inherits(part_result, "artemis_part")) {
    stop("'part_result' must be an artemis_part object")
  }

  if (verbose) {
    cat("[ARTEMIS] PART Clustering Summary \n")
    cat("    Clusters:", part_result$n_clusters, "\n")
    cat("    Outliers:", part_result$n_outliers, "(C0)\n")
    cat("    Total features:", nrow(part_result$data), "\n")
    cat("    Samples:", ncol(part_result$data), "\n")
    cat("    Time:", round(part_result$computation_time, 1), "seconds\n")

    cat("\n    Parameters:\n")
    p <- part_result$parameters
    cat("        q =", p$q, "| min_size =", p$min_size,
        "| B =", p$B, "| Kmax =", p$Kmax, "\n")
    cat("        dist =", p$dist_method, "| linkage =", p$linkage,
        "| scale =", p$scale, "\n")
    cat("        threshold =", round(p$threshold, 4), "\n")

    cat("\n    Cluster sizes:\n")
    sizes <- part_result$cluster_sizes
    for (cl in names(sizes)) {
      cat("        ", cl, ":", sizes[cl], "features\n")
    }
  }

  invisible(NULL)
}


#' Print method for artemis_part
#' @param x An artemis_part object
#' @param ... Additional arguments (ignored)
#' @method print artemis_part
#' @export
print.artemis_part <- function(x, ...) {
  cat("PART clustering:", x$n_clusters, "clusters,",
      nrow(x$data), "features,", ncol(x$data), "samples\n")
  cat("------------------------------\n")
  if (x$n_outliers > 0) cat("    Outliers:", x$n_outliers, "(C0)\n")
  cat("Sizes:", paste(x$cluster_sizes, collapse = ", "), "\n")
  invisible(x)
}


# ==============================================================================
# INTERNAL: Core PART Algorithm
# ==============================================================================

#' Compute distance matrix
#' @noRd
.part_get_dist <- function(X, method, cor_method = "pearson") {
  if (method == "sq.euclidean") {
    dX <- stats::dist(X, method = "euclidean")^2
  } else if (method == "correlation") {
    dX <- stats::as.dist(1 - stats::cor(t(X), method = cor_method))
  } else {
    dX <- stats::dist(X, method = method)
  }
  return(dX)
}


#' Hierarchical clustering wrapper
#' @noRd
.part_do_hclust <- function(d, k, linkage) {
  cl <- stats::hclust(d, method = linkage)
  lab <- stats::cutree(cl, k)
  list(cl = cl, lab = lab)
}


#' K-means wrapper
#' @noRd
.part_do_kmeans <- function(X, k, nstart = 10) {
  cl <- stats::kmeans(X, k, nstart = nstart)
  lab <- cl$cluster
  list(cl = cl, lab = lab)
}


#' Generate cluster partitions for k = 1..Kmax
#' @noRd
.part_find_partition <- function(X, Kmax, dX = NULL, params) {
  if (is.null(dX) && params$cl_method != "kmeans") {
    dX <- .part_get_dist(X, method = params$dist_method, cor_method = params$cor_method)
  }

  cl_lab <- vector("list", Kmax)
  for (k in seq_len(Kmax)) {
    if (params$cl_method == "kmeans") {
      cl_lab[[k]] <- .part_do_kmeans(X, k, nstart = params$nstart)$lab
    } else {
      cl_lab[[k]] <- .part_do_hclust(dX, k = k, linkage = params$linkage)$lab
    }
  }
  return(cl_lab)
}


#' Compute within-cluster dispersion for k = 1..K
#' @noRd
.part_find_w <- function(dX, K, cl_lab) {
  W <- rep(0, K)
  n <- nrow(as.matrix(dX))
  dist_mat <- as.matrix(dX)

  # W for k = 1
  W[1] <- sum(dist_mat) / (2 * n)

  # W for k = 2..K
  for (k in 2:K) {
    lab <- cl_lab[[k]]
    for (i in seq_len(k)) {
      idx <- which(lab == i)
      if (length(idx) < 2) next
      d_k <- dist_mat[idx, idx]
      D_k <- sum(d_k)
      nk <- length(idx)
      W[k] <- W[k] + D_k / (2 * nk)
    }
  }
  return(W)
}


#' Generate uniform reference sample for a single variable
#' @noRd
.part_sim <- function(Xcol) {
  stats::runif(length(Xcol), min = min(Xcol), max = max(Xcol))
}


#' Generate reference within-cluster dispersions
#' @noRd
.part_get_reference_w <- function(X, Kmax, B, ref_gen, params) {
  Wb <- matrix(0, nrow = Kmax, ncol = B)

  # PCA-based reference: transform to PC space
  if (ref_gen == "PC") {
    m <- colMeans(X, na.rm = TRUE)
    Xc <- sweep(X, 2, m)
    s <- svd(Xc)
    newX <- Xc %*% s$v
  }

  for (b in seq_len(B)) {
    if (ref_gen == "PC") {
      U <- apply(newX, 2, .part_sim)
      Z1 <- U %*% t(s$v)
      Z <- sweep(Z1, 2, m, FUN = "+")
    } else {
      Z <- apply(X, 2, .part_sim)
    }

    dZ <- .part_get_dist(Z, method = params$dist_method, cor_method = params$cor_method)
    clW_lab <- .part_find_partition(X = Z, Kmax = Kmax, dX = dZ, params = params)
    Wb[, b] <- .part_find_w(dX = dZ, K = Kmax, cl_lab = clW_lab)
  }

  return(Wb)
}


#' Compute Gap statistic
#' @noRd
.part_gap <- function(X, Kmax, B, ref_gen, cl_lab = NULL, params) {
  n <- nrow(X)

  # Adjust Kmax if needed
  if (n <= Kmax) {
    if (params$cl_method == "hclust") {
      Kmax <- n
    } else {
      Kmax <- n - 1
    }
  }

  dX <- .part_get_dist(X, method = params$dist_method, cor_method = params$cor_method)

  if (is.null(cl_lab)) {
    cl_lab <- .part_find_partition(X = X, Kmax = Kmax, dX = dX, params = params)
  }

  # Observed W
  W <- .part_find_w(dX = dX, K = Kmax, cl_lab = cl_lab)

  # Reference W
  Wb <- .part_get_reference_w(X, Kmax, B, ref_gen, params)

  # Gap = mean(log W*) - log(W)
  L <- apply(log(Wb), 1, mean)
  gap_vals <- L - log(W)

  # Standard error
  sdk <- apply(log(Wb), 1, stats::sd) * sqrt((B - 1) / B)
  sk <- sqrt(1 + (1 / B)) * sdk

  # Gap criterion: smallest k where gap_k >= gap_{k+1} - sk_{k+1}
  if (Kmax < 2) {
    hatK <- 1
  } else {
    diff_vals <- gap_vals[1:(Kmax - 1)] - (gap_vals[2:Kmax] - sk[2:Kmax])
    pos_diff <- which(diff_vals >= 0)

    if (length(pos_diff) == 0 || all(is.na(diff_vals))) {
      hatK <- 1
    } else {
      hatK <- pos_diff[1]
    }
  }

  lab_hatK <- cl_lab[[hatK]]
  list(hatK = hatK, lab_hatK = lab_hatK, gap = gap_vals, sk = sk, W = W)
}


#' Compute stopping threshold from dendrogram heights
#' @noRd
.part_get_threshold <- function(X, q, params) {
  dX <- .part_get_dist(X, method = params$dist_method, cor_method = params$cor_method)
  cl <- stats::hclust(dX, method = params$linkage)
  h <- cl$height
  as.numeric(stats::quantile(h, probs = 1 - q))
}


#' Recursive PART splitting
#' @noRd
.part_recursive <- function(X, Kmax, ind, cl_lab = NULL, params,
                             verbose = FALSE, depth = 0) {

  # Check minimum size for splitting
  if (sum(ind) < (2 * params$min_size)) {
    return(ind)
  }

  n <- sum(ind)
  if (n <= Kmax) {
    if (params$cl_method == "hclust") {
      Kmax <- n
    } else {
      Kmax <- n - 1
    }
  }

  if (is.null(cl_lab)) {
    cl_lab <- .part_find_partition(X = X, Kmax = Kmax, dX = NULL, params = params)
  }

  # Gap statistic to find optimal k
  gap_res <- .part_gap(X = X, Kmax = Kmax, B = params$B,
                        ref_gen = params$ref_gen, cl_lab = cl_lab, params = params)
  hatK <- gap_res$hatK
  lab_hatK <- gap_res$lab_hatK

  # If hatK > 1, check that at least two clusters meet min_size
  if (hatK > 1) {
    if (sum(table(lab_hatK) >= params$min_size) < 2) {
      hatK <- 1
    }
  }

  # --- Case A: hatK == 1 ---
  if (hatK == 1 && length(cl_lab) == 1) {
    return(ind)
  }

  if (hatK == 1) {
    # Try splitting into 2 tentative clusters
    obs1 <- cl_lab[[2]] == 1
    obs2 <- cl_lab[[2]] == 2
    T1 <- X[obs1, , drop = FALSE]
    T2 <- X[obs2, , drop = FALSE]

    # Check dendrogram height
    dX <- .part_get_dist(X, method = params$dist_method, cor_method = params$cor_method)
    hc_res <- .part_do_hclust(dX, k = 1, linkage = params$linkage)$cl
    T_height <- max(hc_res$height)

    if (T_height > params$min_dist) {
      # Create index vectors for tentative clusters
      t1_ind <- ind
      t2_ind <- ind
      active <- which(ind == 1)
      t1_ind[active[!obs1]] <- 0
      t2_ind[active[!obs2]] <- 0

      # Recursive calls
      t1 <- .part_recursive(X = T1, Kmax = params$Kmax_rec, ind = t1_ind,
                             cl_lab = NULL, params = params,
                             verbose = verbose, depth = depth + 1)
      t2 <- .part_recursive(X = T2, Kmax = params$Kmax_rec, ind = t2_ind,
                             cl_lab = NULL, params = params,
                             verbose = verbose, depth = depth + 1)

      if (!is.matrix(t1) && !is.matrix(t2)) {
        return(ind)
      } else {
        return(cbind(t1, t2))
      }
    } else {
      return(ind)
    }
  }

  # --- Case B: hatK > 1 ---
  if (verbose && depth < 2) {
    cat("    Depth", depth, ": splitting into", hatK, "clusters (n =", n, ")\n")
  }

  res <- matrix(NA, nrow = length(ind), ncol = 0)
  for (k in seq_len(hatK)) {
    obs <- lab_hatK == k
    S <- X[obs, , drop = FALSE]
    ind_S <- ind
    active <- which(ind == 1)
    ind_S[active[!obs]] <- 0

    res <- cbind(res, .part_recursive(
      X = S, Kmax = params$Kmax_rec, ind = ind_S,
      cl_lab = NULL, params = params,
      verbose = verbose, depth = depth + 1
    ))
  }

  return(res)
}


#' Assign final cluster labels from recursive output
#' @noRd
.part_get_labels <- function(clusters, min_size) {
  if (!is.matrix(clusters)) {
    clusters <- as.matrix(clusters)
  }

  label <- rep(NA_integer_, nrow(clusters))
  id <- 1L

  for (j in seq_len(ncol(clusters))) {
    members <- which(clusters[, j] == 1)
    if (length(members) < min_size) {
      label[members] <- 0L  # outlier
    } else {
      label[members] <- id
      id <- id + 1L
    }
  }

  return(label)
}


#' Reorder clusters hierarchically and assign colors
#' @noRd
.part_reorder_clusters <- function(labels, mat, dist_method, cor_method = "pearson", linkage = "average") {
  names(labels) <- rownames(mat)

  # Hierarchical clustering of the full matrix for ordering
  dX <- .part_get_dist(mat, method = dist_method, cor_method = cor_method)
  row_clust <- stats::hclust(dX, method = linkage)

  # Get cluster order from dendrogram traversal
  ordered_labels <- labels[row_clust$order]
  # Unique clusters in dendrogram order (excluding outliers initially)
  unique_ordered <- unique(as.character(ordered_labels))

  # Build rename map: original ID → C1, C2, ... (C0 for outliers)
  rename_map <- character()
  cluster_idx <- 1
  for (val in unique_ordered) {
    if (val == "0") {
      rename_map[val] <- "C0"
    } else {
      rename_map[val] <- paste0("C", cluster_idx)
      cluster_idx <- cluster_idx + 1
    }
  }

  # Apply renaming
  new_labels <- rename_map[as.character(labels)]
  names(new_labels) <- names(labels)

  # Count clusters (excluding outliers)
  all_clusters <- sort(unique(new_labels))
  non_outlier <- all_clusters[all_clusters != "C0"]
  n_clusters <- length(non_outlier)
  n_outliers <- sum(new_labels == "C0")

  # Generate colors
  all_display <- if (n_outliers > 0) c(non_outlier, "C0") else non_outlier
  n_colors <- length(all_display)

  # Color palette (matching Set3 from original TiSA, extended via colorRampPalette)
  set3_base <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                 "#FDB462", "#B3DE69", "#FCCDE5")
  cols <- grDevices::colorRampPalette(set3_base)(n_colors)
  names(cols) <- all_display
  if (n_outliers > 0) cols["C0"] <- "#D9D9D9"  # grey for outliers

  # Build cluster_map ordered by cluster
  cluster_map <- data.frame(
    gene = names(new_labels),
    cluster = new_labels,
    cluster_color = cols[new_labels],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # Order: C1 genes, C2 genes, ..., C0 genes
  # Within each cluster, maintain dendrogram order
  dendro_order <- row_clust$order
  dendro_rank <- seq_along(dendro_order)
  names(dendro_rank) <- rownames(mat)[dendro_order]
  cluster_map$dendro_rank <- dendro_rank[cluster_map$gene]
  cluster_map <- cluster_map[order(cluster_map$cluster, cluster_map$dendro_rank), ]
  cluster_map$dendro_rank <- NULL
  rownames(cluster_map) <- cluster_map$gene

  # Cluster sizes (named integer vector)
  sizes_tbl <- table(new_labels)
  sizes_tbl <- sizes_tbl[order(names(sizes_tbl))]
  cluster_sizes <- setNames(as.integer(sizes_tbl), names(sizes_tbl))

  list(
    clusters = new_labels,
    cluster_map = cluster_map,
    n_clusters = n_clusters,
    n_outliers = n_outliers,
    cluster_sizes = cluster_sizes,
    cluster_colors = cols
  )
}
