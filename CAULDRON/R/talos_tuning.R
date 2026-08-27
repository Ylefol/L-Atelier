# ==============================================================================
# TALOS - Embedding Parameter Tuning
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Data-driven parameter sweeps for UMAP (n_neighbors × min_dist) and
# tSNE (perplexity).  Both functions return a talos_embedding_sweep object.
#
# Metrics used:
#   - KNN overlap  : mean Jaccard similarity between PCA and embedding
#                    neighbourhoods (local structure preservation)
#   - Silhouette   : mean silhouette width in 2D embedding space
#   - Composite    : average of both metrics after [0,1] normalisation
#
# Graph construction, clustering, and resolution tuning live in talos_graph.R.
# Core run_ functions live in talos_embedding.R.
#
# All functions prefixed: TALOS_
# ==============================================================================


#' Sweep number of PCA components
#'
#' Choosing the number of PCA components (\code{n_pcs}) affects every
#' downstream step: too few discard biological signal; too many introduce
#' technical noise into the graph and embedding.  This function sweeps a range
#' of \code{n_pcs} values and scores each choice on two complementary metrics:
#'
#' \describe{
#'   \item{Cumulative variance explained}{The percentage of total
#'     gene-expression variance captured by the first \code{n_pcs} components.
#'     Read directly from the stored PCA result — no recomputation required.
#'     Provides context for how much information each cutoff retains, but does
#'     not on its own indicate when additional components stop helping
#'     downstream tasks.}
#'   \item{Graph modularity}{The SNN graph is built from the first \code{n_pcs}
#'     PCA components (at fixed \code{k}, \code{type}, \code{method}, and
#'     \code{resolution}) and community detection is applied.  Graph modularity
#'     measures how cleanly the resulting partition separates relative to a
#'     random null.  This is the primary optimisation target: the \code{n_pcs}
#'     at the elbow of the modularity curve is the point of diminishing returns
#'     — adding more components no longer meaningfully improves cluster
#'     separation.  K, method, and resolution are held fixed; use
#'     \code{\link{TALOS_tune_k}} and \code{\link{TALOS_tune_resolution}}
#'     afterwards to optimise those.}
#' }
#'
#' The elbow in the modularity curve is detected via the kneedle method (point
#' of maximum perpendicular distance from the line connecting the first and
#' last values).
#'
#' @param sce A \code{SingleCellExperiment} with \code{"PCA"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_pca}} first).  The stored
#'   \code{percentVar} attribute is used for cumulative variance — the PCA
#'   result is not recomputed.
#' @param n_pcs_range Integer vector of component counts to test.  Values
#'   exceeding the number of computed PCs are silently dropped.
#'   Default \code{seq(5, 50, by = 5)}.
#' @param k Integer.  Number of nearest neighbours for SNN graph construction
#'   (held fixed).  Default \code{20}.
#' @param type Character.  Edge weight scheme for \code{scran::buildSNNGraph}:
#'   \code{"rank"} (default), \code{"number"}, or \code{"jaccard"}.
#' @param method Character.  Clustering algorithm: \code{"leiden"} (default)
#'   or \code{"louvain"}.
#' @param resolution Numeric.  Clustering resolution (held fixed).
#'   Default \code{1.0}.
#' @param n_runs Integer.  Number of community-detection runs per \code{n_pcs}
#'   value.  Modularity is averaged over all runs to dampen Leiden/Louvain
#'   stochasticity; the representative run (closest to mean modularity) is used
#'   for \code{n_clusters}.  Default \code{5L}.
#' @param seed Integer.  Random seed for the first run; subsequent runs use
#'   \code{seed + 1}, \code{seed + 2}, … Default \code{42L}.
#' @param verbose Logical.  Print per-step progress and a final summary.
#'   Default \code{TRUE}.
#'
#' @return A \code{talos_pc_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per \code{n_pcs} value:
#'       \code{n_pcs}, \code{cum_var} (cumulative \% variance explained),
#'       \code{modularity} (mean over \code{n_runs}), \code{n_clusters}.}
#'     \item{\code{plot}}{Two-panel line plot (cumulative variance explained
#'       and graph modularity vs \code{n_pcs}), with a red dashed vertical
#'       line marking the suggested value.}
#'     \item{\code{best_n_pcs}}{Suggested number of components — the elbow of
#'       the modularity curve.}
#'     \item{\code{params}}{List of sweep parameters for reproducibility.}
#'   }
#' @export
TALOS_tune_pc <- function(sce,
                           n_pcs_range = seq(5L, 50L, by = 5L),
                           k           = 20L,
                           type        = c("rank", "number", "jaccard"),
                           method      = c("leiden", "louvain"),
                           resolution  = 1.0,
                           n_runs      = 5L,
                           seed        = 42L,
                           verbose     = TRUE) {

  type   <- match.arg(type)
  method <- match.arg(method)
  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  pct_var <- attr(reducedDim(sce, "PCA"), "percentVar")
  if (is.null(pct_var))
    stop("No percentVar attribute found on PCA. Re-run TALOS_run_pca().",
         call. = FALSE)

  max_pcs     <- ncol(reducedDim(sce, "PCA"))
  n_pcs_range <- sort(unique(as.integer(n_pcs_range)))
  n_pcs_range <- n_pcs_range[n_pcs_range >= 2L & n_pcs_range <= max_pcs]

  if (length(n_pcs_range) == 0L)
    stop("No valid n_pcs values after clipping to available PCs (max: ",
         max_pcs, ").", call. = FALSE)

  pca_full <- reducedDim(sce, "PCA")
  n_runs   <- max(1L, as.integer(n_runs))

  if (verbose)
    cat(sprintf(
      "\u2500\u2500 TALOS: PC sweep %s\n  n_pcs range : %s\n  k          : %s  |  type: %s\n  method     : %s  |  resolution: %s\n  n_runs     : %s\n%s\n",
      strrep("\u2500", 37),
      paste(n_pcs_range, collapse = ", "),
      k, type, method, resolution, n_runs,
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results <- vector("list", length(n_pcs_range))

  for (i in seq_along(n_pcs_range)) {
    np <- n_pcs_range[i]
    if (verbose)
      cat(sprintf("  (%d/%d) n_pcs = %d ...\n", i, length(n_pcs_range), np))

    pca_i <- pca_full[, seq_len(np), drop = FALSE]
    cum_v <- sum(pct_var[seq_len(np)])

    g   <- scran::buildSNNGraph(t(pca_i), k = as.integer(k), type = type)
    run <- .talos_multi_run_communities(g, method, resolution, seed, n_runs)

    results[[i]] <- data.frame(
      n_pcs      = np,
      cum_var    = round(cum_v,       2),
      modularity = round(run$modularity, 4),
      n_clusters = run$n_clusters
    )
  }

  res_df <- do.call(rbind, results)

  # ── Elbow in modularity curve ─────────────────────────────────────────────────
  best_n_pcs <- .talos_pc_elbow(res_df$n_pcs, res_df$modularity)

  p <- .talos_tune_pc_plot(res_df, best_n_pcs)

  if (verbose)
    cat(sprintf(
      "%s\n  Suggested n_pcs: %d (elbow in graph modularity)\n%s\n",
      strrep("\u2500", 56), best_n_pcs, strrep("\u2500", 56)))

  structure(
    list(results    = res_df,
         plot       = p,
         best_n_pcs = best_n_pcs,
         params     = list(n_pcs_range = n_pcs_range,
                           k           = as.integer(k),
                           type        = type,
                           method      = method,
                           resolution  = resolution,
                           n_runs      = n_runs,
                           seed        = seed)),
    class = "talos_pc_sweep"
  )
}


#' @export
print.talos_pc_sweep <- function(x, ...) {
  cat("\u2500\u2500 TALOS PC sweep \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  n_pcs tested   : %s\n", paste(x$params$n_pcs_range, collapse = ", ")))
  cat(sprintf("  Suggested n_pcs: %d  (elbow in graph modularity)\n", x$best_n_pcs))
  cat(sprintf("  k: %d  |  method: %s  |  resolution: %s\n\n",
              x$params$k, x$params$method, x$params$resolution))
  print(x$results, row.names = FALSE)
  invisible(x)
}


#' Sweep UMAP parameters
#'
#' UMAP has two parameters that most strongly shape the embedding.
#' \code{n_neighbors} governs how many local neighbours each cell considers:
#' small values emphasise fine local structure; large values preserve more
#' global topology at the cost of local detail.  \code{min_dist} controls how
#' tightly points are packed: small values produce compact, well-separated
#' clusters; large values spread points more evenly, which can reveal
#' continuous structure.  This function tests all combinations of both
#' parameters and scores each embedding on two complementary metrics:
#'
#' \describe{
#'   \item{KNN overlap}{Mean Jaccard similarity between each cell's \code{k}
#'     nearest neighbours in PCA space and in the UMAP embedding.  Measures
#'     how faithfully local neighbourhood structure from PCA is preserved in
#'     the 2D layout.  Computed once from a shared KNN reference in PCA space
#'     (\code{BiocNeighbors::findKNN}), which is reused across all embedding
#'     combinations for efficiency.  Ranges from 0 (no overlap) to 1
#'     (identical neighbourhoods).}
#'   \item{Mean silhouette width}{Mean silhouette width computed in 2D UMAP
#'     space using existing cluster labels from \code{cluster_col}.  Measures
#'     how visually separated the clusters are in the final embedding.  Ranges
#'     from \eqn{-1} (misassigned) to \eqn{+1} (well-separated).  Skipped
#'     when \code{cluster_col} is not present in \code{colData} — only KNN
#'     overlap is used for selection in that case.}
#' }
#'
#' Both metrics are normalised to \eqn{[0, 1]} and averaged into a composite
#' score.  The parameter combination with the highest composite score is
#' suggested as the optimum.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"PCA"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_pca}} first).
#' @param n_neighbors_range Integer vector of n_neighbors values to test.
#'   Default \code{c(10, 15, 20, 30, 50)}.
#' @param min_dist_range Numeric vector of min_dist values to test.
#'   Default \code{c(0.05, 0.1, 0.3, 0.5)}.
#' @param n_pcs Integer. Number of PCA components used as UMAP input and for
#'   KNN reference computation. Default \code{30}.
#' @param knn_k Integer. Number of nearest neighbours used for the KNN overlap
#'   metric. Default \code{15}.
#' @param cluster_col Character. \code{colData} column containing cluster
#'   labels used for silhouette computation. Default \code{"cluster"}.
#' @param seed Integer. Random seed for UMAP reproducibility. Default \code{42L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{talos_embedding_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per n_neighbors × min_dist
#'       combination: \code{n_neighbors}, \code{min_dist}, \code{knn_overlap},
#'       \code{mean_sil} (\code{NA} if cluster labels absent), \code{composite}.}
#'     \item{\code{embeddings}}{Named list of cells × 2 matrices — one UMAP
#'       embedding per parameter combination, in the same order as
#'       \code{results}.  Passed directly to
#'       \code{ASPIS_plot_embedding_grid()} for visual comparison.}
#'     \item{\code{plot}}{Tile heatmap of composite score (or KNN overlap if
#'       no clusters) across the n_neighbors × min_dist grid, with the best
#'       combination outlined in red.}
#'     \item{\code{best_params}}{Named list (\code{n_neighbors},
#'       \code{min_dist}) for the combination with the highest composite score.}
#'     \item{\code{params}}{List of sweep parameters for reproducibility.}
#'   }
#' @export
TALOS_tune_umap <- function(sce,
                              n_neighbors_range = c(10L, 15L, 20L, 30L, 50L),
                              min_dist_range    = c(0.05, 0.1, 0.3, 0.5),
                              use_rep           = "PCA",
                              n_pcs             = 30L,
                              knn_k             = 15L,
                              cluster_col       = "cluster",
                              use_graph         = TRUE,
                              seed              = 42L,
                              verbose           = TRUE) {

  run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
              paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
  .talos_check_dimred(sce, use_rep, run_fn)

  rep_mat   <- reducedDim(sce, use_rep)
  n_dims    <- if (use_rep == "PCA") min(as.integer(n_pcs), ncol(rep_mat)) else ncol(rep_mat)
  input_mat <- rep_mat[, seq_len(n_dims), drop = FALSE]
  knn_k     <- as.integer(knn_k)

  has_clusters <- cluster_col %in% names(colData(sce))
  if (!has_clusters && verbose)
    cat(sprintf(
      "  Note: '%s' not found in colData — silhouette will be skipped.\n",
      cluster_col))

  # ── Pre-compute KNN reference in representation space (once) ─────────────────
  if (verbose) cat(sprintf("  Computing KNN reference in %s space ...\n", use_rep))
  knn_ref <- BiocNeighbors::findKNN(input_mat, k = knn_k)$index

  # ── Graph-aligned mode ───────────────────────────────────────────────────────
  # n_neighbors is fixed by the SNN graph; only min_dist is swept.
  if (isTRUE(use_graph)) {
    if (!requireNamespace("uwot", quietly = TRUE))
      stop("Package 'uwot' is required for use_graph = TRUE. ",
           "Install via: install.packages('uwot')", call. = FALSE)
    g <- metadata(sce)$snn_graph
    if (is.null(g))
      stop("use_graph = TRUE requires metadata(sce)$snn_graph. ",
           "Run TALOS_build_graph() first.", call. = FALSE)

    adj   <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE)
    max_w <- max(adj@x)
    if (max_w > 1) adj@x <- adj@x / max_w

    # uwot::umap() treats a sparse X as a DISTANCE matrix, not a similarity
    # matrix — invert the normalized SNN affinity before passing it in.
    # See TALOS_run_umap() for the same fix and the empirical verification.
    adj@x <- 1 - adj@x

    k_graph  <- metadata(sce)$snn_k %||%
                 as.integer(round(mean(igraph::degree(g))))
    n_combos <- length(min_dist_range)

    if (verbose)
      cat(sprintf(
        "── TALOS: UMAP sweep (graph-aligned) %s\n  Input      : SNN graph  |  k = %s\n  min_dist   : %s\n  Combos     : %s  (n_neighbors fixed by graph)\n%s\n",
        strrep("─", 19),
        k_graph,
        paste(min_dist_range, collapse = ", "),
        n_combos,
        strrep("─", 56)))

    results    <- vector("list", n_combos)
    embeddings <- vector("list", n_combos)

    for (i in seq_along(min_dist_range)) {
      md <- min_dist_range[i]
      if (verbose)
        cat(sprintf("  (%d/%d) min_dist = %.2f ...\n", i, n_combos, md))

      set.seed(seed)
      umap_coords <- uwot::umap(
        X            = adj,
        n_components = 2L,
        n_neighbors  = k_graph,
        min_dist     = md,
        verbose      = FALSE
      )
      rownames(umap_coords) <- colnames(sce)

      embeddings[[i]] <- umap_coords
      knn_ov          <- .talos_knn_overlap(knn_ref, umap_coords, knn_k)

      mean_sil <- if (has_clusters) {
        labels <- as.integer(factor(colData(sce)[[cluster_col]]))
        if (length(unique(labels)) > 1L) {
          sil <- cluster::silhouette(labels, dist(umap_coords))
          mean(sil[, "sil_width"])
        } else NA_real_
      } else NA_real_

      results[[i]] <- data.frame(
        n_neighbors = k_graph,
        min_dist    = md,
        knn_overlap = round(knn_ov, 4),
        mean_sil    = if (has_clusters) round(mean_sil, 4) else NA_real_
      )
    }

    res_df <- do.call(rbind, results)
    res_df$composite <- {
      s1 <- .talos_safe_norm(res_df$knn_overlap)
      s2 <- if (has_clusters) .talos_safe_norm(res_df$mean_sil) else rep(0, nrow(res_df))
      if (has_clusters) (s1 + s2) / 2 else s1
    }

    best_i      <- which.max(res_df$composite)
    best_params <- list(n_neighbors = k_graph,
                        min_dist    = res_df$min_dist[best_i])

    p <- .talos_tune_umap_graph_plot(res_df, best_params, has_clusters)

    if (verbose)
      cat(sprintf(
        "%s\n  Best min_dist: %.2f  (composite = %.3f)\n%s\n",
        strrep("─", 56),
        best_params$min_dist, res_df$composite[best_i],
        strrep("─", 56)))

    return(structure(
      list(results     = res_df,
           embeddings  = embeddings,
           plot        = p,
           best_params = best_params,
           params      = list(type              = "umap",
                              use_rep           = use_rep,
                              n_dims            = n_dims,
                              n_neighbors_range = k_graph,
                              min_dist_range    = min_dist_range,
                              knn_k             = knn_k,
                              cluster_col       = cluster_col,
                              use_graph         = TRUE,
                              seed              = seed,
                              has_clusters      = has_clusters)),
      class = "talos_embedding_sweep"
    ))
  }

  # ── Default mode: independent KNN graph ──────────────────────────────────────
  n_neighbors_range <- as.integer(n_neighbors_range)
  n_combos <- length(n_neighbors_range) * length(min_dist_range)

  if (verbose)
    cat(sprintf(
      "── TALOS: UMAP sweep %s\n  n_neighbors: %s\n  min_dist   : %s\n  Input      : %s (%s dims)  |  knn_k: %s\n  Combos     : %s\n%s\n",
      strrep("─", 37),
      paste(n_neighbors_range, collapse = ", "),
      paste(min_dist_range,    collapse = ", "),
      use_rep, n_dims, knn_k, n_combos,
      strrep("─", 56)))

  results    <- vector("list", n_combos)
  embeddings <- vector("list", n_combos)
  run_i      <- 0L

  for (nn in n_neighbors_range) {
    for (md in min_dist_range) {
      run_i <- run_i + 1L
      if (verbose)
        cat(sprintf("  (%d/%d) n_neighbors = %d, min_dist = %.2f ...\n",
                        run_i, n_combos, nn, md))

      set.seed(seed)
      sce_tmp <- scater::runUMAP(sce,
                                  dimred      = use_rep,
                                  n_dimred    = n_dims,
                                  n_neighbors = nn,
                                  min_dist    = md,
                                  name        = "UMAP_sweep")

      emb_mat             <- reducedDim(sce_tmp, "UMAP_sweep")
      embeddings[[run_i]] <- emb_mat
      knn_ov              <- .talos_knn_overlap(knn_ref, emb_mat, knn_k)

      mean_sil <- if (has_clusters) {
        labels <- as.integer(factor(colData(sce_tmp)[[cluster_col]]))
        if (length(unique(labels)) > 1L) {
          sil <- cluster::silhouette(labels, dist(emb_mat))
          mean(sil[, "sil_width"])
        } else NA_real_
      } else NA_real_

      results[[run_i]] <- data.frame(
        n_neighbors = nn,
        min_dist    = md,
        knn_overlap = round(knn_ov, 4),
        mean_sil    = if (has_clusters) round(mean_sil, 4) else NA_real_
      )
    }
  }

  res_df <- do.call(rbind, results)
  res_df$composite <- {
    s1 <- .talos_safe_norm(res_df$knn_overlap)
    s2 <- if (has_clusters) .talos_safe_norm(res_df$mean_sil) else rep(0, nrow(res_df))
    if (has_clusters) (s1 + s2) / 2 else s1
  }

  best_i      <- which.max(res_df$composite)
  best_params <- list(n_neighbors = res_df$n_neighbors[best_i],
                      min_dist    = res_df$min_dist[best_i])

  p <- .talos_tune_umap_plot(res_df, best_params, has_clusters)

  if (verbose)
    cat(sprintf(
      "%s\n  Best combo : n_neighbors = %d, min_dist = %.2f  (composite = %.3f)\n%s\n",
      strrep("─", 56),
      best_params$n_neighbors, best_params$min_dist,
      res_df$composite[best_i],
      strrep("─", 56)))

  structure(
    list(results     = res_df,
         embeddings  = embeddings,
         plot        = p,
         best_params = best_params,
         params      = list(type              = "umap",
                            use_rep           = use_rep,
                            n_dims            = n_dims,
                            n_neighbors_range = n_neighbors_range,
                            min_dist_range    = min_dist_range,
                            knn_k             = knn_k,
                            cluster_col       = cluster_col,
                            use_graph         = FALSE,
                            seed              = seed,
                            has_clusters      = has_clusters)),
    class = "talos_embedding_sweep"
  )
}


#' Sweep tSNE perplexity
#'
#' Perplexity is tSNE's primary parameter: it loosely controls how many
#' neighbours each cell considers when constructing its local probability
#' distribution.  Low perplexity values emphasise very local structure and can
#' fragment continuous populations into disconnected islands; high values
#' integrate broader context but may merge distinct clusters.  Typical
#' effective ranges are 5–50 for most single-cell datasets.  This function
#' tests a range of perplexity values and scores each embedding on two
#' complementary metrics:
#'
#' \describe{
#'   \item{KNN overlap}{Mean Jaccard similarity between each cell's \code{k}
#'     nearest neighbours in PCA space and in the tSNE embedding.  Measures
#'     how faithfully local neighbourhood structure is preserved in the 2D
#'     layout.  A shared KNN reference in PCA space is computed once and
#'     reused across all perplexity values.  Ranges from 0 (no overlap) to
#'     1 (identical neighbourhoods).}
#'   \item{Mean silhouette width}{Mean silhouette width in 2D tSNE space
#'     using existing cluster labels from \code{cluster_col}.  Measures how
#'     visually well-separated the clusters are in the final embedding.  Ranges
#'     from \eqn{-1} (misassigned) to \eqn{+1} (well-separated).  Skipped
#'     when \code{cluster_col} is absent from \code{colData}.}
#' }
#'
#' Both metrics are normalised to \eqn{[0, 1]} and averaged into a composite
#' score.  The perplexity with the highest composite score is suggested.
#' Perplexity values exceeding \eqn{N / 3} (where \eqn{N} is the number of
#' cells) are invalid for tSNE and are silently removed before sweeping.
#'
#' This function returns the same \code{talos_embedding_sweep} class as
#' \code{\link{TALOS_tune_umap}}, so \code{print} and downstream helpers work
#' identically for both.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"PCA"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_pca}} first).
#' @param perplexity_range Numeric vector of perplexity values to test.
#'   Values exceeding \eqn{N/3} are removed automatically.
#'   Default \code{c(10, 20, 30, 50, 100)}.
#' @param n_pcs Integer. PCA components used for tSNE input and KNN reference.
#'   Default \code{30}.
#' @param max_iter Integer. tSNE iterations per run (fixed; not a tuning
#'   target). Default \code{1000}.
#' @param knn_k Integer. Nearest neighbours for the KNN overlap metric.
#'   Default \code{15}.
#' @param cluster_col Character. \code{colData} column for silhouette labels.
#'   Default \code{"cluster"}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{talos_embedding_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per perplexity value:
#'       \code{perplexity}, \code{knn_overlap}, \code{mean_sil} (\code{NA} if
#'       cluster labels absent), \code{composite}.}
#'     \item{\code{embeddings}}{List of cells × 2 matrices — one tSNE
#'       embedding per perplexity value, in the same order as \code{results}.
#'       Passed directly to \code{ASPIS_plot_embedding_grid()} for visual
#'       comparison.}
#'     \item{\code{plot}}{Multi-panel line plot (KNN overlap, mean silhouette
#'       if available, composite vs perplexity), with a red dashed vertical
#'       line at the suggested value.}
#'     \item{\code{best_params}}{Named list (\code{perplexity}) for the value
#'       with the highest composite score.}
#'     \item{\code{params}}{List of sweep parameters for reproducibility.}
#'   }
#' @export
TALOS_tune_tsne <- function(sce,
                              perplexity_range = c(10, 20, 30, 50, 100),
                              use_rep          = "PCA",
                              n_pcs            = 30L,
                              max_iter         = 1000L,
                              knn_k            = 15L,
                              cluster_col      = "cluster",
                              seed             = 42L,
                              verbose          = TRUE) {

  run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
              paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
  .talos_check_dimred(sce, use_rep, run_fn)

  rep_mat <- reducedDim(sce, use_rep)
  n_dims  <- if (use_rep == "PCA") min(as.integer(n_pcs), ncol(rep_mat)) else ncol(rep_mat)
  input_mat <- rep_mat[, seq_len(n_dims), drop = FALSE]
  n_cells <- ncol(sce)
  knn_k   <- as.integer(knn_k)

  # ── Remove invalid perplexity values ─────────────────────────────────────────
  max_perp         <- floor(n_cells / 3)
  perplexity_range <- perplexity_range[perplexity_range <= max_perp]
  if (length(perplexity_range) == 0L)
    stop("All perplexity values exceed N/3 = ", max_perp,
         ". Provide smaller perplexity values.", call. = FALSE)

  has_clusters <- cluster_col %in% names(colData(sce))
  if (!has_clusters && verbose)
    cat(sprintf(
      "  Note: '%s' not found in colData — silhouette will be skipped.\n",
      cluster_col))

  # ── Pre-compute KNN reference (once) ─────────────────────────────────────────
  if (verbose) cat(sprintf("  Computing KNN reference in %s space ...\n", use_rep))
  knn_ref <- BiocNeighbors::findKNN(input_mat, k = knn_k)$index

  n_combos <- length(perplexity_range)
  if (verbose)
    cat(sprintf(
      "\u2500\u2500 TALOS: tSNE sweep %s\n  Perplexity : %s\n  Input      : %s (%s dims)  |  max_iter: %s\n  knn_k      : %s\n%s\n",
      strrep("\u2500", 37),
      paste(perplexity_range, collapse = ", "),
      use_rep, n_dims, max_iter, knn_k,
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results    <- vector("list", n_combos)
  embeddings <- vector("list", n_combos)

  for (i in seq_along(perplexity_range)) {
    perp <- perplexity_range[i]
    if (verbose)
      cat(sprintf("  (%d/%d) perplexity = %g ...\n", i, n_combos, perp))

    set.seed(seed)
    sce_tmp <- scater::runTSNE(sce,
                                dimred     = use_rep,
                                n_dimred   = n_dims,
                                perplexity = perp,
                                max_iter   = as.integer(max_iter),
                                name       = "tSNE_sweep")

    emb_mat      <- reducedDim(sce_tmp, "tSNE_sweep")
    embeddings[[i]] <- emb_mat
    knn_ov       <- .talos_knn_overlap(knn_ref, emb_mat, knn_k)

    mean_sil <- if (has_clusters) {
      labels <- as.integer(factor(colData(sce_tmp)[[cluster_col]]))
      if (length(unique(labels)) > 1L) {
        sil <- cluster::silhouette(labels, dist(emb_mat))
        mean(sil[, "sil_width"])
      } else NA_real_
    } else NA_real_

    results[[i]] <- data.frame(
      perplexity  = perp,
      knn_overlap = round(knn_ov, 4),
      mean_sil    = if (has_clusters) round(mean_sil, 4) else NA_real_
    )
  }

  res_df <- do.call(rbind, results)

  # ── Composite score and best params ──────────────────────────────────────────
  res_df$composite <- {
    s1 <- .talos_safe_norm(res_df$knn_overlap)
    s2 <- if (has_clusters) .talos_safe_norm(res_df$mean_sil) else rep(0, nrow(res_df))
    if (has_clusters) (s1 + s2) / 2 else s1
  }

  best_i      <- which.max(res_df$composite)
  best_params <- list(perplexity = res_df$perplexity[best_i])

  p <- .talos_tune_tsne_plot(res_df, best_params)

  if (verbose)
    cat(sprintf(
      "%s\n  Best perplexity: %g  (composite = %.3f)\n%s\n",
      strrep("\u2500", 56),
      best_params$perplexity, res_df$composite[best_i],
      strrep("\u2500", 56)))

  structure(
    list(results     = res_df,
         embeddings  = embeddings,
         plot        = p,
         best_params = best_params,
         params      = list(type             = "tsne",
                            use_rep          = use_rep,
                            n_dims           = n_dims,
                            perplexity_range = perplexity_range,
                            max_iter         = as.integer(max_iter),
                            knn_k            = knn_k,
                            cluster_col      = cluster_col,
                            seed             = seed,
                            has_clusters     = has_clusters)),
    class = "talos_embedding_sweep"
  )
}


#' @export
print.talos_embedding_sweep <- function(x, ...) {
  type <- x$params$type %||% "embedding"
  cat(sprintf("\u2500\u2500 TALOS %s sweep \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n",
              toupper(type)))

  if (type == "umap") {
    cat(sprintf("  n_neighbors: %s\n",
                paste(x$params$n_neighbors_range, collapse = ", ")))
    cat(sprintf("  min_dist   : %s\n",
                paste(x$params$min_dist_range, collapse = ", ")))
    cat(sprintf("  Best       : n_neighbors = %d, min_dist = %.2f\n\n",
                x$best_params$n_neighbors, x$best_params$min_dist))
  } else {
    cat(sprintf("  Perplexity : %s\n",
                paste(x$params$perplexity_range, collapse = ", ")))
    cat(sprintf("  Best       : perplexity = %g\n\n",
                x$best_params$perplexity))
  }

  print(x$results, row.names = FALSE)
  invisible(x)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Kneedle elbow for TALOS_tune_pc: point of maximum perpendicular distance
# from the line connecting the first and last (n_pcs, modularity) point.
.talos_pc_elbow <- function(x, y) {
  n <- length(x)
  if (n <= 2L) return(x[which.max(y)])
  x_n <- (x - min(x)) / (max(x) - min(x))
  y_n <- (y - min(y)) / (max(y) - min(y) + .Machine$double.eps)
  dx  <- x_n[n] - x_n[1L]
  dy  <- y_n[n] - y_n[1L]
  dists <- abs(dy * x_n - dx * y_n + x_n[n] * y_n[1L] - y_n[n] * x_n[1L]) /
           sqrt(dx^2 + dy^2)
  x[which.max(dists)]
}


# Two-panel line plot for TALOS_tune_pc (cum_var + modularity vs n_pcs).
.talos_tune_pc_plot <- function(res_df, best_n_pcs) {
  metrics <- c("cum_var", "modularity")
  lab_map <- c(cum_var    = "Cumulative variance\nexplained (%)",
               modularity = "Graph modularity")

  df_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(n_pcs  = res_df$n_pcs,
               metric = lab_map[m],
               value  = res_df[[m]])
  }))
  df_long$metric <- factor(df_long$metric, levels = lab_map[metrics])

  ggplot(df_long, aes(x = n_pcs, y = value)) +
    geom_line(colour = "#4E79A7", linewidth = 0.8) +
    geom_point(colour = "#4E79A7", size = 2.5) +
    geom_vline(xintercept = best_n_pcs,
               linetype = "dashed", colour = "#E15759", linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(x       = "Number of PCs",
         y       = NULL,
         title   = "PCA component sweep",
         caption = sprintf("Red dashed: suggested n_pcs = %d", best_n_pcs)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


# Mean Jaccard overlap between PCA KNN and embedding KNN.
# knn_ref_idx: integer matrix (cells × k) from BiocNeighbors::findKNN
# emb_mat    : numeric matrix (cells × 2)
# k          : integer, number of neighbours to use in embedding
.talos_knn_overlap <- function(knn_ref_idx, emb_mat, k) {
  knn_emb <- BiocNeighbors::findKNN(emb_mat, k = k)$index  # cells × k
  n_cells <- nrow(knn_ref_idx)

  overlaps <- vapply(seq_len(n_cells), function(i) {
    ref_set <- knn_ref_idx[i, ]
    emb_set <- knn_emb[i, ]
    length(intersect(ref_set, emb_set)) / length(union(ref_set, emb_set))
  }, numeric(1L))

  mean(overlaps, na.rm = TRUE)
}


# Normalise a numeric vector to [0, 1].  Returns all-0 if range is 0.
.talos_safe_norm <- function(x) {
  x[is.na(x)] <- 0
  rng <- range(x, na.rm = TRUE)
  if (rng[2] - rng[1] < .Machine$double.eps) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}


# Line plot for UMAP sweep in graph-aligned mode (min_dist only).
.talos_tune_umap_graph_plot <- function(res_df, best_params, has_sil) {
  metrics    <- if (has_sil) c("knn_overlap", "mean_sil", "composite") else c("knn_overlap")
  labels_map <- c(knn_overlap = "KNN overlap",
                  mean_sil    = "Mean silhouette (2D)",
                  composite   = "Composite score")

  df_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(min_dist = res_df$min_dist,
               metric   = labels_map[m],
               value    = res_df[[m]])
  }))
  df_long$metric <- factor(df_long$metric, levels = labels_map[metrics])

  ggplot(df_long, aes(x = min_dist, y = value)) +
    geom_line(colour = "#4E79A7", linewidth = 0.8) +
    geom_point(colour = "#4E79A7", size = 2.5) +
    geom_vline(xintercept = best_params$min_dist,
               linetype = "dashed", colour = "#E15759", linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(x       = "min_dist",
         y       = NULL,
         title   = "UMAP min_dist sweep (graph-aligned)",
         caption = sprintf("Red dashed: suggested min_dist = %.2f  |  n_neighbors fixed at %d",
                           best_params$min_dist, best_params$n_neighbors)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


# Tile heatmap for UMAP sweep (n_neighbors × min_dist).
.talos_tune_umap_plot <- function(res_df, best_params, has_sil) {

  score_col <- if (has_sil) "composite" else "knn_overlap"
  score_lab <- if (has_sil) "Composite score\n(KNN overlap + silhouette)"
               else "KNN overlap"

  res_df$n_neighbors <- factor(res_df$n_neighbors)
  res_df$min_dist    <- factor(res_df$min_dist)
  res_df$raw_val     <- res_df[[score_col]]

  best_nn <- factor(best_params$n_neighbors, levels = levels(res_df$n_neighbors))
  best_md <- factor(best_params$min_dist,    levels = levels(res_df$min_dist))

  ggplot(res_df, aes(x = min_dist, y = n_neighbors, fill = raw_val)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = res_df[res_df$n_neighbors == best_nn &
                              res_df$min_dist    == best_md, ],
              aes(x = min_dist, y = n_neighbors),
              fill  = NA, colour = "#E15759", linewidth = 1.5) +
    scale_fill_gradient(low = "#EFF3FF", high = "#2171B5",
                        name = score_lab, na.value = "grey90") +
    labs(x     = "min_dist",
         y     = "n_neighbors",
         title = "UMAP parameter sweep") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
}


# Line plot for tSNE sweep (perplexity).
.talos_tune_tsne_plot <- function(res_df, best_params) {

  metrics    <- c("knn_overlap", "composite")
  labels_map <- c(knn_overlap = "KNN overlap",
                  mean_sil    = "Mean silhouette (2D)",
                  composite   = "Composite score")
  if (!all(is.na(res_df$mean_sil))) metrics <- c("knn_overlap", "mean_sil", "composite")

  df_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(perplexity  = res_df$perplexity,
               metric      = labels_map[m],
               value       = res_df[[m]])
  }))
  df_long$metric <- factor(df_long$metric, levels = labels_map[metrics])

  ggplot(df_long, aes(x = perplexity, y = value)) +
    geom_line(colour = "#4E79A7", linewidth = 0.8) +
    geom_point(colour = "#4E79A7", size = 2.5) +
    geom_vline(xintercept = best_params$perplexity,
               linetype = "dashed", colour = "#E15759", linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(x       = "Perplexity",
         y       = NULL,
         title   = "tSNE perplexity sweep",
         caption = sprintf("Red dashed: suggested perplexity = %g",
                           best_params$perplexity)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}
