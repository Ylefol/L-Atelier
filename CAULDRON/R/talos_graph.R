# ==============================================================================
# TALOS - Graph Construction, Clustering & Graph Parameter Tuning
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# SNN graph construction, community detection, and sweeps over the graph
# parameters k and clustering resolution.
#
# All functions prefixed: TALOS_
# ==============================================================================


#' Build a shared nearest-neighbour graph
#'
#' Constructs a SNN graph from a reduced-dimension representation using
#' \code{scran::buildSNNGraph}.  Defaults to PCA; pass
#' \code{use_rep = "scVI"} to build the graph in the scVI latent space
#' produced by \code{\link{TALOS_run_scvi}}.  The graph is stored in
#' \code{metadata(sce)$snn_graph} for use by \code{\link{TALOS_cluster}}.
#'
#' @param sce A \code{SingleCellExperiment} with the chosen reduced dimension
#'   already populated.
#' @param use_rep Character. Name of the \code{reducedDims} slot to use.
#'   Default \code{"PCA"}.  Use \code{"scVI"} for the batch-corrected latent
#'   embedding from \code{\link{TALOS_run_scvi}}.
#' @param n_pcs Integer. Number of PCA components to use when
#'   \code{use_rep = "PCA"}.  Ignored for other representations (all dims are
#'   used).  Default \code{30}.
#' @param k Integer. Number of nearest neighbours. Higher values produce
#'   broader, more robust clusters. Default \code{20}.
#' @param type Character. Edge weight scheme: \code{"rank"} (default),
#'   \code{"number"}, or \code{"jaccard"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{metadata(sce)$snn_graph} populated.
#' @export
TALOS_build_graph <- function(sce,
                               use_rep = "PCA",
                               n_pcs   = 30L,
                               k       = 20L,
                               type    = c("rank", "number", "jaccard"),
                               verbose = TRUE) {

  type   <- match.arg(type)
  run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
              paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
  .talos_check_dimred(sce, use_rep, run_fn)

  rep_mat <- reducedDim(sce, use_rep)
  n_dims  <- if (use_rep == "PCA") {
    min(as.integer(n_pcs), ncol(rep_mat))
  } else {
    ncol(rep_mat)
  }
  input_mat <- rep_mat[, seq_len(n_dims), drop = FALSE]

  g <- scran::buildSNNGraph(t(input_mat), k = as.integer(k), type = type)

  metadata(sce)$snn_graph <- g

  if (isTRUE(verbose))
    cat(sprintf(
      "\u2500\u2500 TALOS: SNN graph %s\n  Input      : %s (%s dims)\n  k          : %s  |  type: %s\n  Vertices   : %s  |  Edges: %s\n  Stored in  : metadata(sce)$snn_graph\n%s\n",
      strrep("\u2500", 38),
      use_rep, n_dims, k, type,
      format(igraph::vcount(g), big.mark = ","),
      format(igraph::ecount(g), big.mark = ","),
      strrep("\u2500", 56)
    ))

  sce
}


#' Graph-based clustering
#'
#' Applies Leiden or Louvain community detection to the SNN graph built by
#' \code{\link{TALOS_build_graph}}.  Cluster labels are stored in
#' \code{colData(sce)[[cluster_col]]}.
#'
#' @param sce A \code{SingleCellExperiment} with \code{metadata(sce)$snn_graph}
#'   populated (i.e. \code{\link{TALOS_build_graph}} must have been run first).
#' @param method Character. \code{"leiden"} (default) or \code{"louvain"}.
#'   Leiden generally produces better-connected, more reproducible partitions.
#' @param resolution Numeric. Resolution parameter controlling cluster
#'   granularity.  Higher values yield more, smaller clusters. Default
#'   \code{1.0}. Applied to both methods.
#' @param cluster_col Character. \code{colData} column name for storing labels.
#'   Default \code{"cluster"}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with cluster labels in \code{colData(sce)[[cluster_col]]}.
#' @export
TALOS_cluster <- function(sce,
                           method      = c("leiden", "louvain"),
                           resolution  = 1.0,
                           cluster_col = "cluster",
                           seed        = 42L,
                           verbose     = TRUE) {

  method <- match.arg(method)

  g <- metadata(sce)$snn_graph
  if (is.null(g))
    stop("SNN graph not found. Run TALOS_build_graph() first.", call. = FALSE)

  set.seed(seed)
  communities <- if (method == "leiden") {
    igraph::cluster_leiden(g, resolution_parameter = resolution)
  } else {
    igraph::cluster_louvain(g, resolution = resolution)
  }

  colData(sce)[[cluster_col]] <- factor(igraph::membership(communities))

  if (isTRUE(verbose)) {
    cl       <- colData(sce)[[cluster_col]]
    n_cl     <- nlevels(cl)
    cl_sizes <- sort(table(cl), decreasing = TRUE)
    cat(sprintf(
      "\u2500\u2500 TALOS: Clustering (%s, resolution = %s) %s\n  Clusters   : %s\n  Sizes      : %s \u2013 %s (median %s)\n  Stored in  : colData(sce)[[\"%s\"]]\n%s\n",
      method, resolution, strrep("\u2500", 18),
      n_cl,
      format(min(cl_sizes), big.mark = ","),
      format(max(cl_sizes), big.mark = ","),
      format(median(as.integer(cl_sizes)), big.mark = ","),
      cluster_col,
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' Sweep k values for SNN graph construction
#'
#' The \code{k} parameter in \code{\link{TALOS_build_graph}} controls how many
#' nearest neighbours each cell shares edges with.  Small \code{k} values
#' produce fine-grained, loosely connected graphs that tend to over-fragment;
#' large values produce coarser graphs that may merge distinct populations.
#' This function sweeps a range of \code{k} values, builds an SNN graph and
#' applies community detection at each, and scores the resulting partitions on
#' two complementary metrics:
#'
#' \describe{
#'   \item{Modularity}{Graph modularity of the partition — measures how cleanly
#'     communities separate relative to a random null graph.  Higher values
#'     indicate more distinct cluster boundaries.  Always computed.}
#'   \item{Mean silhouette width}{Average silhouette width in the representation
#'     space — for each cell, how much more similar it is to its own cluster
#'     than to the nearest other cluster.  Ranges from \eqn{-1} (misassigned)
#'     to \eqn{+1} (well-separated).  The \code{k} with the highest mean
#'     silhouette is suggested.  Automatically skipped when
#'     \code{ncol(sce) > max_cells_sil} (pairwise distances are \eqn{O(N^2)});
#'     modularity is used instead.}
#' }
#'
#' Resolution and clustering method are held fixed during the sweep — they are
#' not targets here.  Use \code{\link{TALOS_tune_resolution}} afterwards to
#' optimise those parameters with the chosen \code{k}.
#'
#' @param sce A \code{SingleCellExperiment} with the chosen reduced dimension
#'   already populated.
#' @param k_range Integer vector of k values to test.
#'   Default \code{seq(5, 50, by = 5)}.
#' @param use_rep Character. Name of the \code{reducedDims} slot to use for
#'   graph construction and silhouette distances.  Default \code{"PCA"}.  Use
#'   \code{"scVI"} for the batch-corrected latent embedding from
#'   \code{\link{TALOS_run_scvi}}.
#' @param n_pcs Integer. Number of PCA components to use when
#'   \code{use_rep = "PCA"}.  Ignored for other representations (all dims are
#'   used).  Default \code{30}.
#' @param type Character. Edge weight scheme for \code{scran::buildSNNGraph}:
#'   \code{"rank"} (default), \code{"number"}, or \code{"jaccard"}.
#' @param method Character. Clustering algorithm: \code{"leiden"} (default)
#'   or \code{"louvain"}.
#' @param resolution Numeric. Clustering resolution (held fixed during sweep).
#'   Default \code{1.0}.
#' @param compute_silhouette Logical. Compute mean silhouette width.
#'   Automatically disabled when \code{ncol(sce) > max_cells_sil}.
#'   Default \code{TRUE}.
#' @param max_cells_sil Integer. Cell count above which silhouette is skipped.
#'   Default \code{10000}.
#' @param n_runs Integer. Number of community-detection runs per k value.
#'   Results are averaged over runs (modularity mean; representative membership
#'   closest to the mean used for silhouette).  Increasing \code{n_runs}
#'   dampens Leiden/Louvain stochasticity at the cost of proportionally longer
#'   runtime.  Default \code{5L}.
#' @param seed Integer. Random seed for the first run; subsequent runs use
#'   \code{seed + 1}, \code{seed + 2}, … Default \code{42L}.
#' @param verbose Logical. Print per-k progress and a final summary.
#'   Default \code{TRUE}.
#'
#' @return A \code{talos_k_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per k: \code{k},
#'       \code{n_clusters}, \code{modularity} (mean over \code{n_runs}),
#'       \code{mean_sil} (\code{NA} if silhouette was skipped).}
#'     \item{\code{plot}}{Multi-panel line plot (one panel per metric) with a
#'       red dashed vertical line marking the suggested \code{k}.}
#'     \item{\code{best_k}}{Suggested k — the value with the highest mean
#'       silhouette (or highest modularity when silhouette is skipped).}
#'     \item{\code{params}}{List of sweep parameters for reproducibility.}
#'   }
#' @export
TALOS_tune_k <- function(sce,
                          k_range            = seq(5L, 50L, by = 5L),
                          use_rep            = "PCA",
                          n_pcs              = 30L,
                          type               = c("rank", "number", "jaccard"),
                          method             = c("leiden", "louvain"),
                          resolution         = 1.0,
                          compute_silhouette = TRUE,
                          max_cells_sil      = 10000L,
                          n_runs             = 5L,
                          seed               = 42L,
                          verbose            = TRUE) {

  type   <- match.arg(type)
  method <- match.arg(method)
  run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
              paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
  .talos_check_dimred(sce, use_rep, run_fn)

  n_cells   <- ncol(sce)
  rep_mat   <- reducedDim(sce, use_rep)
  n_dims    <- if (use_rep == "PCA") min(as.integer(n_pcs), ncol(rep_mat)) else ncol(rep_mat)
  input_mat <- rep_mat[, seq_len(n_dims), drop = FALSE]
  k_range   <- sort(unique(as.integer(k_range)))

  # ── Silhouette feasibility ───────────────────────────────────────────────────
  do_sil <- isTRUE(compute_silhouette) && n_cells <= as.integer(max_cells_sil)
  if (isTRUE(compute_silhouette) && !do_sil)
    cat(sprintf(
      "Silhouette skipped: %s cells exceeds max_cells_sil (%s). Using modularity.\n",
      format(n_cells, big.mark = ","),
      format(as.integer(max_cells_sil), big.mark = ",")))

  # ── Pre-compute PCA submatrix and distance matrix (once, outside loop) ───────
  dist_mat <- if (do_sil) dist(input_mat) else NULL

  n_runs <- max(1L, as.integer(n_runs))

  if (verbose)
    cat(sprintf(
      "\u2500\u2500 TALOS: k sweep %s\n  k values   : %s\n  Input      : %s (%s dims)  |  type: %s\n  method     : %s  |  resolution: %s\n  n_runs     : %s\n  Silhouette : %s\n%s\n",
      strrep("\u2500", 40),
      paste(k_range, collapse = ", "),
      use_rep, n_dims, type, method, resolution, n_runs,
      if (do_sil) "yes" else "no (modularity used for selection)",
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results <- vector("list", length(k_range))

  for (i in seq_along(k_range)) {
    ki <- k_range[i]
    if (verbose) cat(sprintf("  Testing k = %d ...\n", ki))

    g    <- scran::buildSNNGraph(t(input_mat), k = ki, type = type)
    run  <- .talos_multi_run_communities(g, method, resolution, seed, n_runs)

    mean_sil <- if (do_sil) {
      sil <- cluster::silhouette(run$labels, dist_mat)
      mean(sil[, "sil_width"])
    } else NA_real_

    results[[i]] <- data.frame(
      k          = ki,
      n_clusters = run$n_clusters,
      modularity = round(run$modularity, 4),
      mean_sil   = if (do_sil) round(mean_sil, 4) else NA_real_
    )
  }

  res_df <- do.call(rbind, results)

  best_k <- if (do_sil) res_df$k[which.max(res_df$mean_sil)] else
              res_df$k[which.max(res_df$modularity)]

  p <- .talos_tune_k_plot(res_df, best_k, do_sil)

  if (verbose)
    cat(sprintf(
      "%s\n  Suggested k: %d (%s)\n%s\n",
      strrep("\u2500", 56),
      best_k,
      if (do_sil) "highest mean silhouette width" else "highest modularity",
      strrep("\u2500", 56)))

  structure(
    list(results = res_df,
         plot    = p,
         best_k  = best_k,
         params  = list(k_range    = k_range,
                        use_rep    = use_rep,
                        n_dims     = n_dims,
                        type       = type,
                        method     = method,
                        resolution = resolution,
                        n_runs     = n_runs,
                        silhouette = do_sil)),
    class = "talos_k_sweep"
  )
}


#' @export
print.talos_k_sweep <- function(x, ...) {
  cat("\u2500\u2500 TALOS k sweep \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  k tested   : %s\n", paste(x$params$k_range, collapse = ", ")))
  cat(sprintf("  Suggested k: %d\n", x$best_k))
  cat(sprintf("  Method     : %s  |  resolution: %s\n\n",
              x$params$method, x$params$resolution))
  print(x$results, row.names = FALSE)
  invisible(x)
}


#' Sweep clustering resolution
#'
#' Resolution controls the granularity of Leiden and Louvain community
#' detection: higher values push the algorithm to accept smaller, more numerous
#' clusters; lower values favour larger, coarser partitions.  This function
#' sweeps a range of resolution values on the pre-built SNN graph and scores
#' each partition on four complementary metrics:
#'
#' \describe{
#'   \item{n_clusters}{Number of distinct clusters produced.  Rises
#'     monotonically with resolution; useful for anchoring biological
#'     expectations (e.g. known number of major cell types).}
#'   \item{modularity}{Graph modularity of the partition — how cleanly
#'     community boundaries separate relative to a random null.  Higher is
#'     better, but very high resolution can inflate modularity by fragmenting
#'     real populations.}
#'   \item{mean_sil}{Mean silhouette width in PCA space — how confidently each
#'     cell is assigned to its cluster vs the nearest alternative.  Ranges from
#'     \eqn{-1} (misassigned) to \eqn{+1} (well-separated).  The resolution
#'     per method with the highest mean silhouette is suggested.  Skipped when
#'     \code{ncol(sce) > max_cells_sil}; modularity is used instead.}
#'   \item{min_cl_size}{Size of the smallest cluster.  Near-zero values flag
#'     pathological resolutions that fragment cells into singletons or
#'     micro-clusters, typically a sign the resolution is too high.}
#' }
#'
#' When \code{method = "both"}, Leiden and Louvain are run at every resolution
#' for direct comparison.  The PCA distance matrix for silhouette is computed
#' once before the loop, so runtime scales with the number of resolution values,
#' not the number of cells per iteration.
#'
#' @param sce A \code{SingleCellExperiment} with \code{metadata(sce)$snn_graph}
#'   populated (run \code{\link{TALOS_build_graph}} first).  PCA must also be
#'   present for silhouette computation.
#' @param resolution_range Numeric vector of resolution values to test.
#'   Default \code{seq(0.1, 2.0, by = 0.1)}.
#' @param method Character.  \code{"leiden"}, \code{"louvain"}, or
#'   \code{"both"} (default).
#' @param n_pcs Integer.  PCA components used for the silhouette distance
#'   matrix.  Default \code{30}.
#' @param compute_silhouette Logical.  Compute mean silhouette.  Automatically
#'   disabled when \code{ncol(sce) > max_cells_sil}.  Default \code{TRUE}.
#' @param max_cells_sil Integer.  Cell count above which silhouette is skipped.
#'   Default \code{10000}.
#' @param n_runs Integer. Number of community-detection runs per
#'   method × resolution combination.  Modularity is averaged over all runs;
#'   the representative membership (closest to mean modularity) is used for
#'   silhouette and cluster-size metrics.  Default \code{5L}.
#' @param seed Integer.  Random seed for the first run; subsequent runs use
#'   \code{seed + 1}, \code{seed + 2}, … Default \code{42L}.
#' @param verbose Logical.  Print per-combination progress.  Default \code{TRUE}.
#'
#' @return A \code{talos_resolution_sweep} list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per method × resolution
#'       combination: \code{method}, \code{resolution}, \code{n_clusters},
#'       \code{modularity} (mean over \code{n_runs}), \code{mean_sil}
#'       (\code{NA} if skipped), \code{min_cl_size}.}
#'     \item{\code{plot}}{Four-panel line plot (n_clusters, modularity,
#'       mean_sil, min_cl_size vs resolution), one coloured line per method,
#'       with dashed vertical lines at the best resolution per method.}
#'     \item{\code{best_params}}{Named list (\code{method}, \code{resolution})
#'       for the single combination with the highest score across all methods.}
#'     \item{\code{best_per_method}}{Data frame with the best-scoring row per
#'       clustering algorithm — useful when \code{method = "both"} to compare
#'       Leiden and Louvain optima side by side.}
#'     \item{\code{params}}{List of sweep parameters for reproducibility.}
#'   }
#' @export
TALOS_tune_resolution <- function(sce,
                                   resolution_range   = seq(0.1, 2.0, by = 0.1),
                                   method             = c("both", "leiden", "louvain"),
                                   use_rep            = "PCA",
                                   n_pcs              = 30L,
                                   compute_silhouette = TRUE,
                                   max_cells_sil      = 10000L,
                                   n_runs             = 5L,
                                   seed               = 42L,
                                   verbose            = TRUE) {

  method         <- match.arg(method)
  methods_to_run <- if (method == "both") c("leiden", "louvain") else method

  g <- metadata(sce)$snn_graph
  if (is.null(g))
    stop("SNN graph not found. Run TALOS_build_graph() first.", call. = FALSE)

  resolution_range <- as.numeric(resolution_range)
  n_combos         <- length(resolution_range) * length(methods_to_run)
  n_runs           <- max(1L, as.integer(n_runs))

  # ── Silhouette setup ─────────────────────────────────────────────────────────
  n_cells <- ncol(sce)
  do_sil  <- isTRUE(compute_silhouette) && n_cells <= as.integer(max_cells_sil)

  if (isTRUE(compute_silhouette) && !do_sil)
    cat(sprintf(
      "Silhouette skipped: %s cells exceeds max_cells_sil (%s). Using modularity.\n",
      format(n_cells, big.mark = ","),
      format(as.integer(max_cells_sil), big.mark = ",")))

  dist_mat <- NULL
  if (do_sil) {
    run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
                paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
    .talos_check_dimred(sce, use_rep, run_fn)
    rep_mat  <- reducedDim(sce, use_rep)
    n_dims   <- if (use_rep == "PCA") min(as.integer(n_pcs), ncol(rep_mat)) else ncol(rep_mat)
    sil_mat  <- rep_mat[, seq_len(n_dims), drop = FALSE]
    if (verbose)
      cat(sprintf("  Pre-computing distance matrix in %s space ...\n", use_rep))
    dist_mat <- dist(sil_mat)
  }

  if (verbose)
    cat(sprintf(
      "\u2500\u2500 TALOS: resolution sweep %s\n  Resolutions: %s\n  Method(s)  : %s\n  n_runs     : %s\n  Silhouette : %s\n%s\n",
      strrep("\u2500", 31),
      paste(round(resolution_range, 2), collapse = ", "),
      paste(methods_to_run, collapse = " + "),
      n_runs,
      if (do_sil) "yes" else "no",
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results <- vector("list", n_combos)
  run_i   <- 0L

  for (m in methods_to_run) {
    for (res in resolution_range) {
      run_i <- run_i + 1L
      if (verbose)
        cat(sprintf("  (%d/%d) method = %s, resolution = %.2f ...\n",
                        run_i, n_combos, m, res))

      run      <- .talos_multi_run_communities(g, m, res, seed, n_runs)
      min_size <- min(as.integer(table(run$labels)))

      mean_sil <- if (do_sil && run$n_clusters > 1L) {
        sil <- cluster::silhouette(run$labels, dist_mat)
        mean(sil[, "sil_width"])
      } else NA_real_

      results[[run_i]] <- data.frame(
        method      = m,
        resolution  = res,
        n_clusters  = run$n_clusters,
        modularity  = round(run$modularity, 4),
        mean_sil    = if (do_sil) round(mean_sil, 4) else NA_real_,
        min_cl_size = min_size
      )
    }
  }

  res_df <- do.call(rbind, results)

  # ── Best per method ──────────────────────────────────────────────────────────
  criterion <- if (do_sil) "mean_sil" else "modularity"

  best_per_method <- do.call(rbind, lapply(methods_to_run, function(m) {
    sub      <- res_df[res_df$method == m, ]
    best_row <- sub[which.max(sub[[criterion]]), ]
    best_row$criterion <- criterion
    best_row
  }))

  best_i      <- which.max(res_df[[criterion]])
  best_params <- list(method     = res_df$method[best_i],
                      resolution = res_df$resolution[best_i])

  p <- .talos_tune_resolution_plot(res_df, best_per_method, do_sil,
                                    methods_to_run)

  if (verbose)
    cat(sprintf(
      "%s\n  Best overall: %s, resolution = %.2f (%s)\n%s\n",
      strrep("\u2500", 56),
      best_params$method, best_params$resolution, criterion,
      strrep("\u2500", 56)))

  structure(
    list(results         = res_df,
         plot            = p,
         best_params     = best_params,
         best_per_method = best_per_method,
         params          = list(resolution_range   = resolution_range,
                                method             = method,
                                n_pcs              = n_pcs,
                                n_runs             = n_runs,
                                seed               = seed,
                                silhouette         = do_sil)),
    class = "talos_resolution_sweep"
  )
}


#' @export
print.talos_resolution_sweep <- function(x, ...) {
  cat("\u2500\u2500 TALOS resolution sweep \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  Best overall: %s, resolution = %.2f\n",
              x$best_params$method, x$best_params$resolution))
  cat(sprintf("  Silhouette  : %s\n\n  Best per method:\n",
              if (x$params$silhouette) "yes" else "no"))
  for (i in seq_len(nrow(x$best_per_method))) {
    cat(sprintf("    %-8s resolution = %.2f  n_clusters = %d  (%s)\n",
                x$best_per_method$method[i],
                x$best_per_method$resolution[i],
                x$best_per_method$n_clusters[i],
                x$best_per_method$criterion[i]))
  }
  cat("\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Run community detection n_runs times (incrementing the seed) and return the
# run whose modularity is closest to the mean — a simple way to dampen the
# stochasticity of Leiden/Louvain without committing to a single seed.
#
# Returns a list with:
#   modularity : mean modularity across runs
#   n_clusters : n_clusters of the representative run
#   labels     : membership of the representative run
.talos_multi_run_communities <- function(g, method, resolution, seed, n_runs) {
  runs <- lapply(seq_len(n_runs), function(r) {
    set.seed(seed + r - 1L)
    communities <- if (method == "leiden") {
      igraph::cluster_leiden(g, resolution_parameter = resolution)
    } else {
      igraph::cluster_louvain(g, resolution = resolution)
    }
    labels <- as.integer(igraph::membership(communities))
    list(modularity = igraph::modularity(g, labels),
         n_clusters = length(unique(labels)),
         labels     = labels)
  })

  modularities <- vapply(runs, `[[`, numeric(1L), "modularity")
  mean_modul   <- mean(modularities)
  rep_run      <- runs[[which.min(abs(modularities - mean_modul))]]

  list(modularity = mean_modul,
       n_clusters = rep_run$n_clusters,
       labels     = rep_run$labels)
}


.talos_tune_k_plot <- function(res_df, best_k, do_sil) {

  metrics    <- c("n_clusters", "modularity")
  labels_map <- c(n_clusters = "Number of clusters",
                  modularity = "Modularity",
                  mean_sil   = "Mean silhouette width")
  if (do_sil) metrics <- c(metrics, "mean_sil")

  df_long <- data.frame(
    k            = rep(res_df$k, length(metrics)),
    metric       = rep(metrics, each = nrow(res_df)),
    value        = unlist(res_df[, metrics, drop = FALSE], use.names = FALSE)
  )
  df_long$metric_label <- factor(labels_map[df_long$metric],
                                  levels = labels_map[metrics])

  ggplot(df_long, aes(x = k, y = value)) +
    geom_line(colour = "#4E79A7", linewidth = 0.8) +
    geom_point(colour = "#4E79A7", size = 2.5) +
    geom_vline(xintercept = best_k, linetype = "dashed",
               colour = "#E15759", linewidth = 0.7) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
    labs(x       = "k (nearest neighbours)",
         y       = NULL,
         title   = "SNN graph k sweep",
         caption = sprintf("Red dashed: suggested k = %d", best_k)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


.talos_tune_resolution_plot <- function(res_df, best_per_method, do_sil,
                                         methods_to_run) {

  metrics <- c("n_clusters", "modularity")
  lab_map <- c(n_clusters  = "Number of clusters",
               modularity  = "Modularity",
               mean_sil    = "Mean silhouette (PCA)",
               min_cl_size = "Min. cluster size")
  if (do_sil) metrics <- c(metrics, "mean_sil")
  metrics <- c(metrics, "min_cl_size")

  df_long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(method     = res_df$method,
               resolution = res_df$resolution,
               metric     = lab_map[m],
               value      = res_df[[m]])
  }))
  df_long$metric <- factor(df_long$metric, levels = lab_map[metrics])

  method_colours <- c(leiden = "#4E79A7", louvain = "#F28E2B")
  colours_used   <- method_colours[methods_to_run]

  best_lines <- best_per_method[, c("method", "resolution")]

  cap <- paste(
    vapply(seq_len(nrow(best_per_method)), function(i) {
      sprintf("%s: res = %.2f  (n_cl = %d)",
              best_per_method$method[i],
              best_per_method$resolution[i],
              best_per_method$n_clusters[i])
    }, character(1L)),
    collapse = "   |   ")

  ggplot(df_long, aes(x = resolution, y = value, colour = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(data     = best_lines,
               aes(xintercept = resolution, colour = method),
               linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_color_manual(values = colours_used, name = "Method") +
    labs(x       = "Resolution",
         y       = NULL,
         title   = "Clustering resolution sweep",
         caption = paste("Dashed lines \u2014 best:", cap)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}
