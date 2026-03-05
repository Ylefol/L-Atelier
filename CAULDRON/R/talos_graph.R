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
#' Constructs a SNN graph from the PCA embedding using
#' \code{scran::buildSNNGraph}.  The graph is stored in
#' \code{metadata(sce)$snn_graph} for use by \code{\link{TALOS_cluster}}.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"PCA"} reduced
#'   dimension.
#' @param n_pcs Integer. Number of PCA components to use. Default \code{30}.
#' @param k Integer. Number of nearest neighbours. Higher values produce
#'   broader, more robust clusters. Default \code{20}.
#' @param type Character. Edge weight scheme: \code{"rank"} (default),
#'   \code{"number"}, or \code{"jaccard"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{metadata(sce)$snn_graph} populated.
#' @export
TALOS_build_graph <- function(sce,
                               n_pcs   = 30L,
                               k       = 20L,
                               type    = c("rank", "number", "jaccard"),
                               verbose = TRUE) {

  type <- match.arg(type)
  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_pcs   <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
  pca_mat <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]

  g <- scran::buildSNNGraph(t(pca_mat), k = as.integer(k), type = type)

  metadata(sce)$snn_graph <- g

  if (isTRUE(verbose))
    message(sprintf(
      "\u2500\u2500 TALOS: SNN graph %s\n  Input      : PCA (%s components)\n  k          : %s  |  type: %s\n  Vertices   : %s  |  Edges: %s\n  Stored in  : metadata(sce)$snn_graph\n%s",
      strrep("\u2500", 38),
      n_pcs, k, type,
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
    message(sprintf(
      "\u2500\u2500 TALOS: Clustering (%s, resolution = %s) %s\n  Clusters   : %s\n  Sizes      : %s \u2013 %s (median %s)\n  Stored in  : colData(sce)[[\"%s\"]]\n%s",
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
#' Builds an SNN graph and clusters cells across a range of \code{k} values,
#' returning cluster count, modularity, and (for smaller datasets) mean
#' silhouette width per k.  Use this to choose an appropriate \code{k} before
#' running \code{\link{TALOS_build_graph}} and \code{\link{TALOS_cluster}}
#' in your main workflow.
#'
#' Silhouette width is computed in PCA space using
#' \code{cluster::silhouette()}.  Because pairwise distances are O(N²),
#' silhouette is automatically skipped when \code{ncol(sce) > max_cells_sil};
#' in that case modularity alone is used to select the suggested k.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"PCA"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_pca}} first).
#' @param k_range Integer vector of k values to test.
#'   Default \code{seq(5, 50, by = 5)}.
#' @param n_pcs Integer. Number of PCA components used for graph construction
#'   and silhouette distances. Default \code{30}.
#' @param type Character. Edge weight scheme for \code{scran::buildSNNGraph}:
#'   \code{"rank"} (default), \code{"number"}, or \code{"jaccard"}.
#' @param method Character. Clustering algorithm: \code{"leiden"} (default)
#'   or \code{"louvain"}.
#' @param resolution Numeric. Clustering resolution. Default \code{1.0}.
#' @param compute_silhouette Logical. Compute mean silhouette width.
#'   Automatically disabled when \code{ncol(sce) > max_cells_sil}.
#'   Default \code{TRUE}.
#' @param max_cells_sil Integer. Cell count above which silhouette is skipped.
#'   Default \code{10000}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print per-k progress and a final summary.
#'   Default \code{TRUE}.
#'
#' @return A \code{talos_k_sweep} list with \code{results}, \code{plot},
#'   \code{best_k}, and \code{params}.
#' @export
TALOS_tune_k <- function(sce,
                          k_range            = seq(5L, 50L, by = 5L),
                          n_pcs              = 30L,
                          type               = c("rank", "number", "jaccard"),
                          method             = c("leiden", "louvain"),
                          resolution         = 1.0,
                          compute_silhouette = TRUE,
                          max_cells_sil      = 10000L,
                          seed               = 42L,
                          verbose            = TRUE) {

  type   <- match.arg(type)
  method <- match.arg(method)
  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_cells <- ncol(sce)
  n_pcs   <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
  k_range <- sort(unique(as.integer(k_range)))

  # ── Silhouette feasibility ───────────────────────────────────────────────────
  do_sil <- isTRUE(compute_silhouette) && n_cells <= as.integer(max_cells_sil)
  if (isTRUE(compute_silhouette) && !do_sil)
    message(sprintf(
      "Silhouette skipped: %s cells exceeds max_cells_sil (%s). Using modularity.",
      format(n_cells, big.mark = ","),
      format(as.integer(max_cells_sil), big.mark = ",")))

  # ── Pre-compute PCA submatrix and distance matrix (once, outside loop) ───────
  pca_mat  <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]
  dist_mat <- if (do_sil) dist(pca_mat) else NULL

  if (verbose)
    message(sprintf(
      "\u2500\u2500 TALOS: k sweep %s\n  k values   : %s\n  n_pcs      : %s  |  type: %s\n  method     : %s  |  resolution: %s\n  Silhouette : %s\n%s",
      strrep("\u2500", 40),
      paste(k_range, collapse = ", "),
      n_pcs, type, method, resolution,
      if (do_sil) "yes" else "no (modularity used for selection)",
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results <- vector("list", length(k_range))

  for (i in seq_along(k_range)) {
    ki <- k_range[i]
    if (verbose) message(sprintf("  Testing k = %d ...", ki))

    g <- scran::buildSNNGraph(t(pca_mat), k = ki, type = type)

    set.seed(seed)
    communities <- if (method == "leiden") {
      igraph::cluster_leiden(g, resolution_parameter = resolution)
    } else {
      igraph::cluster_louvain(g, resolution = resolution)
    }

    labels <- as.integer(igraph::membership(communities))
    n_cl   <- length(unique(labels))
    modul  <- igraph::modularity(g, labels)

    mean_sil <- if (do_sil) {
      sil <- cluster::silhouette(labels, dist_mat)
      mean(sil[, "sil_width"])
    } else NA_real_

    results[[i]] <- data.frame(
      k          = ki,
      n_clusters = n_cl,
      modularity = round(modul, 4),
      mean_sil   = if (do_sil) round(mean_sil, 4) else NA_real_
    )
  }

  res_df <- do.call(rbind, results)

  best_k <- if (do_sil) res_df$k[which.max(res_df$mean_sil)] else
              res_df$k[which.max(res_df$modularity)]

  p <- .talos_tune_k_plot(res_df, best_k, do_sil)

  if (verbose)
    message(sprintf(
      "%s\n  Suggested k: %d (%s)\n%s",
      strrep("\u2500", 56),
      best_k,
      if (do_sil) "highest mean silhouette width" else "highest modularity",
      strrep("\u2500", 56)))

  structure(
    list(results = res_df,
         plot    = p,
         best_k  = best_k,
         params  = list(k_range    = k_range,
                        n_pcs      = n_pcs,
                        type       = type,
                        method     = method,
                        resolution = resolution,
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
#' Applies Leiden and/or Louvain community detection across a range of
#' resolution values on the pre-built SNN graph, scoring each combination on
#' four complementary metrics:
#'
#' \describe{
#'   \item{n_clusters}{Number of clusters produced.}
#'   \item{modularity}{Graph modularity of the partition.}
#'   \item{mean_sil}{Mean silhouette width in PCA space.  Skipped when cell
#'     count exceeds \code{max_cells_sil}.}
#'   \item{min_cl_size}{Minimum cluster size — flags pathological resolutions
#'     producing singleton or near-singleton clusters.}
#' }
#'
#' The PCA distance matrix for silhouette is computed once before the sweep
#' loop, so runtime scales with the number of resolution values rather than
#' the number of cells per iteration.
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
#' @param seed Integer.  Random seed.  Default \code{42L}.
#' @param verbose Logical.  Print per-combination progress.  Default \code{TRUE}.
#'
#' @return A \code{talos_resolution_sweep} list with \code{results},
#'   \code{plot}, \code{best_params}, \code{best_per_method}, and \code{params}.
#' @export
TALOS_tune_resolution <- function(sce,
                                   resolution_range   = seq(0.1, 2.0, by = 0.1),
                                   method             = c("both", "leiden", "louvain"),
                                   n_pcs              = 30L,
                                   compute_silhouette = TRUE,
                                   max_cells_sil      = 10000L,
                                   seed               = 42L,
                                   verbose            = TRUE) {

  method         <- match.arg(method)
  methods_to_run <- if (method == "both") c("leiden", "louvain") else method

  g <- metadata(sce)$snn_graph
  if (is.null(g))
    stop("SNN graph not found. Run TALOS_build_graph() first.", call. = FALSE)

  resolution_range <- as.numeric(resolution_range)
  n_runs           <- length(resolution_range) * length(methods_to_run)

  # ── Silhouette setup ─────────────────────────────────────────────────────────
  n_cells <- ncol(sce)
  do_sil  <- isTRUE(compute_silhouette) && n_cells <= as.integer(max_cells_sil)

  if (isTRUE(compute_silhouette) && !do_sil)
    message(sprintf(
      "Silhouette skipped: %s cells exceeds max_cells_sil (%s). Using modularity.",
      format(n_cells, big.mark = ","),
      format(as.integer(max_cells_sil), big.mark = ",")))

  dist_mat <- NULL
  if (do_sil) {
    if (!"PCA" %in% reducedDimNames(sce))
      stop("'PCA' not found. Run TALOS_run_pca() first.", call. = FALSE)
    n_pcs   <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
    pca_mat <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]
    if (verbose) message("  Pre-computing distance matrix in PCA space ...")
    dist_mat <- dist(pca_mat)
  }

  if (verbose)
    message(sprintf(
      "\u2500\u2500 TALOS: resolution sweep %s\n  Resolutions: %s\n  Method(s)  : %s\n  Silhouette : %s\n%s",
      strrep("\u2500", 31),
      paste(round(resolution_range, 2), collapse = ", "),
      paste(methods_to_run, collapse = " + "),
      if (do_sil) "yes" else "no",
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results <- vector("list", n_runs)
  run_i   <- 0L

  for (m in methods_to_run) {
    for (res in resolution_range) {
      run_i <- run_i + 1L
      if (verbose)
        message(sprintf("  (%d/%d) method = %s, resolution = %.2f ...",
                        run_i, n_runs, m, res))

      set.seed(seed)
      communities <- if (m == "leiden") {
        igraph::cluster_leiden(g, resolution_parameter = res)
      } else {
        igraph::cluster_louvain(g, resolution = res)
      }

      labels   <- as.integer(igraph::membership(communities))
      n_cl     <- length(unique(labels))
      modul    <- igraph::modularity(g, labels)
      min_size <- min(as.integer(table(labels)))

      mean_sil <- if (do_sil && n_cl > 1L) {
        sil <- cluster::silhouette(labels, dist_mat)
        mean(sil[, "sil_width"])
      } else NA_real_

      results[[run_i]] <- data.frame(
        method      = m,
        resolution  = res,
        n_clusters  = n_cl,
        modularity  = round(modul, 4),
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
    message(sprintf(
      "%s\n  Best overall: %s, resolution = %.2f (%s)\n%s",
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
