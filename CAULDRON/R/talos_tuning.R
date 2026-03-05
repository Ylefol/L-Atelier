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


#' Sweep UMAP parameters
#'
#' Tests all combinations of \code{n_neighbors} and \code{min_dist} and
#' scores each embedding on two complementary metrics:
#' \enumerate{
#'   \item \strong{KNN overlap} — mean Jaccard similarity between each cell's
#'     k nearest neighbours in PCA space and in the UMAP embedding.  Measures
#'     how faithfully local neighbourhood structure is preserved.
#'   \item \strong{Silhouette width} — mean silhouette width computed in 2D
#'     UMAP space using existing cluster labels.  Skipped when
#'     \code{cluster_col} is not in \code{colData}.
#' }
#' Both metrics are normalised to \eqn{[0, 1]} and averaged into a composite
#' score used to select the best parameter combination.
#'
#' The KNN reference is computed once in PCA space using
#' \code{BiocNeighbors::findKNN}.
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
#' @return A \code{talos_embedding_sweep} list with \code{results},
#'   \code{plot}, \code{best_params}, and \code{params}.
#' @export
TALOS_tune_umap <- function(sce,
                              n_neighbors_range = c(10L, 15L, 20L, 30L, 50L),
                              min_dist_range    = c(0.05, 0.1, 0.3, 0.5),
                              n_pcs             = 30L,
                              knn_k             = 15L,
                              cluster_col       = "cluster",
                              seed              = 42L,
                              verbose           = TRUE) {

  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_pcs             <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
  n_neighbors_range <- as.integer(n_neighbors_range)
  knn_k             <- as.integer(knn_k)

  has_clusters <- cluster_col %in% names(colData(sce))
  if (!has_clusters && verbose)
    message(sprintf(
      "  Note: '%s' not found in colData — silhouette will be skipped.",
      cluster_col))

  pca_mat <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]

  # ── Pre-compute KNN reference in PCA space (once) ────────────────────────────
  if (verbose) message("  Computing KNN reference in PCA space ...")
  knn_ref <- BiocNeighbors::findKNN(pca_mat, k = knn_k)$index   # cells × k

  n_combos <- length(n_neighbors_range) * length(min_dist_range)
  if (verbose)
    message(sprintf(
      "\u2500\u2500 TALOS: UMAP sweep %s\n  n_neighbors: %s\n  min_dist   : %s\n  n_pcs      : %s  |  knn_k: %s\n  Combos     : %s\n%s",
      strrep("\u2500", 37),
      paste(n_neighbors_range, collapse = ", "),
      paste(min_dist_range,    collapse = ", "),
      n_pcs, knn_k, n_combos,
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results    <- vector("list", n_combos)
  embeddings <- vector("list", n_combos)
  run_i      <- 0L

  for (nn in n_neighbors_range) {
    for (md in min_dist_range) {
      run_i <- run_i + 1L
      if (verbose)
        message(sprintf("  (%d/%d) n_neighbors = %d, min_dist = %.2f ...",
                        run_i, n_combos, nn, md))

      set.seed(seed)
      sce_tmp <- scater::runUMAP(sce,
                                  dimred      = "PCA",
                                  n_dimred    = n_pcs,
                                  n_neighbors = nn,
                                  min_dist    = md,
                                  name        = "UMAP_sweep")

      emb_mat            <- reducedDim(sce_tmp, "UMAP_sweep")
      embeddings[[run_i]] <- emb_mat
      knn_ov             <- .talos_knn_overlap(knn_ref, emb_mat, knn_k)

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

  # ── Composite score and best params ──────────────────────────────────────────
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
    message(sprintf(
      "%s\n  Best combo : n_neighbors = %d, min_dist = %.2f  (composite = %.3f)\n%s",
      strrep("\u2500", 56),
      best_params$n_neighbors, best_params$min_dist,
      res_df$composite[best_i],
      strrep("\u2500", 56)))

  structure(
    list(results     = res_df,
         embeddings  = embeddings,
         plot        = p,
         best_params = best_params,
         params      = list(type              = "umap",
                            n_neighbors_range = n_neighbors_range,
                            min_dist_range    = min_dist_range,
                            n_pcs             = n_pcs,
                            knn_k             = knn_k,
                            cluster_col       = cluster_col,
                            seed              = seed,
                            has_clusters      = has_clusters)),
    class = "talos_embedding_sweep"
  )
}


#' Sweep tSNE perplexity
#'
#' Tests a range of perplexity values and scores each tSNE embedding using
#' the same two metrics as \code{\link{TALOS_tune_umap}}: KNN overlap (local
#' structure preservation vs PCA space) and mean silhouette width in 2D
#' embedding space.  Both are normalised to \eqn{[0, 1]} and averaged into a
#' composite score.
#'
#' Perplexity values that exceed \eqn{N / 3} (where \eqn{N} is the number of
#' cells) are automatically removed before sweeping, as they are invalid for
#' tSNE.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"PCA"} in
#'   \code{reducedDims}.
#' @param perplexity_range Numeric vector of perplexity values to test.
#'   Invalid values (> N/3) are removed automatically.
#'   Default \code{c(10, 20, 30, 50, 100)}.
#' @param n_pcs Integer. PCA components used for tSNE and KNN reference.
#'   Default \code{30}.
#' @param max_iter Integer. tSNE iterations per run (fixed; not swept).
#'   Default \code{1000}.
#' @param knn_k Integer. Nearest neighbours for KNN overlap metric.
#'   Default \code{15}.
#' @param cluster_col Character. \code{colData} column for silhouette labels.
#'   Default \code{"cluster"}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{talos_embedding_sweep} list with \code{results},
#'   \code{plot}, \code{best_params}, and \code{params}.
#' @export
TALOS_tune_tsne <- function(sce,
                              perplexity_range = c(10, 20, 30, 50, 100),
                              n_pcs            = 30L,
                              max_iter         = 1000L,
                              knn_k            = 15L,
                              cluster_col      = "cluster",
                              seed             = 42L,
                              verbose          = TRUE) {

  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_pcs  <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
  n_cells <- ncol(sce)
  knn_k  <- as.integer(knn_k)

  # ── Remove invalid perplexity values ─────────────────────────────────────────
  max_perp         <- floor(n_cells / 3)
  perplexity_range <- perplexity_range[perplexity_range <= max_perp]
  if (length(perplexity_range) == 0L)
    stop("All perplexity values exceed N/3 = ", max_perp,
         ". Provide smaller perplexity values.", call. = FALSE)

  has_clusters <- cluster_col %in% names(colData(sce))
  if (!has_clusters && verbose)
    message(sprintf(
      "  Note: '%s' not found in colData — silhouette will be skipped.",
      cluster_col))

  pca_mat <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]

  # ── Pre-compute KNN reference (once) ─────────────────────────────────────────
  if (verbose) message("  Computing KNN reference in PCA space ...")
  knn_ref <- BiocNeighbors::findKNN(pca_mat, k = knn_k)$index

  n_combos <- length(perplexity_range)
  if (verbose)
    message(sprintf(
      "\u2500\u2500 TALOS: tSNE sweep %s\n  Perplexity : %s\n  n_pcs      : %s  |  max_iter: %s\n  knn_k      : %s\n%s",
      strrep("\u2500", 37),
      paste(perplexity_range, collapse = ", "),
      n_pcs, max_iter, knn_k,
      strrep("\u2500", 56)))

  # ── Sweep ────────────────────────────────────────────────────────────────────
  results    <- vector("list", n_combos)
  embeddings <- vector("list", n_combos)

  for (i in seq_along(perplexity_range)) {
    perp <- perplexity_range[i]
    if (verbose)
      message(sprintf("  (%d/%d) perplexity = %g ...", i, n_combos, perp))

    set.seed(seed)
    sce_tmp <- scater::runTSNE(sce,
                                dimred     = "PCA",
                                n_dimred   = n_pcs,
                                perplexity = perp,
                                max_iter   = as.integer(max_iter),
                                name       = "TSNE_sweep")

    emb_mat      <- reducedDim(sce_tmp, "TSNE_sweep")
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
    message(sprintf(
      "%s\n  Best perplexity: %g  (composite = %.3f)\n%s",
      strrep("\u2500", 56),
      best_params$perplexity, res_df$composite[best_i],
      strrep("\u2500", 56)))

  structure(
    list(results     = res_df,
         embeddings  = embeddings,
         plot        = p,
         best_params = best_params,
         params      = list(type             = "tsne",
                            perplexity_range = perplexity_range,
                            n_pcs            = n_pcs,
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
