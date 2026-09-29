###############################################################################
########### Consensus Clustering Embedding Plots ###########
###############################################################################

#' PCA Plot Coloured by Consensus Cluster Assignment
#'
#' @description Projects patients onto the first two principal components of the
#' mean preprocessed data (log1p-transformed, z-scored, averaged across all
#' imputations) and colours points by their consensus cluster assignment at the
#' chosen k.
#'
#' @details
#' PCA is computed on \code{cc_result$scaled_mean_data}, which is already
#' centred and scaled, so \code{prcomp} is called with \code{center=FALSE} and
#' \code{scale.=FALSE}. Cluster separation in PC space is an indirect view of
#' the consensus structure: clusters defined by the consensus matrix (co-
#' clustering frequency) need not be linearly separable in PC space, so partial
#' overlap is expected for clinical data.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param k Integer. Cluster solution to display.
#' @param point_size Numeric. Point size. Default = 1.5.
#' @param alpha Numeric. Point transparency (0–1). Default = 0.7.
#' @param show_ellipse Logical. Draw a 95% confidence ellipse per cluster. Default = TRUE.
#' @param ellipse_level Numeric. Confidence level for ellipses (0–1). Default = 0.95.
#' @param show_centroids Logical. Mark the per-cluster centroid with a filled symbol. Default = TRUE.
#' @param centroid_size Numeric. Size of the centroid symbol. Default = 4.
#' @param title Character. Plot title. If NULL, auto-generated. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_cluster_pca <- function(cc_result, k,
                                     point_size = 1.5, alpha = 0.7,
                                     show_ellipse   = TRUE, ellipse_level = 0.95,
                                     show_centroids = TRUE, centroid_size = 4,
                                     title = NULL) {

  mat  <- cc_result$scaled_mean_data
  asgn <- cc_result$cluster_assignments[[k]]

  common <- intersect(rownames(mat), names(asgn))
  mat    <- mat[common, , drop = FALSE]
  asgn   <- asgn[common]

  pca     <- prcomp(mat, center = FALSE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  clust_levels <- paste0("C", seq_len(k))
  df <- data.frame(
    PC1     = pca$x[, 1],
    PC2     = pca$x[, 2],
    cluster = factor(asgn, levels = seq_len(k), labels = clust_levels)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = cluster)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_brewer(palette = "Set2", name = "Cluster") +
    ggplot2::labs(
      title = if (!is.null(title)) title else paste0("PCA — Consensus Clusters (k=", k, ")"),
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )

  if (show_ellipse)
    p <- p + ggplot2::stat_ellipse(level = ellipse_level, linewidth = 0.7, alpha = 0.8)

  if (show_centroids) {
    centroids <- do.call(rbind, lapply(clust_levels, function(cl) {
      sub <- df[df$cluster == cl, ]
      data.frame(PC1     = mean(sub$PC1),
                 PC2     = mean(sub$PC2),
                 cluster = factor(cl, levels = clust_levels))
    }))
    p <- p +
      ggplot2::geom_point(data = centroids,
                          ggplot2::aes(x = PC1, y = PC2, fill = cluster),
                          shape = 21, size = centroid_size,
                          color = "black", stroke = 1) +
      ggplot2::scale_fill_brewer(palette = "Set2", guide = "none")
  }

  return(p)
}


#' UMAP Plot Coloured by Consensus Cluster Assignment or Item Consensus
#'
#' @description Sibling to \code{AETHER_plot_cluster_pca()}: embeds Train
#' samples (\code{cc_result$scaled_mean_data}) with UMAP instead of PCA, and
#' optionally embeds already Train-standardized new samples (e.g. Test
#' patients) into that SAME fitted embedding via \code{uwot::umap_transform()},
#' overlaid on top. Purely a visualization/reporting layer — the actual
#' out-of-sample confidence number comes from \code{ARTEMIS_project_confidence()}
#' (computed in the real factor space, not this 2D embedding); this plot does
#' not feed back into that computation.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param k Integer. Cluster solution to display.
#' @param item_consensus Output of \code{ARTEMIS_item_consensus()} for the
#'   same k. Required when \code{color_by} is \code{"margin"} or
#'   \code{"assigned_consensus"}; ignored for \code{color_by = "cluster"}.
#' @param color_by Character. \code{"cluster"} (categorical, default),
#'   \code{"margin"}, or \code{"assigned_consensus"} (both continuous, from
#'   \code{item_consensus}).
#' @param new_points Optional matrix/data frame of already Train-standardized
#'   new sample factor scores (same column space as
#'   \code{cc_result$scaled_mean_data}) to overlay on the fitted embedding.
#'   Default = NULL (Train only).
#' @param project_confidence Output of \code{ARTEMIS_project_confidence()} for
#'   the same \code{new_points}. Required when \code{new_color_by != "fixed"}.
#' @param new_color_by Character. How to colour the \code{new_points} overlay:
#'   \code{"fixed"} (default — flat yellow/black marker, no data mapping),
#'   \code{"purity"} (\code{neighbor_cluster_purity} — do a sample's nearest
#'   Train neighbors agree on a cluster), or \code{"margin"}
#'   (\code{neighbor_margin} — how individually solid those neighbors are).
#'   Uses a SEPARATE colour scale (\code{fill}, not \code{color}) from
#'   whatever \code{color_by} maps the Train points to, so Train points stay
#'   coloured by item-consensus margin (their own real quantity) even when
#'   the overlay is coloured by a different, projected-sample-only quantity —
#'   these are deliberately NOT forced to share one scale or one meaning.
#' @param supervised Logical. If TRUE, passes the known cluster assignment at
#'   k as \code{y} to \code{uwot::umap()}, which pulls each cluster's points
#'   together in the embedding. IMPORTANT: this makes clusters look more
#'   separated BY CONSTRUCTION, using labels the embedding wouldn't otherwise
#'   have — it is not evidence of unsupervised structure in the factor space,
#'   and should not be read as such. Useful as a clean illustrative/reference
#'   layout once cluster structure is already trusted from other diagnostics
#'   (consensus matrix, item consensus), not as a diagnostic itself.
#'   Default = FALSE.
#' @param target_weight Numeric (0-1). Only used when \code{supervised=TRUE}:
#'   how much \code{uwot::umap()} weights the label information vs. the
#'   underlying data structure (0 = ignore labels/fully unsupervised, 1 =
#'   ignore data/embed purely by label). Passed through to \code{uwot::umap()}
#'   as-is. Default = 0.5 (uwot's own default).
#' @param n_neighbors Integer. UMAP's own neighborhood-graph size (unrelated
#'   to \code{ARTEMIS_project_confidence()}'s neighbor count). Default = 15.
#' @param min_dist Numeric. UMAP \code{min_dist} parameter. Default = 0.1.
#' @param seed Integer. Random seed for the UMAP fit. Default = 42.
#' @param point_size Numeric. Train point size. Default = 1.5.
#' @param alpha Numeric. Train point transparency. Default = 0.7.
#' @param new_point_size Numeric. Overlaid new-sample point size. Default = 3.
#' @param new_point_shape Integer. Overlaid new-sample point shape. Default =
#'   21 (filled circle with a separate outline/fill — a dark border so
#'   overlaid points read as dots distinct from the borderless Train points,
#'   while still rendering whatever \code{new_color_by} maps to \code{fill}).
#' @param label_new Logical. Label overlaid new points with their sample ID
#'   (uses \code{ggrepel} if installed, falls back to plain text). Default = TRUE.
#' @param title Character. Plot title. If NULL, auto-generated. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @details
#' Requires package \code{uwot} (\code{install.packages("uwot")}).
#'
#' @export
AETHER_plot_cluster_umap <- function(cc_result, k, item_consensus = NULL,
                                      color_by = c("cluster", "margin", "assigned_consensus"),
                                      new_points = NULL,
                                      project_confidence = NULL,
                                      new_color_by = c("fixed", "purity", "margin"),
                                      supervised = FALSE, target_weight = 0.5,
                                      n_neighbors = 15, min_dist = 0.1, seed = 42,
                                      point_size = 1.5, alpha = 0.7,
                                      new_point_size = 3, new_point_shape = 21,
                                      label_new = TRUE,
                                      title = NULL) {

  if (!requireNamespace("uwot", quietly = TRUE))
    stop("Package 'uwot' is required. Install with: install.packages('uwot')")

  color_by     <- match.arg(color_by)
  new_color_by <- match.arg(new_color_by)
  if (color_by != "cluster" && is.null(item_consensus))
    stop("item_consensus (output of ARTEMIS_item_consensus()) is required when color_by = '", color_by, "'")
  if (new_color_by != "fixed" && is.null(project_confidence))
    stop("project_confidence (output of ARTEMIS_project_confidence()) is required when new_color_by = '", new_color_by, "'")

  mat  <- cc_result$scaled_mean_data
  asgn <- cc_result$cluster_assignments[[k]]
  common <- intersect(rownames(mat), names(asgn))
  mat  <- mat[common, , drop = FALSE]
  asgn <- asgn[common]

  if (supervised)
    cat("[AETHER] Supervised UMAP: embedding is shaped by the known cluster",
        "labels (target_weight=", target_weight, ") -- visual separation here",
        "reflects that supervision, not unsupervised structure on its own.\n")

  set.seed(seed)
  fit <- uwot::umap(mat, n_neighbors = n_neighbors, min_dist = min_dist,
                     y = if (supervised) factor(asgn) else NULL,
                     target_weight = target_weight,
                     ret_model = TRUE)
  emb <- fit$embedding
  rownames(emb) <- common

  clust_levels <- paste0("C", seq_len(k))
  df <- data.frame(
    UMAP1   = emb[, 1],
    UMAP2   = emb[, 2],
    cluster = factor(asgn, levels = seq_len(k), labels = clust_levels)
  )

  if (color_by != "cluster") {
    ic_summary <- item_consensus$summary
    rownames(ic_summary) <- ic_summary$sample
    ic_col <- if (color_by == "margin") "margin" else "assigned_consensus"
    df[[color_by]] <- ic_summary[common, ic_col]
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = UMAP1, y = UMAP2))

  if (color_by == "cluster") {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = cluster), size = point_size, alpha = alpha) +
      ggplot2::scale_color_brewer(palette = "Set2", name = "Cluster")
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]), size = point_size, alpha = alpha) +
      ggplot2::scale_color_gradient2(
        low = "firebrick", mid = "grey85", high = "steelblue",
        midpoint = if (color_by == "margin") 0 else 0.5,
        name = if (color_by == "margin") "Item-consensus\nmargin" else "Assigned\nconsensus")
  }

  if (!is.null(new_points)) {
    new_points <- as.matrix(new_points)[, colnames(mat), drop = FALSE]
    new_emb    <- uwot::umap_transform(new_points, fit)
    new_ids    <- rownames(new_points)
    if (is.null(new_ids)) new_ids <- paste0("new", seq_len(nrow(new_points)))
    new_df <- data.frame(UMAP1 = new_emb[, 1], UMAP2 = new_emb[, 2], sample = new_ids)

    if (new_color_by == "fixed") {
      p <- p + ggplot2::geom_point(
        data  = new_df, ggplot2::aes(x = UMAP1, y = UMAP2),
        shape = new_point_shape, size = new_point_size,
        color = "black", fill = "yellow", stroke = 1.2)
    } else {
      pc_summary <- project_confidence$summary
      rownames(pc_summary) <- pc_summary$sample
      pc_col <- if (new_color_by == "purity") "neighbor_cluster_purity" else "neighbor_margin"
      miss   <- setdiff(new_df$sample, pc_summary$sample)
      if (length(miss) > 0)
        stop("new_color_by = '", new_color_by, "' but project_confidence has no entry for: ",
             paste(miss, collapse = ", "))
      new_df[[pc_col]] <- pc_summary[new_df$sample, pc_col]

      p <- p +
        ggplot2::geom_point(
          data  = new_df, ggplot2::aes(x = UMAP1, y = UMAP2, fill = .data[[pc_col]]),
          shape = new_point_shape, size = new_point_size, color = "black", stroke = 1.2) +
        ggplot2::scale_fill_gradient2(
          low = "firebrick", mid = "grey85", high = "steelblue",
          midpoint = if (new_color_by == "purity") 0.5 else 0,
          name = if (new_color_by == "purity") "Neighbor cluster\npurity" else "Neighbor\nmargin")
    }

    if (label_new) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = new_df, ggplot2::aes(x = UMAP1, y = UMAP2, label = sample),
          size = 3, max.overlaps = 20)
      } else {
        p <- p + ggplot2::geom_text(
          data = new_df, ggplot2::aes(x = UMAP1, y = UMAP2, label = sample),
          size = 3, vjust = -1)
      }
    }
  }

  p <- p +
    ggplot2::labs(title = if (!is.null(title)) title else paste0(
      "UMAP", if (supervised) " (label-supervised)" else "",
      " — Consensus Clusters (k=", k, ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right")

  return(p)
}
