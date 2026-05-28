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
