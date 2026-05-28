###############################################################################
########### Consensus Clustering Visualization Functions ###########
###############################################################################

#' Plot Consensus CDF Curves
#'
#' @description Plots empirical cumulative distribution functions (CDFs) of
#' pairwise consensus values for each k. A well-separated clustering shows a
#' CDF with values concentrated near 0 and 1 (flat/horizontal line). The
#' optimal k is where the CDF stops improving substantially.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param k_range Integer vector. Which k values to include. If NULL, uses all. Default = NULL.
#' @param title Character. Plot title. Default = "Consensus CDF".
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_consensus_cdf <- function(cc_result, k_range = NULL, title = "Consensus CDF") {

  if (is.null(k_range)) k_range <- cc_result$k_range

  df_list <- lapply(k_range, function(k) {
    vals <- cc_result$cdf_values[[k]]
    n    <- length(vals)
    data.frame(
      consensus = vals,
      cdf       = seq_len(n) / n,
      k         = factor(k)
    )
  })
  df <- do.call(rbind, df_list)

  pal <- RColorBrewer::brewer.pal(max(3, length(k_range)), "Set2")[seq_len(length(k_range))]

  ggplot2::ggplot(df, ggplot2::aes(x = consensus, y = cdf, color = k)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = pal, name = "k") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "Consensus Index", y = "CDF") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )
}


#' Plot Consensus Matrix Heatmap
#'
#' @description Plots the averaged consensus matrix for a chosen k as a heatmap.
#' Samples are reordered by cluster assignment, making block structure visible.
#' Dark blue = always co-clustered; white = never co-clustered.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param k Integer. Number of clusters to display.
#' @param title Character. Plot title. If NULL, auto-generated. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_consensus_matrix <- function(cc_result, k, title = NULL) {

  cm   <- cc_result$consensus_matrices[[k]]
  asgn <- cc_result$cluster_assignments[[k]]

  # Reorder samples by cluster assignment
  sample_order <- order(asgn)
  cm_ord       <- cm[sample_order, sample_order]
  asgn_ord     <- asgn[sample_order]
  n            <- nrow(cm_ord)

  df <- data.frame(
    row   = rep(seq_len(n), times = n),
    col   = rep(seq_len(n), each  = n),
    value = as.vector(cm_ord)
  )

  # Cluster boundary positions for annotation lines
  boundaries <- cumsum(table(factor(asgn_ord, levels = seq_len(k)))) + 0.5
  boundaries <- boundaries[-length(boundaries)]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#c6dbef", "#2166ac"),
      limits = c(0, 1),
      name   = "Consensus"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = if (!is.null(title)) title else paste0("Consensus Matrix (k=", k, ")"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5),
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank()
    )

  if (length(boundaries) > 0) {
    p <- p +
      ggplot2::geom_hline(yintercept = n - boundaries + 0.5,
                          color = "black", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = boundaries,
                          color = "black", linewidth = 0.4)
  }

  return(p)
}


#' Plot Within-Cluster Sum-of-Squares (WCSS) Elbow
#'
#' @description Plots WCSS against k. The "elbow" — where adding another cluster
#' gives diminishing returns — suggests the optimal k.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param title Character. Plot title. Default = "WCSS Elbow Plot".
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_wcss <- function(cc_result, title = "WCSS Elbow Plot") {

  df <- data.frame(
    k    = cc_result$k_range,
    wcss = cc_result$wcss
  )

  ggplot2::ggplot(df, ggplot2::aes(x = k, y = wcss)) +
    ggplot2::geom_line(color = "#4E79A7", linewidth = 0.9) +
    ggplot2::geom_point(color = "#4E79A7", size = 3) +
    ggplot2::scale_x_continuous(breaks = cc_result$k_range) +
    ggplot2::labs(title = title, x = "Number of Clusters (k)",
                  y = "Within-Cluster Sum of Squares") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}


#' Plot Delta Area (Relative CDF Area Change)
#'
#' @description Plots the relative change in CDF area between consecutive k values.
#' A large positive delta suggests that k explains substantially more structure
#' than k-1. Values plateau or drop near the optimal k.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param title Character. Plot title. Default = "Delta Area".
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_delta_area <- function(cc_result, title = "Delta Area") {

  df <- data.frame(
    k          = cc_result$k_range,
    delta_area = cc_result$delta_area
  )

  ggplot2::ggplot(df, ggplot2::aes(x = factor(k), y = delta_area)) +
    ggplot2::geom_col(fill = "#4E79A7", width = 0.6) +
    ggplot2::labs(title = title, x = "Number of Clusters (k)",
                  y = "Relative Change in CDF Area") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}


#' Multi-Panel Diagnostic Summary for Consensus Clustering
#'
#' @description Combines CDF, WCSS elbow, and delta area into a single
#' three-panel diagnostic figure.
#'
#' @param cc_result Output of \code{ARTEMIS_consensus_cluster()}.
#' @param title Character. Overall figure title. Default = "Consensus Clustering Diagnostics".
#'
#' @return A ggplot object (patchwork layout).
#'
#' @details Requires the \code{patchwork} package.
#'
#' @export
AETHER_plot_consensus_summary <- function(cc_result,
                                           title = "Consensus Clustering Diagnostics") {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }

  p_cdf   <- AETHER_plot_consensus_cdf(cc_result, title = "CDF")
  p_wcss  <- AETHER_plot_wcss(cc_result, title = "WCSS")
  p_delta <- AETHER_plot_delta_area(cc_result, title = "Delta Area")

  (p_cdf | p_wcss | p_delta) +
    patchwork::plot_annotation(title = title,
                               theme = ggplot2::theme(
                                 plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
                               ))
}


#' Stability Histogram for Consensus Clustering
#'
#' @description For each patient, plots how many of the \code{n_runs} independent
#' runs assigned them to the same (modal) cluster. A bar concentrated at \code{n_runs}
#' indicates high stability.
#'
#' @param stability_result Output of \code{ARTEMIS_consensus_stability()}.
#' @param title Character. Plot title. If NULL, auto-generated. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_stability_histogram <- function(stability_result, title = NULL) {

  k      <- stability_result$k
  n_runs <- stability_result$n_runs
  counts <- stability_result$stability_counts

  df <- as.data.frame(table(counts))
  colnames(df) <- c("n_consistent", "n_patients")
  df$n_consistent <- as.integer(as.character(df$n_consistent))

  # Fill full range so bars at missing counts show as 0
  all_counts <- data.frame(n_consistent = seq_len(n_runs))
  df <- merge(all_counts, df, by = "n_consistent", all.x = TRUE)
  df$n_patients[is.na(df$n_patients)] <- 0

  pct_stable <- stability_result$pct_stable

  ggplot2::ggplot(df, ggplot2::aes(x = factor(n_consistent), y = n_patients)) +
    ggplot2::geom_col(fill = "#4E79A7", width = 0.7) +
    ggplot2::labs(
      title = if (!is.null(title)) title else
        paste0("Stability Histogram (k=", k, ")"),
      subtitle = paste0(pct_stable, "% of patients consistently assigned across all ",
                        n_runs, " runs"),
      x = paste0("Times assigned to modal cluster (of ", n_runs, " runs)"),
      y = "Number of Patients"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
    )
}


#' Alluvial Plot of Cluster Stability Across Runs
#'
#' @description Shows how patients flow between cluster assignments across
#' \code{n_runs} independent runs. Stable patients remain within the same
#' stratum; unstable patients cross between strata.
#'
#' @param stability_result Output of \code{ARTEMIS_consensus_stability()}.
#' @param title Character. Plot title. If NULL, auto-generated. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @details Requires the \code{ggalluvial} package. If not available, falls back
#' to \code{AETHER_plot_stability_histogram()}.
#'
#' @export
AETHER_plot_stability_alluvial <- function(stability_result, title = NULL) {

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    cat("[AETHER] ggalluvial not available — falling back to stability histogram.\n")
    cat("    Install with: install.packages('ggalluvial')\n")
    return(AETHER_plot_stability_histogram(stability_result, title = title))
  }

  k      <- stability_result$k
  n_runs <- stability_result$n_runs
  am     <- stability_result$assignment_matrix

  # Long format: patient x run x cluster
  df <- data.frame(
    patient = rep(seq_len(nrow(am)), times = n_runs),
    run     = rep(seq_len(n_runs), each = nrow(am)),
    cluster = as.integer(as.vector(am))
  )
  df$run     <- factor(df$run, labels = paste0("Run ", seq_len(n_runs)))
  df$cluster <- factor(df$cluster)

  # Subsample patients for readability if cohort is large
  n_display <- min(nrow(am), 300)
  if (nrow(am) > n_display) {
    keep_idx  <- sample(nrow(am), n_display)
    df        <- df[df$patient %in% keep_idx, ]
  }

  pal <- RColorBrewer::brewer.pal(max(3, k), "Set2")[seq_len(k)]

  ggplot2::ggplot(df,
    ggplot2::aes(x = run, stratum = cluster,
                 alluvium = patient, fill = cluster)) +
    ggalluvial::geom_flow(alpha = 0.4) +
    ggalluvial::geom_stratum(width = 0.4) +
    ggplot2::scale_fill_manual(values = pal, name = "Cluster") +
    ggplot2::labs(
      title = if (!is.null(title)) title else
        paste0("Cluster Assignment Stability (k=", k, ")"),
      subtitle = if (nrow(stability_result$assignment_matrix) > n_display)
        paste0("Showing ", n_display, " randomly sampled patients") else NULL,
      x = NULL, y = "Number of Patients"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
