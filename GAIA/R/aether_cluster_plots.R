###############################################################################
########### Clustering Visualization Functions ###########
###############################################################################

#' Plot Silhouette for Cluster Quality Assessment
#'
#' @description Visualize silhouette widths for each sample, grouped by cluster.
#' Silhouette width measures how similar a sample is to its own cluster compared
#' to other clusters. Values range from -1 to 1, where higher is better.
#'
#' @param cluster_result An artemis_cluster object from ARTEMIS_cluster_mixed().
#' @param colors Named character vector. Custom colors for clusters.
#'   If NULL (default), uses a default palette.
#' @param show_avg Logical. Show average silhouette lines per cluster.
#'   Default = TRUE.
#' @param title Character. Plot title. Default = "Silhouette Plot".
#'
#' @return A ggplot object
#'
#' @details
#' Interpretation of silhouette width:
#' \itemize{
#'   \item > 0.7: Strong structure
#'   \item 0.5 - 0.7: Reasonable structure
#'   \item 0.25 - 0.5: Weak structure, could be artificial
#'   \item < 0.25: No substantial structure
#' }
#'
#' Negative values indicate samples that may be misclassified (closer to
#' another cluster than their assigned one).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clust <- ARTEMIS_cluster_mixed(data)
#' AETHER_plot_silhouette(clust)
#'
#' }
AETHER_plot_silhouette <- function(cluster_result,
                                    colors = NULL,
                                    show_avg = TRUE,
                                    title = "Silhouette Plot") {

  # Validate input

  if (!inherits(cluster_result, "artemis_cluster")) {
    stop("cluster_result must be an artemis_cluster object from ARTEMIS_cluster_mixed()")
  }

  # Extract silhouette info
  sil <- cluster_result$silhouette
  if (is.null(sil)) {
    stop("No silhouette information in cluster_result. ",
         "This may happen if k=1 or clustering failed.")
  }

  # Build data frame for plotting
  sil_df <- data.frame(
    cluster = factor(sil[, "cluster"]),
    neighbor = sil[, "neighbor"],
    sil_width = sil[, "sil_width"],
    stringsAsFactors = FALSE
  )

  # Sort within each cluster by silhouette width (descending)
  sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil_width), ]

  # Create sample order for plotting (grouped by cluster, sorted by width)
  sil_df$order <- seq_len(nrow(sil_df))

  # Calculate cluster averages
  cluster_avgs <- tapply(sil_df$sil_width, sil_df$cluster, mean)
  overall_avg <- mean(sil_df$sil_width)

  # Default colors
  if (is.null(colors)) {
    n_clusters <- length(unique(sil_df$cluster))
    colors <- setNames(
      scales::hue_pal()(n_clusters),
      levels(sil_df$cluster)
    )
  }

  # Build plot
  p <- ggplot(sil_df, aes(x = order, y = sil_width, fill = cluster)) +
    geom_col(width = 1) +
    coord_flip() +
    scale_fill_manual(values = colors, name = "Cluster") +
    labs(
      title = title,
      subtitle = paste0("Average silhouette width: ", round(overall_avg, 3)),
      x = "Samples (ordered by cluster and silhouette width)",
      y = "Silhouette Width"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Add reference lines
  p <- p +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = overall_avg, linetype = "dashed",
               color = "red", linewidth = 0.8)

  # Add cluster average lines if requested
  if (show_avg) {
    # Calculate positions for cluster average lines
    cluster_positions <- tapply(sil_df$order, sil_df$cluster, function(x) {
      c(min = min(x), max = max(x), mid = mean(x))
    })

    for (cl in names(cluster_positions)) {
      pos <- cluster_positions[[cl]]
      avg <- cluster_avgs[cl]

      # Add segment for cluster average
      p <- p + annotate(
        "segment",
        x = pos["min"] - 0.5, xend = pos["max"] + 0.5,
        y = avg, yend = avg,
        color = "darkgray", linewidth = 0.8, linetype = "dotted"
      )
    }
  }

  # Add interpretation guide
  p <- p +
    annotate("text", x = -Inf, y = 0.75, label = "Strong",
             hjust = -0.1, vjust = -0.5, size = 3, color = "gray40") +
    annotate("text", x = -Inf, y = 0.5, label = "Reasonable",
             hjust = -0.1, vjust = -0.5, size = 3, color = "gray40") +
    annotate("text", x = -Inf, y = 0.25, label = "Weak",
             hjust = -0.1, vjust = -0.5, size = 3, color = "gray40")

  return(p)
}


#' Plot Radar/Spider Chart for Cluster Profiles
#'
#' @description Visualize how clusters differ across multiple variables using
#' a radar (spider) chart. Each axis represents a variable, and each polygon
#' represents a cluster's profile.
#'
#' @param cluster_char An artemis_characterization object from
#'   ARTEMIS_characterize_clusters(), OR a named list of cluster means.
#' @param variables Character vector. Which variables to include in the plot.
#'   If NULL (default), uses top variables from characterization (max 10).
#' @param data Optional data frame. Required if cluster_char doesn't contain
#'   the original data or if you want to use different variables.
#' @param clusters Optional. Cluster assignments vector. Required if providing
#'   data instead of cluster_char.
#' @param normalize Logical. Normalize variables to 0-1 range for comparability.
#'   Default = TRUE.
#' @param colors Named character vector. Custom colors for clusters.
#'   If NULL (default), uses a default palette.
#' @param alpha Numeric. Fill transparency (0-1). Default = 0.2.
#' @param line_width Numeric. Width of polygon borders. Default = 1.
#' @param title Character. Plot title. Default = "Cluster Profiles".
#'
#' @return A ggplot object
#'
#' @details
#' For quantitative variables, the plot shows the mean value per cluster.
#' Variables are normalized to 0-1 range by default so they're comparable
#' on the same axes.
#'
#' Qualitative variables are excluded as they don't have a natural ordering
#' for radar plots. Use the cluster profiles table from characterization
#' to examine qualitative variable distributions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' char <- ARTEMIS_characterize_clusters(clust)
#' AETHER_plot_cluster_radar(char)
#'
#' # With specific variables
#' AETHER_plot_cluster_radar(char, variables = c("age", "score", "measure1"))
#'
#' }
AETHER_plot_cluster_radar <- function(cluster_char,
                                       variables = NULL,
                                       data = NULL,
                                       clusters = NULL,
                                       normalize = TRUE,
                                       colors = NULL,
                                       alpha = 0.2,
                                       line_width = 1,
                                       title = "Cluster Profiles") {

  # ---------------------------------------------------------------------------
  # Extract data and clusters
  # ---------------------------------------------------------------------------
  if (inherits(cluster_char, "artemis_characterization")) {
    # Extract from characterization object
    if (is.null(data)) {
      # Try to get data from the cluster profiles
      # We need the original data for computing means
      stop("Please provide 'data' argument with the original data frame used for clustering")
    }
    if (is.null(clusters)) {
      clusters <- cluster_char$clusters
    }

    # Get top quantitative variables if not specified
    if (is.null(variables)) {
      quanti_vars <- cluster_char$results$variable[cluster_char$results$type == "quantitative"]
      variables <- head(quanti_vars, 10)
    }
  } else {
    stop("cluster_char must be an artemis_characterization object")
  }

  # Validate we have what we need
  if (is.null(data) || is.null(clusters)) {
    stop("Both 'data' and 'clusters' are required")
  }

  # Filter to quantitative variables only
  quanti_mask <- sapply(data[, variables, drop = FALSE], is.numeric)
  if (!all(quanti_mask)) {
    excluded <- variables[!quanti_mask]
    warning("Excluding non-numeric variables: ", paste(excluded, collapse = ", "))
    variables <- variables[quanti_mask]
  }

  if (length(variables) < 3) {
    stop("Radar plot requires at least 3 numeric variables. Found: ", length(variables))
  }

  # ---------------------------------------------------------------------------
  # Compute cluster means
  # ---------------------------------------------------------------------------
  cluster_means <- aggregate(
    data[, variables, drop = FALSE],
    by = list(cluster = clusters),
    FUN = mean,
    na.rm = TRUE
  )

  # ---------------------------------------------------------------------------
  # Normalize to 0-1 if requested
  # ---------------------------------------------------------------------------
  if (normalize) {
    for (v in variables) {
      col_data <- data[[v]]
      min_val <- min(col_data, na.rm = TRUE)
      max_val <- max(col_data, na.rm = TRUE)
      range_val <- max_val - min_val

      if (range_val > 0) {
        cluster_means[[v]] <- (cluster_means[[v]] - min_val) / range_val
      } else {
        cluster_means[[v]] <- 0.5  # Constant variable
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Prepare data for radar plot
  # ---------------------------------------------------------------------------
  n_vars <- length(variables)
  n_clusters <- nrow(cluster_means)

  # Calculate angles for each variable
  angles <- seq(0, 2 * pi, length.out = n_vars + 1)[1:n_vars]

  # Build polygon data for each cluster
  radar_data <- list()

  for (i in seq_len(n_clusters)) {
    cl <- cluster_means$cluster[i]
    values <- as.numeric(cluster_means[i, variables])

    # Close the polygon by repeating first point
    values_closed <- c(values, values[1])
    angles_closed <- c(angles, angles[1])

    # Convert to Cartesian coordinates
    x <- values_closed * cos(angles_closed)
    y <- values_closed * sin(angles_closed)

    radar_data[[i]] <- data.frame(
      cluster = as.character(cl),
      variable = c(variables, variables[1]),
      value = values_closed,
      x = x,
      y = y,
      stringsAsFactors = FALSE
    )
  }

  plot_data <- do.call(rbind, radar_data)
  plot_data$cluster <- factor(plot_data$cluster)

  # ---------------------------------------------------------------------------
  # Create axis/grid data
  # ---------------------------------------------------------------------------
  # Axis lines from center to edge
  axis_data <- data.frame(
    variable = variables,
    x = cos(angles),
    y = sin(angles),
    xend = 0,
    yend = 0,
    stringsAsFactors = FALSE
  )

  # Axis labels
  label_data <- data.frame(
    variable = variables,
    x = 1.15 * cos(angles),
    y = 1.15 * sin(angles),
    stringsAsFactors = FALSE
  )

  # Concentric grid circles
  grid_levels <- c(0.25, 0.5, 0.75, 1.0)
  grid_data <- do.call(rbind, lapply(grid_levels, function(r) {
    theta <- seq(0, 2 * pi, length.out = 100)
    data.frame(
      level = r,
      x = r * cos(theta),
      y = r * sin(theta)
    )
  }))

  # ---------------------------------------------------------------------------
  # Default colors
  # ---------------------------------------------------------------------------
  if (is.null(colors)) {
    colors <- setNames(
      scales::hue_pal()(n_clusters),
      levels(plot_data$cluster)
    )
  }

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot() +
    # Grid circles
    geom_path(data = grid_data, aes(x = x, y = y, group = level),
              color = "gray80", linewidth = 0.3) +
    # Axis lines
    geom_segment(data = axis_data, aes(x = xend, y = yend, xend = x, yend = y),
                 color = "gray60", linewidth = 0.5) +
    # Cluster polygons (filled)
    geom_polygon(data = plot_data,
                 aes(x = x, y = y, fill = cluster, group = cluster),
                 alpha = alpha) +
    # Cluster polygon borders
    geom_path(data = plot_data,
              aes(x = x, y = y, color = cluster, group = cluster),
              linewidth = line_width) +
    # Points at vertices
    geom_point(data = plot_data[!duplicated(paste(plot_data$cluster, plot_data$variable)), ],
               aes(x = x, y = y, color = cluster),
               size = 2) +
    # Axis labels
    geom_text(data = label_data, aes(x = x, y = y, label = variable),
              size = 3, hjust = 0.5, vjust = 0.5) +
    # Styling
    scale_fill_manual(values = colors, name = "Cluster") +
    scale_color_manual(values = colors, name = "Cluster") +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  # Add grid level labels
  if (normalize) {
    p <- p + annotate("text", x = 0.02, y = grid_levels + 0.05,
                      label = grid_levels, size = 2.5, color = "gray50")
  }

  return(p)
}


#' Plot Sankey Diagram for Cluster Transitions Across k Values
#'
#' @description Visualize how samples flow between clusters as k increases.
#' This helps understand cluster stability and identify natural substructure.
#'
#' @param k_eval An artemis_k_evaluation object from ARTEMIS_evaluate_k_range().
#' @param colors Named character vector or NULL. Custom colors for clusters.
#'   If NULL (default), uses a generated palette.
#' @param node_width Numeric. Width of cluster nodes. Default = 30.
#' @param font_size Numeric. Font size for labels. Default = 12.
#' @param title Character. Plot title. Default = "Cluster Transitions Across k".
#' @param save_html Character. If provided, saves plot as standalone HTML file.
#'   Default = NULL (no save).
#'
#' @return A plotly/networkD3 Sankey diagram object
#'
#' @details
#' The Sankey diagram shows:
#' \itemize{
#'   \item Columns: different k values (left to right)
#'   \item Nodes: clusters at each k value
#'   \item Flows: how samples move from clusters at k to clusters at k+1
#'   \item Width: number of samples in each flow
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item Clean splits (one cluster splits into two) suggest natural substructure
#'   \item Messy redistributions suggest artificial divisions
#'   \item Stable cores (samples staying together) indicate robust clusters
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' k_eval <- ARTEMIS_evaluate_k_range(my_data, k_range = 2:5)
#' AETHER_plot_cluster_sankey(k_eval)
#'
#' # Save as HTML for sharing
#' AETHER_plot_cluster_sankey(k_eval, save_html = "cluster_transitions.html")
#'
#' }
AETHER_plot_cluster_sankey <- function(k_eval,
                                        colors = NULL,
                                        node_width = 30,
                                        font_size = 12,
                                        title = "Cluster Transitions Across k",
                                        save_html = NULL) {

  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  if (!inherits(k_eval, "artemis_k_evaluation")) {
    stop("k_eval must be an artemis_k_evaluation object from ARTEMIS_evaluate_k_range()")
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }

  transitions <- k_eval$transitions

  if (is.null(transitions) || nrow(transitions) == 0) {
    stop("No transition data found. Need at least 2 k values to show transitions.")
  }

  # ---------------------------------------------------------------------------
  # Build node list
  # ---------------------------------------------------------------------------
  # Get all unique nodes
  all_nodes <- unique(c(transitions$source, transitions$target))

  # Sort nodes by k value and cluster number for consistent ordering
  node_info <- data.frame(
    node = all_nodes,
    stringsAsFactors = FALSE
  )

  # Extract k value and cluster number
  node_info$k <- as.integer(gsub("k([0-9]+)_C.*", "\\1", node_info$node))
  node_info$cluster <- as.integer(gsub("k[0-9]+_C([0-9]+)", "\\1", node_info$node))

  # Sort by k, then cluster
  node_info <- node_info[order(node_info$k, node_info$cluster), ]

  # Create node index (0-based for plotly)
  node_info$idx <- seq_len(nrow(node_info)) - 1

  # Create lookup
  node_lookup <- setNames(node_info$idx, node_info$node)

  # ---------------------------------------------------------------------------
  # Map transitions to node indices
  # ---------------------------------------------------------------------------
  links <- data.frame(
    source = node_lookup[transitions$source],
    target = node_lookup[transitions$target],
    value = transitions$value,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Generate colors
  # ---------------------------------------------------------------------------
  if (is.null(colors)) {
    # Generate colors - one per cluster across all k
    unique_k <- unique(node_info$k)
    max_clusters <- max(node_info$cluster)

    # Use a palette with enough colors
    base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                     "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5")

    if (max_clusters > length(base_colors)) {
      base_colors <- grDevices::colorRampPalette(base_colors)(max_clusters)
    }

    # Assign same color to same cluster number across k values
    node_info$color <- base_colors[node_info$cluster]
  } else {
    # Use provided colors
    node_info$color <- colors[node_info$node]
    node_info$color[is.na(node_info$color)] <- "#888888"
  }

  # Create link colors (lighter version of source node color)
  link_colors <- sapply(links$source + 1, function(i) {
    col <- node_info$color[i]
    # Make semi-transparent
    rgb_col <- grDevices::col2rgb(col)
    grDevices::rgb(rgb_col[1], rgb_col[2], rgb_col[3], alpha = 100, maxColorValue = 255)
  })

  # ---------------------------------------------------------------------------
  # Create labels with cluster sizes
  # ---------------------------------------------------------------------------
  # Calculate total samples flowing through each node
  node_sizes <- sapply(node_info$node, function(n) {
    # Sum of outgoing flows (for non-terminal) or incoming flows (for terminal k)
    outgoing <- sum(transitions$value[transitions$source == n])
    incoming <- sum(transitions$value[transitions$target == n])
    max(outgoing, incoming)
  })

  node_labels <- paste0(
    "k=", node_info$k, " C", node_info$cluster,
    " (n=", node_sizes, ")"
  )

  # ---------------------------------------------------------------------------
  # Build Sankey diagram
  # ---------------------------------------------------------------------------
  fig <- plotly::plot_ly(
    type = "sankey",
    orientation = "h",

    node = list(
      label = node_labels,
      color = node_info$color,
      pad = 15,
      thickness = node_width,
      line = list(color = "black", width = 0.5)
    ),

    link = list(
      source = links$source,
      target = links$target,
      value = links$value,
      color = link_colors
    )
  )

  # Add layout
  fig <- fig |>
    plotly::layout(
      title = list(text = title, x = 0.5),
      font = list(size = font_size)
    )

  # ---------------------------------------------------------------------------
  # Save HTML if requested
  # ---------------------------------------------------------------------------
  if (!is.null(save_html)) {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
      warning("Package 'htmlwidgets' required to save HTML. ",
              "Install with: install.packages('htmlwidgets')")
    } else {
      htmlwidgets::saveWidget(
        fig,
        file = save_html,
        selfcontained = TRUE,
        title = title
      )
      message("Sankey diagram saved to: ", save_html)
    }
  }

  return(fig)
}


#' Plot k Selection Metrics
#'
#' @description Simple line/bar plot showing silhouette scores across k values.
#' Useful for determining optimal number of clusters.
#'
#' @param k_eval An artemis_k_evaluation object from ARTEMIS_evaluate_k_range().
#' @param highlight_optimal Logical. Highlight the optimal k. Default = TRUE.
#' @param show_sizes Logical. Show cluster sizes as secondary plot. Default = FALSE.
#' @param title Character. Plot title. Default = "Silhouette Score by k".
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' k_eval <- ARTEMIS_evaluate_k_range(my_data)
#' AETHER_plot_k_selection(k_eval)
#'
#' }
AETHER_plot_k_selection <- function(k_eval,
                                     highlight_optimal = TRUE,
                                     show_sizes = FALSE,
                                     title = "Silhouette Score by k") {

  if (!inherits(k_eval, "artemis_k_evaluation")) {
    stop("k_eval must be an artemis_k_evaluation object from ARTEMIS_evaluate_k_range()")
  }

  metrics <- k_eval$metrics
  metrics$k <- factor(metrics$k)

  # Main plot - silhouette scores
  p <- ggplot(metrics, aes(x = k, y = silhouette_avg, group = 1)) +
    geom_line(color = "#4575b4", linewidth = 1) +
    geom_point(size = 3, color = "#4575b4")

  # Highlight optimal
  if (highlight_optimal) {
    optimal_row <- metrics[metrics$k == k_eval$optimal_k, ]
    p <- p +
      geom_point(data = optimal_row, aes(x = k, y = silhouette_avg),
                 size = 5, color = "#d73027", shape = 18) +
      annotate("text",
               x = as.numeric(optimal_row$k),
               y = optimal_row$silhouette_avg,
               label = paste0("Optimal k=", k_eval$optimal_k),
               vjust = -1.5, hjust = 0.5, color = "#d73027", fontface = "bold")
  }

  # Add reference lines for interpretation
  p <- p +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    annotate("text", x = 0.6, y = 0.72, label = "Strong", size = 3, color = "gray40", hjust = 0) +
    annotate("text", x = 0.6, y = 0.52, label = "Reasonable", size = 3, color = "gray40", hjust = 0) +
    annotate("text", x = 0.6, y = 0.27, label = "Weak", size = 3, color = "gray40", hjust = 0)

  p <- p +
    labs(
      title = title,
      x = "Number of Clusters (k)",
      y = "Average Silhouette Width"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, max(metrics$silhouette_avg) * 1.15))

  return(p)
}
