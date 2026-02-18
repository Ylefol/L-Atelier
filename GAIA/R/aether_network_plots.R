# GAIA/Aether/network_plots.R
# Network visualization for PPI and other interaction networks
#
# Provides static (ggraph) and interactive (visNetwork) network plots.
# Works with igraph objects from APOLLO_get_ppi() or any igraph graph.



#' Static PPI Network Plot
#'
#' @description Creates a publication-quality network plot using ggraph.
#' Nodes can be colored by cluster membership, fold change, module, or
#' any named vector. Node size reflects degree or custom metric.
#'
#' @param graph An igraph object (e.g., from APOLLO_get_ppi()).
#' @param color_by Node coloring. One of:
#'   \itemize{
#'     \item Named vector: gene -> value (numeric for gradient, character/factor
#'       for discrete). E.g., cluster assignments or log2FC values.
#'     \item "degree": Color by node degree (connectivity).
#'     \item "betweenness": Color by betweenness centrality.
#'     \item Single color string: uniform color (e.g., "steelblue").
#'     \item NULL: default gray.
#'   }
#' @param size_by Node sizing. One of:
#'   \itemize{
#'     \item "degree" (default): Size by connectivity.
#'     \item "betweenness": Size by betweenness centrality.
#'     \item Named numeric vector: gene -> value.
#'     \item Single number: uniform size.
#'   }
#' @param layout Character. igraph layout algorithm. Default: "fr"
#'   (Fruchterman-Reingold). Other options: "kk" (Kamada-Kawai), "circle",
#'   "tree", "grid", "stress", "dh", "gem", "graphopt".
#' @param show_labels Logical or character. TRUE shows all labels, FALSE hides
#'   all, "hubs" shows labels only for top hub nodes. Default: "hubs".
#' @param hub_n Integer. Number of top hubs to label when show_labels="hubs".
#'   Default: 15.
#' @param colors Color specification for nodes. For discrete: named vector of
#'   colors. For continuous: vector of 2-3 colors for gradient. Default: NULL
#'   (auto-selected based on color_by type).
#' @param edge_alpha Numeric. Edge transparency (0-1). Default: 0.3.
#' @param edge_color Character. Edge color. Default: "gray70".
#' @param size_range Numeric vector of length 2. Min and max node sizes.
#'   Default: c(2, 10).
#' @param title Character. Plot title. Default: NULL (auto-generated).
#' @param seed Integer. Random seed for reproducible layouts. Default: 42.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Basic network colored by degree
#' p <- AETHER_plot_ppi_network(graph, color_by = "degree")
#'
#' # Color by cluster membership
#' clusters <- c(TP53 = "C1", BRCA1 = "C1", EGFR = "C2", MYC = "C2")
#' p <- AETHER_plot_ppi_network(graph, color_by = clusters)
#'
#' # Color by fold change
#' fc <- c(TP53 = 2.1, BRCA1 = -1.5, EGFR = 0.3, MYC = 3.2)
#' p <- AETHER_plot_ppi_network(graph, color_by = fc, colors = c("blue", "white", "red"))
#'
#' }
#' @export
AETHER_plot_ppi_network <- function(graph,
                                     color_by = "degree",
                                     size_by = "degree",
                                     layout = "fr",
                                     show_labels = "hubs",
                                     hub_n = 15,
                                     colors = NULL,
                                     edge_alpha = 0.3,
                                     edge_color = "gray70",
                                     size_range = c(2, 10),
                                     title = NULL,
                                     seed = 42) {

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Package 'ggraph' is required. Install with:\n",
         "  install.packages('ggraph')", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }

  node_names <- igraph::V(graph)$name
  n_nodes <- igraph::vcount(graph)
  n_edges <- igraph::ecount(graph)

  # --- Size variable ---
  if (is.character(size_by) && length(size_by) == 1) {
    if (size_by == "degree") {
      igraph::V(graph)$size_var <- igraph::degree(graph)
    } else if (size_by == "betweenness") {
      igraph::V(graph)$size_var <- igraph::betweenness(graph)
    } else {
      igraph::V(graph)$size_var <- rep(as.numeric(size_by), n_nodes)
    }
  } else if (is.numeric(size_by) && length(size_by) == 1) {
    igraph::V(graph)$size_var <- rep(size_by, n_nodes)
    size_range <- c(size_by, size_by)
  } else if (is.numeric(size_by) && !is.null(names(size_by))) {
    igraph::V(graph)$size_var <- size_by[node_names]
    igraph::V(graph)$size_var[is.na(igraph::V(graph)$size_var)] <- 0
  } else {
    igraph::V(graph)$size_var <- igraph::degree(graph)
  }

  # --- Color variable ---
  is_continuous <- FALSE
  is_discrete <- FALSE

  if (is.null(color_by)) {
    igraph::V(graph)$color_var <- "node"
    is_discrete <- TRUE
  } else if (is.character(color_by) && length(color_by) == 1 && !color_by %in% node_names) {
    if (color_by == "degree") {
      igraph::V(graph)$color_var <- igraph::degree(graph)
      is_continuous <- TRUE
    } else if (color_by == "betweenness") {
      igraph::V(graph)$color_var <- igraph::betweenness(graph)
      is_continuous <- TRUE
    } else {
      # Single color string — handled in aes
      igraph::V(graph)$color_var <- color_by
      is_discrete <- TRUE
    }
  } else if (!is.null(names(color_by))) {
    # Named vector
    mapped <- color_by[node_names]
    if (is.numeric(color_by)) {
      igraph::V(graph)$color_var <- as.numeric(mapped)
      igraph::V(graph)$color_var[is.na(igraph::V(graph)$color_var)] <- 0
      is_continuous <- TRUE
    } else {
      igraph::V(graph)$color_var <- as.character(mapped)
      igraph::V(graph)$color_var[is.na(igraph::V(graph)$color_var)] <- "other"
      is_discrete <- TRUE
    }
  } else {
    igraph::V(graph)$color_var <- "node"
    is_discrete <- TRUE
  }

  # --- Labels ---
  if (is.logical(show_labels) && show_labels) {
    igraph::V(graph)$label_var <- node_names
  } else if (is.character(show_labels) && show_labels == "hubs") {
    degrees <- igraph::degree(graph)
    top_idx <- order(degrees, decreasing = TRUE)[seq_len(min(hub_n, n_nodes))]
    # Only label nodes with at least 1 connection
    top_idx <- top_idx[degrees[top_idx] > 0]
    labels <- rep(NA_character_, n_nodes)
    labels[top_idx] <- node_names[top_idx]
    igraph::V(graph)$label_var <- labels
  } else {
    igraph::V(graph)$label_var <- rep(NA_character_, n_nodes)
  }

  # --- Title ---
  if (is.null(title)) {
    title <- paste0("PPI Network (", n_nodes, " nodes, ", n_edges, " edges)")
  }

  # --- Build plot ---
  set.seed(seed)

  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(alpha = edge_alpha, color = edge_color) +
    ggraph::geom_node_point(aes(size = size_var, color = color_var)) +
    scale_size_continuous(range = size_range, name = "Size") +
    labs(title = title) +
    ggraph::theme_graph(base_family = "") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  # --- Color scale ---
  if (is_continuous) {
    if (is.null(colors)) colors <- c("#2166ac", "#f7f7f7", "#b2182b")
    if (length(colors) == 2) {
      p <- p + scale_color_gradient(low = colors[1], high = colors[2], name = "Value")
    } else {
      mid_val <- median(igraph::V(graph)$color_var, na.rm = TRUE)
      p <- p + scale_color_gradient2(low = colors[1], mid = colors[2],
                                      high = colors[3], midpoint = mid_val,
                                      name = "Value")
    }
  } else if (is_discrete) {
    unique_vals <- unique(igraph::V(graph)$color_var)
    if (length(unique_vals) == 1 && !unique_vals[1] %in% c("node", "other")) {
      # Single color
      p <- p + scale_color_manual(values = setNames(unique_vals, unique_vals), guide = "none")
    } else if (!is.null(colors) && is.character(colors) && !is.null(names(colors))) {
      p <- p + scale_color_manual(values = colors, name = "Group")
    } else {
      p <- p + scale_color_discrete(name = "Group")
    }
  }

  # --- Labels ---
  if (any(!is.na(igraph::V(graph)$label_var))) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggraph::geom_node_text(aes(label = label_var),
                                        repel = TRUE, size = 3,
                                        max.overlaps = 20)
    } else {
      p <- p + ggraph::geom_node_text(aes(label = label_var),
                                        size = 3, vjust = -1)
    }
  }

  return(p)
}


#' Interactive PPI Network Plot
#'
#' @description Creates an interactive network visualization using visNetwork.
#' Nodes can be dragged, zoomed, and hovered for information. Useful for
#' exploring large networks.
#'
#' @param graph An igraph object (e.g., from APOLLO_get_ppi()).
#' @param color_by Node coloring. Same options as AETHER_plot_ppi_network():
#'   named vector, "degree", "betweenness", single color, or NULL.
#' @param size_by Node sizing. Same options as AETHER_plot_ppi_network().
#' @param colors Color specification. Named vector for discrete, or gradient
#'   vector for continuous. Default: NULL (auto-selected).
#' @param size_range Numeric vector of length 2. Min and max node sizes.
#'   Default: c(10, 40).
#' @param show_labels Logical. Show node labels. Default: TRUE.
#' @param physics Logical. Enable physics simulation for layout. Default: TRUE.
#'   Set FALSE for static layout (faster for large networks).
#' @param layout Character. Layout algorithm when physics=FALSE.
#'   Options: "layout_with_fr", "layout_nicely", "layout_in_circle".
#'   Default: "layout_with_fr".
#' @param title Character. Plot title. Default: NULL.
#' @param save_html Character. Path to save as standalone HTML file. Default: NULL.
#' @param seed Integer. Random seed for layout. Default: 42.
#'
#' @return A visNetwork object (rendered in viewer or saved to HTML).
#'
#' @examples
#' \dontrun{
#' # Interactive exploration
#' AETHER_plot_ppi_network_interactive(graph, color_by = "degree")
#'
#' # Save to HTML
#' AETHER_plot_ppi_network_interactive(graph, save_html = "network.html")
#'
#' }
#' @export
AETHER_plot_ppi_network_interactive <- function(graph,
                                                 color_by = "degree",
                                                 size_by = "degree",
                                                 colors = NULL,
                                                 size_range = c(10, 40),
                                                 show_labels = TRUE,
                                                 physics = TRUE,
                                                 layout = "layout_with_fr",
                                                 title = NULL,
                                                 save_html = NULL,
                                                 seed = 42) {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("Package 'visNetwork' is required. Install with:\n",
         "  install.packages('visNetwork')", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }

  node_names <- igraph::V(graph)$name
  n_nodes <- igraph::vcount(graph)

  # --- Size ---
  if (is.character(size_by) && length(size_by) == 1) {
    if (size_by == "degree") {
      size_vals <- igraph::degree(graph)
    } else if (size_by == "betweenness") {
      size_vals <- igraph::betweenness(graph)
    } else {
      size_vals <- rep(mean(size_range), n_nodes)
    }
  } else if (is.numeric(size_by) && length(size_by) == 1) {
    size_vals <- rep(size_by, n_nodes)
  } else if (is.numeric(size_by) && !is.null(names(size_by))) {
    size_vals <- size_by[node_names]
    size_vals[is.na(size_vals)] <- 0
  } else {
    size_vals <- igraph::degree(graph)
  }

  # Rescale to size_range
  if (max(size_vals, na.rm = TRUE) > min(size_vals, na.rm = TRUE)) {
    size_scaled <- size_range[1] + (size_vals - min(size_vals, na.rm = TRUE)) /
      (max(size_vals, na.rm = TRUE) - min(size_vals, na.rm = TRUE)) *
      (size_range[2] - size_range[1])
  } else {
    size_scaled <- rep(mean(size_range), n_nodes)
  }

  # --- Color ---
  if (is.null(color_by)) {
    node_colors <- rep("#97C2FC", n_nodes)
  } else if (is.character(color_by) && length(color_by) == 1 && !color_by %in% node_names) {
    if (color_by == "degree") {
      vals <- igraph::degree(graph)
      node_colors <- .vis_continuous_colors(vals, colors)
    } else if (color_by == "betweenness") {
      vals <- igraph::betweenness(graph)
      node_colors <- .vis_continuous_colors(vals, colors)
    } else {
      node_colors <- rep(color_by, n_nodes)
    }
  } else if (!is.null(names(color_by))) {
    if (is.numeric(color_by)) {
      vals <- color_by[node_names]
      vals[is.na(vals)] <- 0
      node_colors <- .vis_continuous_colors(vals, colors)
    } else {
      groups <- as.character(color_by[node_names])
      groups[is.na(groups)] <- "other"

      if (!is.null(colors) && !is.null(names(colors))) {
        color_map <- colors
      } else {
        unique_groups <- unique(groups)
        pal <- grDevices::hcl.colors(length(unique_groups), "Set2")
        color_map <- setNames(pal, unique_groups)
      }
      node_colors <- unname(color_map[groups])
      node_colors[is.na(node_colors)] <- "#CCCCCC"
    }
  } else {
    node_colors <- rep("#97C2FC", n_nodes)
  }

  # --- Build nodes data.frame ---
  nodes <- data.frame(
    id    = node_names,
    label = if (show_labels) node_names else rep("", n_nodes),
    size  = size_scaled,
    color = node_colors,
    title = paste0("<b>", node_names, "</b><br>",
                   "Degree: ", igraph::degree(graph), "<br>",
                   "Betweenness: ", round(igraph::betweenness(graph), 1)),
    stringsAsFactors = FALSE
  )

  # --- Build edges data.frame ---
  edge_list <- igraph::as_data_frame(graph, what = "edges")
  edges <- data.frame(
    from  = edge_list$from,
    to    = edge_list$to,
    color = "gray80",
    stringsAsFactors = FALSE
  )

  if ("n_resources" %in% colnames(edge_list)) {
    edges$title <- paste0("Resources: ", edge_list$n_resources)
    # Scale edge width by resources
    edges$width <- pmin(edge_list$n_resources, 5)
  }

  # --- Title ---
  if (is.null(title)) {
    title <- paste0("PPI Network (", n_nodes, " nodes, ", nrow(edges), " edges)")
  }

  # --- Build visNetwork ---
  set.seed(seed)

  vis <- visNetwork::visNetwork(nodes, edges, main = title) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE
    ) |>
    visNetwork::visInteraction(
      navigationButtons = TRUE,
      zoomView = TRUE
    )

  if (physics) {
    vis <- vis |>
      visNetwork::visPhysics(
        solver = "forceAtlas2Based",
        stabilization = list(iterations = 200)
      )
  } else {
    vis <- vis |>
      visNetwork::visIgraphLayout(layout = layout, randomSeed = seed) |>
      visNetwork::visPhysics(enabled = FALSE)
  }

  # --- Save HTML ---
  if (!is.null(save_html)) {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
      warning("htmlwidgets package needed for HTML export.")
    } else {
      htmlwidgets::saveWidget(vis, file = save_html, selfcontained = TRUE)
      cat("Saved interactive network to:", save_html, "\n")
    }
  }

  return(vis)
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Map continuous values to a color gradient for visNetwork
#' @keywords internal
.vis_continuous_colors <- function(vals, colors = NULL) {
  if (is.null(colors)) colors <- c("#2166ac", "#f7f7f7", "#b2182b")

  rng <- range(vals, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(colors[ceiling(length(colors) / 2)], length(vals)))
  }

  # Normalize to 0-1
  scaled <- (vals - rng[1]) / (rng[2] - rng[1])
  color_ramp <- grDevices::colorRamp(colors)
  rgb_vals <- color_ramp(scaled)
  grDevices::rgb(rgb_vals[, 1], rgb_vals[, 2], rgb_vals[, 3], maxColorValue = 255)
}
