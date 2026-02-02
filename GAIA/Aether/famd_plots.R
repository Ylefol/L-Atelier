library(ggplot2)
library(dplyr)

###############################################################################
########### FAMD Visualization Functions ###########
###############################################################################

#' Plot FAMD Individuals (Samples)
#'
#' @description Scatter plot of individuals/samples in FAMD space, colored by
#' a grouping variable. Primary visualization for exploring sample clustering.
#'
#' @param famd_result An artemis_famd object from ARTEMIS_famd()
#' @param dims Integer vector of length 2. Which dimensions to plot.
#'   Default = c(1, 2).
#' @param color_by Character, vector, or factor. Can be:
#'   \itemize{
#'     \item Column name from the original data
#'     \item Named vector (names must match FAMD row names - safest for external data)
#'     \item Unnamed vector (must match FAMD sample order exactly)
#'   }
#'   Should be categorical. If NULL, all points are same color.
#' @param point_size Numeric. Size of points. Default = 3.
#' @param alpha Numeric. Point transparency (0-1). Default = 0.8.
#' @param show_ellipse Logical. Draw confidence ellipses around groups.
#'   Default = TRUE if color_by is specified.
#' @param ellipse_level Numeric. Confidence level for ellipses. Default = 0.95.
#' @param ellipse_alpha Numeric. Ellipse fill transparency. Default = 0.1.
#' @param show_labels Logical. Label individual points. Default = FALSE.
#' @param label_size Numeric. Size of point labels. Default = 3.
#' @param colors Named character vector. Custom colors for groups.
#'   If NULL, uses default palette.
#' @param title Character. Plot title. Default = "FAMD - Individuals".
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' result <- ARTEMIS_famd(my_data)
#' AETHER_plot_famd_individuals(result, color_by = "treatment")
#' AETHER_plot_famd_individuals(result, color_by = "batch", dims = c(1, 3))
#'
AETHER_plot_famd_individuals <- function(famd_result,
                                          dims = c(1, 2),
                                          color_by = NULL,
                                          point_size = 3,
                                          alpha = 0.8,
                                          show_ellipse = NULL,
                                          ellipse_level = 0.95,
                                          ellipse_alpha = 0.1,
                                          show_labels = FALSE,
                                          label_size = 3,
                                          colors = NULL,
                                          title = "FAMD - Individuals") {

 # Validate input
  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  # Check dimensions
  max_dim <- ncol(famd_result$ind$coord)
  if (any(dims > max_dim)) {
    stop("Requested dimensions exceed available (max = ", max_dim, ")")
  }

  # Extract coordinates
  plot_data <- famd_result$ind$coord[, dims, drop = FALSE]
  colnames(plot_data) <- c("Dim1", "Dim2")
  plot_data$sample <- rownames(famd_result$ind$coord)

  # Add color variable if specified
  if (!is.null(color_by)) {
    n_samples <- nrow(famd_result$ind$coord)
    famd_rownames <- rownames(famd_result$ind$coord)

    if (is.character(color_by) && length(color_by) == 1) {
      # Case 1: Column name from data
      if (!color_by %in% colnames(famd_result$data_used)) {
        stop("color_by variable '", color_by, "' not found in data. ",
             "Available: ", paste(colnames(famd_result$data_used), collapse = ", "))
      }
      plot_data$group <- famd_result$data_used[[color_by]]

    } else if (length(color_by) == n_samples) {
      # Case 2: Vector of values
      if (!is.null(names(color_by))) {
        # Named vector - match by names (safest)
        if (!all(famd_rownames %in% names(color_by))) {
          missing <- setdiff(famd_rownames, names(color_by))
          stop("Named color_by vector is missing names for samples: ",
               paste(head(missing, 5), collapse = ", "),
               if (length(missing) > 5) paste0(" ... and ", length(missing) - 5, " more"))
        }
        plot_data$group <- color_by[famd_rownames]
      } else {
        # Unnamed vector - trust user on order, but warn
        warning("Using unnamed vector for color_by. Ensure it matches FAMD sample order exactly.")
        plot_data$group <- color_by
      }

    } else {
      stop("color_by must be a column name, or a vector of length ", n_samples,
           " (got length ", length(color_by), ")")
    }

    # Default show_ellipse to TRUE if color_by is specified
    if (is.null(show_ellipse)) show_ellipse <- TRUE
  } else {
    plot_data$group <- "All"
    if (is.null(show_ellipse)) show_ellipse <- FALSE
  }

  # Get variance explained for axis labels
  var_exp <- famd_result$eigenvalues$variance_percent[dims]
  x_lab <- paste0("Dim ", dims[1], " (", round(var_exp[1], 1), "%)")
  y_lab <- paste0("Dim ", dims[2], " (", round(var_exp[2], 1), "%)")

  # Build plot
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2))

  # Add ellipses first (behind points)
  if (show_ellipse && !is.null(color_by)) {
    p <- p + stat_ellipse(
      aes(color = group, fill = group),
      geom = "polygon",
      level = ellipse_level,
      alpha = ellipse_alpha,
      show.legend = FALSE
    )
  }

  # Add points
  if (!is.null(color_by)) {
    p <- p + geom_point(aes(color = group), size = point_size, alpha = alpha)
  } else {
    p <- p + geom_point(size = point_size, alpha = alpha, color = "#4575b4")
  }

  # Add labels if requested
  if (show_labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = sample),
        size = label_size,
        max.overlaps = 20
      )
    } else {
      p <- p + geom_text(aes(label = sample), size = label_size, vjust = -0.5)
    }
  }

  # Apply custom colors
  if (!is.null(colors) && !is.null(color_by)) {
    p <- p + scale_color_manual(values = colors) +
      scale_fill_manual(values = colors)
  }

  # Theme and labels
  p <- p +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = color_by,
      fill = color_by
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    coord_fixed()  # Equal scaling for both axes

  # Add reference lines at origin
  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3)

  return(p)
}


#' Plot FAMD Variables
#'
#' @description Visualize variable contributions to FAMD dimensions. Shows
#' quantitative variables as arrows and qualitative variable categories as points.
#'
#' @param famd_result An artemis_famd object from ARTEMIS_famd()
#' @param dims Integer vector of length 2. Which dimensions to plot.
#'   Default = c(1, 2).
#' @param var_type Character. Which variables to show:
#'   "all" (default), "quanti", or "quali".
#' @param aggregate_quali Logical. If TRUE, aggregate qualitative variables to
#'   variable level (one point per variable) instead of showing each category
#'   separately. Uses contribution-weighted centroid for positioning.
#'   Reduces clutter when many categories exist. Default = FALSE.
#' @param color_by Character. How to color variables:
#'   "type" (quanti vs quali), "contrib" (contribution strength).
#'   Default = "contrib".
#' @param show_labels Logical. Label variables. Default = TRUE.
#' @param label_size Numeric. Size of labels. Default = 3.
#' @param arrow_size Numeric. Size of arrow heads for quanti vars. Default = 0.2.
#' @param contrib_threshold Numeric. Only show variables with contribution above
#'   this threshold (as percentage). Default = 0 (show all).
#' @param title Character. Plot title. Default = "FAMD - Variables".
#'
#' @return A ggplot object
#'
#' @export
#'
AETHER_plot_famd_variables <- function(famd_result,
                                        dims = c(1, 2),
                                        var_type = "all",
                                        aggregate_quali = FALSE,
                                        color_by = "contrib",
                                        show_labels = TRUE,
                                        label_size = 3,
                                        arrow_size = 0.2,
                                        contrib_threshold = 0,
                                        title = "FAMD - Variables") {

  # Validate input
  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  if (!var_type %in% c("all", "quanti", "quali")) {
    stop("var_type must be one of: 'all', 'quanti', 'quali'")
  }

  # Check dimensions
  max_dim <- ncol(famd_result$var$quanti$coord)
  if (any(dims > max_dim)) {
    stop("Requested dimensions exceed available (max = ", max_dim, ")")
  }

  # Build data for quantitative variables
  quanti_data <- NULL
  if (var_type %in% c("all", "quanti")) {
    quanti_coord <- famd_result$var$quanti$coord[, dims, drop = FALSE]
    quanti_contrib <- rowMeans(famd_result$var$quanti$contrib[, dims, drop = FALSE])

    quanti_data <- data.frame(
      variable = rownames(quanti_coord),
      x = quanti_coord[, 1],
      y = quanti_coord[, 2],
      contrib = quanti_contrib,
      type = "Quantitative",
      stringsAsFactors = FALSE
    )
  }

  # Build data for qualitative variables (category levels)
  quali_data <- NULL
  if (var_type %in% c("all", "quali")) {
    quali_coord <- famd_result$var$quali$coord[, dims, drop = FALSE]
    quali_contrib <- famd_result$var$quali$contrib[, dims, drop = FALSE]

    if (aggregate_quali) {
      # Aggregate to variable level using contribution-weighted centroid
      # Category names are typically "VarName_Category" - extract variable name
      category_names <- rownames(quali_coord)

      # Get the original qualitative variable names from the data
      quali_vars <- names(famd_result$data_used)[sapply(famd_result$data_used, function(x) {
        is.factor(x) || is.character(x)
      })]

      # Match each category to its parent variable
      # Try to match by prefix (VarName_Category pattern)
      parent_var <- sapply(category_names, function(cat) {
        # Check which variable this category belongs to
        for (v in quali_vars) {
          # Check if category starts with "varname_"
          if (grepl(paste0("^", v, "_"), cat) || grepl(paste0("^", v, "\\."), cat)) {
            return(v)
          }
        }
        # Fallback: try to extract prefix before underscore/dot
        parts <- strsplit(cat, "[_.]")[[1]]
        if (length(parts) > 1) {
          return(parts[1])
        }
        return(cat)  # Last resort: use category name as-is
      })

      # Compute contribution-weighted centroid per variable
      unique_vars <- unique(parent_var)
      agg_data <- lapply(unique_vars, function(v) {
        idx <- which(parent_var == v)
        if (length(idx) == 1) {
          # Single category - use directly
          # Note: as.numeric() needed because data.frame row isn't a numeric vector
          return(data.frame(
            variable = v,
            x = quali_coord[idx, 1],
            y = quali_coord[idx, 2],
            contrib = mean(as.numeric(quali_contrib[idx, ])),
            stringsAsFactors = FALSE
          ))
        }

        # Multiple categories - weighted centroid
        cat_contrib <- rowMeans(quali_contrib[idx, , drop = FALSE])
        weights <- cat_contrib / sum(cat_contrib)

        # Handle case where all contributions are 0
        if (sum(cat_contrib) == 0) {
          weights <- rep(1 / length(idx), length(idx))
        }

        data.frame(
          variable = v,
          x = sum(quali_coord[idx, 1] * weights),
          y = sum(quali_coord[idx, 2] * weights),
          contrib = sum(cat_contrib),  # Total contribution for the variable
          stringsAsFactors = FALSE
        )
      })

      quali_data <- do.call(rbind, agg_data)
      quali_data$type <- "Qualitative"

    } else {
      # Original behavior: show each category
      quali_contrib_avg <- rowMeans(quali_contrib)

      quali_data <- data.frame(
        variable = rownames(quali_coord),
        x = quali_coord[, 1],
        y = quali_coord[, 2],
        contrib = quali_contrib_avg,
        type = "Qualitative",
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine data
  plot_data <- rbind(quanti_data, quali_data)

  # Apply contribution threshold
  if (contrib_threshold > 0) {
    plot_data <- plot_data[plot_data$contrib >= contrib_threshold, ]
    if (nrow(plot_data) == 0) {
      warning("No variables exceed contribution threshold of ", contrib_threshold, "%")
    }
  }

  # Get variance explained for axis labels
  var_exp <- famd_result$eigenvalues$variance_percent[dims]
  x_lab <- paste0("Dim ", dims[1], " (", round(var_exp[1], 1), "%)")
  y_lab <- paste0("Dim ", dims[2], " (", round(var_exp[2], 1), "%)")

  # Build plot
  p <- ggplot(plot_data, aes(x = x, y = y))

  # Add reference circle for quantitative variables (correlation circle)
  if (var_type %in% c("all", "quanti")) {
    circle_data <- data.frame(
      x = cos(seq(0, 2 * pi, length.out = 100)),
      y = sin(seq(0, 2 * pi, length.out = 100))
    )
    p <- p + geom_path(data = circle_data, aes(x = x, y = y),
                       color = "gray70", linetype = "dashed", inherit.aes = FALSE)
  }

  # Color setup
  if (color_by == "contrib") {
    # Quantitative as arrows
    if (!is.null(quanti_data) && nrow(quanti_data[quanti_data$variable %in% plot_data$variable, ]) > 0) {
      quanti_plot <- plot_data[plot_data$type == "Quantitative", ]
      p <- p + geom_segment(
        data = quanti_plot,
        aes(x = 0, y = 0, xend = x, yend = y, color = contrib),
        arrow = arrow(length = unit(arrow_size, "cm")),
        linewidth = 0.8
      )
    }

    # Qualitative as points
    if (!is.null(quali_data) && nrow(quali_data[quali_data$variable %in% plot_data$variable, ]) > 0) {
      quali_plot <- plot_data[plot_data$type == "Qualitative", ]
      p <- p + geom_point(
        data = quali_plot,
        aes(color = contrib),
        size = 3,
        shape = 17  # Triangle for quali
      )
    }

    p <- p + scale_color_gradient2(
      low = "#3288bd", mid = "#fee08b", high = "#d53e4f",
      midpoint = median(plot_data$contrib),
      name = "Contribution (%)"
    )

  } else if (color_by == "type") {
    type_colors <- c("Quantitative" = "#4575b4", "Qualitative" = "#d73027")

    # Quantitative as arrows
    if (!is.null(quanti_data) && nrow(quanti_data[quanti_data$variable %in% plot_data$variable, ]) > 0) {
      quanti_plot <- plot_data[plot_data$type == "Quantitative", ]
      p <- p + geom_segment(
        data = quanti_plot,
        aes(x = 0, y = 0, xend = x, yend = y, color = type),
        arrow = arrow(length = unit(arrow_size, "cm")),
        linewidth = 0.8
      )
    }

    # Qualitative as points
    if (!is.null(quali_data) && nrow(quali_data[quali_data$variable %in% plot_data$variable, ]) > 0) {
      quali_plot <- plot_data[plot_data$type == "Qualitative", ]
      p <- p + geom_point(
        data = quali_plot,
        aes(color = type),
        size = 3,
        shape = 17
      )
    }

    p <- p + scale_color_manual(values = type_colors, name = "Variable Type")
  }

  # Add labels
  if (show_labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = variable),
        size = label_size,
        max.overlaps = 30
      )
    } else {
      p <- p + geom_text(aes(label = variable), size = label_size, vjust = -0.5)
    }
  }

  # Theme and labels
  p <- p +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    coord_fixed() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3)

  return(p)
}


#' Plot FAMD Contribution Heatmap
#'
#' @description Heatmap showing variable contributions to each FAMD dimension.
#' Cleaner alternative to bar plots for understanding which variables drive
#' which dimensions.
#'
#' @param famd_result An artemis_famd object from ARTEMIS_famd()
#' @param n_dims Integer. Number of dimensions to show. Default = 5.
#' @param var_type Character. Which variables to show:
#'   "all" (default), "quanti", or "quali".
#' @param top_n Integer. Only show top N contributing variables per dimension.
#'   If NULL (default), show all variables.
#' @param title Character. Plot title. Default = "FAMD - Variable Contributions".
#'
#' @return A ggplot object
#'
#' @export
#'
AETHER_plot_famd_contrib <- function(famd_result,
                                      n_dims = 5,
                                      var_type = "all",
                                      top_n = NULL,
                                      title = "FAMD - Variable Contributions") {

  # Validate input
  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  # Limit dimensions
  max_dim <- ncol(famd_result$var$quanti$contrib)
  n_dims <- min(n_dims, max_dim)

  # Build contribution data
  contrib_list <- list()

  if (var_type %in% c("all", "quanti")) {
    quanti_contrib <- famd_result$var$quanti$contrib[, 1:n_dims, drop = FALSE]
    quanti_df <- as.data.frame(quanti_contrib)
    quanti_df$variable <- rownames(quanti_contrib)
    quanti_df$type <- "Quantitative"
    contrib_list$quanti <- quanti_df
  }

  if (var_type %in% c("all", "quali")) {
    # For qualitative, we have category-level contributions
    # Aggregate to variable level by taking max contribution per variable
    quali_contrib <- famd_result$var$quali$contrib[, 1:n_dims, drop = FALSE]
    quali_df <- as.data.frame(quali_contrib)
    quali_df$variable <- rownames(quali_contrib)
    quali_df$type <- "Qualitative"
    contrib_list$quali <- quali_df
  }

  contrib_data <- do.call(rbind, contrib_list)

  # Reshape to long format
  dim_cols <- paste0("Dim.", 1:n_dims)
  # Handle different column naming from FactoMineR
  actual_cols <- colnames(contrib_data)[1:n_dims]

  contrib_long <- reshape(
    contrib_data,
    direction = "long",
    varying = list(actual_cols),
    v.names = "contribution",
    timevar = "dimension",
    times = paste0("Dim ", 1:n_dims),
    idvar = c("variable", "type")
  )

  # Filter to top N if specified
  if (!is.null(top_n)) {
    # Get top contributors across all dimensions
    avg_contrib <- tapply(contrib_long$contribution, contrib_long$variable, mean)
    top_vars <- names(sort(avg_contrib, decreasing = TRUE))[1:min(top_n, length(avg_contrib))]
    contrib_long <- contrib_long[contrib_long$variable %in% top_vars, ]
  }

  # Order variables by total contribution
  var_order <- tapply(contrib_long$contribution, contrib_long$variable, sum)
  contrib_long$variable <- factor(contrib_long$variable,
                                   levels = names(sort(var_order, decreasing = FALSE)))

  # Build heatmap
  p <- ggplot(contrib_long, aes(x = dimension, y = variable, fill = contribution)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#f7f7f7", mid = "#fee08b", high = "#d53e4f",
      midpoint = median(contrib_long$contribution),
      name = "Contribution (%)"
    ) +
    geom_text(aes(label = round(contribution, 1)),
              size = 3, color = "black") +
    labs(
      title = title,
      x = "Dimension",
      y = "Variable"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  # Add type annotation if showing all
  if (var_type == "all") {
    p <- p + facet_grid(type ~ ., scales = "free_y", space = "free_y")
  }

  return(p)
}


#' Plot FAMD Scree Plot
#'
#' @description Bar plot showing variance explained by each FAMD dimension,
#' with cumulative variance line.
#'
#' @param famd_result An artemis_famd object from ARTEMIS_famd()
#' @param n_dims Integer. Number of dimensions to show. Default = 10.
#' @param show_cumulative Logical. Show cumulative variance line. Default = TRUE.
#' @param bar_fill Character. Fill color for bars. Default = "#4575b4".
#' @param line_color Character. Color for cumulative line. Default = "#d73027".
#' @param title Character. Plot title. Default = "FAMD - Scree Plot".
#'
#' @return A ggplot object
#'
#' @export
#'
AETHER_plot_famd_scree <- function(famd_result,
                                    n_dims = 10,
                                    show_cumulative = TRUE,
                                    bar_fill = "#4575b4",
                                    line_color = "#d73027",
                                    title = "FAMD - Scree Plot") {

  # Validate input
  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  # Limit dimensions
  max_dim <- nrow(famd_result$eigenvalues)
  n_dims <- min(n_dims, max_dim)

  # Prepare data
  plot_data <- famd_result$eigenvalues[1:n_dims, ]
  plot_data$dimension <- factor(paste0("Dim ", plot_data$dimension),
                                 levels = paste0("Dim ", 1:n_dims))

  # Build plot
  p <- ggplot(plot_data, aes(x = dimension, y = variance_percent)) +
    geom_col(fill = bar_fill, alpha = 0.8, width = 0.7) +
    geom_text(aes(label = paste0(round(variance_percent, 1), "%")),
              vjust = -0.5, size = 3)

  # Add cumulative line
  if (show_cumulative) {
    # Scale cumulative to fit on same axis
    max_var <- max(plot_data$variance_percent)
    scale_factor <- max_var / 100

    p <- p +
      geom_line(aes(y = cumulative_percent * scale_factor, group = 1),
                color = line_color, linewidth = 1) +
      geom_point(aes(y = cumulative_percent * scale_factor),
                 color = line_color, size = 2) +
      scale_y_continuous(
        name = "Variance Explained (%)",
        limits = c(0, max_var * 1.15),
        sec.axis = sec_axis(~ . / scale_factor,
                            name = "Cumulative Variance (%)")
      )
  } else {
    p <- p +
      scale_y_continuous(
        name = "Variance Explained (%)",
        limits = c(0, max(plot_data$variance_percent) * 1.15)
      )
  }

  p <- p +
    labs(
      title = title,
      x = "Dimension"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  return(p)
}


#' Interactive FAMD Variables Plot
#'
#' @description Interactive plotly version of the FAMD variables plot. Useful
#' for exploring high-cardinality categorical variables where static plots
#' become unreadable.
#'
#' @param famd_result An artemis_famd object from ARTEMIS_famd()
#' @param dims Integer vector of length 2. Which dimensions to plot.
#'   Default = c(1, 2).
#' @param var_type Character. Which variables to show:
#'   "all" (default), "quanti", or "quali".
#' @param point_size Numeric. Size of points. Default = 8.
#' @param show_arrows Logical. Show arrows for quantitative variables.
#'   Default = TRUE.
#' @param title Character. Plot title. Default = "FAMD - Variables (Interactive)".
#' @param save_html Character. If provided, saves plot as standalone HTML file
#'   at this path. Default = NULL (no save).
#'
#' @return A plotly object
#'
#' @details
#' This function creates an interactive scatter plot where you can:
#' \itemize{
#'   \item Hover over points to see variable/category names and contributions
#'   \item Zoom into regions of interest
#'   \item Pan around the plot
#'   \item Toggle variable types on/off via legend
#' }
#'
#' The HTML output is self-contained and can be shared with others who have
#' a modern web browser (no R required to view).
#'
#' @export
#'
#' @examples
#' result <- ARTEMIS_famd(my_data)
#' AETHER_plot_famd_variables_interactive(result)
#'
#' # Save as HTML for sharing
#' AETHER_plot_famd_variables_interactive(result, save_html = "famd_explore.html")
#'
AETHER_plot_famd_variables_interactive <- function(famd_result,
                                                    dims = c(1, 2),
                                                    var_type = "all",
                                                    point_size = 8,
                                                    show_arrows = TRUE,
                                                    title = "FAMD - Variables (Interactive)",
                                                    save_html = NULL) {

  # Check for plotly

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots. ",
         "Install with: install.packages('plotly')")
  }

  # Validate input
  if (!inherits(famd_result, "artemis_famd")) {
    stop("famd_result must be an artemis_famd object from ARTEMIS_famd()")
  }

  if (!var_type %in% c("all", "quanti", "quali")) {
    stop("var_type must be one of: 'all', 'quanti', 'quali'")
  }

  # Check dimensions
  max_dim <- ncol(famd_result$var$quanti$coord)
  if (any(dims > max_dim)) {
    stop("Requested dimensions exceed available (max = ", max_dim, ")")
  }

  # Get variance explained for axis labels
  var_exp <- famd_result$eigenvalues$variance_percent[dims]
  x_lab <- paste0("Dim ", dims[1], " (", round(var_exp[1], 1), "%)")
  y_lab <- paste0("Dim ", dims[2], " (", round(var_exp[2], 1), "%)")

  # ---------------------------------------------------------------------------
  # Build data for quantitative variables
  # ---------------------------------------------------------------------------
  quanti_data <- NULL
  if (var_type %in% c("all", "quanti")) {
    quanti_coord <- famd_result$var$quanti$coord[, dims, drop = FALSE]
    quanti_contrib <- rowMeans(famd_result$var$quanti$contrib[, dims, drop = FALSE])

    quanti_data <- data.frame(
      name = rownames(quanti_coord),
      x = quanti_coord[, 1],
      y = quanti_coord[, 2],
      contrib = quanti_contrib,
      type = "Quantitative",
      parent_var = rownames(quanti_coord),  # Same as name for quanti
      stringsAsFactors = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Build data for qualitative variables (categories)
  # ---------------------------------------------------------------------------
  quali_data <- NULL
  if (var_type %in% c("all", "quali")) {
    quali_coord <- famd_result$var$quali$coord[, dims, drop = FALSE]
    quali_contrib <- rowMeans(famd_result$var$quali$contrib[, dims, drop = FALSE])

    # Try to extract parent variable name from category name
    category_names <- rownames(quali_coord)

    # Get original variable names
    quali_vars <- names(famd_result$data_used)[sapply(famd_result$data_used, function(x) {
      is.factor(x) || is.character(x)
    })]

    # Map categories to parent variables
    parent_var <- sapply(category_names, function(cat) {
      for (v in quali_vars) {
        if (grepl(paste0("^", v, "_"), cat) || grepl(paste0("^", v, "\\."), cat)) {
          return(v)
        }
      }
      # Fallback
      parts <- strsplit(cat, "[_.]")[[1]]
      if (length(parts) > 1) return(parts[1])
      return("Unknown")
    })

    quali_data <- data.frame(
      name = category_names,
      x = quali_coord[, 1],
      y = quali_coord[, 2],
      contrib = quali_contrib,
      type = "Qualitative",
      parent_var = parent_var,
      stringsAsFactors = FALSE
    )
  }

  # Combine data (for reference, though we plot separately)
  plot_data <- rbind(quanti_data, quali_data)

  # Create hover text for each dataset
  if (!is.null(quanti_data) && nrow(quanti_data) > 0) {
    quanti_data$hover_text <- paste0(
      "<b>", quanti_data$name, "</b><br>",
      "Type: Quantitative<br>",
      "Contribution: ", round(quanti_data$contrib, 2), "%<br>",
      "Dim ", dims[1], ": ", round(quanti_data$x, 3), "<br>",
      "Dim ", dims[2], ": ", round(quanti_data$y, 3)
    )
  }

  if (!is.null(quali_data) && nrow(quali_data) > 0) {
    quali_data$hover_text <- paste0(
      "<b>", quali_data$name, "</b><br>",
      "Variable: ", quali_data$parent_var, "<br>",
      "Type: Qualitative<br>",
      "Contribution: ", round(quali_data$contrib, 2), "%<br>",
      "Dim ", dims[1], ": ", round(quali_data$x, 3), "<br>",
      "Dim ", dims[2], ": ", round(quali_data$y, 3)
    )
  }

  # ---------------------------------------------------------------------------
  # Build plotly figure
  # ---------------------------------------------------------------------------

  # Start with empty plot
  fig <- plotly::plot_ly()

  # Add reference circle
  circle_theta <- seq(0, 2 * pi, length.out = 100)
  fig <- fig %>%
    plotly::add_trace(
      x = cos(circle_theta),
      y = sin(circle_theta),
      type = "scatter",
      mode = "lines",
      line = list(color = "gray", dash = "dash", width = 1),
      hoverinfo = "none",
      showlegend = FALSE
    )

  # Add quantitative variables
  if (!is.null(quanti_data) && nrow(quanti_data) > 0) {
    # Add arrows if requested
    if (show_arrows) {
      for (i in seq_len(nrow(quanti_data))) {
        fig <- fig %>%
          plotly::add_trace(
            x = c(0, quanti_data$x[i]),
            y = c(0, quanti_data$y[i]),
            type = "scatter",
            mode = "lines",
            line = list(color = "#4575b4", width = 1.5),
            hoverinfo = "none",
            showlegend = FALSE
          )
      }
    }

    # Add points at arrow tips
    fig <- fig %>%
      plotly::add_trace(
        x = quanti_data$x,
        y = quanti_data$y,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = point_size,
          color = "#4575b4",
          symbol = "circle"
        ),
        text = quanti_data$hover_text,
        hoverinfo = "text",
        name = "Quantitative"
      )
  }

  # Add qualitative variables (categories) - one trace per parent variable for legend control
  if (!is.null(quali_data) && nrow(quali_data) > 0) {
    # Color by parent variable for easier identification
    unique_vars <- unique(quali_data$parent_var)
    n_vars <- length(unique_vars)

    # Generate colors for each parent variable
    if (n_vars <= 12) {
      var_colors <- setNames(
        RColorBrewer::brewer.pal(max(3, n_vars), "Set3")[1:n_vars],
        unique_vars
      )
    } else {
      # For many variables, use a gradient
      var_colors <- setNames(
        grDevices::colorRampPalette(c("#e41a1c", "#377eb8", "#4daf4a",
                                       "#984ea3", "#ff7f00", "#ffff33"))(n_vars),
        unique_vars
      )
    }

    # Add each parent variable as a separate trace for legend toggling
    for (var_name in unique_vars) {
      var_subset <- quali_data[quali_data$parent_var == var_name, ]

      fig <- fig %>%
        plotly::add_trace(
          x = var_subset$x,
          y = var_subset$y,
          type = "scatter",
          mode = "markers",
          marker = list(
            size = point_size,
            color = var_colors[var_name],
            symbol = "triangle-up"
          ),
          text = var_subset$hover_text,
          hoverinfo = "text",
          name = var_name,
          legendgroup = "quali",  # Group all quali together
          legendgrouptitle = list(text = "Qualitative")
        )
    }
  }

  # Add reference lines
  fig <- fig %>%
    plotly::layout(
      title = list(text = title, x = 0.5),
      xaxis = list(
        title = x_lab,
        zeroline = TRUE,
        zerolinecolor = "gray",
        zerolinewidth = 1
      ),
      yaxis = list(
        title = y_lab,
        zeroline = TRUE,
        zerolinecolor = "gray",
        zerolinewidth = 1,
        scaleanchor = "x",
        scaleratio = 1
      ),
      hovermode = "closest",
      legend = list(
        orientation = "v",
        yanchor = "top",
        y = 1,
        xanchor = "left",
        x = 1.02,
        groupclick = "toggleitem"  # Click toggles individual items, not whole group
      )
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
      message("Interactive plot saved to: ", save_html)
    }
  }

  return(fig)
}
