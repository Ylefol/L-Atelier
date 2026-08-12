# GAIA/Aether/cibersort_plots.R
# CIBERSORT deconvolution visualization
#
# Plotting functions for visualizing immune cell type proportions
# from CIBERSORT deconvolution results (artemis_cibersort objects).



# ==============================================================================
# HELPERS
# ==============================================================================

#' Extract proportions matrix from various input types
#' @return Matrix (samples x cell types)
#' @keywords internal
.cibersort_extract_proportions <- function(result) {
  if (inherits(result, "artemis_cibersort")) {
    mat <- result$proportions
  } else if (is.matrix(result) || is.data.frame(result)) {
    mat <- as.matrix(result)
  } else {
    stop("'result' must be an artemis_cibersort object, matrix, or data.frame",
         call. = FALSE)
  }

  # Remove cell types with all-zero proportions
  nonzero <- colSums(mat, na.rm = TRUE) > 0
  mat <- mat[, nonzero, drop = FALSE]

  return(mat)
}


#' Default qualitative palette for cell types (up to 22 for LM22)
#' @keywords internal
.cibersort_cell_palette <- function(n) {
  # Extended qualitative palette for many categories
  base <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#637939"
  )
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  colorRampPalette(base)(n)
}


# ==============================================================================
# STACKED BAR PLOT
# ==============================================================================

#' Plot cell type proportions as stacked bar chart
#'
#' Creates a stacked bar plot showing the estimated cell type proportions
#' for each sample. Useful for comparing immune composition across samples.
#'
#' @param result An artemis_cibersort object or a proportions matrix
#'   (samples x cell types).
#' @param group_by Optional factor or character vector of group labels per
#'   sample. Used to reorder samples by group and add faceting. Length must
#'   match number of samples. Default: NULL.
#' @param colors Named character vector of colors for cell types. Names should
#'   match cell type column names. If NULL, uses a default qualitative palette.
#' @param title Character. Plot title. Default: "Cell Type Proportions".
#' @param text_size Numeric. Base font size. Default: 11.
#' @param legend_position Character. Legend position. Default: "right".
#' @param bar_width Numeric. Width of bars. Default: 0.8.
#' @param margins Numeric vector of length 4. Plot margins in points
#'   (bottom, left, top, right). Default: c(80, 5, 5, 5). Increase bottom
#'   margin for long sample names.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Basic stacked bar
#' AETHER_plot_cibersort_proportions(cibersort_result)
#'
#' # Grouped by condition
#' AETHER_plot_cibersort_proportions(cibersort_result, group_by = sample_groups)
#'
#' }
#' @export
AETHER_plot_cibersort_proportions <- function(result,
                                          group_by = NULL,
                                          colors = NULL,
                                          title = "Cell Type Proportions",
                                          text_size = 11,
                                          legend_position = "right",
                                          bar_width = 0.8,
                                          margins = c(5, 60, 5, 5)) {

  mat <- .cibersort_extract_proportions(result)
  cell_types <- colnames(mat)
  samples <- rownames(mat)

  # Build long-format data
  df <- data.frame(
    sample = rep(samples, times = length(cell_types)),
    cell_type = rep(cell_types, each = length(samples)),
    proportion = as.vector(mat),
    stringsAsFactors = FALSE
  )

  # Handle grouping
  if (!is.null(group_by)) {
    if (length(group_by) != length(samples)) {
      stop("'group_by' length (", length(group_by), ") must match number of samples (",
           length(samples), ")", call. = FALSE)
    }
    group_map <- data.frame(
      sample = samples,
      group = group_by,
      stringsAsFactors = FALSE
    )
    df <- merge(df, group_map, by = "sample")

    # Order samples by group
    sample_order <- samples[order(group_by)]
    df$sample <- factor(df$sample, levels = sample_order)
  } else {
    df$sample <- factor(df$sample, levels = samples)
  }

  # Colors
  if (is.null(colors)) {
    pal <- .cibersort_cell_palette(length(cell_types))
    names(pal) <- cell_types
  } else {
    pal <- colors
  }

  # Plot
  p <- ggplot(df, aes(x = sample, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", width = bar_width) +
    scale_fill_manual(values = pal, name = "Cell Type") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(
      title = title,
      x = NULL,
      y = "Proportion"
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend_position,
      plot.margin = margin(margins[3], margins[4], margins[1], margins[2], unit = "pt")
    )

  # Facet by group if provided
  if (!is.null(group_by)) {
    p <- p + facet_grid(~ group, scales = "free_x", space = "free_x")
  }

  return(p)
}


# ==============================================================================
# HEATMAP
# ==============================================================================

#' Plot cell type proportions as heatmap
#'
#' Creates a heatmap of cell type proportions (cell types as rows, samples
#' as columns) using pheatmap.
#'
#' @param result An artemis_cibersort object or a proportions matrix
#'   (samples x cell types).
#' @param annotation_col Data.frame of sample annotations for column
#'   annotation. Rownames must match sample names. Default: NULL.
#' @param cluster_rows Logical. Cluster cell types. Default: TRUE.
#' @param cluster_cols Logical. Cluster samples. Default: TRUE.
#' @param scale Character. Scale by "row", "column", or "none".
#'   Default: "none" (proportions are already comparable).
#' @param color_palette Character vector of colors for heatmap gradient.
#'   Default: white-to-dark-red sequential palette.
#' @param title Character. Plot title. Default: "Cell Type Proportions".
#' @param legend_title Character. Title for the color legend. Default:
#'   "Proportion".
#' @param show_rownames Logical. Show cell type names. Default: TRUE.
#' @param show_colnames Logical. Show sample names. Default: auto (TRUE if
#'   <= 30 samples).
#' @param ... Additional arguments passed to pheatmap::pheatmap().
#'
#' @return A pheatmap object (invisibly).
#'
#' @examples
#' \dontrun{
#' AETHER_plot_cibersort_heatmap(cibersort_result)
#'
#' # With sample annotations
#' annot <- data.frame(condition = groups, row.names = sample_names)
#' AETHER_plot_cibersort_heatmap(cibersort_result, annotation_col = annot)
#'
#' }
#' @export
AETHER_plot_cibersort_heatmap <- function(result,
                                      annotation_col = NULL,
                                      cluster_rows = TRUE,
                                      cluster_cols = TRUE,
                                      scale = "none",
                                      color_palette = NULL,
                                      title = "Cell Type Proportions",
                                      legend_title = "Proportion",
                                      show_rownames = TRUE,
                                      show_colnames = TRUE,
                                      ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')")
  }

  mat <- .cibersort_extract_proportions(result)

  # Transpose: cell types as rows, samples as columns
  mat <- t(mat)

  # Auto-detect column name display
  if (is.null(show_colnames)) {
    show_colnames <- ncol(mat) <= 30
  }

  # Default color palette: white to dark red
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("#FFFFFF", "#FEE0D2", "#FC9272", "#DE2D26"))(100)
  }

  # Build legend labels with title
  n_colors <- length(color_palette)
  legend_breaks <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 5)
  legend_labels <- round(legend_breaks, 3)
  legend_labels[length(legend_labels)] <- paste0(legend_labels[length(legend_labels)],
                                                  "\n", legend_title)

  p <- pheatmap::pheatmap(
    mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale = scale,
    color = color_palette,
    annotation_col = annotation_col,
    main = title,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    legend_breaks = legend_breaks,
    legend_labels = legend_labels,
    ...
  )

  invisible(p)
}


# ==============================================================================
# BOXPLOT BY CONDITION
# ==============================================================================

#' Plot cell type proportions as boxplots by condition
#'
#' Creates faceted boxplots showing cell type proportions grouped by
#' condition. Each facet represents a cell type with boxplots per group.
#'
#' @param result An artemis_cibersort object or a proportions matrix
#'   (samples x cell types).
#' @param group_by Character vector of condition labels per sample. Required.
#'   Length must match number of samples. If the vector is **named** (names =
#'   group labels, values = color hex codes), the names are used as group
#'   labels and the values as group colors, overriding the \code{colors}
#'   parameter.
#' @param cell_types Character vector of cell types to include. Default: NULL
#'   (all non-zero cell types).
#' @param top_n Integer. Show only the top N cell types by mean proportion.
#'   Default: NULL (show all).
#' @param colors Named character vector of colors per group. If NULL, uses
#'   a default palette. Ignored if \code{group_by} is a named vector.
#' @param title Character. Plot title. Default: "Cell Type Proportions by Group".
#' @param text_size Numeric. Base font size. Default: 11.
#' @param point_size Numeric. Size of jittered points. Default: 2.
#' @param show_points Logical. Overlay jittered data points. Default: TRUE.
#' @param ncol Integer. Number of facet columns. Default: NULL (auto).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' AETHER_plot_cibersort_boxplot(cibersort_result, group_by = sample_conditions)
#'
#' # Top 10 cell types only
#' AETHER_plot_cibersort_boxplot(cibersort_result, group_by = groups, top_n = 10)
#'
#' # Named vector: names = groups, values = colors
#' group_colors <- c(CTRL = "#1f77b4", UVC = "#d62728")
#' groups_named <- group_colors[sample_sheet$condition]
#' AETHER_plot_cibersort_boxplot(cibersort_result, group_by = groups_named)
#'
#' }
#' @export
AETHER_plot_cibersort_boxplot <- function(result,
                                      group_by,
                                      cell_types = NULL,
                                      top_n = NULL,
                                      colors = NULL,
                                      title = "Cell Type Proportions by Group",
                                      text_size = 11,
                                      point_size = 2,
                                      show_points = TRUE,
                                      ncol = NULL) {

  if (missing(group_by)) {
    stop("'group_by' is required for boxplot visualization.", call. = FALSE)
  }

  mat <- .cibersort_extract_proportions(result)
  samples <- rownames(mat)

  # Detect named vector: names = group labels, values = colors
  group_by_colors <- NULL
  if (!is.null(names(group_by))) {
    group_by_colors <- group_by           # values are colors
    group_by <- names(group_by)           # names are group labels
  }

  if (length(group_by) != length(samples)) {
    stop("'group_by' length (", length(group_by), ") must match number of samples (",
         length(samples), ")", call. = FALSE)
  }

  # Subset to requested cell types
  if (!is.null(cell_types)) {
    valid <- cell_types[cell_types %in% colnames(mat)]
    if (length(valid) == 0) {
      stop("None of the specified cell_types found in results.", call. = FALSE)
    }
    mat <- mat[, valid, drop = FALSE]
  }

  # Filter to top N by mean proportion
  if (!is.null(top_n) && ncol(mat) > top_n) {
    mean_props <- colMeans(mat, na.rm = TRUE)
    top_types <- names(sort(mean_props, decreasing = TRUE))[seq_len(top_n)]
    mat <- mat[, top_types, drop = FALSE]
  }

  all_cell_types <- colnames(mat)

  # Build long-format data
  df <- data.frame(
    sample = rep(samples, times = length(all_cell_types)),
    cell_type = rep(all_cell_types, each = length(samples)),
    proportion = as.vector(mat),
    group = rep(group_by, times = length(all_cell_types)),
    stringsAsFactors = FALSE
  )

  # Order cell types by mean proportion (descending)
  mean_order <- names(sort(colMeans(mat, na.rm = TRUE), decreasing = TRUE))
  df$cell_type <- factor(df$cell_type, levels = mean_order)

  # Colors: named group_by > explicit colors param > default palette
  groups <- unique(group_by)
  if (!is.null(group_by_colors)) {
    # Colors extracted from named group_by vector
    color_map <- unique(data.frame(
      group = group_by,
      color = group_by_colors,
      stringsAsFactors = FALSE
    ))
    pal <- setNames(color_map$color, color_map$group)
  } else if (!is.null(colors)) {
    pal <- colors
  } else {
    default_pal <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                     "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#aec7e8",
                     "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                     "#f7b6d2", "#dbdb8d", "#9edae5", "#393b79", "#637939")
    if (length(groups) > length(default_pal)) {
      pal <- colorRampPalette(default_pal)(length(groups))
    } else {
      pal <- default_pal[seq_along(groups)]
    }
    names(pal) <- groups
  }

  # Plot
  p <- ggplot(df, aes(x = group, y = proportion, fill = group)) +
    geom_boxplot(outlier.shape = if (show_points) NA else 19,
                 alpha = 0.7) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = ncol) +
    labs(
      title = title,
      x = NULL,
      y = "Proportion",
      fill = "Group"
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = text_size - 1),
      legend.position = "bottom"
    )

  if (show_points) {
    p <- p + geom_jitter(aes(color = group), width = 0.15,
                         size = point_size, alpha = 0.7,
                         show.legend = FALSE) +
      scale_color_manual(values = pal)
  }

  return(p)
}
