#' 2D Density Scatter Plot
#'
#' Plots a 2D histogram where color encodes the number of data points in each
#' bin. Designed to accept output from
#' \code{\link{ARTEMIS_prepare_bw_comparison}} but works with any data frame
#' containing two numeric columns.
#'
#' @param data data.frame with at least two numeric columns.
#' @param x_col Character. Column name for x-axis values. Default \code{"x"}.
#' @param y_col Character. Column name for y-axis values. Default \code{"y"}.
#' @param x_label Character or \code{NULL}. X-axis label. Falls back to
#'   \code{attr(data, "x_label")} then \code{x_col} if not supplied.
#' @param y_label Character or \code{NULL}. Y-axis label. Falls back to
#'   \code{attr(data, "y_label")} then \code{y_col} if not supplied.
#' @param title Character or \code{NULL}. Plot title. Default \code{NULL}.
#' @param nbins Integer. Number of bins in each dimension. Default 100.
#' @param log_axes Logical. If \code{TRUE}, apply \code{log1p} to x and y
#'   values before plotting. Useful for direct signal comparisons where a few
#'   high-signal outliers compress the distribution into a corner. Axis labels
#'   are automatically wrapped in \code{log1p(...)}. Default \code{FALSE}.
#' @param log_density Logical. If \code{TRUE}, apply \code{log1p} to the bin
#'   count colour scale. Useful when a few bins dominate the colour range.
#'   Default \code{FALSE}.
#' @param fill_cap_quantile Numeric (0–1) or \code{NULL}. Caps the fill colour
#'   scale at this quantile of non-empty bin counts, preventing a single
#'   dominant bin (e.g. the zero-signal spike at the origin) from compressing
#'   the rest of the palette. Bins above the cap are painted the maximum
#'   colour (\code{scales::squish}). Default \code{0.99}. Set to \code{NULL}
#'   to disable capping.
#' @param palette Character or \code{NULL}. Named RColorBrewer palette for the
#'   fill scale. Default \code{"PiYG"}. Set to \code{NULL} to use
#'   \code{color_low}/\code{color_high} instead.
#' @param palette_direction Integer. \code{-1} reverses the palette order
#'   (default); \code{1} uses the palette as-is. For \code{"PiYG"} the
#'   RColorBrewer order is pink→green, so \code{-1} gives green (low) → pink
#'   (high). For sequential palettes like \code{"Blues"} (light→dark),
#'   \code{1} is typically preferred.
#' @param color_low Character or \code{NULL}. Low-density fill colour. Only
#'   used when \code{palette = NULL}. Default \code{NULL} (falls back to
#'   \code{"white"}).
#' @param color_high Character or \code{NULL}. High-density fill colour. Only
#'   used when \code{palette = NULL}. Default \code{NULL} (falls back to
#'   \code{"darkblue"}).
#' @param show_zero_lines Logical. If \code{TRUE}, draw dashed reference lines
#'   at x = 0 and y = 0. Useful for log2 ratio plots. Default \code{FALSE}.
#' @param show_diagonal Logical. If \code{TRUE}, draw a dotted y = x reference
#'   line. Useful for direct signal comparisons where deviation from the
#'   diagonal indicates change. Default \code{FALSE}.
#' @param quadrant_labels \code{NULL}, \code{"auto"}, or a 4-element character
#'   vector. When \code{"auto"}, labels the four quadrants created by the zero
#'   reference lines using subject names parsed from the axis labels, e.g.
#'   \code{"K27ac(+) ATAC(+)"} (Q1 top-right), \code{"K27ac(-) ATAC(+)"}
#'   (Q2 top-left), \code{"K27ac(-) ATAC(-)"} (Q3 bottom-left),
#'   \code{"K27ac(+) ATAC(-)"} (Q4 bottom-right). Subject names are extracted
#'   by stripping \code{" log2(...)"} from the axis labels. Supply a 4-element
#'   vector to use fully custom labels in the same corner order. Intended for
#'   use with \code{show_zero_lines = TRUE}. Default \code{NULL} (no labels).
#'
#' @return A ggplot2 object.
#'
#' @export
AETHER_plot_density_scatter <- function(
  data,
  x_col              = "x",
  y_col              = "y",
  x_label            = NULL,
  y_label            = NULL,
  title              = NULL,
  nbins              = 100L,
  log_axes           = FALSE,
  log_density        = FALSE,
  fill_cap_quantile  = 0.99,
  palette            = "PiYG",
  palette_direction  = -1L,
  color_low          = NULL,
  color_high         = NULL,
  show_zero_lines    = FALSE,
  show_diagonal      = FALSE,
  quadrant_labels    = NULL
) {
  if (!x_col %in% colnames(data))
    stop("[AETHER] x_col '", x_col, "' not found in data", call. = FALSE)
  if (!y_col %in% colnames(data))
    stop("[AETHER] y_col '", y_col, "' not found in data", call. = FALSE)

  x_lab <- if (!is.null(x_label))                    x_label
            else if (!is.null(attr(data, "x_label"))) attr(data, "x_label")
            else x_col
  y_lab <- if (!is.null(y_label))                    y_label
            else if (!is.null(attr(data, "y_label"))) attr(data, "y_label")
            else y_col

  plot_df <- data.frame(
    x = data[[x_col]],
    y = data[[y_col]],
    stringsAsFactors = FALSE
  )

  if (log_axes) {
    plot_df$x <- log1p(plot_df$x)
    plot_df$y <- log1p(plot_df$y)
    x_lab     <- paste0("log1p(", x_lab, ")")
    y_lab     <- paste0("log1p(", y_lab, ")")
  }

  if (log_density) {
    fill_mapping <- aes(fill = after_stat(log1p(count)))
    legend_name  <- "log1p(count)"
  } else {
    fill_mapping <- aes(fill = after_stat(count))
    legend_name  <- "count"
  }

  # Pre-bin to derive the fill scale cap before plotting
  count_cap <- NULL
  if (!is.null(fill_cap_quantile) &&
      diff(range(plot_df$x, na.rm = TRUE)) > 0 &&
      diff(range(plot_df$y, na.rm = TRUE)) > 0) {
    x_breaks  <- seq(min(plot_df$x, na.rm = TRUE), max(plot_df$x, na.rm = TRUE),
                     length.out = nbins + 1L)
    y_breaks  <- seq(min(plot_df$y, na.rm = TRUE), max(plot_df$y, na.rm = TRUE),
                     length.out = nbins + 1L)
    bin_counts <- as.integer(table(
      cut(plot_df$x, breaks = x_breaks, include.lowest = TRUE),
      cut(plot_df$y, breaks = y_breaks, include.lowest = TRUE)
    ))
    nonzero <- bin_counts[bin_counts > 0L]
    if (length(nonzero)) {
      count_cap <- as.numeric(quantile(nonzero, fill_cap_quantile))
      if (log_density) count_cap <- log1p(count_cap)
    }
  }

  fill_scale <- if (!is.null(palette)) {
    pal_colors <- RColorBrewer::brewer.pal(11, palette)
    if (palette_direction == -1L) pal_colors <- rev(pal_colors)
    scale_fill_gradientn(
      colors = pal_colors,
      name   = legend_name,
      limits = if (!is.null(count_cap)) c(0, count_cap) else NULL,
      oob    = scales::squish
    )
  } else {
    scale_fill_gradient(
      low    = if (!is.null(color_low))  color_low  else "white",
      high   = if (!is.null(color_high)) color_high else "darkblue",
      name   = legend_name,
      limits = if (!is.null(count_cap)) c(0, count_cap) else NULL,
      oob    = scales::squish
    )
  }

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_bin2d(mapping = fill_mapping, bins = nbins) +
    fill_scale +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_bw()

  if (show_zero_lines)
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
                 color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
                 color = "grey50")

  if (show_diagonal)
    p <- p +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                  linewidth = 0.4, color = "grey50")

  if (!is.null(quadrant_labels)) {
    labels <- if (identical(quadrant_labels, "auto")) {
      .subj <- function(lab) sub(" log2\\(.*", "", lab)
      xn <- .subj(x_lab)
      yn <- .subj(y_lab)
      c(paste0(xn, "(+) ", yn, "(+)"),
        paste0(xn, "(-) ", yn, "(+)"),
        paste0(xn, "(-) ", yn, "(-)"),
        paste0(xn, "(+) ", yn, "(-)"))
    } else {
      if (length(quadrant_labels) != 4L)
        stop("[AETHER] quadrant_labels must be 'auto' or a 4-element character vector ",
             "(order: Q1 top-right, Q2 top-left, Q3 bottom-left, Q4 bottom-right)",
             call. = FALSE)
      quadrant_labels
    }
    x_rng <- range(plot_df$x, na.rm = TRUE)
    y_rng <- range(plot_df$y, na.rm = TRUE)
    x_pad <- diff(x_rng) * 0.04
    y_pad <- diff(y_rng) * 0.04
    q_x   <- c(x_rng[2] - x_pad, x_rng[1] + x_pad, x_rng[1] + x_pad, x_rng[2] - x_pad)
    q_y   <- c(y_rng[2] - y_pad, y_rng[2] - y_pad, y_rng[1] + y_pad, y_rng[1] + y_pad)
    q_hj  <- c(1, 0, 0, 1)
    q_vj  <- c(1, 1, 0, 0)
    for (i in seq_along(labels))
      p <- p + annotate("text", x = q_x[i], y = q_y[i], label = labels[i],
                        hjust = q_hj[i], vjust = q_vj[i], size = 3,
                        color = "grey30", fontface = "italic")
  }

  p
}
