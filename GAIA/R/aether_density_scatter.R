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
#' @param log_density Logical. If \code{TRUE}, apply \code{log1p} to the bin
#'   count colour scale. Useful when a few bins dominate the colour range.
#'   Default \code{FALSE}.
#' @param color_low Character. Colour for low-density bins. Default
#'   \code{"white"}.
#' @param color_high Character. Colour for high-density bins. Default
#'   \code{"darkblue"}.
#' @param show_zero_lines Logical. If \code{TRUE}, draw dashed reference lines
#'   at x = 0 and y = 0. Useful for log2 ratio plots. Default \code{FALSE}.
#'
#' @return A ggplot2 object.
#'
#' @export
AETHER_plot_density_scatter <- function(
  data,
  x_col           = "x",
  y_col           = "y",
  x_label         = NULL,
  y_label         = NULL,
  title           = NULL,
  nbins           = 100L,
  log_density     = FALSE,
  color_low       = "white",
  color_high      = "darkblue",
  show_zero_lines = FALSE
) {
  if (!x_col %in% colnames(data))
    stop("[AETHER] x_col '", x_col, "' not found in data", call. = FALSE)
  if (!y_col %in% colnames(data))
    stop("[AETHER] y_col '", y_col, "' not found in data", call. = FALSE)

  x_lab <- if (!is.null(x_label))            x_label
            else if (!is.null(attr(data, "x_label"))) attr(data, "x_label")
            else x_col
  y_lab <- if (!is.null(y_label))            y_label
            else if (!is.null(attr(data, "y_label"))) attr(data, "y_label")
            else y_col

  plot_df <- data.frame(
    x = data[[x_col]],
    y = data[[y_col]],
    stringsAsFactors = FALSE
  )

  if (log_density) {
    fill_mapping <- aes(fill = after_stat(log1p(count)))
    legend_name  <- "log1p(count)"
  } else {
    fill_mapping <- aes(fill = after_stat(count))
    legend_name  <- "count"
  }

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_bin2d(mapping = fill_mapping, bins = nbins) +
    scale_fill_gradient(low = color_low, high = color_high, name = legend_name) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_bw()

  if (show_zero_lines)
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
                 color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
                 color = "grey50")

  p
}
