###############################################################################
########### Coefficient Visualization for Regression Models ###########
###############################################################################
#
# Forest plots and related visualizations for regression coefficients,
# particularly suited for displaying odds ratios from logistic regression
# in clinical/epidemiological research.
#
###############################################################################


#' Forest Plot of Regression Coefficients
#'
#' @description Creates a forest plot displaying coefficients (or odds ratios)
#' with confidence intervals. Standard visualization for presenting regression
#' results in clinical research.
#'
#' @param coef_table A data.frame from \code{ARTEMIS_extract_coefficients()},
#'   or any data.frame with columns: variable, estimate, ci_lower, ci_upper.
#' @param null_line Numeric. Reference line value (default = 1 for odds ratios,
#'   use 0 for raw coefficients). Set to NULL to omit.
#' @param log_scale Logical. Use log scale for x-axis (default = TRUE for odds
#'   ratios, FALSE for coefficients).
#' @param sort_by Character. How to sort variables:
#'   \itemize{
#'     \item "estimate": By effect size (default)
#'     \item "significance": By p-value
#'     \item "name": Alphabetically
#'     \item "none": Keep original order
#'   }
#' @param decreasing Logical. Sort direction (default = TRUE, largest first).
#' @param colors Named vector or character. Colors for significant/non-significant:
#'   \itemize{
#'     \item Default: c(significant = "#E41A1C", nonsignificant = "#377EB8")
#'     \item Or single color for all points
#'   }
#' @param point_size Numeric. Size of point estimates (default = 3).
#' @param line_width Numeric. Width of confidence interval lines (default = 0.8).
#' @param text_size Numeric. Base text size (default = 11).
#' @param title Character. Plot title (default = "Forest Plot").
#' @param xlab Character. X-axis label (default based on log_scale).
#' @param show_values Logical. Show estimate values on right side (default = TRUE).
#' @param sig_level Numeric. Significance threshold for coloring (default = 0.05).
#'
#' @return A ggplot object that can be further customized or saved.
#'
#' @details
#' Forest plots are the standard way to present regression coefficients in
#' medical/clinical research. Key features:
#' \itemize{
#'   \item Each row = one variable
#'   \item Point = effect estimate (odds ratio or coefficient)
#'   \item Horizontal line = confidence interval
#'   \item Vertical reference line at null value (1 for OR, 0 for coef)
#' }
#'
#' Interpretation for odds ratios:
#' \itemize{
#'   \item OR > 1 and CI doesn't cross 1: Significant positive association
#'   \item OR < 1 and CI doesn't cross 1: Significant negative association
#'   \item CI crosses 1: Not statistically significant
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic forest plot of odds ratios
#' coef_table <- ARTEMIS_extract_coefficients(final_model, verbose = FALSE)
#' AETHER_plot_coefficients(coef_table)
#'
#' # For raw coefficients
#' coef_table <- ARTEMIS_extract_coefficients(final_model, format = "coefficient",
#'                                             verbose = FALSE)
#' AETHER_plot_coefficients(coef_table, null_line = 0, log_scale = FALSE)
#'
#' # Sort by p-value
#' AETHER_plot_coefficients(coef_table, sort_by = "significance")
#'
#' }
AETHER_plot_coefficients <- function(coef_table,
                                      null_line = 1,
                                      log_scale = TRUE,
                                      sort_by = "estimate",
                                      decreasing = TRUE,
                                      colors = NULL,
                                      point_size = 3,
                                      line_width = 0.8,
                                      text_size = 11,
                                      title = "Forest Plot",
                                      xlab = NULL,
                                      show_values = TRUE,
                                      sig_level = 0.05) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  # Validate input
  required_cols <- c("variable", "estimate", "ci_lower", "ci_upper")
  missing_cols <- setdiff(required_cols, colnames(coef_table))
  if (length(missing_cols) > 0) {
    stop("coef_table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Work with a copy
  plot_data <- coef_table

  # Ensure numeric
  plot_data$estimate <- as.numeric(plot_data$estimate)
  plot_data$ci_lower <- as.numeric(plot_data$ci_lower)
  plot_data$ci_upper <- as.numeric(plot_data$ci_upper)

  # Handle p-values (may be character if formatted)
  if ("p_value" %in% colnames(plot_data)) {
    plot_data$p_numeric <- as.numeric(plot_data$p_value)
    plot_data$significant <- plot_data$p_numeric < sig_level
  } else {
    plot_data$significant <- TRUE  # Assume all significant if no p-value
  }

  # Sort data
  sort_by <- tolower(sort_by)
  if (sort_by == "estimate") {
    order_idx <- order(plot_data$estimate, decreasing = decreasing)
  } else if (sort_by == "significance" && "p_numeric" %in% colnames(plot_data)) {
    order_idx <- order(plot_data$p_numeric, decreasing = !decreasing)
  } else if (sort_by == "name") {
    order_idx <- order(plot_data$variable, decreasing = decreasing)
  } else {
    order_idx <- seq_len(nrow(plot_data))
  }

  plot_data <- plot_data[order_idx, ]
  plot_data$variable <- factor(plot_data$variable, levels = rev(plot_data$variable))

  # Set up colors
  if (is.null(colors)) {
    colors <- c(significant = "#E41A1C", nonsignificant = "#377EB8")
  } else if (length(colors) == 1) {
    colors <- c(significant = colors, nonsignificant = colors)
  }

  plot_data$color_group <- ifelse(plot_data$significant, "significant", "nonsignificant")

  # Create label text for annotations (done before ggplot so it's in the data)
  plot_data$label <- sprintf("%.2f [%.2f, %.2f]",
                             plot_data$estimate,
                             plot_data$ci_lower,
                             plot_data$ci_upper)

  # Set default x-axis label
  if (is.null(xlab)) {
    xlab <- if (log_scale) "Odds Ratio (log scale)" else "Coefficient"
  }

  # Build the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = variable)) +
    # Confidence interval lines
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lower, xend = ci_upper, y = variable, yend = variable,
                   color = color_group),
      linewidth = line_width
    ) +
    # Point estimates
    ggplot2::geom_point(
      ggplot2::aes(color = color_group),
      size = point_size
    ) +
    # Color scale
    ggplot2::scale_color_manual(
      values = colors,
      labels = c(significant = paste0("p < ", sig_level),
                 nonsignificant = paste0("p >= ", sig_level)),
      name = "Significance"
    ) +
    # Labels
    ggplot2::labs(
      title = title,
      x = xlab,
      y = NULL
    ) +
    # Theme
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Add null reference line
  if (!is.null(null_line)) {
    p <- p + ggplot2::geom_vline(
      xintercept = null_line,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5
    )
  }

  # Log scale if requested
  if (log_scale) {
    # Calculate nice breaks for log scale
    all_vals <- c(plot_data$ci_lower, plot_data$ci_upper, plot_data$estimate)
    min_val <- min(all_vals[all_vals > 0], na.rm = TRUE)
    max_val <- max(all_vals, na.rm = TRUE)

    # Generate log-spaced breaks
    log_min <- floor(log10(min_val))
    log_max <- ceiling(log10(max_val))
    breaks <- 10^seq(log_min, log_max, by = 1)

    # Add 0.5, 2, 5 type breaks if range is small
    if (log_max - log_min <= 2) {
      minor_breaks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
      breaks <- minor_breaks[minor_breaks >= min_val * 0.8 & minor_breaks <= max_val * 1.2]
    }

    p <- p + ggplot2::scale_x_log10(breaks = breaks)
  }

  # Add value annotations - to the right of bar, or left if bar extends too far right
  if (show_values) {
    # Determine threshold: if ci_upper is in the top 30% of the range, put label on left
    max_upper <- max(plot_data$ci_upper, na.rm = TRUE)
    min_lower <- min(plot_data$ci_lower, na.rm = TRUE)

    if (log_scale) {
      # For log scale, work in log space
      log_range <- log10(max_upper) - log10(min_lower)
      threshold <- 10^(log10(min_lower) + 0.7 * log_range)

      # Position and alignment for each bar
      plot_data$label_x <- ifelse(
        plot_data$ci_upper > threshold,
        plot_data$ci_lower / 1.1,  # Left of bar
        plot_data$ci_upper * 1.1   # Right of bar
      )
    } else {
      # For linear scale
      data_range <- max_upper - min_lower
      threshold <- min_lower + 0.7 * data_range
      offset <- data_range * 0.02

      plot_data$label_x <- ifelse(
        plot_data$ci_upper > threshold,
        plot_data$ci_lower - offset,  # Left of bar
        plot_data$ci_upper + offset   # Right of bar
      )
    }

    # Set hjust based on position (1 = right-align for left labels, 0 = left-align for right labels)
    plot_data$label_hjust <- ifelse(
      plot_data$ci_upper > threshold,
      1,  # Right-align (text extends left from position)
      0   # Left-align (text extends right from position)
    )

    p <- p +
      ggplot2::geom_text(
        data = plot_data,
        ggplot2::aes(x = label_x, y = variable, label = label, hjust = label_hjust),
        size = text_size / 4,
        color = "gray30",
        inherit.aes = FALSE
      )
  }

  return(p)
}


#' Plot Lambda Selection Path
#'
#' @description Visualizes the cross-validation error across lambda values,
#' highlighting lambda.min and lambda.1se choices.
#'
#' @param cv_result Output from \code{ARTEMIS_select_lambda()}.
#' @param show_1se Logical. Show lambda.1se in addition to lambda.min (default = TRUE).
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#'
#' @export
#'
AETHER_plot_lambda_path <- function(cv_result,
                                     show_1se = TRUE,
                                     title = "Cross-Validation for Lambda Selection") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Build plot data
  plot_data <- data.frame(
    lambda = cv_result$cv_fit$lambda,
    cve = cv_result$cve,
    cvse = cv_result$cvse
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log10(lambda), y = cve)) +
    # Error ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cve - cvse, ymax = cve + cvse),
      alpha = 0.2,
      fill = "steelblue"
    ) +
    # CV error line
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 1.5) +
    # Lambda.min line
    ggplot2::geom_vline(
      xintercept = log10(cv_result$lambda.min),
      linetype = "dashed",
      color = "#E41A1C",
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = log10(cv_result$lambda.min),
      y = max(plot_data$cve),
      label = paste0("lambda.min\n(", cv_result$n_selected_min, " groups)"),
      hjust = -0.1,
      vjust = 1,
      size = 3,
      color = "#E41A1C"
    ) +
    # Labels
    ggplot2::labs(
      title = title,
      x = expression(log[10](lambda)),
      y = "Cross-Validation Error"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Add lambda.1se if requested
  if (show_1se) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = log10(cv_result$lambda.1se),
        linetype = "dashed",
        color = "#377EB8",
        linewidth = 0.8
      ) +
      ggplot2::annotate(
        "text",
        x = log10(cv_result$lambda.1se),
        y = max(plot_data$cve) * 0.9,
        label = paste0("lambda.1se\n(", cv_result$n_selected_1se, " groups)"),
        hjust = 1.1,
        vjust = 1,
        size = 3,
        color = "#377EB8"
      )
  }

  return(p)
}
