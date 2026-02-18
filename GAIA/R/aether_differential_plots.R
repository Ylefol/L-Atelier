###############################################################################
########### Differential / Fold-Change Visualizations ###########
###############################################################################

#' Waterfall Plot of Fold Changes
#'
#' @description Creates a waterfall plot where regions are ranked by their
#' log2 fold change. Each bar represents one region, sorted from most negative
#' to most positive change. Useful for identifying outliers and seeing the
#' distribution of effects across regions.
#'
#' @param shift_result Result from ARTEMIS_global_shift_test(), or a data.frame
#'   with columns 'region_id' and 'log2FC'.
#' @param n_label Integer. Number of top/bottom regions to label (default = 5).
#'   Set to 0 to disable labels.
#' @param fc_threshold Numeric. log2FC threshold for coloring bars as
#'   "significant" (default = 0.5, i.e., ~1.4-fold). Bars beyond this threshold
#'   are colored differently.
#' @param title Character. Plot title (default = "Waterfall Plot").
#' @param colors Named character vector with colors for "down", "neutral", "up"
#'   (default uses blue/gray/red scheme).
#'
#' @return A ggplot object.
#'
#' @details
#' The waterfall plot shows:
#' - X-axis: regions ranked by log2FC (most decreased on left, most increased on right)
#' - Y-axis: log2 fold change value
#' - Bar color: blue for decreased (< -threshold), red for increased (> threshold),
#'   gray for neutral
#' - Labels: top N most extreme regions on each end
#'
#' This visualization quickly answers:
#' - How many regions go up vs down?
#' - Are there outliers driving the effect?
#' - Is the effect uniform or heterogeneous?
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From shift test result
#' p <- AETHER_plot_waterfall(shift_result)
#'
#' # Label more regions
#' p <- AETHER_plot_waterfall(shift_result, n_label = 10)
#'
#' # Stricter threshold for coloring
#' p <- AETHER_plot_waterfall(shift_result, fc_threshold = 1.0)
#'
#' }
AETHER_plot_waterfall <- function(shift_result,
                                   n_label = 5,
                                   fc_threshold = 0.5,
                                   title = "Waterfall Plot of Fold Changes",
                                   colors = NULL) {

  # Extract per_region data

if (is.list(shift_result) && "per_region" %in% names(shift_result)) {
    df <- shift_result$per_region
    groups <- if ("groups" %in% names(shift_result)) shift_result$groups else c("A", "B")
  } else if (is.data.frame(shift_result)) {
    df <- shift_result
    groups <- c("A", "B")
  } else {
    stop("shift_result must be output from ARTEMIS_global_shift_test() or a data.frame")
  }

  # Check required columns
  if (!"log2FC" %in% colnames(df)) {
    stop("Data must contain 'log2FC' column")
  }
  if (!"region_id" %in% colnames(df)) {
    df$region_id <- paste0("region_", seq_len(nrow(df)))
  }

  # Default colors
  if (is.null(colors)) {
    colors <- c(
      down = "#2166ac",     # Blue
      neutral = "#878787",  # Gray
      up = "#b2182b"        # Red
    )
  }

  # Sort by log2FC
  df <- df[order(df$log2FC), ]
  df$rank <- seq_len(nrow(df))

  # Assign direction category
  df$direction <- ifelse(df$log2FC < -fc_threshold, "down",
                         ifelse(df$log2FC > fc_threshold, "up", "neutral"))
  df$direction <- factor(df$direction, levels = c("down", "neutral", "up"))

  # Identify regions to label
  df$label <- ""
  if (n_label > 0) {
    n_label <- min(n_label, floor(nrow(df) / 2))
    # Label most negative
    df$label[1:n_label] <- df$region_id[1:n_label]
    # Label most positive
    df$label[(nrow(df) - n_label + 1):nrow(df)] <- df$region_id[(nrow(df) - n_label + 1):nrow(df)]
  }

  # Count regions in each category
  n_down <- sum(df$direction == "down")
  n_up <- sum(df$direction == "up")
  n_neutral <- sum(df$direction == "neutral")

  # Create subtitle with counts
  subtitle <- sprintf("Decreased: %d | Unchanged: %d | Increased: %d (threshold: |log2FC| > %.2f)",
                      n_down, n_neutral, n_up, fc_threshold)

  # Build plot
  p <- ggplot(df, aes(x = rank, y = log2FC, fill = direction)) +
    geom_bar(stat = "identity", width = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-fc_threshold, fc_threshold),
               linetype = "dashed", color = "gray40", linewidth = 0.3) +
    scale_fill_manual(values = colors,
                      labels = c(down = paste0("Down (n=", n_down, ")"),
                                 neutral = paste0("Neutral (n=", n_neutral, ")"),
                                 up = paste0("Up (n=", n_up, ")"))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("Regions ranked by fold change (n = ", nrow(df), ")"),
      y = paste0("log2 Fold Change (", groups[2], " / ", groups[1], ")"),
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 10),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    )

  # Add labels for extreme regions with angled text
  if (n_label > 0) {
    label_df <- df[df$label != "", ]

    # Calculate y-axis range to set appropriate margins
    y_range <- max(abs(df$log2FC))
    margin_expand <- y_range * 0.3  # 30% extra space for labels

    p <- p +
      geom_text(
        data = label_df,
        aes(label = label),
        angle = 90,
        hjust = ifelse(label_df$log2FC < 0, 1.1, -0.1),
        vjust = 0.5,
        size = 2.5,
        color = "black"
      ) +
      coord_cartesian(
        clip = "off",
        ylim = c(-y_range - margin_expand, y_range + margin_expand)
      ) +
      theme(
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt")
      )
  }

  return(p)
}


#' Volcano-Style Plot for Anchor Analysis
#'
#' @description Creates a plot showing log2 fold change vs mean signal strength.
#' Useful for identifying whether high-signal or low-signal regions show
#' different behavior.
#'
#' @param shift_result Result from ARTEMIS_global_shift_test().
#' @param fc_threshold Numeric. log2FC threshold for highlighting (default = 0.5).
#' @param n_label Integer. Number of extreme points to label (default = 5).
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#'
#' @export
#'
AETHER_plot_fc_vs_signal <- function(shift_result,
                                      fc_threshold = 0.5,
                                      n_label = 5,
                                      title = "Fold Change vs Signal Strength") {

  if (!is.list(shift_result) || !"per_region" %in% names(shift_result)) {
    stop("shift_result must be output from ARTEMIS_global_shift_test()")
  }

  df <- shift_result$per_region
  groups <- if ("groups" %in% names(shift_result)) shift_result$groups else c("A", "B")

  # Get mean signal columns
  mean_cols <- grep("^mean_", colnames(df), value = TRUE)
  if (length(mean_cols) < 2) {
    stop("Could not find mean signal columns in per_region data")
  }

  # Calculate average signal across both conditions
  df$avg_signal <- (df[[mean_cols[1]]] + df[[mean_cols[2]]]) / 2
  df$log_avg_signal <- log2(df$avg_signal + 1)

  # Assign direction
  df$direction <- ifelse(df$log2FC < -fc_threshold, "down",
                         ifelse(df$log2FC > fc_threshold, "up", "neutral"))
  df$direction <- factor(df$direction, levels = c("down", "neutral", "up"))

  colors <- c(down = "#2166ac", neutral = "#878787", up = "#b2182b")

  # Identify extreme points to label
  df$label <- ""
  if (n_label > 0) {
    # Most extreme by fold change
    ordered_idx <- order(abs(df$log2FC), decreasing = TRUE)
    top_idx <- ordered_idx[1:min(n_label, length(ordered_idx))]
    df$label[top_idx] <- df$region_id[top_idx]
  }

  p <- ggplot(df, aes(x = log_avg_signal, y = log2FC, color = direction)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-fc_threshold, fc_threshold),
               linetype = "dashed", color = "gray40", linewidth = 0.3) +
    scale_color_manual(values = colors) +
    labs(
      title = title,
      x = "log2(mean signal + 1)",
      y = paste0("log2 Fold Change (", groups[2], " / ", groups[1], ")"),
      color = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (n_label > 0 && any(df$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = df[df$label != "", ],
      aes(label = label),
      size = 3,
      max.overlaps = 20
    )
  }

  return(p)
}
