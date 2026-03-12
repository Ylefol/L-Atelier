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


# ==============================================================================
# DEA plot helpers (internal)
# ==============================================================================

#' Extract data.frame from a DEA result object or data.frame
#' @keywords internal
.dea_extract_df <- function(dea_result, required) {
  if (is.list(dea_result) && !is.data.frame(dea_result) && "results" %in% names(dea_result)) {
    df <- dea_result$results
  } else if (is.data.frame(dea_result)) {
    df <- dea_result
  } else {
    stop("dea_result must be a data.frame or a list with a $results element (from ARTEMIS_perform_dea()).",
         call. = FALSE)
  }
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Required column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  df
}


#' Resolve label column for DEA plots
#' @keywords internal
.dea_resolve_label_col <- function(df, label_col) {
  if (!is.null(label_col)) {
    if (!label_col %in% colnames(df)) {
      stop("label_col '", label_col, "' not found. ",
           "Available columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
    }
    return(label_col)
  }
  if ("gene_name" %in% colnames(df)) return("gene_name")
  if ("feature_id" %in% colnames(df)) return("feature_id")
  stop("Cannot auto-detect label column. Provide label_col explicitly.", call. = FALSE)
}


#' Assign four-category significance labels to a DEA data.frame
#'
#' Adds a .cat column: "up", "down", "low_reg", "non_sig".
#' NA filter values are treated as non-significant.
#' @keywords internal
.dea_assign_categories <- function(df, filter_choice, p_thresh, l2fc_thresh) {
  fval <- df[[filter_choice]]
  lfc  <- df$log2FoldChange

  sig  <- !is.na(fval) & fval < p_thresh
  up   <- sig & !is.na(lfc) & lfc >  l2fc_thresh
  down <- sig & !is.na(lfc) & lfc < -l2fc_thresh
  low  <- sig & !is.na(lfc) & abs(lfc) <= l2fc_thresh

  df$.cat <- "non_sig"
  df$.cat[up]   <- "up"
  df$.cat[down] <- "down"
  df$.cat[low]  <- "low_reg"
  df
}


# ==============================================================================
# Volcano and MA plots
# ==============================================================================

#' Volcano Plot for Differential Expression/Accessibility Results
#'
#' Creates a four-category volcano plot (-log10 p-value vs log2 fold change)
#' with embedded counts in the legend and optional gene labeling. Categories:
#' up-regulated, down-regulated, low-regulation (significant but below l2FC
#' threshold), and non-significant.
#'
#' @param dea_result A data.frame of DEA results, or a list with a
#'   \code{$results} element (as returned by \code{ARTEMIS_perform_dea()}).
#'   Required columns: \code{log2FoldChange}, \code{pvalue}, and the column
#'   named by \code{filter_choice}.
#' @param label_col Character. Column to use for point labels. If NULL,
#'   auto-selects \code{gene_name} if present, otherwise \code{feature_id}.
#' @param genes_of_interest Character vector. Values in \code{label_col} to
#'   always label. Default: NULL.
#' @param show_non_sig_interest Logical. If FALSE, genes of interest that do
#'   not meet both thresholds are not labeled. Default: TRUE.
#' @param label_top_n Integer. Label the top N features by \code{filter_choice}
#'   regardless of direction. Default: 0 (disabled).
#' @param filter_choice Column for significance filtering: \code{"padj"} or
#'   \code{"pvalue"}. Default: \code{"padj"}.
#' @param l2fc_thresh Numeric. log2 fold change threshold. Default: 1.
#' @param p_thresh Numeric. Significance threshold. Default: 0.05.
#' @param title Character. Plot title. Default: "Volcano Plot".
#' @param colors Named character vector with colors for \code{"up"},
#'   \code{"down"}, \code{"low_reg"}, \code{"non_sig"}. Default uses
#'   red/blue/green/gray.
#' @param point_size Numeric. Point size. Default: 0.8.
#' @param point_alpha Numeric. Point transparency. Default: 0.7.
#'
#' @return A ggplot object.
#'
#' @details
#' The four categories are:
#' \itemize{
#'   \item \strong{up-reg}: \code{log2FoldChange > l2fc_thresh} AND
#'     \code{filter_choice < p_thresh}
#'   \item \strong{down-reg}: \code{log2FoldChange < -l2fc_thresh} AND
#'     \code{filter_choice < p_thresh}
#'   \item \strong{low-regulation}: \code{|log2FoldChange| <= l2fc_thresh} AND
#'     \code{filter_choice < p_thresh}
#'   \item \strong{non-significant}: \code{filter_choice >= p_thresh} or NA
#' }
#' A horizontal dashed line marks the pvalue of the least-significant gene
#' that still passes the filter threshold, showing the actual significance
#' boundary in the plotted space. Vertical dashed lines mark
#' \code{±l2fc_thresh}.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_volcano(dea_result, title = "KO vs WT")
#'
#' p <- AETHER_plot_volcano(dea_result, label_top_n = 10,
#'                          genes_of_interest = c("SMUG1", "OGG1"))
#'
#' p <- AETHER_plot_volcano(dea_result, filter_choice = "pvalue", p_thresh = 0.01)
#' }
#' @export
AETHER_plot_volcano <- function(dea_result,
                                 label_col = NULL,
                                 genes_of_interest = NULL,
                                 show_non_sig_interest = TRUE,
                                 label_top_n = 0,
                                 filter_choice = "padj",
                                 l2fc_thresh = 1,
                                 p_thresh = 0.05,
                                 title = "Volcano Plot",
                                 colors = NULL,
                                 point_size = 0.8,
                                 point_alpha = 0.7) {

  df <- .dea_extract_df(dea_result, required = c("log2FoldChange", "pvalue", filter_choice))
  label_col <- .dea_resolve_label_col(df, label_col)
  df <- .dea_assign_categories(df, filter_choice, p_thresh, l2fc_thresh)
  df <- df[!is.na(df$pvalue), , drop = FALSE]

  if (is.null(colors)) {
    colors <- c(up = "#B31B21", down = "#1465AC", low_reg = "green", non_sig = "darkgray")
  }

  n_up      <- sum(df$.cat == "up")
  n_down    <- sum(df$.cat == "down")
  n_low_reg <- sum(df$.cat == "low_reg")
  n_non_sig <- sum(df$.cat == "non_sig")

  lbl <- c(
    up      = paste0("up-reg | ", filter_choice, "<", p_thresh, " (n=", n_up, ")"),
    down    = paste0("down-reg | ", filter_choice, "<", p_thresh, " (n=", n_down, ")"),
    low_reg = paste0("low-regulation | ", filter_choice, "<", p_thresh, " (n=", n_low_reg, ")"),
    non_sig = paste0("non-significant | ", filter_choice, "\u2265", p_thresh, " (n=", n_non_sig, ")")
  )

  df$Significance <- factor(lbl[df$.cat],
                             levels = c(lbl["up"], lbl["down"], lbl["low_reg"], lbl["non_sig"]))
  df$.cat <- NULL
  # Non-sig drawn first (background), significant on top
  df <- df[order(df$Significance, decreasing = TRUE), ]

  # Horizontal line at pvalue of the least-significant passing gene
  sig_mask   <- df[[filter_choice]] < p_thresh & !is.na(df[[filter_choice]])
  sig_line_y <- if (any(sig_mask)) -log10(max(df$pvalue[sig_mask], na.rm = TRUE)) else NA_real_

  # Label data: genes of interest
  labs_interest <- df[0, ]
  if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
    labs_interest <- df[df[[label_col]] %in% genes_of_interest, , drop = FALSE]
    if (!show_non_sig_interest) {
      labs_interest <- labs_interest[labs_interest$Significance != lbl["non_sig"], , drop = FALSE]
    }
  }
  labs_interest$.label <- labs_interest[[label_col]]

  # Label data: top N
  labs_top <- df[0, ]
  if (label_top_n > 0) {
    labs_top <- head(df[order(df[[filter_choice]]), ], label_top_n)
  }
  labs_top$.label <- labs_top[[label_col]]

  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_vline(xintercept = c(-l2fc_thresh, l2fc_thresh),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(
      values = setNames(unname(colors[c("up", "down", "low_reg", "non_sig")]),
                        unname(lbl[c("up", "down", "low_reg", "non_sig")])),
      breaks = unname(lbl),
      name   = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab(expression("Log"[2]*"Fold Change")) +
    ylab(expression("-log"[10]*"(p-value)")) +
    ggtitle(title) +
    theme_light() +
    theme(
      text             = element_text(size = 10),
      plot.title       = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.position  = "bottom",
      legend.direction = "vertical"
    )

  if (!is.na(sig_line_y)) {
    p <- p + geom_hline(yintercept = sig_line_y, linetype = "dashed",
                        color = "black", linewidth = 0.4)
  }

  if (nrow(labs_interest) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_interest, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 1, segment.colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  if (nrow(labs_top) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_top, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 0.5, colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  return(p)
}


#' MA Plot for Differential Expression/Accessibility Results
#'
#' Creates a four-category MA plot (log2 fold change vs log2 mean expression)
#' with embedded counts in the legend and optional gene labeling.
#'
#' @param dea_result A data.frame of DEA results, or a list with a
#'   \code{$results} element (as returned by \code{ARTEMIS_perform_dea()}).
#'   Required columns: \code{baseMean}, \code{log2FoldChange}, and the column
#'   named by \code{filter_choice}.
#' @param label_col Character. Column to use for point labels. If NULL,
#'   auto-selects \code{gene_name} if present, otherwise \code{feature_id}.
#' @param genes_of_interest Character vector. Values in \code{label_col} to
#'   label. Default: NULL.
#' @param label_top_n Integer. Label the top N features by \code{filter_choice}.
#'   Default: 0 (disabled).
#' @param filter_choice Column for significance filtering: \code{"padj"} or
#'   \code{"pvalue"}. Default: \code{"padj"}.
#' @param l2fc_thresh Numeric. log2 fold change threshold. Default: 1.
#' @param p_thresh Numeric. Significance threshold. Default: 0.05.
#' @param title Character. Plot title. Default: "MA Plot".
#' @param colors Named character vector with colors for \code{"up"},
#'   \code{"down"}, \code{"low_reg"}, \code{"non_sig"}. Default uses
#'   red/blue/green/gray.
#' @param point_size Numeric. Point size. Default: 0.8.
#' @param point_alpha Numeric. Point transparency. Default: 0.7.
#'
#' @return A ggplot object.
#'
#' @details
#' X-axis is \code{log2(baseMean + 1)}. Y-axis is \code{log2FoldChange}.
#' A solid red line marks y = 0. Dashed black lines mark \code{±l2fc_thresh}.
#' Non-significant points are drawn first (background) with significant points
#' on top.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_ma(dea_result, title = "KO vs WT")
#'
#' p <- AETHER_plot_ma(dea_result, label_top_n = 10,
#'                     genes_of_interest = c("SMUG1", "OGG1"))
#' }
#' @export
AETHER_plot_ma <- function(dea_result,
                            label_col = NULL,
                            genes_of_interest = NULL,
                            label_top_n = 0,
                            filter_choice = "padj",
                            l2fc_thresh = 1,
                            p_thresh = 0.05,
                            title = "MA Plot",
                            colors = NULL,
                            point_size = 0.8,
                            point_alpha = 0.7) {

  df <- .dea_extract_df(dea_result, required = c("baseMean", "log2FoldChange", filter_choice))
  label_col <- .dea_resolve_label_col(df, label_col)
  df <- .dea_assign_categories(df, filter_choice, p_thresh, l2fc_thresh)
  df <- df[!is.na(df$baseMean) & !is.na(df$log2FoldChange), , drop = FALSE]
  df$baseMean_log2 <- log2(df$baseMean + 1)

  if (is.null(colors)) {
    colors <- c(up = "#B31B21", down = "#1465AC", low_reg = "green", non_sig = "darkgray")
  }

  n_up      <- sum(df$.cat == "up")
  n_down    <- sum(df$.cat == "down")
  n_low_reg <- sum(df$.cat == "low_reg")
  n_non_sig <- sum(df$.cat == "non_sig")

  lbl <- c(
    up      = paste0("up-reg | ", filter_choice, "<", p_thresh, " (n=", n_up, ")"),
    down    = paste0("down-reg | ", filter_choice, "<", p_thresh, " (n=", n_down, ")"),
    low_reg = paste0("low-regulation | ", filter_choice, "<", p_thresh, " (n=", n_low_reg, ")"),
    non_sig = paste0("non-significant | ", filter_choice, "\u2265", p_thresh, " (n=", n_non_sig, ")")
  )

  df$Significance <- factor(lbl[df$.cat],
                             levels = c(lbl["up"], lbl["down"], lbl["low_reg"], lbl["non_sig"]))
  df$.cat <- NULL
  # Non-sig drawn first (background), significant on top
  df <- df[order(df$Significance, decreasing = TRUE), ]

  # Label data: genes of interest
  labs_interest <- df[0, ]
  if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
    labs_interest <- df[df[[label_col]] %in% genes_of_interest, , drop = FALSE]
  }
  labs_interest$.label <- labs_interest[[label_col]]

  # Label data: top N
  labs_top <- df[0, ]
  if (label_top_n > 0) {
    labs_top <- head(df[order(df[[filter_choice]]), ], label_top_n)
  }
  labs_top$.label <- labs_top[[label_col]]

  p <- ggplot(df, aes(x = baseMean_log2, y = log2FoldChange, color = Significance)) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", linewidth = 0.5) +
    geom_hline(yintercept = c(-l2fc_thresh, l2fc_thresh),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(
      values = setNames(unname(colors[c("up", "down", "low_reg", "non_sig")]),
                        unname(lbl[c("up", "down", "low_reg", "non_sig")])),
      breaks = unname(lbl),
      name   = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab(expression("log"[2]*"(baseMean + 1)")) +
    ylab(expression("log"[2]*"Fold Change")) +
    ggtitle(title) +
    theme_light() +
    theme(
      text             = element_text(size = 10),
      plot.title       = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.position  = "bottom",
      legend.direction = "vertical"
    )

  if (nrow(labs_interest) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_interest, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 1, segment.colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  if (nrow(labs_top) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_top, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 0.5, colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  return(p)
}


# ==============================================================================
# DEG-Peak Scatter Plot
# ==============================================================================

#' Scatter Plot of DEG LFC vs Distance to Nearest ATAC Peak
#'
#' @description Visualises the output of \code{APOLLO_link_degs_to_peaks()} as
#' a scatter plot of log2 fold change against signed distance to the TSS
#' (negative = upstream, positive = downstream). Each point is one peak-DEG
#' pair. Point size encodes statistical significance (-log10 padj); colour
#' encodes expression direction.
#'
#' @param deg_peaks An \code{apollo_deg_peaks} object from
#'   \code{APOLLO_link_degs_to_peaks()}.
#' @param color_up Character. Colour for upregulated DEGs.
#'   Default \code{"#D6604D"} (red).
#' @param color_down Character. Colour for downregulated DEGs.
#'   Default \code{"#2166AC"} (blue).
#' @param alpha Numeric. Point transparency. Default \code{0.7}.
#' @param size_range Numeric vector of length 2. Min and max point sizes
#'   (mapped to -log10 padj). Default \code{c(1, 5)}.
#' @param label_top Integer. Number of top genes (by padj) to label with
#'   \code{ggrepel}. Default \code{10}. Set to \code{0} to suppress labels.
#' @param label_col Character. Column in \code{$linked} used for labels.
#'   Default \code{NULL} (uses the \code{gene_col} from the linking step).
#' @param show_tss_line Logical. Draw a vertical dashed line at distance = 0
#'   (the TSS). Default \code{TRUE}.
#' @param show_lfc_line Logical. Draw a horizontal dashed line at LFC = 0.
#'   Default \code{TRUE}.
#' @param title Character. Plot title. Default \code{NULL} (auto-generated).
#' @param x_limits Numeric vector of length 2 or \code{NULL}. X-axis limits
#'   in bp. Default \code{NULL} (auto from data).
#'
#' @return A \code{ggplot} object.
#'
#' @export
AETHER_plot_deg_peak_scatter <- function(deg_peaks,
                                          color_up      = "#D6604D",
                                          color_down    = "#2166AC",
                                          alpha         = 0.7,
                                          size_range    = c(1, 5),
                                          label_top     = 10,
                                          label_col     = NULL,
                                          show_tss_line = TRUE,
                                          show_lfc_line = TRUE,
                                          title         = NULL,
                                          x_limits      = NULL) {

  if (!inherits(deg_peaks, "apollo_deg_peaks"))
    stop("deg_peaks must be an apollo_deg_peaks object from ",
         "APOLLO_link_degs_to_peaks()", call. = FALSE)

  df <- deg_peaks$linked

  if (nrow(df) == 0)
    stop("No peak-DEG pairs in deg_peaks$linked. Nothing to plot.", call. = FALSE)

  gene_col <- deg_peaks$params$gene_col
  if (is.null(label_col)) label_col <- gene_col
  if (!label_col %in% colnames(df))
    stop("label_col '", label_col, "' not found in deg_peaks$linked", call. = FALSE)

  if (!"log2FoldChange" %in% colnames(df))
    stop("deg_peaks$linked must contain 'log2FoldChange'", call. = FALSE)
  if (!"padj" %in% colnames(df))
    stop("deg_peaks$linked must contain 'padj'", call. = FALSE)

  df$.direction <- ifelse(df$log2FoldChange >= 0, "Up", "Down")
  df$.neg_log_p <- -log10(pmax(df$padj, 1e-300))
  df$.label     <- df[[label_col]]

  if (is.null(title)) {
    dist_kb <- round(deg_peaks$params$distance / 1000)
    title <- paste0("DEG LFC vs distance to TSS  (\u00B1", dist_kb,
                    " kb)  \u2014  ",
                    deg_peaks$stats$n_pairs, " peak-DEG pairs")
  }

  p <- ggplot(df, aes(x = .data$distance_to_tss,
                      y = .data$log2FoldChange)) +
    geom_point(
      aes(color = .data$.direction, size = .data$.neg_log_p),
      alpha = alpha
    ) +
    scale_color_manual(values = c("Up" = color_up, "Down" = color_down),
                       name   = "Direction") +
    scale_size_continuous(range = size_range,
                          name  = expression(-log[10](padj))) +
    scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      limits = x_limits
    ) +
    labs(
      x     = paste0("Signed distance to TSS (bp)\n",
                     "\u25C4 upstream  |  downstream \u25BA"),
      y     = expression(log[2]~"Fold Change"),
      title = title
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.title       = element_text(size = 10)
    )

  if (show_tss_line)
    p <- p + geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.4)
  if (show_lfc_line)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.4)

  # Label top genes by padj — one label per gene, placed at its closest peak
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    df_genes <- df[order(df$abs_distance), ]
    df_genes <- df_genes[!duplicated(df_genes[[gene_col]]), ]
    df_genes <- df_genes[order(df_genes$padj), ]
    labs_df  <- head(df_genes, label_top)

    if (nrow(labs_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data        = labs_df,
        mapping     = aes(label = .data$.label),
        size        = 3,
        color       = "black",
        show.legend = FALSE,
        box.padding = 0.3,
        max.overlaps = 20
      )
    }
  }

  p
}
