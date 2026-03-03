# GAIA/Aether/activity_plots.R
# Visualization functions for activity inference results
#
# Functions for plotting TF and pathway activity scores from decoupleR analysis.


# ==============================================================================
# ACTIVITY HEATMAPS
# ==============================================================================

#' Plot activity heatmap
#'
#' Creates a heatmap of activity scores (TFs or pathways vs samples).
#'
#' @param result A decoupler_result or decoupler_comparison object, OR a matrix
#'   of activity scores (sources as rows, samples as columns).
#' @param top_n Integer. Number of top sources to display (by absolute mean
#'   activity). Default: 50. Use NULL to show all.
#' @param cluster_rows Logical. Cluster rows (sources). Default: TRUE.
#' @param cluster_cols Logical. Cluster columns (samples). Default: TRUE.
#' @param scale Character. Scale data by "row", "column", or "none". Default: "row".
#' @param annotation_col Data.frame of sample annotations (rownames = sample names).
#'   Default: NULL.
#' @param color_palette Character vector of colors for heatmap. Default: blue-white-red.
#' @param title Character. Plot title. Default: auto-generated.
#' @param show_rownames Logical. Show source names. Default: TRUE if <= 50 sources.
#' @param show_colnames Logical. Show sample names. Default: TRUE if <= 30 samples.
#' @param ... Additional arguments passed to pheatmap::pheatmap().
#'
#' @return A pheatmap object (invisibly).
#'
#' @details
#' If a decoupler_comparison object is provided, uses the consensus scores
#' by default. To plot a specific method's results, extract it first:
#' comparison$results$ulm$activities
#'
#' @examples
#' \dontrun{
#' # From decoupler result
#' AETHER_plot_activity_heatmap(tf_result)
#'
#' # With sample annotations
#' annot <- data.frame(
#'   condition = c("ctrl", "ctrl", "treat", "treat"),
#'   row.names = colnames(tf_result$activities)
#' )
#' AETHER_plot_activity_heatmap(tf_result, annotation_col = annot)
#'
#' # From comparison (uses consensus)
#' AETHER_plot_activity_heatmap(comparison)
#'
#' }
#' @export
AETHER_plot_activity_heatmap <- function(result,
                                          top_n = 50,
                                          cluster_rows = TRUE,
                                          cluster_cols = TRUE,
                                          scale = "row",
                                          annotation_col = NULL,
                                          color_palette = NULL,
                                          title = NULL,
                                          show_rownames = NULL,
                                          show_colnames = NULL,
                                          ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')")
  }

  # Extract activity matrix
  if (inherits(result, "decoupler_comparison")) {
    mat <- result$consensus
    if (is.null(title)) title <- "Consensus Activity Scores"
  } else if (inherits(result, "decoupler_result")) {
    mat <- result$activities
    if (is.null(title)) {
      title <- paste(result$method_name, "Activity Scores")
    }
  } else if (is.matrix(result)) {
    mat <- result
    if (is.null(title)) title <- "Activity Scores"
  } else {
    stop("'result' must be a decoupler_result, decoupler_comparison, or matrix")
  }

  # Remove rows with all NA
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]

  # Select top sources by absolute mean activity
  if (!is.null(top_n) && nrow(mat) > top_n) {
    mean_activity <- rowMeans(abs(mat), na.rm = TRUE)
    top_sources <- names(sort(mean_activity, decreasing = TRUE))[1:top_n]
    mat <- mat[top_sources, , drop = FALSE]
  }

  # Default display options
  if (is.null(show_rownames)) {
    show_rownames <- nrow(mat) <= 50
  }
  if (is.null(show_colnames)) {
    show_colnames <- ncol(mat) <= 30
  }

  # Default color palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  }

  # Create heatmap
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
    ...
  )

  invisible(p)
}


#' Plot top activities as bar chart
#'
#' Creates a bar chart showing top positive and negative activities.
#' Useful for visualizing which TFs/pathways are most active.
#'
#' @param result A decoupler_result object or activity matrix.
#' @param sample Character or integer. Which sample to plot. If NULL and multiple
#'   samples, uses mean across samples. Default: NULL.
#' @param top_n Integer. Number of top positive and negative to show. Default: 15.
#' @param title Character. Plot title. Default: auto-generated.
#' @param x_label Character. X-axis label. Default: "Activity Score".
#' @param fill_colors Character vector of length 2. Colors for positive and
#'   negative activities. Default: c("#B2182B", "#2166AC").
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Mean across all samples
#' AETHER_plot_top_activities(tf_result)
#'
#' # Specific sample
#' AETHER_plot_top_activities(tf_result, sample = "sample1")
#'
#' }
#' @export
AETHER_plot_top_activities <- function(result,
                                        sample = NULL,
                                        top_n = 15,
                                        title = NULL,
                                        x_label = "Activity Score",
                                        fill_colors = c("#B2182B", "#2166AC")) {

  # Extract activity matrix
  if (inherits(result, "decoupler_result")) {
    mat <- result$activities
    if (is.null(title)) {
      title <- paste("Top", result$method_name, "Activities")
    }
  } else if (is.matrix(result)) {
    mat <- result
    if (is.null(title)) title <- "Top Activities"
  } else {
    stop("'result' must be a decoupler_result or matrix")
  }

  # Get activity vector
  if (is.null(sample)) {
    # Mean across samples
    activity <- rowMeans(mat, na.rm = TRUE)
    subtitle <- "Mean across samples"
  } else if (is.character(sample) && sample %in% colnames(mat)) {
    activity <- mat[, sample]
    subtitle <- paste("Sample:", sample)
  } else if (is.numeric(sample) && sample <= ncol(mat)) {
    activity <- mat[, sample]
    subtitle <- paste("Sample:", colnames(mat)[sample])
  } else {
    stop("Invalid sample specification")
  }

  # Get top positive and negative
  activity <- sort(activity, decreasing = TRUE)
  top_pos <- head(activity[activity > 0], top_n)
  top_neg <- tail(activity[activity < 0], top_n)

  # Combine and create data frame
  selected <- c(top_pos, top_neg)
  df <- data.frame(
    source = factor(names(selected), levels = names(selected)),
    activity = selected,
    direction = ifelse(selected > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )

  # Reorder for plotting (most positive at top)
  df$source <- factor(df$source, levels = rev(names(selected)))

  # Create plot
  p <- ggplot(df, aes(x = activity, y = source, fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Positive" = fill_colors[1], "Negative" = fill_colors[2])) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = NULL,
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank()
    )

  return(p)
}


# ==============================================================================
# METHOD COMPARISON PLOTS
# ==============================================================================

#' Plot method comparison correlation heatmap
#'
#' Visualizes the correlation between different methods' activity scores.
#'
#' @param comparison A decoupler_comparison object.
#' @param title Character. Plot title. Default: "Method Correlation".
#' @param color_palette Character vector of colors. Default: white to blue.
#' @param show_values Logical. Display correlation values. Default: TRUE.
#' @param ... Additional arguments passed to pheatmap::pheatmap().
#'
#' @return A pheatmap object (invisibly).
#'
#' @examples
#' \dontrun{
#' comparison <- ARTEMIS_decoupler_compare_methods(mat, network, methods = c("ulm", "mlm", "wsum"))
#' AETHER_plot_method_correlation(comparison)
#'
#' }
#' @export
AETHER_plot_method_correlation <- function(comparison,
                                            title = "Method Correlation",
                                            color_palette = NULL,
                                            show_values = TRUE,
                                            ...) {

  if (!inherits(comparison, "decoupler_comparison")) {
    stop("'comparison' must be a decoupler_comparison object")
  }

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.")
  }

  cor_mat <- comparison$correlations

  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("#FFFFFF", "#2166AC"))(100)
  }

  p <- pheatmap::pheatmap(
    cor_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    main = title,
    display_numbers = show_values,
    number_format = "%.2f",
    ...
  )

  invisible(p)
}


#' Plot method agreement / variability
#'
#' Shows mean activity with min-max range across methods for each source.
#' Visualizes which TFs/pathways have consistent results across methods vs
#' those where methods disagree.
#'
#' @param comparison A decoupler_comparison object.
#' @param top_n Integer. Number of sources to show. Default: 30.
#' @param order_by Character. How to order sources: "activity" (absolute mean),
#'   "variance" (SD across methods), or "name". Default: "activity".
#' @param title Character. Plot title. Default: "Activity Across Methods".
#' @param point_size Numeric. Size of mean activity points. Default: 2.
#' @param colors Character vector of length 2. Colors for positive and negative
#'   mean activities. Default: c("#B2182B", "#2166AC").
#'
#' @return A ggplot object.
#'
#' @details
#' The plot shows:
#' - Point: mean activity across methods
#' - Error bars: range (min to max) across methods
#' - Color: whether mean activity is positive (red) or negative (blue)
#'
#' Sources where methods agree will have small error bars. Sources where
#' methods disagree (or even disagree on direction) will have large bars
#' crossing zero.
#'
#' @examples
#' \dontrun{
#' AETHER_plot_method_agreement(comparison)
#' AETHER_plot_method_agreement(comparison, order_by = "variance")
#'
#' }
#' @export
AETHER_plot_method_agreement <- function(comparison,
                                          top_n = 30,
                                          order_by = "activity",
                                          title = "Activity Across Methods",
                                          point_size = 2,
                                          colors = c("#B2182B", "#2166AC")) {

  if (!inherits(comparison, "decoupler_comparison")) {
    stop("'comparison' must be a decoupler_comparison object")
  }

  # Get summary

  df <- comparison$summary

  # Order sources
  if (order_by == "activity") {
    df <- df[order(abs(df$mean_activity), decreasing = TRUE), ]
  } else if (order_by == "variance") {
    df <- df[order(df$sd_activity, decreasing = TRUE), ]
  } else if (order_by == "name") {
    df <- df[order(df$source), ]
  }

  df <- head(df, top_n)

  # Determine direction
  df$direction <- ifelse(df$mean_activity >= 0, "Positive", "Negative")

  # Order for plotting (reversed so top is at top of plot)
  df$source <- factor(df$source, levels = rev(df$source))

  # Create plot
  p <- ggplot(df, aes(x = mean_activity, y = source, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(xmin = min_activity, xmax = max_activity),
      height = 0.3,
      linewidth = 0.5,
      orientation = "y"
    ) +
    geom_point(size = point_size) +
    scale_color_manual(
      values = c("Positive" = colors[1], "Negative" = colors[2]),
      name = "Direction"
    ) +
    labs(
      title = title,
      subtitle = paste("Top", top_n, "sources | Point = mean, bars = range across methods"),
      x = "Activity Score",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )

  return(p)
}


#' Plot activity comparison across methods
#'
#' Creates a faceted plot showing activity distributions for each method.
#'
#' @param comparison A decoupler_comparison object.
#' @param sources Character vector. Specific sources to plot. If NULL, uses
#'   top sources by mean activity. Default: NULL.
#' @param top_n Integer. If sources is NULL, how many top sources to show.
#'   Default: 10.
#' @param plot_type Character. "boxplot" or "violin". Default: "boxplot".
#' @param title Character. Plot title. Default: "Activity by Method".
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Top 10 sources across methods
#' AETHER_plot_activity_by_method(comparison)
#'
#' # Specific sources
#' AETHER_plot_activity_by_method(comparison, sources = c("TP53", "MYC", "STAT3"))
#'
#' }
#' @export
AETHER_plot_activity_by_method <- function(comparison,
                                            sources = NULL,
                                            top_n = 10,
                                            plot_type = "boxplot",
                                            title = "Activity by Method") {

  if (!inherits(comparison, "decoupler_comparison")) {
    stop("'comparison' must be a decoupler_comparison object")
  }

  # Select sources
  if (is.null(sources)) {
    sources <- head(comparison$summary$source, top_n)
  }

  # Build long-format data
  df_list <- lapply(names(comparison$results), function(method) {
    mat <- comparison$results[[method]]$activities
    available_sources <- intersect(sources, rownames(mat))
    if (length(available_sources) == 0) return(NULL)

    long <- data.frame(
      method = method,
      source = rep(available_sources, each = ncol(mat)),
      sample = rep(colnames(mat), times = length(available_sources)),
      activity = as.vector(t(mat[available_sources, , drop = FALSE])),
      stringsAsFactors = FALSE
    )
    return(long)
  })

  df <- do.call(rbind, df_list)

  if (nrow(df) == 0) {
    stop("No data available for specified sources")
  }

  # Create plot
  p <- ggplot(df, aes(x = method, y = activity, fill = method))

  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8)
  } else {
    p <- p + geom_boxplot(alpha = 0.7)
  }

  p <- p +
    facet_wrap(~source, scales = "free_y") +
    labs(
      title = title,
      x = "Method",
      y = "Activity Score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(p)
}


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Plot activity vs trait correlation
#'
#' Scatter plot showing correlation between activity scores and a sample trait.
#'
#' @param result A decoupler_result object.
#' @param trait Numeric vector of trait values (one per sample), OR a data.frame
#'   with a single numeric column.
#' @param source Character. Which source (TF/pathway) to plot. If NULL,
#'   shows the source with highest correlation. Default: NULL.
#' @param method Character. Correlation method: "pearson" or "spearman".
#'   Default: "pearson".
#' @param show_labels Logical. If TRUE, labels points with sample names using
#'   ggrepel for automatic non-overlapping positioning. Default: FALSE.
#' @param label_size Numeric. Size of point labels. Default: 3.
#' @param title Character. Plot title. Default: auto-generated.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Plot most correlated TF
#' AETHER_plot_activity_trait(tf_result, trait = sample_ages)
#'
#' # Specific TF with labels
#' AETHER_plot_activity_trait(tf_result, trait = sample_ages, source = "TP53",
#'                            show_labels = TRUE)
#'
#' # With non-overlapping labels (requires ggrepel)
#' AETHER_plot_activity_trait(tf_result, trait = sample_ages, show_labels = "repel")
#'
#' }
#' @export
AETHER_plot_activity_trait <- function(result,
                                        trait,
                                        source = NULL,
                                        method = "pearson",
                                        show_labels = FALSE,
                                        label_size = 3,
                                        title = NULL) {

  # Extract activities
  if (inherits(result, "decoupler_result")) {
    mat <- result$activities
  } else {
    stop("'result' must be a decoupler_result object")
  }

  # Process trait
  if (is.data.frame(trait)) {
    if (ncol(trait) != 1) {
      stop("If trait is a data.frame, it must have exactly one column")
    }
    trait_name <- colnames(trait)[1]
    trait <- trait[[1]]
  } else {
    trait_name <- "Trait"
  }

  # Align samples
  if (!is.null(names(trait))) {
    common <- intersect(colnames(mat), names(trait))
    if (length(common) == 0) {
      stop("No matching sample names between activities and trait")
    }
    mat <- mat[, common, drop = FALSE]
    trait <- trait[common]
  } else if (length(trait) != ncol(mat)) {
    stop("Length of trait must match number of samples")
  }

  # Select source
  if (is.null(source)) {
    # Find source with highest correlation
    cors <- apply(mat, 1, function(x) cor(x, trait, method = method, use = "complete.obs"))
    source <- names(which.max(abs(cors)))
    if (is.null(title)) {
      title <- paste("Most Correlated:", source)
    }
  } else if (!source %in% rownames(mat)) {
    stop("Source '", source, "' not found in activities")
  }

  activity <- mat[source, ]

  # Compute correlation
  cor_val <- cor(activity, trait, method = method, use = "complete.obs")
  cor_test <- cor.test(activity, trait, method = method)
  pval <- cor_test$p.value

  # Build data frame
  df <- data.frame(
    activity = activity,
    trait = trait,
    sample = colnames(mat)
  )

  # Create plot
  if (is.null(title)) {
    title <- source
  }

  subtitle <- sprintf("r = %.3f, p = %.2e (%s)", cor_val, pval, method)

  p <- ggplot(df, aes(x = activity, y = trait)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "#2166AC", fill = "#92C5DE")

  # Add labels if requested (use ggrepel by default if available)
  if (!isFALSE(show_labels)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = sample),
        size = label_size,
        max.overlaps = 20,
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.size = 0.3
      )
    } else {
      # Fallback to basic labels if ggrepel not installed
      p <- p + geom_text(aes(label = sample), size = label_size,
                         vjust = -0.5, hjust = 0.5, check_overlap = TRUE)
    }
  }

  p <- p +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste(source, "Activity"),
      y = trait_name
    ) +
    theme_minimal()

  return(p)
}
