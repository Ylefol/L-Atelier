library(ggplot2)
###############################################################################
########### Sparse CCA Visualization ###########
###############################################################################

#' Plot Feature Weights from Sparse CCA
#'
#' @description Creates a stem plot visualization of canonical weights from sCCA,
#' highlighting non-zero features and optionally showing true signal regions.
#' Prints summary statistics to console.
#'
#' @param weights Numeric vector of canonical weights (from sCCA result)
#' @param method_name Character string for plot title (default = "sCCA")
#' @param true_features Optional numeric vector c(start, end) indicating true signal region.
#'   If provided, region will be highlighted with blue shading. Useful for simulated data.
#' @param threshold Threshold for considering weights as non-zero (default = 1e-6)
#'
#' @return ggplot2 object. Also prints summary statistics as side effect.
#'
#' @details
#' Creates a stem plot where:
#' - Non-zero weights (|weight| > threshold) are shown in orange
#' - Zero/negligible weights are shown in gray
#' - True signal regions (if specified) are highlighted with blue background
#' - Prints: number of non-zero features, max weight, mean absolute weight
#'
#' @export
#'
#' @examples
#' # Simulate sCCA result
#' weights <- c(rnorm(20, mean = 0.5), rep(0, 80))
#'
#' # Basic plot
#' p <- AETHER_plot_feature_weights(weights, method_name = "ConvCCA")
#'
#' # Plot with true signal region highlighted
#' p <- AETHER_plot_feature_weights(weights, method_name = "ConvCCA",
#'                           true_features = c(1, 20))
#'
AETHER_plot_feature_weights <- function(weights, method_name = "sCCA",
                                  true_features = NULL, threshold = 1e-6) {
  # Create data frame
  df <- data.frame(
    feature_index = 1:length(weights),
    weight = as.numeric(weights)
  )

  # Classify features as zero or non-zero based on threshold
  df$is_nonzero <- abs(df$weight) > threshold

  # Create the plot
  p <- ggplot(df, aes(x = feature_index, y = weight)) +
    # Add horizontal reference line at y = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    # Add vertical segments from 0 to each weight (stem plot style)
    geom_segment(aes(x = feature_index, xend = feature_index,
                     y = 0, yend = weight, color = is_nonzero),
                 linewidth = 0.8, alpha = 0.7) +
    # Add points at the end of each segment
    geom_point(aes(color = is_nonzero), size = 1.5, alpha = 0.8) +
    # Color scheme
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray70"),
                       labels = c("TRUE" = "Non-zero", "FALSE" = "Zero"),
                       name = "Weight Status") +
    # Add true feature region if provided
    {if (!is.null(true_features)) {
      annotate("rect", xmin = true_features[1], xmax = true_features[2],
               ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue")
    }} +
    # Labels and theme
    labs(
      title = paste0("Feature Weights: ", method_name),
      subtitle = if (!is.null(true_features)) {
        paste0("Blue region shows true signal features (",
               true_features[1], "-", true_features[2], ")")
      } else NULL,
      x = "Feature Index",
      y = "Canonical Weight"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.title = element_text(size = 12),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  # Print summary statistics
  n_nonzero <- sum(df$is_nonzero)
  cat(sprintf("\n%s Summary:\n", method_name))
  cat(sprintf("  Non-zero features: %d / %d (%.1f%%)\n",
              n_nonzero, length(weights), 100 * n_nonzero / length(weights)))
  cat(sprintf("  Max weight: %.4f\n", max(abs(weights))))
  cat(sprintf("  Mean absolute weight: %.4f\n", mean(abs(weights))))

  return(p)
}

#' Plot Feature Weights from Multi-Dataset Sparse CCA
#'
#' @description Creates faceted stem plots showing canonical weights from multi-dataset sCCA.
#' Each dataset gets its own panel. Highlights non-zero features and optionally shows
#' true signal regions. Prints summary statistics for each dataset.
#'
#' @param weights_list List of numeric vectors, each containing canonical weights for one dataset.
#'   Length of list = number of datasets integrated.
#' @param method_name Character string for plot title (default = "sCCA")
#' @param dataset_names Optional character vector of dataset names for panel labels.
#'   If NULL, uses "Dataset 1", "Dataset 2", etc.
#' @param true_features_list Optional list of numeric vectors, each c(start, end) indicating
#'   true signal region for that dataset. NULL elements mean no highlighting for that dataset.
#'   Useful for simulated data validation.
#' @param threshold Threshold for considering weights as non-zero (default = 1e-6)
#'
#' @return ggplot2 object with faceted panels. Also prints summary statistics as side effect.
#'
#' @details
#' Creates faceted stem plots (one panel per dataset) where:
#' - Non-zero weights (|weight| > threshold) are shown in orange
#' - Zero/negligible weights are shown in gray
#' - True signal regions (if specified) are highlighted with blue background
#' - Each panel shows one dataset's canonical weights
#' - Prints for each dataset: number of non-zero features, max weight, mean absolute weight
#'
#' @export
#'
#' @examples
#' # Simulate multi-dataset sCCA results
#' weights1 <- c(rnorm(20, mean = 0.5), rep(0, 80))
#' weights2 <- c(rnorm(15, mean = 0.6), rep(0, 85))
#' weights3 <- c(rnorm(10, mean = 0.4), rep(0, 40))
#'
#' # Basic plot
#' p <- AETHER_plot_multi_feature_weights(
#'   weights_list = list(weights1, weights2, weights3),
#'   method_name = "multi.convCCA",
#'   dataset_names = c("RNA-seq", "ATAC-seq", "CUT&TAG")
#' )
#'
#' # Plot with true signal regions highlighted
#' p <- AETHER_plot_multi_feature_weights(
#'   weights_list = list(weights1, weights2, weights3),
#'   method_name = "multi.convCCA",
#'   dataset_names = c("RNA-seq", "ATAC-seq", "CUT&TAG"),
#'   true_features_list = list(c(1, 20), c(1, 15), c(1, 10))
#' )
#'
AETHER_plot_multi_feature_weights <- function(weights_list,
                                        method_name = "sCCA",
                                        dataset_names = NULL,
                                        true_features_list = NULL,
                                        threshold = 1e-6) {

  # Input validation
  n_datasets <- length(weights_list)

  # Generate dataset names if not provided
  if (is.null(dataset_names)) {
    dataset_names <- paste0("Dataset ", 1:n_datasets)
  }

  # Build long-format data frame
  df_list <- list()
  true_regions <- list()

  for (i in 1:n_datasets) {
    weights <- weights_list[[i]]

    # Create data frame for this dataset
    df_list[[i]] <- data.frame(
      dataset = dataset_names[i],
      feature_index = 1:length(weights),
      weight = as.numeric(weights),
      is_nonzero = abs(weights) > threshold
    )

    # Store true feature regions if provided
    if (!is.null(true_features_list) && !is.null(true_features_list[[i]])) {
      true_regions[[i]] <- data.frame(
        dataset = dataset_names[i],
        xmin = true_features_list[[i]][1],
        xmax = true_features_list[[i]][2]
      )
    }
  }

  # Combine all datasets
  df <- do.call(rbind, df_list)
  df$dataset <- factor(df$dataset, levels = dataset_names)

  # Create the plot
  p <- ggplot(df, aes(x = feature_index, y = weight)) +
    # Add horizontal reference line at y = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    # Add vertical segments from 0 to each weight (stem plot style)
    geom_segment(aes(x = feature_index, xend = feature_index,
                     y = 0, yend = weight, color = is_nonzero),
                 linewidth = 0.8, alpha = 0.7) +
    # Add points at the end of each segment
    geom_point(aes(color = is_nonzero), size = 1.5, alpha = 0.8) +
    # Color scheme
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray70"),
                       labels = c("TRUE" = "Non-zero", "FALSE" = "Zero"),
                       name = "Weight Status") +
    # Facet by dataset
    facet_wrap(~dataset, ncol = 1, scales = "free_x") +
    # Labels and theme
    labs(
      title = paste0("Feature Weights: ", method_name),
      subtitle = if (!is.null(true_features_list)) {
        "Blue regions show true signal features"
      } else NULL,
      x = "Feature Index",
      y = "Canonical Weight"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.title = element_text(size = 12),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11)
    )

  # Add true feature regions if provided
  if (length(true_regions) > 0) {
    true_regions_df <- do.call(rbind, true_regions)
    p <- p + geom_rect(data = true_regions_df,
                       aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       alpha = 0.1, fill = "blue", inherit.aes = FALSE)
  }

  # Print summary statistics for each dataset
  cat(sprintf("\n%s Summary:\n", method_name))
  for (i in 1:n_datasets) {
    weights <- weights_list[[i]]
    n_nonzero <- sum(abs(weights) > threshold)
    cat(sprintf("\n%s:\n", dataset_names[i]))
    cat(sprintf("  Non-zero features: %d / %d (%.1f%%)\n",
                n_nonzero, length(weights), 100 * n_nonzero / length(weights)))
    cat(sprintf("  Max weight: %.4f\n", max(abs(weights))))
    cat(sprintf("  Mean absolute weight: %.4f\n", mean(abs(weights))))
  }
  cat("\n")

  return(p)
}


#' Plot CV Performance Comparison Between Coarse and Fine Grid Search
#'
#' @description Creates marginal CV performance plots comparing coarse and fine
#' grid search results. Shows how CV score changes with each tau parameter,
#' demonstrating the refinement gained from fine grid search.
#'
#' @param coarse_cv_results Data frame from MINERVA_cv_scca with coarse grid.
#'   Must contain columns: param1, param2, ..., paramN, mean_score
#' @param fine_cv_results Data frame from MINERVA_cv_scca with fine grid.
#'   Must contain columns: param1, param2, ..., paramN, mean_score
#' @param n_datasets Integer. Number of datasets (tau parameters).
#' @param param_names Optional character vector of parameter names for axis labels.
#'   If NULL, uses "tau1", "tau2", etc. (default = NULL)
#' @param metric_name Character string for y-axis label (default = "CV Score")
#' @param x_axis_type Character string controlling x-axis range. Options:
#'   "full" = fixed 0-1 range showing complete parameter space (default),
#'   "auto" = zoom to data range (fine grid region)
#'
#' @return ggplot2 object with faceted panels showing marginal CV performance
#'   curves for coarse and fine grids.
#'
#' @details
#' Creates a multi-panel plot with one panel per tau parameter. Each panel shows:
#' \itemize{
#'   \item **Coarse grid curve**: CV performance across coarse tau values
#'     (averaged over other tau parameters)
#'   \item **Fine grid curve**: CV performance in refined range
#'     (averaged over other tau parameters)
#'   \item **Markers**: Best value from fine grid highlighted
#' }
#'
#' Marginal plots are computed by:
#' \enumerate{
#'   \item Grouping by each tau parameter
#'   \item Averaging CV score across all other tau values
#'   \item Plotting the marginal relationship
#' }
#'
#' This visualization demonstrates:
#' \itemize{
#'   \item Initial broad exploration (coarse grid)
#'   \item Focused refinement (fine grid in promising region)
#'   \item Performance improvement from fine-tuning
#' }
#'
#' @export
#'
#' @examples
#' # After running coarse and fine CV
#'
#' # Full view showing complete 0-1 parameter space
#' p_full <- AETHER_plot_cv_comparison(
#'   coarse_cv_results = coarse_cv$cv_results,
#'   fine_cv_results = fine_cv$cv_results,
#'   n_datasets = 3,
#'   param_names = c("RNA-seq tau", "ATAC-seq tau", "CUT&TAG tau"),
#'   x_axis_type = "full"
#' )
#' print(p_full)
#'
#' # Zoomed view focusing on fine grid region
#' p_zoom <- AETHER_plot_cv_comparison(
#'   coarse_cv_results = coarse_cv$cv_results,
#'   fine_cv_results = fine_cv$cv_results,
#'   n_datasets = 3,
#'   param_names = c("RNA-seq tau", "ATAC-seq tau", "CUT&TAG tau"),
#'   x_axis_type = "auto"
#' )
#' print(p_zoom)
#'
AETHER_plot_cv_comparison <- function(coarse_cv_results,
                                       fine_cv_results,
                                       n_datasets,
                                       param_names = NULL,
                                       metric_name = "CV Score",
                                       x_axis_type = "full") {

  # Input validation
  if (!is.data.frame(coarse_cv_results) || !is.data.frame(fine_cv_results)) {
    stop("coarse_cv_results and fine_cv_results must be data frames")
  }

  if (!is.numeric(n_datasets) || n_datasets < 1) {
    stop("n_datasets must be a positive integer")
  }

  if (!x_axis_type %in% c("full", "auto")) {
    stop("x_axis_type must be either 'full' or 'auto'")
  }

  # Generate parameter names if not provided
  if (is.null(param_names)) {
    param_names <- paste0("tau", 1:n_datasets)
  }

  # Expected column names
  expected_cols <- c(paste0("param", 1:n_datasets), "mean_score")

  # Check that required columns exist
  if (!all(expected_cols %in% colnames(coarse_cv_results))) {
    stop("coarse_cv_results missing required columns: ",
         paste(setdiff(expected_cols, colnames(coarse_cv_results)), collapse = ", "))
  }
  if (!all(expected_cols %in% colnames(fine_cv_results))) {
    stop("fine_cv_results missing required columns: ",
         paste(setdiff(expected_cols, colnames(fine_cv_results)), collapse = ", "))
  }

  # Compute marginal scores for each parameter
  marginal_data_list <- list()

  for (i in 1:n_datasets) {
    param_col <- paste0("param", i)

    # Coarse marginals: group by param_i, average score (excluding NAs)
    coarse_marginal <- aggregate(
      mean_score ~ get(param_col),
      data = coarse_cv_results,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    colnames(coarse_marginal) <- c("param_value", "cv_score")
    coarse_marginal$grid_type <- "Coarse"
    coarse_marginal$parameter <- param_names[i]

    # Fine marginals: group by param_i, average score (excluding NAs)
    fine_marginal <- aggregate(
      mean_score ~ get(param_col),
      data = fine_cv_results,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    colnames(fine_marginal) <- c("param_value", "cv_score")
    fine_marginal$grid_type <- "Fine"
    fine_marginal$parameter <- param_names[i]

    # Combine coarse and fine for this parameter
    marginal_data_list[[i]] <- rbind(coarse_marginal, fine_marginal)
  }

  # Combine all parameters into one data frame
  plot_data <- do.call(rbind, marginal_data_list)
  plot_data$parameter <- factor(plot_data$parameter, levels = param_names)
  plot_data$grid_type <- factor(plot_data$grid_type, levels = c("Coarse", "Fine"))

  # Find best value from coarse grid for each parameter
  best_coarse_values <- data.frame()
  for (i in 1:n_datasets) {
    coarse_data <- plot_data[plot_data$parameter == param_names[i] &
                              plot_data$grid_type == "Coarse", ]
    if (nrow(coarse_data) > 0) {
      best_idx <- which.max(coarse_data$cv_score)
      best_coarse_values <- rbind(best_coarse_values, coarse_data[best_idx, ])
    }
  }

  # Find best value from fine grid for each parameter
  best_fine_values <- data.frame()
  for (i in 1:n_datasets) {
    fine_data <- plot_data[plot_data$parameter == param_names[i] &
                            plot_data$grid_type == "Fine", ]
    if (nrow(fine_data) > 0) {
      best_idx <- which.max(fine_data$cv_score)
      best_fine_values <- rbind(best_fine_values, fine_data[best_idx, ])
    }
  }

  # Create the plot
  p <- ggplot(plot_data, aes(x = param_value, y = cv_score,
                              color = grid_type, linetype = grid_type)) +
    # Lines connecting points
    geom_line(linewidth = 1, alpha = 0.8) +
    # Points at each measurement
    geom_point(size = 2.5, alpha = 0.9) +
    # Highlight best values from coarse grid
    {if (nrow(best_coarse_values) > 0) {
      geom_point(data = best_coarse_values,
                 aes(x = param_value, y = cv_score),
                 color = "black", size = 4, shape = 21,
                 fill = "#0072B2", stroke = 1.5,
                 inherit.aes = FALSE)
    }} +
    # Highlight best values from fine grid
    {if (nrow(best_fine_values) > 0) {
      geom_point(data = best_fine_values,
                 aes(x = param_value, y = cv_score),
                 color = "black", size = 4, shape = 21,
                 fill = "#E69F00", stroke = 1.5,
                 inherit.aes = FALSE)
    }} +
    # Color and line type scheme
    scale_color_manual(
      values = c("Coarse" = "#0072B2", "Fine" = "#E69F00"),
      name = "Grid Type"
    ) +
    scale_linetype_manual(
      values = c("Coarse" = "dashed", "Fine" = "solid"),
      name = "Grid Type"
    ) +
    # Conditional x-axis scaling based on x_axis_type
    {if (x_axis_type == "full") {
      scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.1)
      )
    }} +
    # Facet by parameter with conditional scales
    facet_wrap(~parameter, ncol = n_datasets,
               scales = if (x_axis_type == "full") "free_y" else "free") +
    # Labels and theme
    labs(
      title = "Cross-Validation Performance: Coarse vs Fine Grid",
      subtitle = "Marginal CV scores averaged across other parameters. Best coarse (blue) and fine (gold) values marked.",
      x = "Parameter Value",
      y = metric_name
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11)
    )

  # Print summary
  cat("\n=== CV Comparison Summary ===\n")
  cat(sprintf("Number of parameters: %d\n", n_datasets))
  cat(sprintf("Coarse grid combinations: %d\n", nrow(coarse_cv_results)))
  cat(sprintf("Fine grid combinations: %d\n", nrow(fine_cv_results)))

  cat("\nBest values from coarse grid:\n")
  for (i in 1:n_datasets) {
    if (i <= nrow(best_coarse_values)) {
      cat(sprintf("  %s: %.3f (CV score: %.4f)\n",
                  param_names[i],
                  best_coarse_values$param_value[i],
                  best_coarse_values$cv_score[i]))
    }
  }

  cat("\nBest values from fine grid:\n")
  for (i in 1:n_datasets) {
    if (i <= nrow(best_fine_values)) {
      cat(sprintf("  %s: %.3f (CV score: %.4f)\n",
                  param_names[i],
                  best_fine_values$param_value[i],
                  best_fine_values$cv_score[i]))
    }
  }
  cat("\n")

  return(p)
}


#' Plot Best CV Score Comparison Between Coarse and Fine Grid
#'
#' @description Creates a simple bar chart comparing the best CV scores from
#' coarse and fine grid search. Designed for clarity in presentations and
#' publications for non-technical audiences.
#'
#' @param coarse_cv_results Data frame from MINERVA_cv_scca with coarse grid.
#'   Must contain column: mean_score
#' @param fine_cv_results Data frame from MINERVA_cv_scca with fine grid.
#'   Must contain column: mean_score
#' @param metric_name Character string for y-axis label (default = "CV Score")
#'
#' @return ggplot2 object showing bar chart of best scores from each grid.
#'
#' @details
#' Creates a simple two-bar comparison showing:
#' \itemize{
#'   \item Best CV score from coarse grid search
#'   \item Best CV score from fine grid search (after refinement)
#' }
#'
#' This visualization clearly demonstrates whether the fine grid refinement
#' improved performance over the initial coarse search, without the complexity
#' of marginal plots.
#'
#' The function also prints:
#' \itemize{
#'   \item Best score from each grid
#'   \item Absolute improvement (fine - coarse)
#'   \item Percentage improvement
#' }
#'
#' @export
#'
#' @examples
#' # After running coarse and fine CV
#' p <- AETHER_plot_best_cv_scores(
#'   coarse_cv_results = coarse_cv$cv_results,
#'   fine_cv_results = fine_cv$cv_results,
#'   metric_name = "Correlation"
#' )
#' print(p)
#'
AETHER_plot_best_cv_scores <- function(coarse_cv_results,
                                        fine_cv_results,
                                        metric_name = "CV Score") {

  # Input validation
  if (!is.data.frame(coarse_cv_results) || !is.data.frame(fine_cv_results)) {
    stop("coarse_cv_results and fine_cv_results must be data frames")
  }

  if (!"mean_score" %in% colnames(coarse_cv_results)) {
    stop("coarse_cv_results must contain 'mean_score' column")
  }

  if (!"mean_score" %in% colnames(fine_cv_results)) {
    stop("fine_cv_results must contain 'mean_score' column")
  }

  # Extract best scores
  best_coarse <- max(coarse_cv_results$mean_score, na.rm = TRUE)
  best_fine <- max(fine_cv_results$mean_score, na.rm = TRUE)

  # Create data frame for plotting
  plot_data <- data.frame(
    grid_type = factor(c("Coarse Grid", "Fine Grid"),
                       levels = c("Coarse Grid", "Fine Grid")),
    best_score = c(best_coarse, best_fine)
  )

  # Calculate improvement
  improvement <- best_fine - best_coarse
  pct_improvement <- (improvement / best_coarse) * 100

  # Create the plot
  p <- ggplot(plot_data, aes(x = grid_type, y = best_score, fill = grid_type)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
    # Add value labels on top of bars
    geom_text(aes(label = sprintf("%.4f", best_score)),
              vjust = -0.5, size = 5, fontface = "bold") +
    # Color scheme matching the marginal plot
    scale_fill_manual(
      values = c("Coarse Grid" = "#0072B2", "Fine Grid" = "#E69F00")
    ) +
    # Labels and theme
    labs(
      title = "Best CV Score: Coarse vs Fine Grid Search",
      subtitle = sprintf("Improvement: %.4f (%.2f%%)", improvement, pct_improvement),
      x = NULL,
      y = metric_name
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # Expand y-axis slightly to make room for labels
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  # Print summary
  cat("\n=== Best CV Score Comparison ===\n")
  cat(sprintf("Coarse grid search: %.4f\n", best_coarse))
  cat(sprintf("Fine grid search:   %.4f\n", best_fine))
  cat(sprintf("Improvement:        %.4f (%.2f%%)\n", improvement, pct_improvement))
  cat("\n")

  return(p)
}