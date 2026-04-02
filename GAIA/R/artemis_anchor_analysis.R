#' Artemis - Anchor-Based Analysis Functions
#'
#' @description Functions for statistical analysis of signal at anchor regions,
#' comparing conditions (e.g., WT vs KO) at predefined genomic loci.


#' Global Shift Analysis at Anchor Regions
#'
#' @description Tests whether signal at anchor regions is systematically
#' different between two conditions (e.g., WT vs KO). Useful for answering
#' "Is there an overall effect at these sites?"
#'
#' @param quant_result A quantification result from ELEUTHIA_quantify_bed(),
#'   containing counts matrix, annotation, and targets with group column.
#' @param group_col Character. Column name in targets containing group labels
#'   (default = "group").
#' @param group_a Character. First group label (default = "WT"). This is the
#'   reference/baseline group.
#' @param group_b Character. Second group label (default = "SMUG1_KO"). This is
#'   the comparison group.
#' @param summary_method Character. How to summarize signal per region per group:
#'   "mean" (default) or "median".
#' @param test Character. Statistical test to use: "wilcoxon" (default, paired
#'   Wilcoxon signed-rank test) or "t.test" (paired t-test).
#' @param pseudocount Numeric. Added to counts before log transformation to
#'   avoid log(0) (default = 1).
#' @param verbose Logical. Print results summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{per_region}{Data.frame with per-region statistics: region_id,
#'     mean_a, mean_b, log2FC, raw_diff}
#'   \item{global_test}{List with test statistic, p.value, and test name}
#'   \item{summary}{List with summary statistics: n_regions, n_up, n_down,
#'     median_log2FC, mean_log2FC}
#'   \item{groups}{Character vector of group names used (c(group_a, group_b))}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates mean/median signal per region for each group
#'   \item Computes log2 fold change (group_b / group_a) per region
#'   \item Tests for systematic shift using paired test across all regions
#'   \item Returns per-region statistics and global test result
#' }
#'
#' A significant p-value indicates that signal at anchor regions is
#' systematically higher or lower in group_b compared to group_a.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare ATAC signal at ChIP anchors between WT and KO
#' result <- ARTEMIS_global_shift_test(atac_at_anchors,
#'                                      group_a = "WT",
#'                                      group_b = "SMUG1_KO")
#'
#' # View summary
#' result$summary
#'
#' # Get per-region fold changes
#' head(result$per_region)
#'
#' }
ARTEMIS_global_shift_test <- function(quant_result,
                                       group_col = "group",
                                       group_a = "WT",
                                       group_b = "SMUG1_KO",
                                       summary_method = "mean",
                                       test = "wilcoxon",
                                       pseudocount = 1,
                                       verbose = TRUE) {

  # Validate inputs
  if (!is.list(quant_result) || !all(c("counts", "targets") %in% names(quant_result))) {
    stop("quant_result must be a list with 'counts' and 'targets' elements")
  }

  counts <- quant_result$counts
  targets <- quant_result$targets

  if (!group_col %in% colnames(targets)) {
    stop("group_col '", group_col, "' not found in targets")
  }

  groups <- targets[[group_col]]

  if (!group_a %in% groups) {
    stop("group_a '", group_a, "' not found in targets$", group_col,
         ". Available: ", paste(unique(groups), collapse = ", "))
  }
  if (!group_b %in% groups) {
    stop("group_b '", group_b, "' not found in targets$", group_col,
         ". Available: ", paste(unique(groups), collapse = ", "))
  }

  summary_method <- match.arg(summary_method, c("mean", "median"))
  test <- match.arg(test, c("wilcoxon", "t.test"))

  # Get sample indices for each group

  idx_a <- which(groups == group_a)
  idx_b <- which(groups == group_b)

  n_a <- length(idx_a)
  n_b <- length(idx_b)
  n_regions <- nrow(counts)

  if (verbose) {
    cat("[ARTEMIS] Global shift analysis:\n")
    cat("    Group A (", group_a, "): ", n_a, " samples\n", sep = "")
    cat("    Group B (", group_b, "): ", n_b, " samples\n", sep = "")
    cat("    Regions: ", n_regions, "\n", sep = "")
    cat("    Summary method: ", summary_method, "\n", sep = "")
    cat("    Test: ", test, "\n\n", sep = "")
  }

  # Calculate per-region summary for each group
  summary_fn <- if (summary_method == "mean") rowMeans else function(x) apply(x, 1, median)

  signal_a <- summary_fn(counts[, idx_a, drop = FALSE])
  signal_b <- summary_fn(counts[, idx_b, drop = FALSE])

  # Calculate log2 fold change (with pseudocount)
  log2FC <- log2((signal_b + pseudocount) / (signal_a + pseudocount))
  raw_diff <- signal_b - signal_a

  # Create per-region results
  region_ids <- if (!is.null(rownames(counts))) {
    rownames(counts)
  } else if (!is.null(quant_result$annotation$peak_id)) {
    quant_result$annotation$peak_id
  } else {
    paste0("region_", seq_len(n_regions))
  }

  per_region <- data.frame(
    region_id = region_ids,
    mean_a = signal_a,
    mean_b = signal_b,
    log2FC = log2FC,
    raw_diff = raw_diff,
    stringsAsFactors = FALSE
  )
  colnames(per_region)[2:3] <- paste0(summary_method, "_", c(group_a, group_b))

  # Global statistical test (paired across regions)
  if (test == "wilcoxon") {
    test_result <- wilcox.test(signal_b, signal_a, paired = TRUE)
    test_name <- "Paired Wilcoxon signed-rank test"
  } else {
    test_result <- t.test(signal_b, signal_a, paired = TRUE)
    test_name <- "Paired t-test"
  }

  # Summary statistics
  n_up <- sum(log2FC > 0, na.rm = TRUE)
  n_down <- sum(log2FC < 0, na.rm = TRUE)
  n_unchanged <- sum(log2FC == 0, na.rm = TRUE)

  summary_stats <- list(
    n_regions = n_regions,
    n_up = n_up,
    n_down = n_down,
    n_unchanged = n_unchanged,
    median_log2FC = median(log2FC, na.rm = TRUE),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    sd_log2FC = sd(log2FC, na.rm = TRUE)
  )

  if (verbose) {
    cat("[ARTEMIS] Results:\n")
    cat("    Regions with increased signal in ", group_b, ": ", n_up,
        " (", round(100 * n_up / n_regions, 1), "%)\n", sep = "")
    cat("    Regions with decreased signal in ", group_b, ": ", n_down,
        " (", round(100 * n_down / n_regions, 1), "%)\n", sep = "")
    cat("    Median log2FC: ", round(summary_stats$median_log2FC, 3), "\n", sep = "")
    cat("    Mean log2FC: ", round(summary_stats$mean_log2FC, 3), "\n", sep = "")
    cat("\n")
    cat(" Global test (", test_name, "):\n", sep = "")
    cat("    p-value: ", format.pval(test_result$p.value, digits = 3), "\n", sep = "")

    if (test_result$p.value < 0.05) {
      direction <- if (summary_stats$median_log2FC > 0) "INCREASED" else "DECREASED"
      cat("    Interpretation: Signal at anchor regions is significantly ",
          direction, " in ", group_b, " vs ", group_a, "\n", sep = "")
    } else {
      cat("    Interpretation: No significant systematic shift detected\n")
    }
  }

  return(list(
    per_region = per_region,
    global_test = list(
      statistic = test_result$statistic,
      p.value = test_result$p.value,
      test = test_name
    ),
    summary = summary_stats,
    groups = c(group_a, group_b)
  ))
}


#' Plot Global Shift Results
#'
#' @description Visualizes the results of ARTEMIS_global_shift_test() with
#' a density plot of log2 fold changes and a paired dot plot.
#'
#' @param shift_result Result from ARTEMIS_global_shift_test().
#' @param plot_type Character. Type of plot: "density" (log2FC distribution),
#'   "paired" (paired dot plot), or "both" (default).
#' @param highlight_threshold Numeric. log2FC threshold for highlighting
#'   strongly changed regions (default = 1, i.e., 2-fold).
#' @param title Character. Plot title (default = "Global Shift Analysis").
#'
#' @return A ggplot object (or list of ggplot objects if plot_type = "both").
#'
#' @export
#'
ARTEMIS_plot_global_shift <- function(shift_result,
                                       plot_type = "both",
                                       highlight_threshold = 1,
                                       title = "Global Shift Analysis") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  plot_type <- match.arg(plot_type, c("density", "paired", "both"))

  per_region <- shift_result$per_region
  groups <- shift_result$groups
  p_value <- shift_result$global_test$p.value
  median_fc <- shift_result$summary$median_log2FC

  plots <- list()

  # Density plot of log2FC
  if (plot_type %in% c("density", "both")) {
    p_density <- ggplot2::ggplot(per_region, ggplot2::aes(x = log2FC)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_vline(xintercept = median_fc, linetype = "solid",
                          color = "red", linewidth = 1) +
      ggplot2::annotate("text",
                        x = median_fc,
                        y = Inf,
                        label = paste0("median = ", round(median_fc, 2)),
                        vjust = 2, hjust = -0.1, color = "red", size = 3.5) +
      ggplot2::labs(
        title = title,
        subtitle = paste0("p = ", format.pval(p_value, digits = 2),
                          " (", shift_result$global_test$test, ")"),
        x = paste0("log2 Fold Change (", groups[2], " / ", groups[1], ")"),
        y = "Density"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(color = "gray40")
      )

    plots$density <- p_density
  }

  # Paired dot plot
  if (plot_type %in% c("paired", "both")) {
    # Reshape for paired plot
    col_a <- grep(groups[1], colnames(per_region), value = TRUE)[1]
    col_b <- grep(groups[2], colnames(per_region), value = TRUE)[1]

    paired_data <- data.frame(
      region = rep(per_region$region_id, 2),
      group = rep(groups, each = nrow(per_region)),
      signal = c(per_region[[col_a]], per_region[[col_b]])
    )
    paired_data$group <- factor(paired_data$group, levels = groups)

    p_paired <- ggplot2::ggplot(paired_data,
                                 ggplot2::aes(x = group, y = log2(signal + 1))) +
      ggplot2::geom_line(ggplot2::aes(group = region),
                         alpha = 0.3, color = "gray50") +
      ggplot2::geom_point(ggplot2::aes(color = group), alpha = 0.6, size = 2) +
      ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.5,
                            color = "red", linewidth = 0.8) +
      ggplot2::scale_color_manual(values = c("steelblue", "darkorange")) +
      ggplot2::labs(
        title = paste0("Signal at Anchor Regions (n = ", nrow(per_region), ")"),
        subtitle = paste0("p = ", format.pval(p_value, digits = 2)),
        x = "Group",
        y = "log2(signal + 1)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )

    plots$paired <- p_paired
  }

  if (plot_type == "both") {
    return(plots)
  } else {
    return(plots[[1]])
  }
}
