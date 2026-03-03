# GAIA/Eleuthia/activity_export.R
# Export decoupleR activity inference results
#
# Saves decoupler_result and decoupler_comparison objects to CSV/RDS files.
# Compatible with DEMETER_load_activity() for reloading.
#
# Requires: export_shared.R (for .save_ggplot, .save_pheatmap, .plot_paths)


#' Export decoupleR activity results
#'
#' Saves activity inference results (single method or multi-method comparison)
#' to CSV and RDS files. Handles both \code{decoupler_result} and
#' \code{decoupler_comparison} objects.
#'
#' @param result A \code{decoupler_result} or \code{decoupler_comparison} object
#'   from ARTEMIS_run_decoupler(), ARTEMIS_infer_tf_activity(),
#'   ARTEMIS_infer_pathway_activity(), or ARTEMIS_decoupler_compare_methods().
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param prefix Character. Prefix for output filenames. Default: "activity".
#' @param save_plots Logical. Save diagnostic plots. Default: TRUE.
#' @param plot_format Character. Plot format: "png", "pdf", or "both".
#'   Default: "png".
#' @param top_n Integer. Number of top sources to show in heatmap and bar
#'   plots. Default: 30.
#' @param save_rds Logical. Save full R objects as RDS for exact reloading.
#'   Default: TRUE.
#' @param save_long Logical. Save full long-format results from decoupleR.
#'   Default: TRUE. Can produce large files for comparisons.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' For a \code{decoupler_result} (single method), exports:
#' \itemize{
#'   \item Activity matrix (sources x samples) as CSV
#'   \item P-value matrix as CSV (if available)
#'   \item Full long-format results as CSV (optional)
#'   \item Metadata text file (method, network info, gene overlap)
#'   \item Plots: activity heatmap, top activities bar chart
#'   \item RDS file of the full object (optional)
#' }
#'
#' For a \code{decoupler_comparison} (multi-method), exports:
#' \itemize{
#'   \item Per-method activity matrices in subdirectory
#'   \item Consensus scores as CSV
#'   \item Summary table (per-source stats across methods) as CSV
#'   \item Method correlation matrix as CSV
#'   \item Method statistics as CSV
#'   \item Metadata text file
#'   \item Plots: activity heatmap, top activities, method correlation, method agreement
#'   \item RDS file of the full object (optional)
#' }
#'
#' Saved RDS files can be reloaded with \code{DEMETER_load_activity()}.
#'
#' @examples
#' \dontrun{
#' # Export single-method result
#' tf_result <- ARTEMIS_infer_tf_activity(expr_matrix, method = "ulm")
#' ELEUTHIA_export_activity_results(tf_result, "results/tf_activity/")
#'
#' # Export multi-method comparison
#' comparison <- ARTEMIS_decoupler_compare_methods(
#'   expr_matrix, network, methods = c("ulm", "mlm", "wsum")
#' )
#' ELEUTHIA_export_activity_results(comparison, "results/tf_comparison/")
#'
#' # Reload later
#' tf_result <- DEMETER_load_activity("results/tf_activity/activity_result.rds")
#'
#' }
#' @export
ELEUTHIA_export_activity_results <- function(result,
                                              output_dir,
                                              prefix = "activity",
                                              save_plots = TRUE,
                                              plot_format = "png",
                                              top_n = 30,
                                              save_rds = TRUE,
                                              save_long = TRUE,
                                              verbose = TRUE) {

  # Validate input
  is_single <- inherits(result, "decoupler_result")
  is_comparison <- inherits(result, "decoupler_comparison")

  if (!is_single && !is_comparison) {
    stop("'result' must be a decoupler_result or decoupler_comparison object.\n",
         "  Use ARTEMIS_run_decoupler(), ARTEMIS_infer_tf_activity(), ",
         "ARTEMIS_infer_pathway_activity(), or ARTEMIS_decoupler_compare_methods().")
  }

  if (verbose) cat("=== Exporting Activity Results ===\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created directory:", output_dir, "\n")
  }

  files_created <- character(0)

  if (is_single) {
    files_created <- .export_single_result(result, output_dir, prefix,
                                            save_plots, plot_format, top_n,
                                            save_rds, save_long, verbose)
  } else {
    files_created <- .export_comparison_result(result, output_dir, prefix,
                                                save_plots, plot_format, top_n,
                                                save_rds, save_long, verbose)
  }

  if (verbose) {
    cat("\n--- Export Summary ---\n")
    cat("Total files created:", length(files_created), "\n")
    cat("Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}


# ==============================================================================
# INTERNAL: Export single decoupler_result
# ==============================================================================

.export_single_result <- function(result, output_dir, prefix,
                                   save_plots, plot_format, top_n,
                                   save_rds, save_long, verbose) {
  files_created <- character(0)

  if (verbose) {
    cat("Type: Single method result\n")
    cat("Method:", result$method_name, "(", result$method, ")\n")
    cat("Sources:", result$n_sources, "| Samples:", result$n_samples, "\n")
  }

  # --------------------------------------------------------------------------
  # Activity matrix
  # --------------------------------------------------------------------------
  act_file <- file.path(output_dir, paste0(prefix, "_activities.csv"))
  act_df <- as.data.frame(result$activities)
  act_df$source <- rownames(act_df)
  act_df <- act_df[, c("source", setdiff(names(act_df), "source"))]
  write.csv(act_df, act_file, row.names = FALSE)
  files_created <- c(files_created, act_file)
  if (verbose) cat("  Activities:", act_file, "\n")

  # --------------------------------------------------------------------------
  # P-value matrix
  # --------------------------------------------------------------------------
  if (!is.null(result$pvalues)) {
    pval_file <- file.path(output_dir, paste0(prefix, "_pvalues.csv"))
    pval_df <- as.data.frame(result$pvalues)
    pval_df$source <- rownames(pval_df)
    pval_df <- pval_df[, c("source", setdiff(names(pval_df), "source"))]
    write.csv(pval_df, pval_file, row.names = FALSE)
    files_created <- c(files_created, pval_file)
    if (verbose) cat("  P-values:", pval_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Full long-format results
  # --------------------------------------------------------------------------
  if (save_long && !is.null(result$results_long)) {
    long_file <- file.path(output_dir, paste0(prefix, "_results_long.csv"))
    write.csv(as.data.frame(result$results_long), long_file, row.names = FALSE)
    files_created <- c(files_created, long_file)
    if (verbose) cat("  Long-format results:", long_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Metadata
  # --------------------------------------------------------------------------
  meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))
  meta_lines <- c(
    "=== DecoupleR Activity Results ===",
    paste("Method:", result$method_name, "(", result$method, ")"),
    paste("Sources:", result$n_sources),
    paste("Samples:", result$n_samples),
    paste("Gene overlap:", result$gene_overlap),
    paste("Min source size:", result$minsize),
    ""
  )

  # Analysis type
  if (!is.null(result$analysis_type)) {
    meta_lines <- c(meta_lines, paste("Analysis type:", result$analysis_type))
  }

  # Network info
  if (!is.null(result$network_info)) {
    ni <- result$network_info
    meta_lines <- c(meta_lines,
      "",
      "--- Network Info ---",
      paste("Database:", ifelse(is.null(ni$database), "custom", ni$database)),
      paste("Organism:", ifelse(is.null(ni$organism), "unknown", ni$organism)),
      paste("Network type:", ifelse(is.null(ni$network_type), "unknown", ni$network_type)),
      paste("Network sources:", ni$n_sources),
      paste("Network targets:", ni$n_targets)
    )
  }

  meta_lines <- c(meta_lines, "", paste("Exported:", Sys.time()))
  writeLines(meta_lines, meta_file)
  files_created <- c(files_created, meta_file)
  if (verbose) cat("  Metadata:", meta_file, "\n")

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  if (save_plots) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }

    if (verbose) cat("\nGenerating plots...\n")

    # Activity heatmap (pheatmap - needs special save)
    tryCatch({
      .save_pheatmap(
        function() AETHER_plot_activity_heatmap(result, top_n = top_n),
        plot_dir, "activity_heatmap", plot_format,
        width = max(8, result$n_samples * 0.3 + 4),
        height = max(6, min(top_n, result$n_sources) * 0.25 + 3)
      )
      files_created <- c(files_created,
                          .plot_paths(plot_dir, "activity_heatmap", plot_format))
      if (verbose) cat("  Activity heatmap\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create activity heatmap -", e$message, "\n")
    })

    # Top activities bar chart (ggplot)
    tryCatch({
      p <- AETHER_plot_top_activities(result, top_n = min(top_n, 20))
      plot_files <- .save_ggplot(p, plot_dir, "top_activities", plot_format,
                                  width = 8, height = max(6, min(top_n, 20) * 0.3 + 2))
      files_created <- c(files_created, plot_files)
      if (verbose) cat("  Top activities\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create top activities plot -", e$message, "\n")
    })

    if (verbose) cat("  Plots saved to:", plot_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, "_result.rds"))
    saveRDS(result, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("  RDS:", rds_file, "\n")
  }

  return(files_created)
}


# ==============================================================================
# INTERNAL: Export decoupler_comparison
# ==============================================================================

.export_comparison_result <- function(result, output_dir, prefix,
                                       save_plots, plot_format, top_n,
                                       save_rds, save_long, verbose) {
  files_created <- character(0)

  if (verbose) {
    cat("Type: Multi-method comparison\n")
    cat("Methods:", paste(result$methods, collapse = ", "), "\n")
    cat("Common sources:", length(result$common_sources), "\n")
    cat("Common samples:", length(result$common_samples), "\n")
  }

  # --------------------------------------------------------------------------
  # Per-method activity matrices (in subdirectory)
  # --------------------------------------------------------------------------
  methods_dir <- file.path(output_dir, paste0(prefix, "_methods"))
  if (!dir.exists(methods_dir)) {
    dir.create(methods_dir)
  }

  for (m in names(result$results)) {
    method_result <- result$results[[m]]

    # Activity matrix
    act_file <- file.path(methods_dir, paste0(m, "_activities.csv"))
    act_df <- as.data.frame(method_result$activities)
    act_df$source <- rownames(act_df)
    act_df <- act_df[, c("source", setdiff(names(act_df), "source"))]
    write.csv(act_df, act_file, row.names = FALSE)
    files_created <- c(files_created, act_file)

    # P-values
    if (!is.null(method_result$pvalues)) {
      pval_file <- file.path(methods_dir, paste0(m, "_pvalues.csv"))
      pval_df <- as.data.frame(method_result$pvalues)
      pval_df$source <- rownames(pval_df)
      pval_df <- pval_df[, c("source", setdiff(names(pval_df), "source"))]
      write.csv(pval_df, pval_file, row.names = FALSE)
      files_created <- c(files_created, pval_file)
    }

    # Long-format results per method
    if (save_long && !is.null(method_result$results_long)) {
      long_file <- file.path(methods_dir, paste0(m, "_results_long.csv"))
      write.csv(as.data.frame(method_result$results_long), long_file, row.names = FALSE)
      files_created <- c(files_created, long_file)
    }
  }

  if (verbose) cat("  Per-method results:", methods_dir, "/\n")

  # --------------------------------------------------------------------------
  # Consensus scores
  # --------------------------------------------------------------------------
  if (!is.null(result$consensus)) {
    cons_file <- file.path(output_dir, paste0(prefix, "_consensus.csv"))
    cons_df <- as.data.frame(result$consensus)
    cons_df$source <- rownames(cons_df)
    cons_df <- cons_df[, c("source", setdiff(names(cons_df), "source"))]
    write.csv(cons_df, cons_file, row.names = FALSE)
    files_created <- c(files_created, cons_file)
    if (verbose) cat("  Consensus scores:", cons_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Summary table (per-source stats across methods)
  # --------------------------------------------------------------------------
  if (!is.null(result$summary)) {
    summary_file <- file.path(output_dir, paste0(prefix, "_summary.csv"))
    write.csv(result$summary, summary_file, row.names = FALSE)
    files_created <- c(files_created, summary_file)
    if (verbose) cat("  Summary:", summary_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Method correlations
  # --------------------------------------------------------------------------
  if (!is.null(result$correlations)) {
    cor_file <- file.path(output_dir, paste0(prefix, "_method_correlations.csv"))
    cor_df <- as.data.frame(result$correlations)
    cor_df$method <- rownames(cor_df)
    cor_df <- cor_df[, c("method", setdiff(names(cor_df), "method"))]
    write.csv(cor_df, cor_file, row.names = FALSE)
    files_created <- c(files_created, cor_file)
    if (verbose) cat("  Method correlations:", cor_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Method statistics
  # --------------------------------------------------------------------------
  if (!is.null(result$method_stats)) {
    stats_file <- file.path(output_dir, paste0(prefix, "_method_stats.csv"))
    write.csv(result$method_stats, stats_file, row.names = FALSE)
    files_created <- c(files_created, stats_file)
    if (verbose) cat("  Method stats:", stats_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Metadata
  # --------------------------------------------------------------------------
  meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))
  meta_lines <- c(
    "=== DecoupleR Method Comparison Results ===",
    paste("Methods:", paste(result$methods, collapse = ", ")),
    paste("Number of methods:", result$n_methods),
    paste("Common sources:", length(result$common_sources)),
    paste("Common samples:", length(result$common_samples)),
    ""
  )

  # Method correlation summary
  if (!is.null(result$correlations)) {
    lower_tri <- result$correlations[lower.tri(result$correlations)]
    meta_lines <- c(meta_lines,
      "--- Method Agreement ---",
      paste("Mean inter-method correlation:", round(mean(lower_tri), 3)),
      paste("Min inter-method correlation:", round(min(lower_tri), 3)),
      paste("Max inter-method correlation:", round(max(lower_tri), 3)),
      ""
    )
  }

  # Network info
  if (!is.null(result$network_info)) {
    ni <- result$network_info
    meta_lines <- c(meta_lines,
      "--- Network Info ---",
      paste("Database:", ifelse(is.null(ni$database), "custom", ni$database)),
      paste("Organism:", ifelse(is.null(ni$organism), "unknown", ni$organism)),
      paste("Network type:", ifelse(is.null(ni$network_type), "unknown", ni$network_type)),
      paste("Network sources:", ni$n_sources),
      paste("Network targets:", ni$n_targets)
    )
  }

  meta_lines <- c(meta_lines, "", paste("Exported:", Sys.time()))
  writeLines(meta_lines, meta_file)
  files_created <- c(files_created, meta_file)
  if (verbose) cat("  Metadata:", meta_file, "\n")

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  if (save_plots) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }

    if (verbose) cat("\nGenerating plots...\n")

    # Activity heatmap (consensus scores, pheatmap)
    tryCatch({
      .save_pheatmap(
        function() AETHER_plot_activity_heatmap(result, top_n = top_n),
        plot_dir, "activity_heatmap", plot_format,
        width = max(8, length(result$common_samples) * 0.3 + 4),
        height = max(6, min(top_n, length(result$common_sources)) * 0.25 + 3)
      )
      files_created <- c(files_created,
                          .plot_paths(plot_dir, "activity_heatmap", plot_format))
      if (verbose) cat("  Activity heatmap (consensus)\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create activity heatmap -", e$message, "\n")
    })

    # Top activities bar chart (from first method as representative, ggplot)
    tryCatch({
      first_method <- result$results[[1]]
      p <- AETHER_plot_top_activities(first_method, top_n = min(top_n, 20))
      plot_files <- .save_ggplot(p, plot_dir, "top_activities", plot_format,
                                  width = 8, height = max(6, min(top_n, 20) * 0.3 + 2))
      files_created <- c(files_created, plot_files)
      if (verbose) cat("  Top activities (", result$methods[1], ")\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create top activities plot -", e$message, "\n")
    })

    # Method correlation heatmap (pheatmap)
    tryCatch({
      n_methods <- result$n_methods
      .save_pheatmap(
        function() AETHER_plot_method_correlation(result),
        plot_dir, "method_correlation", plot_format,
        width = max(5, n_methods * 1.2 + 2),
        height = max(5, n_methods * 1.0 + 2)
      )
      files_created <- c(files_created,
                          .plot_paths(plot_dir, "method_correlation", plot_format))
      if (verbose) cat("  Method correlation\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create method correlation plot -", e$message, "\n")
    })

    # Method agreement (ggplot)
    tryCatch({
      p <- AETHER_plot_method_agreement(result, top_n = top_n)
      plot_files <- .save_ggplot(p, plot_dir, "method_agreement", plot_format,
                                  width = 10, height = max(6, top_n * 0.25 + 2))
      files_created <- c(files_created, plot_files)
      if (verbose) cat("  Method agreement\n")
    }, error = function(e) {
      if (verbose) cat("  Warning: Could not create method agreement plot -", e$message, "\n")
    })

    if (verbose) cat("  Plots saved to:", plot_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, "_comparison.rds"))
    saveRDS(result, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("  RDS:", rds_file, "\n")
  }

  return(files_created)
}


# ==============================================================================
# NOTE: Plot save helpers (.save_ggplot, .save_pheatmap, .plot_paths)
# are defined in export_shared.R — source that file before this one.
# ==============================================================================
