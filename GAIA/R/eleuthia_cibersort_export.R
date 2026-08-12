# GAIA/Eleuthia/cibersort_export.R
# Export CIBERSORT deconvolution results
#
# Saves artemis_cibersort objects to CSV/RDS files with optional plots.
# Follows the same pattern as activity_export.R.
#
# Requires: activity_export.R (for .save_ggplot, .save_pheatmap helpers)
#           cibersort_plots.R (for visualization functions)


#' Export CIBERSORT deconvolution results
#'
#' Saves cell type proportions, quality diagnostics, metadata, and optional
#' plots from an \code{artemis_cibersort} object.
#'
#' @param result An \code{artemis_cibersort} object from
#'   \code{ARTEMIS_cibersort()}.
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param prefix Character. Prefix for output filenames. Default: "cibersort".
#' @param group_by Optional character/factor vector of condition labels per
#'   sample. Used for grouped boxplot. Supports named vectors (names = groups,
#'   values = colors). Default: NULL.
#' @param save_plots Logical. Save diagnostic plots. Default: TRUE.
#' @param plot_format Character. Plot format: "png", "pdf", or "both".
#'   Default: "png".
#' @param top_n Integer. Top N cell types for boxplot (by mean proportion).
#'   Default: NULL (show all non-zero).
#' @param save_rds Logical. Save full R object as RDS for reloading.
#'   Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Exports the following files:
#' \itemize{
#'   \item Proportions matrix (samples x cell types) as CSV
#'   \item Diagnostics (correlation, RMSE, p-values) as CSV
#'   \item Metadata text file (parameters, quality summary)
#'   \item Plots: stacked bar, heatmap, boxplot (if group_by provided)
#'   \item RDS file of the full object (optional)
#' }
#'
#' @examples
#' \dontrun{
#' ELEUTHIA_export_cibersort_results(cib_result, "results/cibersort/")
#'
#' # With grouped boxplot
#' ELEUTHIA_export_cibersort_results(
#'   cib_result, "results/cibersort/",
#'   group_by = sample_conditions
#' )
#'
#' }
#' @export
ELEUTHIA_export_cibersort_results <- function(result,
                                               output_dir,
                                               prefix = "cibersort",
                                               group_by = NULL,
                                               save_plots = TRUE,
                                               plot_format = "png",
                                               top_n = NULL,
                                               save_rds = TRUE,
                                               verbose = TRUE) {

  if (!inherits(result, "artemis_cibersort")) {
    stop("'result' must be an artemis_cibersort object from ARTEMIS_cibersort().",
         call. = FALSE)
  }

  if (verbose) cat("[ELEUTHIA] Exporting CIBERSORT Results \n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("    Created directory:", output_dir, "\n")
  }

  files_created <- character(0)

  if (verbose) {
    cat("    Samples:", result$metadata$n_samples, "\n")
    cat("    Cell types:", result$metadata$n_cell_types, "\n")
    cat("    Permutations:", result$metadata$perm, "\n")
    cat("    QN:", result$metadata$QN, "\n")
  }

  # --------------------------------------------------------------------------
  # Proportions CSV
  # --------------------------------------------------------------------------
  prop_file <- file.path(output_dir, paste0(prefix, "_proportions.csv"))
  prop_df <- as.data.frame(result$proportions)
  prop_df$sample <- rownames(prop_df)
  prop_df <- prop_df[, c("sample", setdiff(names(prop_df), "sample"))]
  write.csv(prop_df, prop_file, row.names = FALSE)
  files_created <- c(files_created, prop_file)
  if (verbose) cat("        Proportions:", prop_file, "\n")

  # --------------------------------------------------------------------------
  # Diagnostics CSV (correlation, RMSE, p-values)
  # --------------------------------------------------------------------------
  diag_df <- data.frame(
    sample = rownames(result$proportions),
    stringsAsFactors = FALSE
  )

  if (!is.null(result$correlations)) {
    diag_df$correlation <- result$correlations
  }
  if (!is.null(result$rmse)) {
    diag_df$rmse <- result$rmse
  }
  if (!is.null(result$p_values)) {
    diag_df$p_value <- result$p_values
  }

  if (ncol(diag_df) > 1) {
    diag_file <- file.path(output_dir, paste0(prefix, "_diagnostics.csv"))
    write.csv(diag_df, diag_file, row.names = FALSE)
    files_created <- c(files_created, diag_file)
    if (verbose) cat("        Diagnostics:", diag_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Metadata
  # --------------------------------------------------------------------------
  meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))
  meta_lines <- c(
    "=== CIBERSORT Deconvolution Results ===",
    paste("Samples:", result$metadata$n_samples),
    paste("Cell types:", result$metadata$n_cell_types),
    paste("Permutations:", result$metadata$perm),
    paste("Quantile normalization:", result$metadata$QN),
    paste("Signature matrix:", result$metadata$sig_matrix),
    ""
  )

  if (!is.null(result$metadata$mixture_file)) {
    meta_lines <- c(meta_lines,
      paste("Mixture file:", result$metadata$mixture_file)
    )
  }

  # Quality summary
  if (!is.null(result$correlations)) {
    meta_lines <- c(meta_lines,
      "",
      "--- Quality Metrics ---",
      paste("Mean correlation:", round(mean(result$correlations, na.rm = TRUE), 4)),
      paste("Min correlation:", round(min(result$correlations, na.rm = TRUE), 4)),
      paste("Max correlation:", round(max(result$correlations, na.rm = TRUE), 4))
    )
  }
  if (!is.null(result$rmse)) {
    meta_lines <- c(meta_lines,
      paste("Mean RMSE:", round(mean(result$rmse, na.rm = TRUE), 4)),
      paste("Min RMSE:", round(min(result$rmse, na.rm = TRUE), 4)),
      paste("Max RMSE:", round(max(result$rmse, na.rm = TRUE), 4))
    )
  }

  # Top cell types by mean proportion
  mean_props <- sort(colMeans(result$proportions), decreasing = TRUE)
  top_cells <- head(mean_props[mean_props > 0.01], 10)
  if (length(top_cells) > 0) {
    meta_lines <- c(meta_lines,
      "",
      "--- Top Cell Types (mean proportion > 1%) ---"
    )
    for (i in seq_along(top_cells)) {
      meta_lines <- c(meta_lines,
        paste0("  ", names(top_cells)[i], ": ", round(top_cells[i] * 100, 1), "%")
      )
    }
  }

  if (!is.null(result$metadata$elapsed_seconds)) {
    meta_lines <- c(meta_lines,
      "",
      paste("Runtime:", round(result$metadata$elapsed_seconds, 1), "seconds")
    )
  }

  meta_lines <- c(meta_lines, "", paste("Exported:", Sys.time()))
  writeLines(meta_lines, meta_file)
  files_created <- c(files_created, meta_file)
  if (verbose) cat("[ELEUTHIA]   Metadata:", meta_file, "\n")

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  if (save_plots) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }

    if (verbose) cat("[ELEUTHIA] Generating plots...\n")

    # Stacked bar plot (ggplot)
    tryCatch({
      p <- AETHER_plot_cibersort_proportions(result, group_by = group_by)
      plot_files <- .save_ggplot(p, plot_dir, "cell_proportions", plot_format,
                                  width = max(8, result$metadata$n_samples * 0.4 + 3),
                                  height = 7)
      files_created <- c(files_created, plot_files)
      if (verbose) cat("[ELEUTHIA]   Cell proportions (stacked bar)\n")
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA]   Warning: Could not create proportions plot -", e$message, "\n")
    })

    # Heatmap (pheatmap)
    tryCatch({
      # Build annotation if group_by provided
      annot_col <- NULL
      if (!is.null(group_by)) {
        group_labels <- if (!is.null(names(group_by))) names(group_by) else group_by
        annot_col <- data.frame(
          condition = group_labels,
          row.names = rownames(result$proportions)
        )
      }

      .save_pheatmap(
        function() AETHER_plot_cibersort_heatmap(result, annotation_col = annot_col),
        plot_dir, "cell_heatmap", plot_format,
        width = max(8, result$metadata$n_samples * 0.3 + 4),
        height = max(6, result$metadata$n_cell_types * 0.25 + 3)
      )
      files_created <- c(files_created,
                          .plot_paths(plot_dir, "cell_heatmap", plot_format))
      if (verbose) cat("[ELEUTHIA]   Cell heatmap\n")
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA]   Warning: Could not create heatmap -", e$message, "\n")
    })

    # Boxplot by group (only if group_by provided, ggplot)
    if (!is.null(group_by)) {
      tryCatch({
        p <- AETHER_plot_cibersort_boxplot(result, group_by = group_by, top_n = top_n)
        n_cell_types <- if (!is.null(top_n)) min(top_n, result$metadata$n_cell_types) else result$metadata$n_cell_types
        plot_files <- .save_ggplot(p, plot_dir, "cell_boxplot", plot_format,
                                    width = 16,
                                    height = max(6, ceiling(n_cell_types / 4) * 3))
        files_created <- c(files_created, plot_files)
        if (verbose) cat("[ELEUTHIA]   Cell boxplot (by group)\n")
      }, error = function(e) {
        if (verbose) cat("[ELEUTHIA]   Warning: Could not create boxplot -", e$message, "\n")
      })
    }

    if (verbose) cat("[ELEUTHIA]   Plots saved to:", plot_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, "_result.rds"))
    saveRDS(result, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("[ELEUTHIA]   RDS:", rds_file, "\n")
  }

  if (verbose) {
    cat("[ELEUTHIA] Export Summary \n")
    cat("    Total files created:", length(files_created), "\n")
    cat("    Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}
