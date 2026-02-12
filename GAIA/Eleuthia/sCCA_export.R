#' Eleuthia - Save Results Functions
#'
#' @description Functions for saving analysis results, plots, and summaries
#' in an organized structure.


#' Save sCCA Analysis Results
#'
#' @description Comprehensive save function for sCCA cross-validation analysis.
#' Saves plots (PNG and RDS), CV results, best model, parameter grids, and
#' generates a plain text summary.
#'
#' @param output_dir Character string. Directory where all results will be saved.
#' @param coarse_cv_results List. Output from MINERVA_cv_scca for coarse grid.
#'   Must contain $cv_results and $best_model.
#' @param fine_cv_results List. Output from MINERVA_cv_scca for fine grid.
#'   Must contain $cv_results and $best_model.
#' @param plots List of ggplot objects with named elements:
#'   \itemize{
#'     \item cv_marginal_full: Marginal plot with full x-axis
#'     \item cv_marginal_zoom: Marginal plot zoomed to data (optional)
#'     \item cv_best_scores: Bar plot of best CV scores
#'     \item feature_weights: Multi-feature weight plot
#'   }
#' @param method Character string. Method name ("ConvCCA" or "RelPMDCCA").
#' @param X List of datasets used in analysis (for metadata).
#' @param pilot_results List. Output from MINERVA_pilot_compare_methods (optional).
#' @param prefix Character string. Prefix for saved files (default = "scca").
#' @param metadata List. Additional metadata to include in summary:
#'   \itemize{
#'     \item k: Number of CV folds
#'     \item lambda: Lambda parameter
#'     \item nIter: Number of iterations
#'     \item penalty: Penalty type
#'     \item metric: CV metric used
#'   }
#'
#' @return Invisibly returns a list with paths to all saved files.
#'
#' @details
#' Creates the following file structure in output_dir:
#' \preformatted{
#' output_dir/
#'   ├── plots/
#'   │   ├── {prefix}_cv_marginal_full.png
#'   │   ├── {prefix}_cv_marginal_zoom.png
#'   │   ├── {prefix}_cv_best_scores.png
#'   │   └── {prefix}_feature_weights.png
#'   ├── plot_objects/
#'   │   └── {prefix}_plots.rds (all ggplot objects)
#'   ├── code_elements/
#'   │   ├── {prefix}_coarse_cv.rds
#'   │   ├── {prefix}_fine_cv.rds
#'   │   ├── {prefix}_best_model.rds
#'   │   ├── {prefix}_coarse_grid.csv
#'   │   └── {prefix}_fine_grid.csv
#'   └── {prefix}_summary.txt
#' }
#'
#' @export
#'
#' @examples
#' # After running full CV workflow
#' saved_files <- ELEUTHIA_save_scca_results(
#'   output_dir = "results/analysis_2026_01_14",
#'   coarse_cv_results = cv_result_coarse,
#'   fine_cv_results = cv_result_fine,
#'   plots = list(
#'     cv_marginal_full = p_marginal_full,
#'     cv_marginal_zoom = p_marginal_zoom,
#'     cv_best_scores = p_best_scores,
#'     feature_weights = p_weights
#'   ),
#'   method = "ConvCCA",
#'   X = sim_data$datasets,
#'   pilot_results = pilot_result,
#'   metadata = list(k = 5, lambda = 10, nIter = 100, penalty = "LASSO")
#' )
#'
ELEUTHIA_save_scca_results <- function(output_dir,
                                        coarse_cv_results,
                                        fine_cv_results,
                                        plots,
                                        method,
                                        X,
                                        pilot_results = NULL,
                                        prefix = "scca",
                                        metadata = list()) {

  # Input validation
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create subdirectories
  plots_dir <- file.path(output_dir, "plots")
  plot_objects_dir <- file.path(output_dir, "plot_objects")
  code_elements_dir <- file.path(output_dir, "code_elements")

  dir.create(plots_dir, showWarnings = FALSE)
  dir.create(plot_objects_dir, showWarnings = FALSE)
  dir.create(code_elements_dir, showWarnings = FALSE)

  # Track all saved files
  saved_files <- list()

  cat("Saving sCCA analysis results...\n")

  # ============================================================================
  # 1. Save plots as PNG
  # ============================================================================
  cat("  - Saving plots (PNG)...\n")

  if (!is.null(plots$cv_marginal_full)) {
    path <- file.path(plots_dir, paste0(prefix, "_cv_marginal_full.png"))
    ggsave(path, plot = plots$cv_marginal_full, width = 12, height = 4, dpi = 300)
    saved_files$cv_marginal_full_png <- path
  }

  if (!is.null(plots$cv_marginal_zoom)) {
    path <- file.path(plots_dir, paste0(prefix, "_cv_marginal_zoom.png"))
    ggsave(path, plot = plots$cv_marginal_zoom, width = 12, height = 4, dpi = 300)
    saved_files$cv_marginal_zoom_png <- path
  }

  if (!is.null(plots$cv_best_scores)) {
    path <- file.path(plots_dir, paste0(prefix, "_cv_best_scores.png"))
    ggsave(path, plot = plots$cv_best_scores, width = 6, height = 4, dpi = 300)
    saved_files$cv_best_scores_png <- path
  }

  if (!is.null(plots$feature_weights)) {
    path <- file.path(plots_dir, paste0(prefix, "_feature_weights.png"))
    ggsave(path, plot = plots$feature_weights, width = 10, height = 8, dpi = 300)
    saved_files$feature_weights_png <- path
  }

  # ============================================================================
  # 2. Save plot objects as RDS
  # ============================================================================
  cat("  - Saving plot objects (RDS)...\n")
  path <- file.path(plot_objects_dir, paste0(prefix, "_plots.rds"))
  saveRDS(plots, path)
  saved_files$plot_objects <- path

  # ============================================================================
  # 3. Save CV results
  # ============================================================================
  cat("  - Saving CV results...\n")

  path <- file.path(code_elements_dir, paste0(prefix, "_coarse_cv.rds"))
  saveRDS(coarse_cv_results, path)
  saved_files$coarse_cv <- path

  path <- file.path(code_elements_dir, paste0(prefix, "_fine_cv.rds"))
  saveRDS(fine_cv_results, path)
  saved_files$fine_cv <- path

  # ============================================================================
  # 4. Save best model separately
  # ============================================================================
  cat("  - Saving best model...\n")
  path <- file.path(code_elements_dir, paste0(prefix, "_best_model.rds"))
  saveRDS(fine_cv_results$best_model, path)
  saved_files$best_model <- path

  # ============================================================================
  # 5. Save parameter grids as CSV
  # ============================================================================
  cat("  - Saving parameter grids...\n")

  path <- file.path(code_elements_dir, paste0(prefix, "_coarse_grid.csv"))
  write.csv(coarse_cv_results$cv_results, path, row.names = FALSE)
  saved_files$coarse_grid <- path

  path <- file.path(code_elements_dir, paste0(prefix, "_fine_grid.csv"))
  write.csv(fine_cv_results$cv_results, path, row.names = FALSE)
  saved_files$fine_grid <- path

  # ============================================================================
  # 6. Generate and save plain text summary
  # ============================================================================
  cat("  - Generating summary...\n")

  summary_txt <- ELEUTHIA_generate_scca_summary(
    coarse_cv_results = coarse_cv_results,
    fine_cv_results = fine_cv_results,
    pilot_results = pilot_results,
    method = method,
    X = X,
    metadata = metadata,
    saved_files = saved_files
  )

  path <- file.path(output_dir, paste0(prefix, "_summary.txt"))
  writeLines(summary_txt, path)
  saved_files$summary <- path

  cat(sprintf("\nAll results saved to: %s\n", output_dir))
  cat(sprintf("Summary available at: %s\n", path))

  invisible(saved_files)
}


#' Generate sCCA Analysis Summary (Plain Text)
#'
#' @description Internal function to generate a plain text summary
#' of sCCA cross-validation analysis.
#'
#' @param coarse_cv_results List from MINERVA_cv_scca
#' @param fine_cv_results List from MINERVA_cv_scca
#' @param pilot_results List from MINERVA_pilot_compare_methods (optional)
#' @param method Character string
#' @param X List of datasets
#' @param metadata List of additional metadata
#' @param saved_files List of saved file paths
#'
#' @return Character vector with plain text content (one element per line)
#'
#' @keywords internal
#'
ELEUTHIA_generate_scca_summary <- function(coarse_cv_results,
                                            fine_cv_results,
                                            pilot_results,
                                            method,
                                            X,
                                            metadata,
                                            saved_files) {

  # Initialize text content
  txt <- c()

  # Header
  txt <- c(txt, "================================================================================")
  txt <- c(txt, "sCCA ANALYSIS SUMMARY")
  txt <- c(txt, "================================================================================")
  txt <- c(txt, "")
  txt <- c(txt, paste0("Date:     ", format(Sys.time(), "%Y-%m-%d %H:%M")))
  txt <- c(txt, paste0("Method:   ", method))
  txt <- c(txt, paste0("Datasets: ", length(X)))
  txt <- c(txt, "")

  # Dataset dimensions
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "DATASET DIMENSIONS")
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "")
  for (i in seq_along(X)) {
    txt <- c(txt, sprintf("  Dataset %d: %d samples × %d features",
                        i, nrow(X[[i]]), ncol(X[[i]])))
  }
  txt <- c(txt, "")

  # Pilot results (if available)
  if (!is.null(pilot_results)) {
    txt <- c(txt, "--------------------------------------------------------------------------------")
    txt <- c(txt, "PILOT METHOD COMPARISON")
    txt <- c(txt, "--------------------------------------------------------------------------------")
    txt <- c(txt, "")
    txt <- c(txt, sprintf("%-15s %12s %10s %16s %s", "Method", "Correlation", "Sparsity", "Composite Score", "Winner"))
    txt <- c(txt, strrep("-", 80))

    for (i in 1:nrow(pilot_results$results)) {
      row <- pilot_results$results[i, ]
      winner <- if (row$method == pilot_results$best_method) "✓" else ""
      txt <- c(txt, sprintf("%-15s %12.4f %9.1f%% %16.4f %s",
                          row$method,
                          row$correlation,
                          row$sparsity_mean * 100,
                          row$composite_score,
                          winner))
    }
    txt <- c(txt, "")
  }

  # Coarse grid search
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "COARSE GRID SEARCH")
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "")

  # Extract unique tau values from coarse grid
  n_datasets <- length(X)
  tau_cols <- paste0("param", 1:n_datasets)
  coarse_grid <- coarse_cv_results$cv_results[, tau_cols, drop = FALSE]
  unique_taus <- lapply(coarse_grid, function(x) sort(unique(x)))

  txt <- c(txt, paste0("  Tau values:  [", paste(unique_taus[[1]], collapse = ", "), "]"))
  txt <- c(txt, paste0("  Grid size:   ", nrow(coarse_cv_results$cv_results), " combinations"))

  # Best coarse tau
  best_coarse_idx <- which.max(coarse_cv_results$cv_results$mean_score)
  best_coarse_tau <- coarse_cv_results$cv_results[best_coarse_idx, tau_cols, drop = FALSE]
  tau_str <- paste(sprintf("Dataset%d=%.2f", 1:n_datasets, as.numeric(best_coarse_tau)),
                   collapse = ", ")
  txt <- c(txt, paste0("  Best tau:    ", tau_str))
  txt <- c(txt, sprintf("  CV Score:    %.4f",
                      coarse_cv_results$cv_results$mean_score[best_coarse_idx]))
  txt <- c(txt, "")

  # Fine grid refinement
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "FINE GRID REFINEMENT")
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "")

  # Extract unique tau values from fine grid
  fine_grid <- fine_cv_results$cv_results[, tau_cols, drop = FALSE]
  unique_fine_taus <- lapply(fine_grid, function(x) sort(unique(x[!is.na(x)])))

  txt <- c(txt, paste0("  Fine tau values: [", paste(unique_fine_taus[[1]], collapse = ", "), "]"))

  # Count non-NA combinations
  n_non_na <- sum(!is.na(fine_cv_results$cv_results$mean_score))
  txt <- c(txt, sprintf("  Grid size:       %d total combinations (%d with valid CV scores)",
                      nrow(fine_cv_results$cv_results), n_non_na))

  # Best fine tau
  best_fine_idx <- which.max(fine_cv_results$cv_results$mean_score)
  best_fine_tau <- fine_cv_results$cv_results[best_fine_idx, tau_cols, drop = FALSE]
  tau_str <- paste(sprintf("Dataset%d=%.2f", 1:n_datasets, as.numeric(best_fine_tau)),
                   collapse = ", ")
  txt <- c(txt, paste0("  Best tau:        ", tau_str))

  best_fine_score <- fine_cv_results$cv_results$mean_score[best_fine_idx]
  best_coarse_score <- coarse_cv_results$cv_results$mean_score[best_coarse_idx]

  txt <- c(txt, sprintf("  CV Score:        %.4f", best_fine_score))

  improvement <- ((best_fine_score - best_coarse_score) / best_coarse_score) * 100
  txt <- c(txt, sprintf("  Improvement:     %+.1f%%", improvement))
  txt <- c(txt, "")

  # Final model sparsity
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "FINAL MODEL")
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "")
  txt <- c(txt, "Sparsity:")

  for (i in seq_along(fine_cv_results$best_model$W)) {
    w <- fine_cv_results$best_model$W[[i]]
    n_nonzero <- sum(abs(w) > 1e-6)
    n_total <- length(w)
    sparsity_pct <- (n_nonzero / n_total) * 100

    txt <- c(txt, sprintf("  Dataset %d: %d/%d features (%.1f%%)",
                        i, n_nonzero, n_total, sparsity_pct))
  }
  txt <- c(txt, "")

  # Analysis parameters
  if (length(metadata) > 0) {
    txt <- c(txt, "--------------------------------------------------------------------------------")
    txt <- c(txt, "ANALYSIS PARAMETERS")
    txt <- c(txt, "--------------------------------------------------------------------------------")
    txt <- c(txt, "")
    if (!is.null(metadata$k)) txt <- c(txt, paste0("  CV folds:    ", metadata$k))
    if (!is.null(metadata$lambda)) txt <- c(txt, paste0("  Lambda:      ", metadata$lambda))
    if (!is.null(metadata$nIter)) txt <- c(txt, paste0("  Iterations:  ", metadata$nIter))
    if (!is.null(metadata$penalty)) txt <- c(txt, paste0("  Penalty:     ", metadata$penalty))
    if (!is.null(metadata$metric)) txt <- c(txt, paste0("  Metric:      ", metadata$metric))
    txt <- c(txt, "")
  }

  # Files saved
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "FILES SAVED")
  txt <- c(txt, "--------------------------------------------------------------------------------")
  txt <- c(txt, "")
  txt <- c(txt, "Plots:")
  if (!is.null(saved_files$cv_marginal_full_png))
    txt <- c(txt, sprintf("  %s", basename(saved_files$cv_marginal_full_png)))
  if (!is.null(saved_files$cv_marginal_zoom_png))
    txt <- c(txt, sprintf("  %s", basename(saved_files$cv_marginal_zoom_png)))
  if (!is.null(saved_files$cv_best_scores_png))
    txt <- c(txt, sprintf("  %s", basename(saved_files$cv_best_scores_png)))
  if (!is.null(saved_files$feature_weights_png))
    txt <- c(txt, sprintf("  %s", basename(saved_files$feature_weights_png)))
  txt <- c(txt, "")

  txt <- c(txt, "Data Objects:")
  if (!is.null(saved_files$coarse_cv))
    txt <- c(txt, sprintf("  %s - Coarse CV results", basename(saved_files$coarse_cv)))
  if (!is.null(saved_files$fine_cv))
    txt <- c(txt, sprintf("  %s - Fine CV results", basename(saved_files$fine_cv)))
  if (!is.null(saved_files$best_model))
    txt <- c(txt, sprintf("  %s - Best final model", basename(saved_files$best_model)))
  if (!is.null(saved_files$plot_objects))
    txt <- c(txt, sprintf("  %s - Plot objects (ggplot)", basename(saved_files$plot_objects)))
  txt <- c(txt, "")

  txt <- c(txt, "Parameter Grids:")
  if (!is.null(saved_files$coarse_grid))
    txt <- c(txt, sprintf("  %s", basename(saved_files$coarse_grid)))
  if (!is.null(saved_files$fine_grid))
    txt <- c(txt, sprintf("  %s", basename(saved_files$fine_grid)))

  txt <- c(txt, "")
  txt <- c(txt, "================================================================================")

  return(txt)
}
