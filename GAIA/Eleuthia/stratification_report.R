###############################################################################
########### Patient Stratification Pipeline Report ###########
###############################################################################

#' Generate Stratification Pipeline Report
#'
#' @description Generate a comprehensive report from the patient stratification
#' pipeline, including clustering results, stability assessment, variable
#' importance, and optionally save tables and plots.
#'
#' @param cluster_result An artemis_cluster object from ARTEMIS_cluster_mixed().
#' @param stability_result Optional. An artemis_stability object from
#'   ARTEMIS_cluster_stability().
#' @param char_result Optional. An artemis_characterization object from
#'   ARTEMIS_characterize_clusters().
#' @param famd_result Optional. An artemis_famd object from ARTEMIS_famd().
#' @param comparison_result Optional. An artemis_variable_comparison object
#'   from ARTEMIS_compare_variable_importance().
#' @param qc_result Optional. A hades_qc object from HADES_qc_mixed_data().
#' @param output_dir Character. Directory to save output files. If NULL (default),
#'   only prints to console without saving files.
#' @param prefix Character. Prefix for output filenames. Default = "stratification".
#' @param save_tables Logical. Save key tables as CSV files. Default = TRUE
#'   (when output_dir is provided).
#' @param save_plots Logical. Save plots as PNG files. Default = TRUE
#'   (when output_dir is provided).
#' @param plot_width Numeric. Width of saved plots in inches. Default = 10.
#' @param plot_height Numeric. Height of saved plots in inches. Default = 8.
#' @param verbose Logical. Print report to console. Default = TRUE.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{summary}{Character vector of the text report}
#'   \item{tables}{List of data frames that were/would be exported}
#'   \item{files_saved}{Character vector of files saved (if output_dir provided)}
#' }
#'
#' @details
#' This function consolidates results from the stratification pipeline into
#' a single report. It can:
#' \itemize{
#'   \item Print a formatted summary to console
#'   \item Save the summary as a text file
#'   \item Export key tables (cluster assignments, variable importance, etc.) as CSV
#'   \item Generate and save standard plots (silhouette, FAMD, radar)
#' }
#'
#' At minimum, cluster_result is required. Other results are optional and
#' will be included if provided.
#'
#' @export
#'
#' @examples
#' # Generate console report only
#' ELEUTHIA_stratification_report(cluster_result)
#'
#' # Full report with file export
#' ELEUTHIA_stratification_report(
#'   cluster_result = clust,
#'   stability_result = stab,
#'   char_result = char,
#'   famd_result = famd,
#'   output_dir = "results/stratification"
#' )
#'
ELEUTHIA_stratification_report <- function(cluster_result,
                                            stability_result = NULL,
                                            char_result = NULL,
                                            famd_result = NULL,
                                            comparison_result = NULL,
                                            qc_result = NULL,
                                            output_dir = NULL,
                                            prefix = "stratification",
                                            save_tables = TRUE,
                                            save_plots = TRUE,
                                            plot_width = 10,
                                            plot_height = 8,
                                            verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(cluster_result, "artemis_cluster")) {
    stop("cluster_result must be an artemis_cluster object from ARTEMIS_cluster_mixed()")
  }

  # Create output directory if specified
  files_saved <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  # Initialize report
  report <- character(0)
  tables <- list()

  # Helper to add lines to report
  add_line <- function(...) {
    report <<- c(report, paste0(...))
  }

  add_separator <- function(char = "=", width = 70) {
    add_line(paste(rep(char, width), collapse = ""))
  }

  # ---------------------------------------------------------------------------
  # Header
  # ---------------------------------------------------------------------------
  add_separator()
  add_line("PATIENT STRATIFICATION PIPELINE REPORT")
  add_separator()
  add_line("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add_line("")

  # ---------------------------------------------------------------------------
  # Section 1: Data Overview
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("1. DATA OVERVIEW")
  add_separator("-", 70)

  data_used <- cluster_result$data_used
  n_samples <- nrow(data_used)
  n_vars <- ncol(data_used)
  n_quanti <- sum(sapply(data_used, is.numeric))
  n_quali <- n_vars - n_quanti

  add_line("Samples: ", n_samples)
  add_line("Variables: ", n_vars, " (", n_quanti, " quantitative, ", n_quali, " qualitative)")
  add_line("")

  # QC summary if provided
  if (!is.null(qc_result) && inherits(qc_result, "hades_qc")) {
    add_line("QC Summary:")
    add_line("  Variables flagged (any issue): ", qc_result$summary$variables_flagged_any)
    add_line("  Samples flagged (high NA): ", qc_result$summary$samples_flagged_na)
    if (!is.null(qc_result$summary$n_high_correlation_pairs)) {
      add_line("  High correlation pairs: ", qc_result$summary$n_high_correlation_pairs)
    }
    add_line("")
  }

  # ---------------------------------------------------------------------------
  # Section 2: Clustering Results
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("2. CLUSTERING RESULTS")
  add_separator("-", 70)

  add_line("Method: ", cluster_result$method)
  add_line("Number of clusters (k): ", cluster_result$k)
  add_line("Average silhouette width: ", round(cluster_result$silhouette_avg, 3))
  add_line("")

  # Cluster sizes
  cluster_sizes <- table(cluster_result$clusters)
  add_line("Cluster sizes:")
  for (cl in names(cluster_sizes)) {
    pct <- round(cluster_sizes[cl] / n_samples * 100, 1)
    add_line("  Cluster ", cl, ": ", cluster_sizes[cl], " samples (", pct, "%)")
  }
  add_line("")

  # Silhouette interpretation
  sil_avg <- cluster_result$silhouette_avg
  if (sil_avg > 0.7) {
    sil_interp <- "Strong structure"
  } else if (sil_avg > 0.5) {
    sil_interp <- "Reasonable structure"
  } else if (sil_avg > 0.25) {
    sil_interp <- "Weak structure"
  } else {
    sil_interp <- "No substantial structure"
  }
  add_line("Silhouette interpretation: ", sil_interp)
  add_line("")

  # K selection info if auto
  if (!is.null(cluster_result$k_selection)) {
    add_line("K selection (silhouette scores):")
    for (i in seq_len(nrow(cluster_result$k_selection))) {
      row <- cluster_result$k_selection[i, ]
      marker <- if (row$k == cluster_result$k) " <- selected" else ""
      add_line("  k=", row$k, ": ", round(row$silhouette_avg, 3), marker)
    }
    add_line("")
  }

  # Create cluster assignments table
  cluster_table <- data.frame(
    sample = if (!is.null(rownames(data_used))) rownames(data_used) else paste0("sample_", seq_len(n_samples)),
    cluster = cluster_result$clusters,
    stringsAsFactors = FALSE
  )

  # Add silhouette width per sample
  if (!is.null(cluster_result$silhouette)) {
    cluster_table$silhouette_width <- cluster_result$silhouette[, "sil_width"]
  }

  tables$cluster_assignments <- cluster_table

  # ---------------------------------------------------------------------------
  # Section 3: Stability Assessment
  # ---------------------------------------------------------------------------
  if (!is.null(stability_result) && inherits(stability_result, "artemis_stability")) {
    add_separator("-", 70)
    add_line("3. CLUSTER STABILITY")
    add_separator("-", 70)

    add_line("Bootstrap iterations: ", stability_result$n_boot)
    add_line("Mean Jaccard similarity: ", round(stability_result$jaccard_mean, 3))
    add_line("Interpretation: ", stability_result$interpretation)
    add_line("")

    add_line("Per-cluster stability:")
    for (cl in names(stability_result$jaccard_per_cluster)) {
      add_line("  Cluster ", cl, ": ", round(stability_result$jaccard_per_cluster[cl], 3))
    }
    add_line("")
  }

  # ---------------------------------------------------------------------------
  # Section 4: Variable Importance
  # ---------------------------------------------------------------------------
  if (!is.null(char_result) && inherits(char_result, "artemis_characterization")) {
    add_separator("-", 70)
    add_line("4. VARIABLE IMPORTANCE (Cluster Discrimination)")
    add_separator("-", 70)

    results <- char_result$results

    n_large <- sum(results$effect_interpretation == "Large", na.rm = TRUE)
    n_medium <- sum(results$effect_interpretation == "Medium", na.rm = TRUE)
    n_small <- sum(results$effect_interpretation == "Small", na.rm = TRUE)

    add_line("Variables with large effect: ", n_large)
    add_line("Variables with medium effect: ", n_medium)
    add_line("Variables with small effect: ", n_small)
    add_line("")

    add_line("Top 15 discriminating variables:")
    add_line(sprintf("  %-3s %-25s %-12s %8s %12s",
                     "#", "Variable", "Type", "Effect", "Interpretation"))
    add_line(paste(rep("-", 65), collapse = ""))

    n_show <- min(15, nrow(results))
    for (i in seq_len(n_show)) {
      row <- results[i, ]
      add_line(sprintf("  %-3d %-25s %-12s %8.3f %12s",
                       i,
                       substr(row$variable, 1, 25),
                       row$type,
                       row$effect_size,
                       row$effect_interpretation))
    }
    add_line("")

    tables$variable_importance <- results
  }

  # ---------------------------------------------------------------------------
  # Section 5: FAMD Results
  # ---------------------------------------------------------------------------
  if (!is.null(famd_result) && inherits(famd_result, "artemis_famd")) {
    add_separator("-", 70)
    add_line("5. FAMD DIMENSIONALITY REDUCTION")
    add_separator("-", 70)

    eigenvalues <- famd_result$eigenvalues

    add_line("Variance explained by top dimensions:")
    n_show <- min(5, nrow(eigenvalues))
    for (i in seq_len(n_show)) {
      add_line(sprintf("  Dim %d: %5.1f%% (cumulative: %5.1f%%)",
                       i,
                       eigenvalues$variance_percent[i],
                       eigenvalues$cumulative_percent[i]))
    }
    add_line("")

    # Dimensions needed for 80% variance
    dims_80 <- which(eigenvalues$cumulative_percent >= 80)[1]
    if (!is.na(dims_80)) {
      add_line("Dimensions for 80% variance: ", dims_80)
    }
    add_line("")

    tables$famd_eigenvalues <- eigenvalues
  }

  # ---------------------------------------------------------------------------
  # Section 6: Variable Comparison
  # ---------------------------------------------------------------------------
  if (!is.null(comparison_result) && inherits(comparison_result, "artemis_variable_comparison")) {
    add_separator("-", 70)
    add_line("6. VARIABLE REFINEMENT ANALYSIS")
    add_separator("-", 70)

    add_line("Low threshold: ", comparison_result$summary$low_threshold * 100, " percentile")
    add_line("Variables low on cluster effect only: ",
             sum(comparison_result$comparison$low_cluster & !comparison_result$comparison$low_famd))
    add_line("Variables low on FAMD contrib only: ",
             sum(comparison_result$comparison$low_famd & !comparison_result$comparison$low_cluster))
    add_line("Variables low on BOTH (removal candidates): ",
             comparison_result$summary$n_low_both)
    add_line("")

    if (length(comparison_result$low_signal) > 0) {
      add_line("Candidates for removal:")
      for (v in comparison_result$low_signal) {
        add_line("  - ", v)
      }
      add_line("")
    }

    tables$variable_comparison <- comparison_result$comparison
  }

  # ---------------------------------------------------------------------------
  # Section 7: Recommendations
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("7. RECOMMENDATIONS")
  add_separator("-", 70)

  recommendations <- character(0)

  # Based on silhouette
  if (sil_avg < 0.25) {
    recommendations <- c(recommendations,
                         "- Cluster structure is weak. Consider different k or variable selection.")
  } else if (sil_avg < 0.5) {
    recommendations <- c(recommendations,
                         "- Cluster structure is moderate. Variable refinement may improve results.")
  }

  # Based on stability
  if (!is.null(stability_result)) {
    if (stability_result$jaccard_mean < 0.6) {
      recommendations <- c(recommendations,
                           "- Clusters are unstable. Consider removing noisy variables.")
    }
  }

  # Based on variable comparison
  if (!is.null(comparison_result) && comparison_result$summary$n_low_both > 0) {
    recommendations <- c(recommendations,
                         paste0("- ", comparison_result$summary$n_low_both,
                                " low-signal variables identified for potential removal."))
  }

  if (length(recommendations) == 0) {
    add_line("Clustering appears robust. Proceed with interpretation.")
  } else {
    for (rec in recommendations) {
      add_line(rec)
    }
  }
  add_line("")

  add_separator()
  add_line("END OF REPORT")
  add_separator()

  # ---------------------------------------------------------------------------
  # Print to console
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat(paste(report, collapse = "\n"), "\n")
  }

  # ---------------------------------------------------------------------------
  # Save files
  # ---------------------------------------------------------------------------
  if (!is.null(output_dir)) {

    # Save text report
    report_file <- file.path(output_dir, paste0(prefix, "_report.txt"))
    writeLines(report, report_file)
    files_saved <- c(files_saved, report_file)

    # Save tables
    if (save_tables) {
      for (table_name in names(tables)) {
        table_file <- file.path(output_dir, paste0(prefix, "_", table_name, ".csv"))
        write.csv(tables[[table_name]], table_file, row.names = FALSE)
        files_saved <- c(files_saved, table_file)
      }
    }

    # Save plots
    if (save_plots) {
      # Silhouette plot
      tryCatch({
        sil_plot <- AETHER_plot_silhouette(cluster_result)
        sil_file <- file.path(output_dir, paste0(prefix, "_silhouette.png"))
        ggplot2::ggsave(sil_file, sil_plot, width = plot_width, height = plot_height)
        files_saved <- c(files_saved, sil_file)
      }, error = function(e) {
        warning("Could not save silhouette plot: ", e$message)
      })

      # K selection and Sankey plots (only if k="auto" was used)
      if (!is.null(cluster_result$k_evaluation)) {
        # K selection plot
        tryCatch({
          k_sel_plot <- AETHER_plot_k_selection(cluster_result$k_evaluation)
          k_sel_file <- file.path(output_dir, paste0(prefix, "_k_selection.png"))
          ggplot2::ggsave(k_sel_file, k_sel_plot, width = plot_width, height = plot_height)
          files_saved <- c(files_saved, k_sel_file)
        }, error = function(e) {
          warning("Could not save k selection plot: ", e$message)
        })

        # Cluster Sankey diagram (saved as HTML since it's interactive)
        tryCatch({
          sankey_file <- file.path(output_dir, paste0(prefix, "_cluster_sankey.html"))
          AETHER_plot_cluster_sankey(cluster_result$k_evaluation, save_html = sankey_file)
          files_saved <- c(files_saved, sankey_file)
        }, error = function(e) {
          warning("Could not save cluster Sankey plot: ", e$message)
        })
      }

      # FAMD plots
      if (!is.null(famd_result)) {
        tryCatch({
          # Individuals plot (colored by cluster)
          cluster_vec <- setNames(as.character(cluster_result$clusters),
                                   rownames(cluster_result$data_used))
          ind_plot <- AETHER_plot_famd_individuals(famd_result, color_by = cluster_vec)
          ind_file <- file.path(output_dir, paste0(prefix, "_famd_individuals.png"))
          ggplot2::ggsave(ind_file, ind_plot, width = plot_width, height = plot_height)
          files_saved <- c(files_saved, ind_file)

          # Scree plot
          scree_plot <- AETHER_plot_famd_scree(famd_result)
          scree_file <- file.path(output_dir, paste0(prefix, "_famd_scree.png"))
          ggplot2::ggsave(scree_file, scree_plot, width = plot_width, height = plot_height)
          files_saved <- c(files_saved, scree_file)

          # Variables plot
          var_plot <- AETHER_plot_famd_variables(famd_result, aggregate_quali = TRUE)
          var_file <- file.path(output_dir, paste0(prefix, "_famd_variables.png"))
          ggplot2::ggsave(var_file, var_plot, width = plot_width, height = plot_height)
          files_saved <- c(files_saved, var_file)
        }, error = function(e) {
          warning("Could not save FAMD plots: ", e$message)
        })
      }

      # Radar plot
      if (!is.null(char_result)) {
        tryCatch({
          radar_plot <- AETHER_plot_cluster_radar(char_result,
                                                   data = cluster_result$data_used,
                                                   clusters = cluster_result$clusters)
          radar_file <- file.path(output_dir, paste0(prefix, "_cluster_radar.png"))
          ggplot2::ggsave(radar_file, radar_plot, width = plot_width, height = plot_height, bg='white')
          files_saved <- c(files_saved, radar_file)
        }, error = function(e) {
          warning("Could not save radar plot: ", e$message)
        })
      }
    }

    if (verbose) {
      cat("\nFiles saved to:", output_dir, "\n")
      for (f in files_saved) {
        cat("  -", basename(f), "\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  invisible(list(
    summary = report,
    tables = tables,
    files_saved = files_saved
  ))
}


#' Quick Summary of Stratification Results
#'
#' @description Print a quick one-page summary of stratification results
#' without file export. Useful for interactive exploration.
#'
#' @param cluster_result An artemis_cluster object.
#' @param stability_result Optional. An artemis_stability object.
#' @param char_result Optional. An artemis_characterization object.
#'
#' @return Invisibly returns NULL.
#'
#' @export
ELEUTHIA_quick_summary <- function(cluster_result,
                                    stability_result = NULL,
                                    char_result = NULL) {

  cat("===== STRATIFICATION QUICK SUMMARY =====\n\n")

  # Clustering
  cat("CLUSTERING:\n")
  cat("  k =", cluster_result$k, "| Method:", cluster_result$method, "\n")
  cat("  Silhouette:", round(cluster_result$silhouette_avg, 3))

  sil <- cluster_result$silhouette_avg
  if (sil > 0.7) cat(" (Strong)\n")
  else if (sil > 0.5) cat(" (Reasonable)\n")
  else if (sil > 0.25) cat(" (Weak)\n")
  else cat(" (Poor)\n")

  cat("  Cluster sizes:", paste(table(cluster_result$clusters), collapse = " / "), "\n\n")

  # Stability
  if (!is.null(stability_result)) {
    cat("STABILITY:\n")
    cat("  Jaccard:", round(stability_result$jaccard_mean, 3),
        "-", stability_result$interpretation, "\n\n")
  }

  # Top variables
  if (!is.null(char_result)) {
    cat("TOP 5 DISCRIMINATING VARIABLES:\n")
    top5 <- head(char_result$results, 5)
    for (i in seq_len(nrow(top5))) {
      cat(sprintf("  %d. %s (effect=%.3f)\n",
                  i, top5$variable[i], top5$effect_size[i]))
    }
    cat("\n")
  }

  cat("=========================================\n")

  invisible(NULL)
}
