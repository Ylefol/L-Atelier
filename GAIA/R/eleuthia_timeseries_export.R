# GAIA/Eleuthia/timeseries_export.R
# Export time series analysis and PART clustering results
#
# Saves PART, time series DEA, gene selection, and enrichment results to
# CSV/RDS files and optionally generates plots.


#' Export Time Series Analysis Results
#'
#' Comprehensive export function for the time series analysis pipeline. Accepts
#' PART clustering results, time series DEA results, gene selection, and
#' enrichment results. Follows the same conventions as
#' \code{ELEUTHIA_export_wgcna_results()} and
#' \code{ELEUTHIA_export_activity_results()}.
#'
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part()}.
#'   NULL to skip PART export. Default: NULL.
#' @param ts_de One or a list of \code{artemis_ts_de} objects from
#'   \code{ARTEMIS_timeseries_conditional()} or
#'   \code{ARTEMIS_timeseries_temporal()}. NULL to skip. Default: NULL.
#' @param gene_selection Character vector of selected gene IDs (e.g., from
#'   \code{ARTEMIS_select_de_genes()}). NULL to skip. Default: NULL.
#' @param enrichment Enrichment result from \code{APOLLO_enrich_gost()}.
#'   NULL to skip. Delegates to \code{ELEUTHIA_export_enrichment()}.
#'   Default: NULL.
#' @param sample_info Data.frame with sample metadata (for trajectory plots).
#'   Must have rownames matching sample columns. Default: NULL.
#' @param time_col Character. Timepoint column in sample_info. Default: "timepoint".
#' @param group_col Character. Group column in sample_info. Default: "group".
#' @param save_plots Logical. Generate and save plots. Default: TRUE.
#' @param plot_format Character. "png", "pdf", or "both". Default: "png".
#' @param prefix Character. Prefix for output filenames. Default: "timeseries".
#' @param save_rds Logical. Save full R objects as RDS. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Depending on which arguments are provided, exports:
#'
#' \strong{PART results} (when \code{part_result} provided):
#' \itemize{
#'   \item Cluster assignments CSV (gene → cluster)
#'   \item Cluster map CSV (gene, cluster, color)
#'   \item Cluster summary CSV (cluster, size, color)
#'   \item Parameters text file
#'   \item Plots: heatmap, cluster trajectories (if sample_info provided),
#'     cluster means (if sample_info provided)
#' }
#'
#' \strong{Time series DEA} (when \code{ts_de} provided):
#' \itemize{
#'   \item Per-comparison full result tables in subdirectory
#'   \item Summary CSV (timepoints, n_up, n_down)
#'   \item Plots: DEA summary bar chart
#' }
#'
#' \strong{Gene selection} (when \code{gene_selection} provided):
#' \itemize{
#'   \item Selected genes CSV
#'   \item Gene summary CSV (per-comparison breakdown, if available)
#' }
#'
#' \strong{Enrichment} (when \code{enrichment} provided):
#' \itemize{
#'   \item Delegates to \code{ELEUTHIA_export_enrichment()} in subdirectory
#' }
#'
#' @examples
#' \dontrun{
#' ELEUTHIA_export_timeseries_results(
#'   output_dir = "results/timeseries/",
#'   part_result = part_res,
#'   ts_de = list(cond_de, temp_de),
#'   gene_selection = selected_genes,
#'   sample_info = sample_info
#' )
#'
#' }
#' @export
ELEUTHIA_export_timeseries_results <- function(output_dir,
                                                part_result = NULL,
                                                ts_de = NULL,
                                                gene_selection = NULL,
                                                enrichment = NULL,
                                                sample_info = NULL,
                                                time_col = "timepoint",
                                                group_col = "group",
                                                save_plots = TRUE,
                                                plot_format = "png",
                                                prefix = "timeseries",
                                                save_rds = TRUE,
                                                verbose = TRUE) {

  # Validate: at least one result object
  if (is.null(part_result) && is.null(ts_de) &&
      is.null(gene_selection) && is.null(enrichment)) {
    stop("At least one result object must be provided ",
         "(part_result, ts_de, gene_selection, or enrichment)")
  }

  plot_format <- match.arg(plot_format, c("png", "pdf", "both"))

  if (verbose) cat("[ELEUTHIA] Exporting Time Series Results \n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("[ELEUTHIA] Created directory:", output_dir, "\n")
  }

  files_created <- character(0)

  # ==========================================================================
  # PART RESULTS
  # ==========================================================================
  if (!is.null(part_result)) {
    if (!inherits(part_result, "artemis_part")) {
      stop("'part_result' must be an artemis_part object")
    }

    if (verbose) {
      cat("[ELEUTHIA] --- PART Clustering ---\n")
      cat("[ELEUTHIA] Clusters:", part_result$n_clusters, "| Features:", nrow(part_result$data),
          "| Outliers:", part_result$n_outliers, "\n")
    }

    # Cluster assignments (gene → cluster)
    assign_file <- file.path(output_dir, paste0(prefix, "_cluster_assignments.csv"))
    assign_df <- data.frame(
      gene = names(part_result$clusters),
      cluster = unname(part_result$clusters),
      stringsAsFactors = FALSE
    )
    write.csv(assign_df, assign_file, row.names = FALSE)
    files_created <- c(files_created, assign_file)
    if (verbose) cat("[ELEUTHIA]   Cluster assignments:", assign_file, "\n")

    # Cluster map (gene, cluster, color — ordered by cluster)
    map_file <- file.path(output_dir, paste0(prefix, "_cluster_map.csv"))
    write.csv(part_result$cluster_map, map_file, row.names = FALSE)
    files_created <- c(files_created, map_file)
    if (verbose) cat("[ELEUTHIA]   Cluster map:", map_file, "\n")

    # Cluster summary
    summary_file <- file.path(output_dir, paste0(prefix, "_cluster_summary.csv"))
    cl_names <- names(part_result$cluster_sizes)
    summary_df <- data.frame(
      cluster = cl_names,
      size = unname(part_result$cluster_sizes),
      color = unname(part_result$cluster_colors[cl_names]),
      stringsAsFactors = FALSE
    )
    write.csv(summary_df, summary_file, row.names = FALSE)
    files_created <- c(files_created, summary_file)
    if (verbose) cat("[ELEUTHIA]   Cluster summary:", summary_file, "\n")

    # Parameters
    param_file <- file.path(output_dir, paste0(prefix, "_part_parameters.txt"))
    p <- part_result$parameters
    param_lines <- c(
      "=== PART Clustering Parameters ===",
      paste("Clusters:", part_result$n_clusters),
      paste("Outliers:", part_result$n_outliers),
      paste("Total features:", nrow(part_result$data)),
      paste("Samples:", ncol(part_result$data)),
      paste("Computation time:", round(part_result$computation_time, 1), "seconds"),
      "",
      "--- Parameters ---",
      paste("q:", p$q),
      paste("min_size:", p$min_size),
      paste("B:", p$B),
      paste("Kmax:", p$Kmax),
      paste("dist_method:", p$dist_method),
      paste("linkage:", p$linkage),
      paste("scale:", p$scale),
      paste("seed:", ifelse(is.null(p$seed), "NULL", p$seed)),
      paste("threshold:", round(p$threshold, 6)),
      "",
      paste("Exported:", Sys.time())
    )
    writeLines(param_lines, param_file)
    files_created <- c(files_created, param_file)
    if (verbose) cat("[ELEUTHIA]   Parameters:", param_file, "\n")

    # PART plots
    if (save_plots) {
      plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
      if (!dir.exists(plot_dir)) dir.create(plot_dir)

      if (verbose) cat("[ELEUTHIA]   Generating plots...\n")

      # Heatmap (ComplexHeatmap — uses save_path directly)
      if (!is.null(sample_info)) {
        tryCatch({
          n_genes <- nrow(part_result$data)
          n_samples <- ncol(part_result$data)
          hm_width <- max(8, n_samples * 0.4 + 4)
          hm_height <- max(6, min(n_genes * 0.02 + 3, 20))

          if (plot_format %in% c("png", "both")) {
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              group_col = group_col, time_col = time_col,
              save_path = file.path(plot_dir, "part_heatmap.png"),
              width = hm_width, height = hm_height
            )
          }
          if (plot_format %in% c("pdf", "both")) {
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              group_col = group_col, time_col = time_col,
              save_path = file.path(plot_dir, "part_heatmap.pdf"),
              width = hm_width, height = hm_height
            )
          }
          files_created <- c(files_created,
                              .ts_plot_paths(plot_dir, "part_heatmap", plot_format))
          if (verbose) cat("[ELEUTHIA]     Heatmap\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]     Warning: Could not create heatmap -", e$message, "\n")
        })
      } else {
        if (verbose) cat("[ELEUTHIA]     Skipping heatmap (requires sample_info)\n")
      }

      # Cluster trajectories (ggplot, requires sample_info)
      if (!is.null(sample_info)) {
        tryCatch({
          p <- AETHER_plot_cluster_trajectories(
            part_result, sample_info,
            time_col = time_col, group_col = group_col
          )
          .ts_save_ggplot(p, plot_dir, "cluster_trajectories", plot_format,
                           width = max(10, part_result$n_clusters * 3),
                           height = max(6, ceiling(part_result$n_clusters / 3) * 3 + 2))
          files_created <- c(files_created,
                              .ts_plot_paths(plot_dir, "cluster_trajectories", plot_format))
          if (verbose) cat("[ELEUTHIA]     Cluster trajectories\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]     Warning: Could not create trajectory plot -", e$message, "\n")
        })

        # Cluster means overlay
        tryCatch({
          p <- AETHER_plot_cluster_means(
            part_result, sample_info,
            time_col = time_col, group_col = group_col
          )
          .ts_save_ggplot(p, plot_dir, "cluster_means", plot_format,
                           width = 10, height = 7)
          files_created <- c(files_created,
                              .ts_plot_paths(plot_dir, "cluster_means", plot_format))
          if (verbose) cat("[ELEUTHIA]     Cluster means overlay\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]     Warning: Could not create cluster means plot -", e$message, "\n")
        })
      }

      if (verbose) cat("[ELEUTHIA]   Plots saved to:", plot_dir, "/\n")
    }

    # RDS
    if (save_rds) {
      rds_file <- file.path(output_dir, paste0(prefix, "_part_result.rds"))
      saveRDS(part_result, rds_file)
      files_created <- c(files_created, rds_file)
      if (verbose) cat("[ELEUTHIA]   RDS:", rds_file, "\n")
    }
  }

  # ==========================================================================
  # TIME SERIES DEA
  # ==========================================================================
  if (!is.null(ts_de)) {
    # Normalize to list
    if (inherits(ts_de, "artemis_ts_de")) {
      ts_de_list <- list(ts_de)
    } else if (is.list(ts_de)) {
      ts_de_list <- ts_de
    } else {
      stop("'ts_de' must be an artemis_ts_de object or a list of them")
    }

    if (verbose) cat("[ELEUTHIA] Time Series DEA \n")

    de_dir <- file.path(output_dir, paste0(prefix, "_dea"))
    if (!dir.exists(de_dir)) dir.create(de_dir)

    for (i in seq_along(ts_de_list)) {
      ts_obj <- ts_de_list[[i]]
      if (!inherits(ts_obj, "artemis_ts_de")) {
        warning("Skipping non-artemis_ts_de element at position ", i)
        next
      }

      type_label <- ts_obj$type
      if (verbose) cat("[ELEUTHIA]   Type:", type_label, "| Comparisons:", length(ts_obj$results), "\n")

      # Summary CSV
      sum_file <- file.path(de_dir, paste0(type_label, "_summary.csv"))
      write.csv(ts_obj$summary, sum_file, row.names = FALSE)
      files_created <- c(files_created, sum_file)
      if (verbose) cat("[ELEUTHIA]     Summary:", sum_file, "\n")

      # Per-comparison results
      comp_dir <- file.path(de_dir, paste0(type_label, "_results"))
      if (!dir.exists(comp_dir)) dir.create(comp_dir)

      for (exp_name in names(ts_obj$results)) {
        de_res <- ts_obj$results[[exp_name]]
        if (!is.null(de_res$results)) {
          res_file <- file.path(comp_dir, paste0(exp_name, ".csv"))
          write.csv(de_res$results, res_file, row.names = FALSE)
          files_created <- c(files_created, res_file)
        }
      }
      if (verbose) cat("[ELEUTHIA]     Per-comparison results:", comp_dir, "/\n")

      # DEA summary plot
      if (save_plots) {
        plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
        if (!dir.exists(plot_dir)) dir.create(plot_dir)

        tryCatch({
          p <- AETHER_plot_timeseries_summary(ts_obj)
          .ts_save_ggplot(p, plot_dir,
                           paste0("dea_summary_", type_label), plot_format,
                           width = max(7, nrow(ts_obj$summary) * 1.2 + 2),
                           height = 6)
          files_created <- c(files_created,
                              .ts_plot_paths(plot_dir, paste0("dea_summary_", type_label),
                                              plot_format))
          if (verbose) cat("[ELEUTHIA]     DEA summary plot\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]     Warning: Could not create DEA summary plot -", e$message, "\n")
        })
      }

      # RDS
      if (save_rds) {
        rds_file <- file.path(de_dir, paste0(type_label, "_result.rds"))
        saveRDS(ts_obj, rds_file)
        files_created <- c(files_created, rds_file)
        if (verbose) cat("[ELEUTHIA]     RDS:", rds_file, "\n")
      }
    }
  }

  # ==========================================================================
  # GENE SELECTION
  # ==========================================================================
  if (!is.null(gene_selection)) {
    if (verbose) cat("[ELEUTHIA] Gene Selection \n")

    # Selected genes
    genes_file <- file.path(output_dir, paste0(prefix, "_selected_genes.csv"))
    genes_df <- data.frame(gene = gene_selection, stringsAsFactors = FALSE)
    write.csv(genes_df, genes_file, row.names = FALSE)
    files_created <- c(files_created, genes_file)
    if (verbose) cat("[ELEUTHIA]   Selected genes:", length(gene_selection), "->", genes_file, "\n")

    # Gene summary (if available as attribute)
    gene_summary <- attr(gene_selection, "gene_summary")
    if (!is.null(gene_summary)) {
      gsummary_file <- file.path(output_dir, paste0(prefix, "_gene_summary.csv"))
      write.csv(gene_summary, gsummary_file, row.names = FALSE)
      files_created <- c(files_created, gsummary_file)
      if (verbose) cat("[ELEUTHIA]   Gene summary:", gsummary_file, "\n")
    }

    # Selection parameters (if available as attribute)
    sel_params <- attr(gene_selection, "parameters")
    if (!is.null(sel_params)) {
      sel_file <- file.path(output_dir, paste0(prefix, "_selection_parameters.txt"))
      param_lines <- c(
        "=== Gene Selection Parameters ===",
        paste("Total selected:", length(gene_selection)),
        paste("log2FC threshold:", sel_params$l2fc_thresh),
        paste("P-value column:", sel_params$p_col),
        paste("P-value threshold:", sel_params$p_thresh),
        paste("Combination method:", ifelse(sel_params$union, "union", "intersection")),
        "",
        paste("Exported:", Sys.time())
      )
      writeLines(param_lines, sel_file)
      files_created <- c(files_created, sel_file)
      if (verbose) cat("[ELEUTHIA]   Selection parameters:", sel_file, "\n")
    }
  }

  # ==========================================================================
  # ENRICHMENT
  # ==========================================================================
  if (!is.null(enrichment)) {
    if (verbose) cat("[ELEUTHIA] Enrichment \n")

    enrich_dir <- file.path(output_dir, paste0(prefix, "_enrichment"))
    tryCatch({
      enrich_files <- ELEUTHIA_export_enrichment(
        enrichment,
        output_dir = enrich_dir,
        prefix = prefix,
        save_plots = save_plots,
        plot_format = plot_format,
        save_rds = save_rds,
        verbose = verbose
      )
      files_created <- c(files_created, enrich_files)
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA]   Warning: Could not export enrichment -", e$message, "\n")
    })
  }

  # ==========================================================================
  # SUMMARY
  # ==========================================================================
  if (verbose) {
    cat("[ELEUTHIA] Export Summary \n")
    cat("    Total files created:", length(files_created), "\n")
    cat("    Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}


# ==============================================================================
# INTERNAL: Plot save helpers (prefixed to avoid collision with activity_export)
# ==============================================================================

#' Save a ggplot to file(s)
#' @noRd
.ts_save_ggplot <- function(p, plot_dir, name, plot_format,
                              width = 10, height = 8) {
  if (plot_format %in% c("png", "both")) {
    png_file <- file.path(plot_dir, paste0(name, ".png"))
    ggplot2::ggsave(png_file, p, width = width, height = height,
                     dpi = 150, bg = "white")
  }
  if (plot_format %in% c("pdf", "both")) {
    pdf_file <- file.path(plot_dir, paste0(name, ".pdf"))
    ggplot2::ggsave(pdf_file, p, width = width, height = height)
  }
}


#' Save a pheatmap to file(s) via device wrapper
#' @noRd
.ts_save_pheatmap <- function(plot_fn, plot_dir, name, plot_format,
                                width = 10, height = 8) {
  if (plot_format %in% c("png", "both")) {
    png_file <- file.path(plot_dir, paste0(name, ".png"))
    grDevices::png(png_file, width = width, height = height,
                    units = "in", res = 150)
    plot_fn()
    grDevices::dev.off()
  }
  if (plot_format %in% c("pdf", "both")) {
    pdf_file <- file.path(plot_dir, paste0(name, ".pdf"))
    grDevices::pdf(pdf_file, width = width, height = height)
    plot_fn()
    grDevices::dev.off()
  }
}


#' Get expected plot file paths
#' @noRd
.ts_plot_paths <- function(plot_dir, name, plot_format) {
  paths <- character(0)
  if (plot_format %in% c("png", "both")) {
    paths <- c(paths, file.path(plot_dir, paste0(name, ".png")))
  }
  if (plot_format %in% c("pdf", "both")) {
    paths <- c(paths, file.path(plot_dir, paste0(name, ".pdf")))
  }
  return(paths)
}
