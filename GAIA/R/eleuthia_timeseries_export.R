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
#' @param part_divergence An \code{artemis_part_divergence} object from
#'   \code{ARTEMIS_part_trajectory_divergence()}, or NULL to skip. Only used
#'   alongside \code{part_result} (same PART plots directory/prefix); exports
#'   the results table and, when \code{save_plots = TRUE}, the
#'   \code{AETHER_plot_cluster_divergence()} ranking plot. Default: NULL.
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
#' @param group_colors Named character vector of colors for groups (names must
#'   match values in \code{sample_info[[group_col]]}). Forwarded to
#'   \code{AETHER_plot_part_heatmap()} and \code{AETHER_plot_cluster_trajectories()}.
#'   NULL = each function's own default (auto-assigned by first-appearance
#'   order of the group column). Default: NULL.
#' @param save_plots Logical. Generate and save plots. Default: TRUE.
#' @param save_full_gostplot Logical. When \code{enrichment} is exported, also
#'   generate the full (unfiltered) interactive g:GOSt Manhattan plot per
#'   module/cluster -- forwarded as-is to \code{ELEUTHIA_export_enrichment()}.
#'   Default: TRUE.
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
                                                part_divergence = NULL,
                                                ts_de = NULL,
                                                gene_selection = NULL,
                                                enrichment = NULL,
                                                sample_info = NULL,
                                                sample_col = NULL,
                                                time_col = "timepoint",
                                                group_col = "group",
                                                group_colors = NULL,
                                                save_plots = TRUE,
                                                save_full_gostplot = TRUE,
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

    # Cluster trajectory divergence (optional -- ARTEMIS_part_trajectory_divergence())
    if (!is.null(part_divergence)) {
      if (!inherits(part_divergence, "artemis_part_divergence")) {
        stop("'part_divergence' must be an artemis_part_divergence object ",
             "from ARTEMIS_part_trajectory_divergence()")
      }
      divergence_file <- file.path(output_dir, paste0(prefix, "_cluster_divergence.csv"))
      write.csv(part_divergence$results, divergence_file, row.names = FALSE)
      files_created <- c(files_created, divergence_file)
      if (verbose) cat("[ELEUTHIA]   Cluster divergence:", divergence_file, "\n")

      if (save_rds) {
        divergence_rds <- file.path(output_dir, paste0(prefix, "_cluster_divergence.rds"))
        saveRDS(part_divergence, divergence_rds)
        files_created <- c(files_created, divergence_rds)
        if (verbose) cat("[ELEUTHIA]   Cluster divergence RDS:", divergence_rds, "\n")
      }
    }

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
          hm_height <- max(6, min(n_genes * 0.02 + 3, 20))
          # Width capped at height (rather than a flat value) so very large
          # sample counts can't produce an excessively wide/short image.
          hm_width  <- min(hm_height, max(8, n_samples * 0.4 + 4))

          if (plot_format %in% c("png", "both")) {
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              sample_col = sample_col,
              group_col = group_col, time_col = time_col,
              group_colors = group_colors,
              save_path = file.path(plot_dir, "part_heatmap.png"),
              width = hm_width, height = hm_height
            )
          }
          if (plot_format %in% c("pdf", "both")) {
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              sample_col = sample_col,
              group_col = group_col, time_col = time_col,
              group_colors = group_colors,
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
            time_col = time_col, group_col = group_col,
            colors = group_colors
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

      # Cluster trajectory divergence ranking (optional, no sample_info needed --
      # already baked into part_divergence)
      if (!is.null(part_divergence)) {
        tryCatch({
          p <- AETHER_plot_cluster_divergence(part_divergence)
          .ts_save_ggplot(p, plot_dir, "cluster_divergence", plot_format,
                           width = 8, height = max(4, part_result$n_clusters * 0.35 + 2))
          files_created <- c(files_created,
                              .ts_plot_paths(plot_dir, "cluster_divergence", plot_format))
          if (verbose) cat("[ELEUTHIA]     Cluster divergence ranking\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]     Warning: Could not create cluster divergence plot -", e$message, "\n")
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
        save_full_gostplot = save_full_gostplot,
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


#' Build a browsable table of time series DE results across comparisons
#'
#' @description Returns a wide data.frame -- one row per feature significant
#' in at least one comparison, with a \code{log2FoldChange}/\code{padj}
#' column pair per comparison -- meant to be rendered with
#' \code{DT::datatable()} directly in a report, so a viewer can search/sort
#' across every comparison's result for one feature at a glance. Mirrors
#' CAULDRON's \code{TALARIA_build_de_comparison_table()} (same shape, same
#' intended \code{DT::datatable()} pairing), adapted to accept ARTEMIS's own
#' \code{artemis_ts_de} result objects directly instead of requiring the
#' caller to pre-shape a \code{gene}-column data.frame per comparison.
#'
#' @param ts_de One \code{artemis_ts_de} object (from
#'   \code{ARTEMIS_timeseries_conditional_limma()} or
#'   \code{ARTEMIS_timeseries_temporal_limma()}), or a list of them -- the
#'   same input shape accepted by
#'   \code{\link{ELEUTHIA_export_timeseries_results}}'s \code{ts_de}
#'   argument. Every comparison across all supplied objects is flattened
#'   into its own pair of output columns, named \code{"<type>_<comparison>
#'   log2FC"/"padj"} when more than one \code{artemis_ts_de} object is
#'   supplied (so e.g. a conditional and a temporal comparison sharing a
#'   name don't collide), or just \code{"<comparison> log2FC"/"padj"} when
#'   only one is supplied.
#' @param id_map Named character vector (\code{feature_id -> gene/protein
#'   symbol}), or \code{NULL} (default) to leave \code{Feature} as the only
#'   row identifier. When supplied, an additional \code{Symbol} column is
#'   added (\code{NA} for any \code{feature_id} absent from \code{id_map}).
#' @param sig_col Character. Column used to decide which features are
#'   included (i.e. significant in at least one comparison). Default
#'   \code{"padj"}. Set to \code{"pvalue"} for pipelines that call
#'   significance on the raw p-value instead of FDR.
#' @param sig_thresh Numeric. Cutoff applied to \code{sig_col}. Default
#'   \code{0.05}.
#' @param l2fc_thresh Numeric. Absolute log2 fold change cutoff used to
#'   decide which features are included. Default \code{0.5}.
#' @param digits Integer. Rounding applied to \code{log2FoldChange} columns
#'   (\code{padj} columns are kept to 3 significant figures regardless).
#'   Default \code{3}.
#'
#' @return A data.frame with a \code{Feature} column (and a \code{Symbol}
#'   column if \code{id_map} given), plus a \code{"<comparison> log2FC"} /
#'   \code{"<comparison> padj"} column pair per comparison.
#'
#' @seealso \code{\link{ELEUTHIA_export_timeseries_results}}
#'
#' @examples
#' \dontrun{
#' tbl <- ELEUTHIA_build_de_comparison_table(list(cond_de, temp_de), id_map = id_map)
#' DT::datatable(tbl, rownames = FALSE)
#' }
#' @export
ELEUTHIA_build_de_comparison_table <- function(ts_de,
                                                id_map      = NULL,
                                                sig_col     = "padj",
                                                sig_thresh  = 0.05,
                                                l2fc_thresh = 0.5,
                                                digits      = 3) {

  if (inherits(ts_de, "artemis_ts_de")) {
    ts_de_list <- list(ts_de)
  } else if (is.list(ts_de)) {
    ts_de_list <- ts_de
  } else {
    stop("[ELEUTHIA] 'ts_de' must be an artemis_ts_de object or a list of them.", call. = FALSE)
  }
  for (x in ts_de_list) {
    if (!inherits(x, "artemis_ts_de"))
      stop("[ELEUTHIA] Every element of 'ts_de' must be an artemis_ts_de object.", call. = FALSE)
  }

  multi_type <- length(ts_de_list) > 1
  de_list <- list()
  for (ts_obj in ts_de_list) {
    for (exp_name in names(ts_obj$results)) {
      label <- if (multi_type) paste0(ts_obj$type, "_", exp_name) else exp_name
      de_list[[label]] <- ts_obj$results[[exp_name]]$results
    }
  }

  req_cols <- unique(c("feature_id", "log2FoldChange", "padj", sig_col))
  for (nm in names(de_list)) {
    missing_cols <- setdiff(req_cols, names(de_list[[nm]]))
    if (length(missing_cols) > 0)
      stop("[ELEUTHIA] Comparison '", nm, "' is missing column(s): ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  comparisons <- names(de_list)

  # Union of features significant in at least one comparison -- every value
  # in a selected/listed feature's row is still shown for every comparison,
  # significant there or not.
  sig_features <- unique(unlist(lapply(de_list, function(df) {
    mask <- !is.na(df[[sig_col]]) & !is.na(df$log2FoldChange) &
      df[[sig_col]] < sig_thresh & abs(df$log2FoldChange) >= l2fc_thresh
    df$feature_id[mask]
  })))
  sig_features <- sort(sig_features)

  if (length(sig_features) == 0)
    stop("[ELEUTHIA] No feature was significant (", sig_col, " < ", sig_thresh,
         ", |log2FoldChange| >= ", l2fc_thresh, ") in any comparison in 'ts_de'.", call. = FALSE)

  # feature x comparison matrices, NA where a feature wasn't tested in that
  # comparison.
  log2fc_mat <- matrix(NA_real_, nrow = length(sig_features), ncol = length(comparisons),
                        dimnames = list(sig_features, comparisons))
  padj_mat <- log2fc_mat

  for (nm in comparisons) {
    df  <- de_list[[nm]]
    idx <- match(sig_features, df$feature_id)
    hit <- !is.na(idx)
    log2fc_mat[hit, nm] <- df$log2FoldChange[idx[hit]]
    padj_mat[hit, nm]   <- df$padj[idx[hit]]
  }

  out <- data.frame(Feature = sig_features, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(id_map)) out$Symbol <- unname(id_map[sig_features])
  for (nm in comparisons) {
    out[[paste0(nm, " log2FC")]] <- round(log2fc_mat[, nm], digits)
    out[[paste0(nm, " padj")]]   <- signif(padj_mat[, nm], 3)
  }
  rownames(out) <- NULL
  out
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
