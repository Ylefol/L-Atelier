# GAIA/Eleuthia/dea_export.R
# Export standard DEA + PART clustering + enrichment results
#
# For non-time-series differential expression analysis workflows.


#' Export DEA + PART + Enrichment Results
#'
#' Comprehensive export function for the standard differential expression
#' analysis pipeline: DEA → gene selection → PART clustering → enrichment.
#'
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param dea_result DEA result from \code{ARTEMIS_differential_counts()}.
#'   NULL to skip DEA export. Default: NULL.
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part()}.
#'   NULL to skip PART export. Default: NULL.
#' @param gene_selection Character vector of selected gene IDs (e.g., from
#'   \code{ARTEMIS_select_de_genes()}). NULL to skip. Default: NULL.
#' @param enrichment Enrichment result from \code{APOLLO_enrich_gost()}.
#'   NULL to skip. Delegates to \code{ELEUTHIA_export_enrichment()}.
#'   Default: NULL.
#' @param l2fc_thresh Numeric. log2 fold change threshold for significant
#'   genes export. Default: 1.0.
#' @param p_thresh Numeric. Adjusted p-value threshold for significant
#'   genes export. Default: 0.05.
#' @param sample_info Data.frame with sample metadata (for heatmap annotations).
#'   Must have rownames matching sample columns. Default: NULL.
#' @param group_col Character. Group column in sample_info. Default: "group".
#' @param save_plots Logical. Generate and save plots. Default: TRUE.
#' @param plot_format Character. "png", "pdf", or "both". Default: "png".
#' @param prefix Character. Prefix for output filenames. Default: "dea".
#' @param save_rds Logical. Save full R objects as RDS. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Depending on which arguments are provided, exports:
#'
#' \strong{DEA results} (when \code{dea_result} provided):
#' \itemize{
#'   \item All results CSV (full DESeq2 output)
#'   \item Significant results CSV (filtered by thresholds)
#'   \item Metadata text file (comparison, thresholds, counts)
#' }
#'
#' \strong{PART results} (when \code{part_result} provided):
#' \itemize{
#'   \item Cluster assignments CSV (gene → cluster)
#'   \item Cluster summary CSV (cluster, n_genes, color)
#'   \item Parameters text file
#'   \item Heatmap plot (if sample_info provided)
#' }
#'
#' \strong{Gene selection} (when \code{gene_selection} provided):
#' \itemize{
#'   \item Selected genes CSV
#' }
#'
#' \strong{Enrichment} (when \code{enrichment} provided):
#' \itemize{
#'   \item Delegates to \code{ELEUTHIA_export_enrichment()} in subdirectory
#' }
#'
#' @examples
#' \dontrun{
#' # Full pipeline export
#' ELEUTHIA_export_dea_results(
#'   output_dir = "results/dea_analysis/",
#'   dea_result = dea_res,
#'   part_result = part_res,
#'   gene_selection = sig_genes,
#'   enrichment = enrich_res,
#'   sample_info = sample_sheet
#' )
#'
#' # DEA only
#' ELEUTHIA_export_dea_results(
#'   output_dir = "results/",
#'   dea_result = dea_res,
#'   prefix = "wt_vs_ko"
#' )
#'
#' }
#' @export
ELEUTHIA_export_dea_results <- function(output_dir,
                                         dea_result = NULL,
                                         part_result = NULL,
                                         gene_selection = NULL,
                                         enrichment = NULL,
                                         l2fc_thresh = 1.0,
                                         p_thresh = 0.05,
                                         sample_info = NULL,
                                         group_col = "group",
                                         save_plots = TRUE,
                                         plot_format = "png",
                                         prefix = "dea",
                                         save_rds = TRUE,
                                         verbose = TRUE) {

  # ==========================================================================
  # Validate inputs
  # ==========================================================================

  if (is.null(dea_result) && is.null(part_result) &&
      is.null(gene_selection) && is.null(enrichment)) {
    stop("At least one result object must be provided ",
         "(dea_result, part_result, gene_selection, or enrichment)")
  }

  # ==========================================================================
  # Setup directories
  # ==========================================================================

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (save_plots) {
    plot_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
  }

  files_created <- character(0)

  if (verbose) {
    cat("[ELEUTHIA] Exporting DEA Analysis Results \n")
    cat("    Output directory:", output_dir, "\n\n")
  }

  # ==========================================================================
  # DEA Results
  # ==========================================================================

  if (!is.null(dea_result)) {
    if (verbose) cat("[ELEUTHIA] DEA Results \n")

    # Check it's a valid DEA result
    if (!is.list(dea_result) || !"results" %in% names(dea_result)) {
      warning("dea_result does not appear to be from ARTEMIS_differential_counts(). ",
              "Expected a list with 'results' element.")
    } else {
      res <- dea_result$results

      # All results
      all_file <- file.path(output_dir, paste0(prefix, "_all_results.csv"))
      write.csv(res, all_file, row.names = FALSE)
      files_created <- c(files_created, all_file)
      if (verbose) cat("[ELEUTHIA]   All results:", basename(all_file), "\n")

      # Significant results
      sig_res <- res[!is.na(res$padj) &
                       res$padj < p_thresh &
                       abs(res$log2FoldChange) >= l2fc_thresh, ]

      sig_file <- file.path(output_dir, paste0(prefix, "_significant.csv"))
      write.csv(sig_res, sig_file, row.names = FALSE)
      files_created <- c(files_created, sig_file)
      if (verbose) cat("[ELEUTHIA]   Significant:", basename(sig_file),
                       "(", nrow(sig_res), "genes )\n")

      # Metadata
      meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))
      meta_lines <- c(
        "DEA Export Metadata",
        paste0("Generated: ", Sys.time()),
        "",
        "Comparison:",
        paste0("  ", dea_result$comparison %||% "Unknown"),
        "",
        "Thresholds:",
        paste0("  |log2FC| >= ", l2fc_thresh),
        paste0("  padj < ", p_thresh),
        "",
        "Summary:",
        paste0("  Total features tested: ", nrow(res)),
        paste0("  Significant (up): ", sum(sig_res$log2FoldChange > 0)),
        paste0("  Significant (down): ", sum(sig_res$log2FoldChange < 0)),
        paste0("  Total significant: ", nrow(sig_res))
      )
      writeLines(meta_lines, meta_file)
      files_created <- c(files_created, meta_file)
      if (verbose) cat("[ELEUTHIA]   Metadata:", basename(meta_file), "\n")

      # RDS
      if (save_rds) {
        rds_file <- file.path(output_dir, paste0(prefix, "_dea_result.rds"))
        saveRDS(dea_result, rds_file)
        files_created <- c(files_created, rds_file)
        if (verbose) cat("[ELEUTHIA]   RDS:", basename(rds_file), "\n")
      }

      # Volcano plot
      if (save_plots) {
        tryCatch({
          p_vol <- AETHER_plot_volcano(
            dea_result,
            l2fc_thresh = l2fc_thresh,
            p_thresh = p_thresh
          )
          vol_files <- .save_ggplot(p_vol, plot_dir,
                                    paste0(prefix, "_volcano"), plot_format,
                                    width = 8, height = 7)
          files_created <- c(files_created, vol_files)
          if (verbose) cat("[ELEUTHIA]   Volcano plot:", basename(vol_files[1]), "\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]   Warning: Could not create volcano plot -", e$message, "\n")
        })

        # MA plot
        tryCatch({
          p_ma <- AETHER_plot_ma(
            dea_result,
            l2fc_thresh = l2fc_thresh,
            p_thresh = p_thresh
          )
          ma_files <- .save_ggplot(p_ma, plot_dir,
                                   paste0(prefix, "_ma"), plot_format,
                                   width = 8, height = 7)
          files_created <- c(files_created, ma_files)
          if (verbose) cat("[ELEUTHIA]   MA plot:", basename(ma_files[1]), "\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]   Warning: Could not create MA plot -", e$message, "\n")
        })
      }
    }
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # PART Results
  # ==========================================================================

  if (!is.null(part_result)) {
    if (verbose) cat("[ELEUTHIA] PART Clustering Results \n")

    if (!inherits(part_result, "artemis_part")) {
      warning("part_result does not appear to be an artemis_part object.")
    } else {
      # Cluster assignments
      cluster_file <- file.path(output_dir, paste0(prefix, "_part_clusters.csv"))
      write.csv(part_result$cluster_map, cluster_file, row.names = FALSE)
      files_created <- c(files_created, cluster_file)
      if (verbose) cat("[ELEUTHIA]   Cluster assignments:", basename(cluster_file), "\n")

      # Cluster summary
      summary_df <- data.frame(
        cluster = names(part_result$cluster_sizes),
        n_genes = unname(part_result$cluster_sizes),
        color = unname(part_result$cluster_colors[names(part_result$cluster_sizes)])
      )
      summary_file <- file.path(output_dir, paste0(prefix, "_part_summary.csv"))
      write.csv(summary_df, summary_file, row.names = FALSE)
      files_created <- c(files_created, summary_file)
      if (verbose) cat("[ELEUTHIA]   Cluster summary:", basename(summary_file),
                       "(", part_result$n_clusters, "clusters )\n")

      # Parameters
      param_file <- file.path(output_dir, paste0(prefix, "_part_parameters.txt"))
      params <- part_result$parameters
      param_lines <- c(
        "PART Clustering Parameters",
        paste0("Generated: ", Sys.time()),
        "",
        paste0("q (splitting threshold): ", params$q),
        paste0("min_size: ", params$min_size),
        paste0("B (bootstrap iterations): ", params$B),
        paste0("Kmax: ", params$Kmax),
        paste0("dist_method: ", params$dist_method),
        paste0("linkage: ", params$linkage),
        paste0("scale: ", params$scale),
        paste0("seed: ", params$seed %||% "NULL"),
        "",
        "Results:",
        paste0("  Clusters found: ", part_result$n_clusters),
        paste0("  Outlier genes: ", part_result$n_outliers),
        paste0("  Total genes: ", nrow(part_result$cluster_map)),
        paste0("  Computation time: ", round(part_result$computation_time, 2), " sec")
      )
      writeLines(param_lines, param_file)
      files_created <- c(files_created, param_file)
      if (verbose) cat("[ELEUTHIA]   Parameters:", basename(param_file), "\n")

      # Heatmap plot
      if (save_plots && !is.null(sample_info)) {
        tryCatch({
          n_genes <- nrow(part_result$data)
          n_samples <- ncol(part_result$data)
          hm_height <- max(6, min(n_genes * 0.02 + 3, 20))
          # Width capped at height (rather than a flat value) so very large
          # sample counts can't produce an excessively wide/short image.
          hm_width  <- min(hm_height, max(8, n_samples * 0.4 + 4))

          if (plot_format %in% c("png", "both")) {
            hm_file <- file.path(plot_dir, paste0(prefix, "_part_heatmap.png"))
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              group_col = group_col,
              time_col = NULL,  # Not time series
              save_path = hm_file,
              width = hm_width,
              height = hm_height
            )
            files_created <- c(files_created, hm_file)
            if (verbose) cat("[ELEUTHIA]   Heatmap:", basename(hm_file), "\n")
          }

          if (plot_format %in% c("pdf", "both")) {
            hm_file <- file.path(plot_dir, paste0(prefix, "_part_heatmap.pdf"))
            AETHER_plot_part_heatmap(
              part_result, sample_info,
              group_col = group_col,
              time_col = NULL,
              save_path = hm_file,
              width = hm_width,
              height = hm_height
            )
            files_created <- c(files_created, hm_file)
          }
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]   Warning: Could not create heatmap -", e$message, "\n")
        })

        # Cluster means by group — compact companion to the heatmap
        tryCatch({
          p_means <- AETHER_plot_cluster_group_means(
            part_result, sample_info,
            group_col = group_col
          )
          means_files <- .save_ggplot(p_means, plot_dir,
                                      paste0(prefix, "_part_cluster_means"), plot_format,
                                      width = max(8, part_result$n_clusters * 1.2), height = 6)
          files_created <- c(files_created, means_files)
          if (verbose) cat("[ELEUTHIA]   Cluster means by group:", basename(means_files[1]), "\n")
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA]   Warning: Could not create cluster means plot -", e$message, "\n")
        })
      } else if (save_plots && is.null(sample_info)) {
        if (verbose) cat("[ELEUTHIA]   Note: sample_info not provided, skipping heatmap\n")
      }

      # RDS
      if (save_rds) {
        rds_file <- file.path(output_dir, paste0(prefix, "_part_result.rds"))
        saveRDS(part_result, rds_file)
        files_created <- c(files_created, rds_file)
        if (verbose) cat("[ELEUTHIA]   RDS:", basename(rds_file), "\n")
      }
    }
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # Gene Selection
  # ==========================================================================

  if (!is.null(gene_selection)) {
    if (verbose) cat("[ELEUTHIA] Gene Selection \n")

    genes_file <- file.path(output_dir, paste0(prefix, "_selected_genes.csv"))
    genes_df <- data.frame(gene_id = gene_selection)
    write.csv(genes_df, genes_file, row.names = FALSE)
    files_created <- c(files_created, genes_file)
    if (verbose) cat("[ELEUTHIA]   Selected genes:", basename(genes_file),
                     "(", length(gene_selection), "genes )\n\n")
  }

  # ==========================================================================
  # Enrichment
  # ==========================================================================

  if (!is.null(enrichment)) {
    if (verbose) cat("[ELEUTHIA] Enrichment Results \n")

    enrich_dir <- file.path(output_dir, "enrichment")

    tryCatch({
      enrich_files <- ELEUTHIA_export_enrichment(
        enrichment = enrichment,
        output_dir = enrich_dir,
        prefix = prefix,
        by_source = FALSE,
        by_module = TRUE,
        save_plots = save_plots,
        plot_format = plot_format,
        save_rds = save_rds,
        verbose = verbose
      )
      files_created <- c(files_created, enrich_files)
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA]   Warning: Could not export enrichment -", e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # Final summary
  # ==========================================================================

  if (verbose) {
    cat("[ELEUTHIA] Export Complete \n")
    cat("    Total files created:", length(files_created), "\n")
    cat("    Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}


# Helper for NULL coalescing (if not already available)
`%||%` <- function(a, b) if (is.null(a)) b else a
