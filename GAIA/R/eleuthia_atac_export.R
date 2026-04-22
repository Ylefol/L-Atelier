# GAIA/Eleuthia/atac_export.R
# Orchestrating export function for ATAC-seq analysis results
#
# Delegates to dedicated exporters for each result type:
#   - Annotated peaks         → CSV per peak set
#   - DEA results             → ELEUTHIA_export_dea_results()
#   - Enrichment              → ELEUTHIA_export_enrichment()
#   - HOMER motif results     → ELEUTHIA_export_homer_results()
#   - Sequence composition    → ELEUTHIA_export_sequence_composition()


#' Export ATAC-seq Analysis Results
#'
#' Orchestrating export function for the standard ATAC-seq analysis pipeline.
#' Delegates to dedicated exporters for each result type; any \code{NULL}
#' argument is silently skipped.
#'
#' @param output_dir Character. Root directory to save all results.
#'   Created if needed.
#' @param annotated_peaks Named list of annotated peak data.frames (e.g.,
#'   one entry per group: \code{list(WT = wt_peaks, KO = ko_peaks)}).
#'   Each data.frame is saved as \code{{name}_annotated_peaks.csv} in the
#'   \code{peaks/} subdirectory. Default: NULL.
#' @param dea_result DEA result from \code{ARTEMIS_differential_counts()}.
#'   Passed to \code{ELEUTHIA_export_dea_results()}. Default: NULL.
#' @param enrichment Enrichment result from \code{APOLLO_enrich_gost()}.
#'   Passed to \code{ELEUTHIA_export_enrichment()}. Default: NULL.
#' @param homer_result A \code{homer_motif} or \code{homer_motif_batch}
#'   object. Passed to \code{ELEUTHIA_export_homer_results()}. Default: NULL.
#' @param tss_enrichment A \code{hades_tss_enrichment} object from
#'   \code{HADES_tss_enrichment()}. Scores are saved as CSV; profile and score
#'   plots are saved via \code{AETHER_plot_tss_enrichment()}. Default: NULL.
#' @param tss_color_by Character. \code{color_by} argument passed to
#'   \code{AETHER_plot_tss_enrichment()}. Any sample sheet metadata column
#'   carried through by \code{HADES_tss_enrichment()} can be used.
#'   Default: \code{"group"}.
#' @param sequence_composition A data.frame from
#'   \code{APOLLO_sequence_composition()}. Passed to
#'   \code{ELEUTHIA_export_sequence_composition()}. Default: NULL.
#' @param prefix Character. Filename prefix for DEA, enrichment, HOMER, and
#'   composition exports. Default: "atac".
#' @param l2fc_thresh Numeric. log2FC threshold for DEA significance.
#'   Default: 1.0.
#' @param p_thresh Numeric. Adjusted p-value threshold for DEA significance.
#'   Default: 0.05.
#' @param sample_info Data.frame with sample metadata for DEA heatmap.
#'   Default: NULL.
#' @param group_col Character. Group column in \code{sample_info}.
#'   Default: \code{"group"}.
#' @param save_plots Logical. Generate and save plots for all result types.
#'   Default: TRUE.
#' @param plot_format Character. "png", "pdf", or "both". Default: "png".
#' @param save_sequences Logical. Save raw sequences from sequence composition
#'   to a separate file. Default: FALSE.
#' @param save_rds Logical. Save full R objects as RDS for all result types.
#'   Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible named list of character vectors, one element per
#'   result type exported (e.g., \code{$peaks}, \code{$dea}, \code{$homer}).
#'   Each element contains the file paths created by the corresponding exporter.
#'
#' @details
#' Output subdirectory layout:
#' \preformatted{
#' output_dir/
#'   peaks/                        # annotated peak CSVs + annotation barplot
#'   dea/                          # DEA results, plots, RDS
#'     plots/                      # volcano, MA, PART heatmap
#'   enrichment/                   # enrichment results (always top-level)
#'   homer/                        # HOMER motifs, plots, RDS
#'   composition/                  # sequence composition CSV + summary
#'   tss_enrichment/               # TSS enrichment scores CSV, plots, RDS
#' }
#'
#' @examples
#' \dontrun{
#' ELEUTHIA_export_atac_results(
#'   output_dir        = "results/atac/",
#'   annotated_peaks   = list(WT = wt_peaks, KO = ko_peaks),
#'   dea_result        = dea_res,
#'   enrichment        = enrich_res,
#'   homer_result      = homer_batch,
#'   sequence_composition = comp,
#'   prefix            = "smug1_atac"
#' )
#' }
#' @export
ELEUTHIA_export_atac_results <- function(output_dir,
                                          annotated_peaks = NULL,
                                          dea_result = NULL,
                                          enrichment = NULL,
                                          homer_result = NULL,
                                          sequence_composition = NULL,
                                          tss_enrichment = NULL,
                                          tss_color_by = "group",
                                          prefix = "atac",
                                          l2fc_thresh = 1.0,
                                          p_thresh = 0.05,
                                          sample_info = NULL,
                                          group_col = "group",
                                          save_plots = TRUE,
                                          plot_format = "png",
                                          save_sequences = FALSE,
                                          save_rds = TRUE,
                                          verbose = TRUE) {

  if (is.null(annotated_peaks) && is.null(dea_result) &&
      is.null(enrichment) && is.null(homer_result) &&
      is.null(sequence_composition) && is.null(tss_enrichment)) {
    stop("At least one result object must be provided.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) {
    cat("[ELEUTHIA] Exporting ATAC-seq Analysis Results \n")
    cat("    Output directory:", output_dir, "\n\n")
  }

  all_files <- list()

  # ==========================================================================
  # Annotated peaks (named list of data.frames)
  # ==========================================================================
  if (!is.null(annotated_peaks)) {
    if (verbose) cat("[ELEUTHIA] Annotated Peaks \n")

    if (!is.list(annotated_peaks)) {
      warning("annotated_peaks must be a named list of data.frames. Skipping.")
    } else {
      peaks_dir <- file.path(output_dir, "peaks")
      if (!dir.exists(peaks_dir)) {
        dir.create(peaks_dir)
      }

      peak_files <- character(0)
      peak_names <- names(annotated_peaks)
      if (is.null(peak_names)) {
        peak_names <- paste0("set_", seq_along(annotated_peaks))
      }

      for (i in seq_along(annotated_peaks)) {
        nm  <- peak_names[i]
        df  <- annotated_peaks[[i]]
        if (!is.data.frame(df)) {
          if (verbose) cat("[ELEUTHIA]   Skipping", nm, "(not a data.frame)\n")
          next
        }
        out_file <- file.path(peaks_dir,
                               paste0(nm, "_annotated_peaks.csv"))
        write.csv(df, out_file, row.names = FALSE)
        peak_files <- c(peak_files, out_file)
        if (verbose) cat("[ELEUTHIA]   ", nm, ":", basename(out_file),
                         "(", nrow(df), "peaks )\n")
      }

      # Annotation barplot (all peak sets combined into one plot)
      if (save_plots) {
        valid_peaks <- Filter(is.data.frame, annotated_peaks)
        if (length(valid_peaks) > 0) {
          tryCatch({
            p <- AETHER_plot_annotation_bar(valid_peaks)
            bar_files <- .save_ggplot(
              p, peaks_dir, paste0(prefix, "_annotation_bar"), plot_format,
              width  = 8,
              height = max(3, length(valid_peaks) * 0.8 + 1.5)
            )
            peak_files <- c(peak_files, bar_files)
            if (verbose) cat("[ELEUTHIA]   Annotation barplot:", basename(bar_files[1]), "\n")
          }, error = function(e) {
            if (verbose) cat("[ELEUTHIA]   Warning: Could not create annotation barplot -",
                             e$message, "\n")
          })
        }
      }

      all_files$peaks <- peak_files
    }
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # DEA Results
  # ==========================================================================
  if (!is.null(dea_result)) {
    dea_dir <- file.path(output_dir, "dea")
    tryCatch({
      dea_files <- ELEUTHIA_export_dea_results(
        output_dir  = dea_dir,
        dea_result  = dea_result,
        enrichment  = NULL,       # enrichment always goes to its own folder
        l2fc_thresh = l2fc_thresh,
        p_thresh    = p_thresh,
        sample_info = sample_info,
        group_col   = group_col,
        save_plots  = save_plots,
        plot_format = plot_format,
        prefix      = prefix,
        save_rds    = save_rds,
        verbose     = verbose
      )
      all_files$dea <- dea_files
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: DEA export failed -", e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # Enrichment (always top-level, independent of DEA)
  # ==========================================================================
  if (!is.null(enrichment)) {
    enrich_dir <- file.path(output_dir, "enrichment")
    tryCatch({
      enrich_files <- ELEUTHIA_export_enrichment(
        enrichment  = enrichment,
        output_dir  = enrich_dir,
        prefix      = prefix,
        by_module   = TRUE,
        save_plots  = save_plots,
        plot_format = plot_format,
        save_rds    = save_rds,
        verbose     = verbose
      )
      all_files$enrichment <- enrich_files
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: Enrichment export failed -", e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # HOMER results
  # ==========================================================================
  if (!is.null(homer_result)) {
    homer_dir <- file.path(output_dir, "homer")
    tryCatch({
      homer_files <- ELEUTHIA_export_homer_results(
        homer_result = homer_result,
        output_dir   = homer_dir,
        prefix       = prefix,
        save_plots   = save_plots,
        plot_format  = plot_format,
        save_rds     = save_rds,
        verbose      = verbose
      )
      all_files$homer <- homer_files
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: HOMER export failed -", e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # Sequence composition
  # ==========================================================================
  if (!is.null(sequence_composition)) {
    comp_dir <- file.path(output_dir, "composition")
    tryCatch({
      comp_files <- ELEUTHIA_export_sequence_composition(
        composition    = sequence_composition,
        output_dir     = comp_dir,
        prefix         = prefix,
        save_sequences = save_sequences,
        save_rds       = save_rds,
        verbose        = verbose
      )
      all_files$composition <- comp_files
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: Composition export failed -", e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # TSS enrichment QC
  # ==========================================================================
  if (!is.null(tss_enrichment)) {
    tss_dir <- file.path(output_dir, "tss_enrichment")
    if (!dir.exists(tss_dir)) dir.create(tss_dir, recursive = TRUE)

    tryCatch({
      if (verbose) cat("[ELEUTHIA] TSS Enrichment QC\n")
      tss_files <- character(0)

      # Scores CSV
      scores_file <- file.path(tss_dir, paste0(prefix, "_tss_enrichment_scores.csv"))
      write.csv(tss_enrichment$scores, scores_file, row.names = FALSE)
      tss_files <- c(tss_files, scores_file)
      if (verbose) cat("    Scores:", basename(scores_file), "\n")

      # Plots
      if (save_plots) {
        plts <- AETHER_plot_tss_enrichment(tss_enrichment, color_by = tss_color_by)

        profile_files <- .save_ggplot(
          plts$profile, tss_dir, paste0(prefix, "_tss_profile"),
          plot_format, width = 8, height = 4
        )
        score_files <- .save_ggplot(
          plts$score, tss_dir, paste0(prefix, "_tss_scores"),
          plot_format, width = 6, height = 4
        )
        tss_files <- c(tss_files, profile_files, score_files)
        if (verbose) {
          cat("    Profile plot:", basename(profile_files[1]), "\n")
          cat("    Score plot:  ", basename(score_files[1]),  "\n")
        }
      }

      # RDS
      if (save_rds) {
        rds_file <- file.path(tss_dir, paste0(prefix, "_tss_enrichment.rds"))
        saveRDS(tss_enrichment, rds_file)
        tss_files <- c(tss_files, rds_file)
        if (verbose) cat("    RDS:", basename(rds_file), "\n")
      }

      all_files$tss_enrichment <- tss_files
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: TSS enrichment export failed -",
                       e$message, "\n")
    })
    if (verbose) cat("\n")
  }

  # ==========================================================================
  # Final summary
  # ==========================================================================
  total_files <- sum(lengths(all_files))
  if (verbose) {
    cat("[ELEUTHIA] Export Complete \n")
    cat("    Total files created:", total_files, "\n")
    cat("    Output directory:", output_dir, "\n")
  }

  invisible(all_files)
}
