# GAIA/Eleuthia/homer_export.R
# Export HOMER motif enrichment results
#
# Handles both homer_motif (single peak set) and homer_motif_batch objects.
# Saves known/denovo CSVs, metadata, and optional plots.


#' Export HOMER Motif Enrichment Results
#'
#' Saves HOMER motif enrichment results to CSV files, a metadata summary,
#' and optional plots. Works with both single-set (\code{homer_motif}) and
#' batch (\code{homer_motif_batch}) objects.
#'
#' @param homer_result A \code{homer_motif} or \code{homer_motif_batch}
#'   object from \code{APOLLO_homer_motif_enrichment()},
#'   \code{APOLLO_homer_motif_enrichment_batch()}, or
#'   \code{APOLLO_load_homer_results()} / \code{APOLLO_load_homer_batch()}.
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param prefix Character. Prefix for output filenames. Default: "homer".
#' @param top_n Integer. Top N motifs for plots (passed to plot functions).
#'   Default: 15.
#' @param q_thresh Numeric. q-value threshold for plot filtering. Default: 0.05.
#' @param save_plots Logical. Save diagnostic plots. Default: TRUE.
#' @param plot_format Character. Plot format: "png", "pdf", or "both".
#'   Default: "png".
#' @param save_rds Logical. Save full R object as RDS. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Exports the following files:
#' \itemize{
#'   \item \strong{Known motifs CSV} — all known motif results
#'     (\code{{prefix}_known_motifs.csv})
#'   \item \strong{De novo motifs CSV} — de novo motifs if available
#'     (\code{{prefix}_denovo_motifs.csv})
#'   \item \strong{Metadata TXT} — run parameters and summary
#'     (\code{{prefix}_metadata.txt})
#'   \item \strong{Dotplot} — top enriched motifs
#'     (via \code{AETHER_plot_motif_enrichment})
#'   \item \strong{Heatmap} — motif enrichment across peak sets, batch only
#'     (via \code{AETHER_plot_motif_comparison})
#'   \item \strong{RDS} — full object for reloading
#' }
#'
#' For batch objects, the known motifs CSV contains the combined table with a
#' \code{peak_set} column identifying each peak set.
#'
#' @examples
#' \dontrun{
#' # Single peak set
#' ELEUTHIA_export_homer_results(homer_res, "results/homer/")
#'
#' # Batch
#' ELEUTHIA_export_homer_results(homer_batch, "results/homer/",
#'                               prefix = "atac_homer")
#' }
#' @export
ELEUTHIA_export_homer_results <- function(homer_result,
                                           output_dir,
                                           prefix = "homer",
                                           top_n = 15,
                                           q_thresh = 0.05,
                                           save_plots = TRUE,
                                           plot_format = "png",
                                           save_rds = TRUE,
                                           verbose = TRUE) {

  is_batch <- inherits(homer_result, "homer_motif_batch")
  is_single <- inherits(homer_result, "homer_motif")

  if (!is_batch && !is_single) {
    stop("'homer_result' must be a homer_motif or homer_motif_batch object.",
         call. = FALSE)
  }

  if (verbose) {
    cat("[ELEUTHIA] Exporting HOMER Motif Results \n")
    cat("[ELEUTHIA] Type:", if (is_batch) "Batch" else "Single peak set", "\n")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("[ELEUTHIA] Created directory:", output_dir, "\n")
  }

  files_created <- character(0)
  meta <- homer_result$metadata

  # --------------------------------------------------------------------------
  # Known motifs CSV
  # --------------------------------------------------------------------------
  if (is_batch) {
    known_df <- homer_result$combined
  } else {
    known_df <- homer_result$known
  }

  if (!is.null(known_df) && nrow(known_df) > 0) {
    known_file <- file.path(output_dir, paste0(prefix, "_known_motifs.csv"))
    write.csv(known_df, known_file, row.names = FALSE)
    files_created <- c(files_created, known_file)
    if (verbose) cat("[ELEUTHIA]   Known motifs:", known_file,
                     "(", nrow(known_df), "rows )\n")
  } else {
    if (verbose) cat("[ELEUTHIA]   No known motif results to export.\n")
  }

  # --------------------------------------------------------------------------
  # De novo motifs CSV (single only; batch stores per-set in $results)
  # --------------------------------------------------------------------------
  denovo_df <- NULL
  if (is_single) {
    denovo_df <- homer_result$denovo
  } else if (is_batch) {
    # Collect de novo from each result in the batch
    denovo_list <- lapply(names(homer_result$results), function(nm) {
      dv <- homer_result$results[[nm]]$denovo
      if (!is.null(dv) && nrow(dv) > 0) {
        dv$peak_set <- nm
        dv
      }
    })
    denovo_list <- Filter(Negate(is.null), denovo_list)
    if (length(denovo_list) > 0) {
      denovo_df <- do.call(rbind, denovo_list)
      rownames(denovo_df) <- NULL
    }
  }

  if (!is.null(denovo_df) && nrow(denovo_df) > 0) {
    denovo_file <- file.path(output_dir, paste0(prefix, "_denovo_motifs.csv"))
    write.csv(denovo_df, denovo_file, row.names = FALSE)
    files_created <- c(files_created, denovo_file)
    if (verbose) cat("[ELEUTHIA]   De novo motifs:", denovo_file,
                     "(", nrow(denovo_df), "rows )\n")
  }

  # --------------------------------------------------------------------------
  # Metadata TXT
  # --------------------------------------------------------------------------
  meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))

  meta_lines <- c(
    "=== HOMER Motif Enrichment Results ===",
    paste("Type:", if (is_batch) "Batch" else "Single"),
    paste("Exported:", Sys.time()),
    ""
  )

  if (is_batch) {
    n_sets <- meta$n_sets %||% length(homer_result$results)
    meta_lines <- c(meta_lines,
      paste("Peak sets:", n_sets),
      paste("Set names:", paste(names(homer_result$results), collapse = ", ")),
      ""
    )
  }

  # Run parameters from metadata
  param_fields <- c("genome", "size", "mask", "denovo", "denovo_n",
                    "nproc", "bed_file", "bg_file")
  for (fld in param_fields) {
    val <- meta[[fld]]
    if (!is.null(val)) {
      meta_lines <- c(meta_lines, paste0(fld, ": ", val))
    }
  }

  # Timing
  if (!is.null(meta$elapsed_seconds)) {
    meta_lines <- c(meta_lines,
      paste("Runtime:", round(meta$elapsed_seconds, 1), "seconds")
    )
  }

  # Summary counts
  meta_lines <- c(meta_lines, "")
  if (!is.null(known_df) && nrow(known_df) > 0) {
    sig_known <- sum(!is.na(known_df$q_value) & known_df$q_value < q_thresh,
                     na.rm = TRUE)
    meta_lines <- c(meta_lines,
      paste("Known motifs total:", nrow(known_df)),
      paste(paste0("Known motifs (q < ", q_thresh, "):"), sig_known)
    )
  }
  if (!is.null(denovo_df) && nrow(denovo_df) > 0) {
    meta_lines <- c(meta_lines, paste("De novo motifs:", nrow(denovo_df)))
  }

  # Sequence counts per set
  if (is_single) {
    if (!is.null(meta$n_target_seqs)) {
      meta_lines <- c(meta_lines, paste("Target sequences:", meta$n_target_seqs))
    }
    if (!is.null(meta$n_background_seqs)) {
      meta_lines <- c(meta_lines, paste("Background sequences:", meta$n_background_seqs))
    }
  } else if (is_batch) {
    meta_lines <- c(meta_lines, "")
    meta_lines <- c(meta_lines, "Per-set sequence counts:")
    for (nm in names(homer_result$results)) {
      m <- homer_result$results[[nm]]$metadata
      n_t <- m$n_target_seqs %||% "?"
      n_b <- m$n_background_seqs %||% "?"
      meta_lines <- c(meta_lines,
        paste0("  ", nm, ": ", n_t, " targets, ", n_b, " background")
      )
    }
  }

  writeLines(meta_lines, meta_file)
  files_created <- c(files_created, meta_file)
  if (verbose) cat("[ELEUTHIA]   Metadata:", meta_file, "\n")

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  if (save_plots && !is.null(known_df) && nrow(known_df) > 0) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
    if (verbose) cat("[ELEUTHIA] Generating plots...\n")

    # Dotplot (works for both single and batch)
    tryCatch({
      p <- AETHER_plot_motif_enrichment(
        homer_result,
        top_n   = top_n,
        q_thresh = q_thresh
      )
      n_sets_approx <- if (is_batch) length(homer_result$results) else 1
      plot_files <- .save_ggplot(
        p, plot_dir, paste0(prefix, "_dotplot"), plot_format,
        width  = max(8, n_sets_approx * 1.5 + 4),
        height = max(6, top_n * 0.35 + 3)
      )
      files_created <- c(files_created, plot_files)
      if (verbose) cat("[ELEUTHIA]   Dotplot\n")
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA]   Warning: Could not create dotplot -", e$message, "\n")
    })

    # Heatmap (batch only)
    if (is_batch) {
      tryCatch({
        n_sets <- length(homer_result$results)
        .save_pheatmap(
          function() AETHER_plot_motif_comparison(
            homer_result,
            top_n    = top_n,
            q_thresh = q_thresh
          ),
          plot_dir, paste0(prefix, "_heatmap"), plot_format,
          width  = max(8, n_sets * 1.2 + 4),
          height = max(6, top_n * 0.35 + 3)
        )
        files_created <- c(files_created,
                            .plot_paths(plot_dir,
                                        paste0(prefix, "_heatmap"),
                                        plot_format))
        if (verbose) cat("[ELEUTHIA]   Heatmap (batch comparison)\n")
      }, error = function(e) {
        if (verbose) cat("[ELEUTHIA]   Warning: Could not create heatmap -", e$message, "\n")
      })
    }

    if (verbose) cat("[ELEUTHIA]   Plots saved to:", plot_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, "_result.rds"))
    saveRDS(homer_result, rds_file)
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
