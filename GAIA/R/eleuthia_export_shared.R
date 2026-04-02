# GAIA/Eleuthia/export_shared.R
# Shared/generic export functions used across multiple pipelines
#
# Contains export utilities that are not specific to a single analysis type
# (WGCNA, TiSA, decoupleR, CIBERSORT, etc.) but can be used by any of them.
#
# Includes shared plot-saving helpers (.save_ggplot, .save_pheatmap, .plot_paths)
# used by activity_export.R, cibersort_export.R, and other export modules.


# ==============================================================================
# SHARED PLOT SAVE HELPERS
# ==============================================================================

#' Save a ggplot to file(s)
#' @param p A ggplot object.
#' @param plot_dir Directory to save to.
#' @param name Base filename (without extension).
#' @param plot_format "png", "pdf", or "both".
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @return Character vector of saved file paths.
#' @noRd
.save_ggplot <- function(p, plot_dir, name, plot_format,
                          width = 10, height = 8) {
  saved <- character(0)
  if (plot_format %in% c("png", "both")) {
    png_file <- file.path(plot_dir, paste0(name, ".png"))
    ggplot2::ggsave(png_file, p, width = width, height = height, dpi = 150, bg = "white")
    saved <- c(saved, png_file)
  }
  if (plot_format %in% c("pdf", "both")) {
    pdf_file <- file.path(plot_dir, paste0(name, ".pdf"))
    ggplot2::ggsave(pdf_file, p, width = width, height = height)
    saved <- c(saved, pdf_file)
  }
  return(saved)
}


#' Save a pheatmap to file(s)
#'
#' pheatmap objects draw on creation, so we call the plot function inside
#' a graphics device to capture the output.
#' @param plot_fn A function that creates the pheatmap (called inside device).
#' @param plot_dir Directory to save to.
#' @param name Base filename (without extension).
#' @param plot_format "png", "pdf", or "both".
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @noRd
.save_pheatmap <- function(plot_fn, plot_dir, name, plot_format,
                            width = 10, height = 8) {
  if (plot_format %in% c("png", "both")) {
    png_file <- file.path(plot_dir, paste0(name, ".png"))
    grDevices::png(png_file, width = width, height = height, units = "in", res = 150)
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


#' Get expected plot file paths for a given name and format
#' @param plot_dir Directory containing plots.
#' @param name Base filename (without extension).
#' @param plot_format "png", "pdf", or "both".
#' @return Character vector of expected file paths.
#' @noRd
.plot_paths <- function(plot_dir, name, plot_format) {
  paths <- character(0)
  if (plot_format %in% c("png", "both")) {
    paths <- c(paths, file.path(plot_dir, paste0(name, ".png")))
  }
  if (plot_format %in% c("pdf", "both")) {
    paths <- c(paths, file.path(plot_dir, paste0(name, ".pdf")))
  }
  return(paths)
}


# ==============================================================================
# SHARED EXPORT FUNCTIONS
# ==============================================================================


#' Export enrichment results to files
#'
#' Exports gprofiler2 enrichment results (from APOLLO_enrich_gost) to CSV files
#' and optionally generates dotplots for each source.
#'
#' @param enrichment A gost_enrichment object from APOLLO_enrich_gost().
#' @param output_dir Output directory path.
#' @param prefix Prefix for output filenames. Default: "enrichment".
#' @param by_source Logical. Create separate files per source (GO:BP, KEGG, etc.). Default: FALSE.
#' @param by_module Logical. Create separate files per module/cluster (all sources
#'   combined per query). Default: TRUE.
#' @param save_plots Logical. Generate dotplots per source. Default: TRUE.
#' @param plot_top_n Integer. Number of top terms per module in dotplots. Default: 10.
#' @param plot_format Character. Plot file format: "png", "pdf", or "both". Default: "png".
#' @param plot_width Numeric. Plot width in inches. Default: 10.
#' @param plot_height_per_module Numeric. Plot height per module in inches. Default: 3.
#' @param save_rds Logical. Save full R object as RDS. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible list of file paths created.
#'
#' @details
#' Creates the following files:
#' - <prefix>_combined.csv: All results in one table
#' - <prefix>_summary.csv: Summary counts per module/source
#' - <prefix>_by_module/<module>.csv: Results split by module/cluster (default)
#' - <prefix>_by_source/<source>.csv: Results split by source (if by_source = TRUE)
#' - <prefix>_plots/dotplot_<source>.png: Dotplot per source (if save_plots = TRUE)
#' - <prefix>.rds: Full gost_enrichment object (if save_rds = TRUE)
#'
#' The "module" column contains whatever labels were used in the original query
#' (e.g., "C1", "C2" for PART clusters, "blue", "turquoise" for WGCNA modules).
#'
#' @examples
#' \dontrun{
#' enrich <- APOLLO_enrich_gost(module_genes)
#' ELEUTHIA_export_enrichment(enrich, "results/enrichment")
#'
#' # Split by source instead of module
#' ELEUTHIA_export_enrichment(enrich, "results/enrichment",
#'                            by_module = FALSE, by_source = TRUE)
#'
#' }
#' @export
ELEUTHIA_export_enrichment <- function(enrichment,
                                        output_dir,
                                        prefix = "enrichment",
                                        by_source = FALSE,
                                        by_module = TRUE,
                                        save_plots = TRUE,
                                        plot_top_n = 10,
                                        plot_format = "png",
                                        plot_width = 10,
                                        plot_height_per_module = 3,
                                        save_rds = TRUE,
                                        verbose = TRUE) {

  if (!inherits(enrichment, "gost_enrichment")) {
    stop("enrichment must be a gost_enrichment object from APOLLO_enrich_gost()")
  }

  if (verbose) cat("[ELEUTHIA] Exporting Enrichment Results \n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("[ELEUTHIA] Created directory:", output_dir, "\n")
  }

  files_created <- character(0)
  combined <- enrichment$combined

  # --------------------------------------------------------------------------
  # Handle list columns (gprofiler2 returns some columns as lists)
  # Convert to comma-separated strings for CSV export
  # --------------------------------------------------------------------------
  .flatten_list_cols <- function(df) {
    if (nrow(df) == 0) return(df)

    list_cols <- sapply(df, is.list)
    for (col in names(list_cols)[list_cols]) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x) || length(x) == 0) {
          NA_character_
        } else {
          paste(x, collapse = ",")
        }
      })
    }
    return(df)
  }

  # --------------------------------------------------------------------------
  # Combined results
  # --------------------------------------------------------------------------
  if (nrow(combined) > 0) {
    combined_file <- file.path(output_dir, paste0(prefix, "_combined.csv"))
    write.csv(.flatten_list_cols(combined), combined_file, row.names = FALSE)
    files_created <- c(files_created, combined_file)
    if (verbose) cat("[ELEUTHIA] Combined results:", combined_file, "\n")
  } else {
    if (verbose) cat("[ELEUTHIA] No significant results to export.\n")
  }

  # --------------------------------------------------------------------------
  # Summary
  # --------------------------------------------------------------------------
  summary_file <- file.path(output_dir, paste0(prefix, "_summary.csv"))
  write.csv(enrichment$summary, summary_file, row.names = FALSE)
  files_created <- c(files_created, summary_file)
  if (verbose) cat("[ELEUTHIA] Summary:", summary_file, "\n")

  # --------------------------------------------------------------------------
  # By source
  # --------------------------------------------------------------------------
  if (by_source && nrow(combined) > 0) {
    source_dir <- file.path(output_dir, paste0(prefix, "_by_source"))
    if (!dir.exists(source_dir)) {
      dir.create(source_dir)
    }

    sources <- unique(combined$source)
    for (src in sources) {
      src_df <- combined[combined$source == src, ]
      if (nrow(src_df) > 0) {
        # Clean source name for filename (GO:BP -> GO_BP)
        src_clean <- gsub(":", "_", src)
        src_file <- file.path(source_dir, paste0(src_clean, ".csv"))
        write.csv(.flatten_list_cols(src_df), src_file, row.names = FALSE)
        files_created <- c(files_created, src_file)
      }
    }
    if (verbose) cat("[ELEUTHIA] By source:", source_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # By module
  # --------------------------------------------------------------------------
  if (by_module && nrow(combined) > 0) {
    module_dir <- file.path(output_dir, paste0(prefix, "_by_module"))
    if (!dir.exists(module_dir)) {
      dir.create(module_dir)
    }

    modules <- unique(combined$module)
    for (mod in modules) {
      mod_df <- combined[combined$module == mod, ]
      if (nrow(mod_df) > 0) {
        mod_file <- file.path(module_dir, paste0(mod, ".csv"))
        write.csv(.flatten_list_cols(mod_df), mod_file, row.names = FALSE)
        files_created <- c(files_created, mod_file)
      }
    }
    if (verbose) cat("    By module:", module_dir, "/\n")
  }

  # --------------------------------------------------------------------------
  # Generate dotplots per source
  # --------------------------------------------------------------------------
  if (save_plots && nrow(combined) > 0) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }

    sources <- unique(combined$source)

    if (verbose) cat("    Generating dotplots...\n")

    for (src in sources) {
      src_df <- combined[combined$source == src, ]
      if (nrow(src_df) == 0) next

      # Clean source name for filename
      src_clean <- gsub(":", "_", src)

      # Calculate height based on wrapped line count of unique terms.
      # Terms are wrapped at width = 50 (matches AETHER_plot_gost_dotplot default).
      # A term wrapping to N lines occupies N times the vertical space of a
      # single-line term, so we sum lines rather than count terms.
      src_df_ordered <- src_df[order(src_df$p_value), ]
      top_terms <- do.call(rbind, lapply(split(src_df_ordered, src_df_ordered$module), head, plot_top_n))
      unique_term_names <- unique(top_terms$term_name)
      total_lines <- sum(sapply(unique_term_names, function(nm) {
        length(strwrap(nm, width = 50))
      }))
      # ~0.3 inches per line, minimum 6 inches, maximum 24 inches
      plot_height <- min(24, max(6, total_lines * 0.3))

      # Generate plot using AETHER function
      p <- tryCatch({
        AETHER_plot_gost_dotplot(
          gost_result = src_df,
          source = src,
          top_n = plot_top_n
        )
      }, error = function(e) {
        if (verbose) cat("[ELEUTHIA] Warning: Could not create plot for", src, "-", e$message, "\n")
        NULL
      })

      if (!is.null(p)) {
        # Save in requested format(s)
        if (plot_format %in% c("png", "both")) {
          png_file <- file.path(plot_dir, paste0("dotplot_", src_clean, ".png"))
          ggplot2::ggsave(png_file, p, width = plot_width, height = plot_height,
                          dpi = 150, bg = "white")
          files_created <- c(files_created, png_file)
        }

        if (plot_format %in% c("pdf", "both")) {
          pdf_file <- file.path(plot_dir, paste0("dotplot_", src_clean, ".pdf"))
          ggplot2::ggsave(pdf_file, p, width = plot_width, height = plot_height)
          files_created <- c(files_created, pdf_file)
        }

        if (verbose) cat("  ", src, "\n", sep = "")
      }
    }

    if (verbose) cat("[ELEUTHIA] Plots saved to:", plot_dir, "/\n")

    # ------------------------------------------------------------------------
    # GO DAG plots (BP, MF, CC) — one per ontology if present in results
    # ------------------------------------------------------------------------
    go_ont_map <- c("GO:BP" = "BP", "GO:MF" = "MF", "GO:CC" = "CC")
    for (go_src in names(go_ont_map)) {
      if (!go_src %in% sources) next
      ont <- go_ont_map[[go_src]]
      if (verbose) cat("[ELEUTHIA] Generating GO DAG (", ont, ")...\n", sep = "")

      p_dag <- tryCatch(
        AETHER_plot_go_dag(enrichment, ont = ont, verbose = FALSE),
        error = function(e) {
          if (verbose) cat("[ELEUTHIA] Warning: GO DAG (", ont, ") failed - ", e$message, "\n", sep = "")
          NULL
        }
      )

      if (!is.null(p_dag)) {
        dag_base <- file.path(plot_dir, paste0("go_dag_", tolower(ont)))
        if (plot_format %in% c("png", "both")) {
          f <- paste0(dag_base, ".png")
          ggplot2::ggsave(f, p_dag, width = 14, height = 10, dpi = 150, bg = "white")
          files_created <- c(files_created, f)
        }
        if (plot_format %in% c("pdf", "both")) {
          f <- paste0(dag_base, ".pdf")
          ggplot2::ggsave(f, p_dag, width = 14, height = 10)
          files_created <- c(files_created, f)
        }
        if (verbose) cat("[ELEUTHIA] GO DAG (", ont, ") saved\n", sep = "")
      }
    }

    # ------------------------------------------------------------------------
    # Enrichment map per non-GO source
    # ------------------------------------------------------------------------
    non_go_sources <- sources[!grepl("^GO:", sources)]
    for (src in non_go_sources) {
      if (nrow(combined[combined$source == src, ]) == 0L) next
      if (verbose) cat("[ELEUTHIA] Generating enrichment map (", src, ")...\n", sep = "")
      src_clean <- gsub(":", "_", src)

      p_emap <- tryCatch(
        AETHER_plot_enrichment_map(enrichment, source = src, verbose = FALSE),
        error = function(e) {
          if (verbose) cat("[ELEUTHIA] Warning: EMAP (", src, ") failed - ", e$message, "\n", sep = "")
          NULL
        }
      )

      if (!is.null(p_emap)) {
        emap_base <- file.path(plot_dir, paste0("emap_", src_clean))
        if (plot_format %in% c("png", "both")) {
          f <- paste0(emap_base, ".png")
          ggplot2::ggsave(f, p_emap, width = 12, height = 10, dpi = 150, bg = "white")
          files_created <- c(files_created, f)
        }
        if (plot_format %in% c("pdf", "both")) {
          f <- paste0(emap_base, ".pdf")
          ggplot2::ggsave(f, p_emap, width = 12, height = 10)
          files_created <- c(files_created, f)
        }
        if (verbose) cat("[ELEUTHIA] EMAP (", src, ") saved\n", sep = "")
      }
    }
  }

  # --------------------------------------------------------------------------
  # Save RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, ".rds"))
    saveRDS(enrichment, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("[ELEUTHIA] RDS:", rds_file, "\n")
  }

  if (verbose) {
    cat("\n[ELEUTHIA]  Export Summary \n")
    cat("    Total files created:", length(files_created), "\n")
    cat("    Total terms exported:", nrow(combined), "\n")
  }

  invisible(files_created)
}
