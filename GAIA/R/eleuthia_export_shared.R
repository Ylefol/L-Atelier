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
#'   combined per query). Unlike \code{combined}/\code{by_source}, these hold
#'   EVERY evaluated term (not just significance-filtered ones) with a
#'   \code{significant} column to filter on, since \code{enrichment$results}
#'   already stores the full unfiltered per-module gost output. Default: TRUE.
#' @param save_plots Logical. Generate dotplots per source. Default: TRUE.
#' @param save_full_gostplot Logical. Also generate one traditional, interactive
#'   g:GOSt Manhattan plot per module (saved as HTML), showing ALL evaluated
#'   terms regardless of significance. Unlike the significance-filtered dotplots,
#'   this is generated even for modules with zero significant terms. Default: TRUE.
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
#' - <prefix>_combined.csv: Significance-filtered results in one table
#' - <prefix>_summary.csv: Summary counts per module/source
#' - <prefix>_by_module/<module>.csv: ALL evaluated terms per module/cluster,
#'   unfiltered for significance (has a `significant` column) (default)
#' - <prefix>_by_source/<source>.csv: Significance-filtered results split by source
#'   (if by_source = TRUE)
#' - <prefix>_plots/dotplot_<source>.png: Dotplot per source (if save_plots = TRUE)
#' - <prefix>_plots/gostplot_full_<module>.html: Interactive, unfiltered g:GOSt
#'   Manhattan plot per module (if save_full_gostplot = TRUE)
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
                                        save_full_gostplot = TRUE,
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
  # Pulled from enrichment$results (the FULL, unfiltered per-module gost
  # objects), not from `combined` — so these CSVs hold every evaluated term
  # (with its own `significant` column to filter on later), not just the
  # significance-filtered subset. Still respects min_term_size/max_term_size.
  # --------------------------------------------------------------------------
  if (by_module && length(enrichment$results) > 0) {
    module_dir <- file.path(output_dir, paste0(prefix, "_by_module"))
    if (!dir.exists(module_dir)) {
      dir.create(module_dir)
    }

    min_ts <- enrichment$metadata$min_term_size %||% 1
    max_ts <- enrichment$metadata$max_term_size %||% Inf

    for (mod in names(enrichment$results)) {
      mod_df <- enrichment$results[[mod]]$result
      if (is.null(mod_df) || nrow(mod_df) == 0) next

      mod_df <- mod_df[mod_df$term_size >= min_ts & mod_df$term_size <= max_ts, ]
      if (nrow(mod_df) > 0) {
        mod_df$module <- mod
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
  # Full (unfiltered) gostplots — one interactive HTML per module
  # Independent of `combined`/significance filtering: enrichment$results holds
  # every evaluated term per module, so modules with zero significant terms
  # still get a plot (the whole point — seeing non-significant clusters' biology).
  # --------------------------------------------------------------------------
  if (save_plots && save_full_gostplot && length(enrichment$results) > 0) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    if (verbose) cat("[ELEUTHIA] Generating full (unfiltered) gostplots...\n")

    for (mod in names(enrichment$results)) {
      html_file <- file.path(plot_dir, paste0("gostplot_full_", mod, ".html"))
      ok <- tryCatch({
        AETHER_plot_gost_full(enrichment$results[[mod]], save_html = html_file)
        TRUE
      }, error = function(e) {
        if (verbose) cat("[ELEUTHIA] Warning: full gostplot failed for", mod, "-", e$message, "\n")
        FALSE
      })
      if (ok) {
        files_created <- c(files_created, html_file)
        if (verbose) cat("  ", mod, "\n", sep = "")
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


#' Export GSEA results to files
#'
#' Exports fgsea results (from \code{APOLLO_gsea()}) to CSV files and
#' optionally generates an NES dotplot and per-pathway running-score
#' enrichment plots for significant pathways.
#'
#' @param gsea_result An \code{apollo_gsea} object from \code{APOLLO_gsea()}.
#' @param output_dir Output directory path.
#' @param prefix Prefix for output filenames. Default: "gsea".
#' @param plot_top_n Integer. Number of top pathways (by padj, split
#'   enriched/depleted) shown in the dotplot. Default: 20.
#' @param save_plots Logical. Generate dotplot and enrichment plots. Default: TRUE.
#' @param do_enrichment_plots Logical. Generate per-pathway running-score
#'   plots for significant pathways. Default: TRUE.
#' @param max_enrichment_plots Integer. Cap on the number of per-pathway
#'   plots generated (ordered by padj), to avoid producing hundreds of files
#'   when many pathways are significant. Default: 20.
#' @param plot_format Character. "png", "pdf", or "both". Default: "png".
#' @param plot_width Numeric. Dotplot width in inches. Default: 10.
#' @param save_rds Logical. Save full apollo_gsea object as RDS. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Creates the following files:
#' - <prefix>_results.csv: full fgsea results (all tested gene sets)
#' - <prefix>_significant.csv: results filtered to the FDR threshold used in
#'   APOLLO_gsea()
#' - <prefix>_plots/dotplot.<ext>: NES dotplot (top plot_top_n pathways)
#' - <prefix>_plots/enrichment_<pathway>.<ext>: running-score plot per
#'   significant pathway (up to max_enrichment_plots), if
#'   do_enrichment_plots = TRUE
#' - <prefix>.rds: full apollo_gsea object (if save_rds = TRUE)
#'
#' @examples
#' \dontrun{
#' ranked <- APOLLO_rank_from_de(de, rank_by = "t")
#' gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))
#' ELEUTHIA_export_gsea(gsea_result, "results/gsea")
#' }
#' @export
ELEUTHIA_export_gsea <- function(gsea_result,
                                  output_dir,
                                  prefix               = "gsea",
                                  plot_top_n           = 20,
                                  save_plots           = TRUE,
                                  do_enrichment_plots  = TRUE,
                                  max_enrichment_plots = 20,
                                  plot_format          = "png",
                                  plot_width           = 10,
                                  save_rds             = TRUE,
                                  verbose              = TRUE) {

  if (!inherits(gsea_result, "apollo_gsea")) {
    stop("gsea_result must be an apollo_gsea object from APOLLO_gsea()")
  }

  if (verbose) cat("[ELEUTHIA] Exporting GSEA Results\n")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("[ELEUTHIA] Created directory:", output_dir, "\n")
  }

  files_created <- character(0)
  results_df    <- gsea_result$results

  # --------------------------------------------------------------------------
  # Flatten any list columns before CSV export
  # --------------------------------------------------------------------------
  .flatten_list_cols <- function(df) {
    if (nrow(df) == 0) return(df)
    list_cols <- sapply(df, is.list)
    for (col in names(list_cols)[list_cols]) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ",")
      })
    }
    df
  }

  # --------------------------------------------------------------------------
  # Full + significant results CSVs
  # --------------------------------------------------------------------------
  results_file <- file.path(output_dir, paste0(prefix, "_results.csv"))
  write.csv(.flatten_list_cols(results_df), results_file, row.names = FALSE)
  files_created <- c(files_created, results_file)
  if (verbose) cat("[ELEUTHIA] Results:", results_file, "\n")

  sig_file <- file.path(output_dir, paste0(prefix, "_significant.csv"))
  write.csv(.flatten_list_cols(gsea_result$significant), sig_file, row.names = FALSE)
  files_created <- c(files_created, sig_file)
  if (verbose) cat("[ELEUTHIA] Significant:", sig_file, "\n")

  # --------------------------------------------------------------------------
  # Plots
  # --------------------------------------------------------------------------
  if (save_plots && nrow(results_df) > 0) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    if (verbose) cat("    Generating dotplot...\n")

    p_dot <- tryCatch({
      AETHER_plot_gsea_dotplot(gsea_result, top_n = plot_top_n)
    }, error = function(e) {
      if (verbose) cat("[ELEUTHIA] Warning: Could not create dotplot -", e$message, "\n")
      NULL
    })

    if (!is.null(p_dot)) {
      dot_height <- min(20, max(6, plot_top_n * 0.3))
      saved <- .save_ggplot(p_dot, plot_dir, "dotplot", plot_format,
                             width = plot_width, height = dot_height)
      files_created <- c(files_created, saved)
      if (verbose) cat("[ELEUTHIA] Dotplot saved\n")
    }

    # ------------------------------------------------------------------------
    # Per-pathway running-score plots (significant pathways, capped)
    # ------------------------------------------------------------------------
    if (do_enrichment_plots && nrow(gsea_result$significant) > 0) {
      sig_ordered   <- gsea_result$significant[order(gsea_result$significant$padj), ]
      plot_pathways <- utils::head(sig_ordered$pathway, max_enrichment_plots)

      if (verbose) {
        cap_note <- if (nrow(sig_ordered) > length(plot_pathways)) {
          paste0(" (capped from ", nrow(sig_ordered), ")")
        } else ""
        cat("    Generating", length(plot_pathways), "enrichment plot(s)", cap_note, "...\n")
      }

      for (pw in plot_pathways) {
        p_enr <- tryCatch({
          AETHER_plot_gsea_enrichment(gsea_result, pathway = pw)
        }, error = function(e) {
          if (verbose) cat("[ELEUTHIA] Warning: Enrichment plot failed for '",
                            pw, "' -", e$message, "\n")
          NULL
        })

        if (!is.null(p_enr)) {
          pw_clean <- gsub("[^A-Za-z0-9_]+", "_", pw)
          saved <- .save_ggplot(p_enr, plot_dir, paste0("enrichment_", pw_clean),
                                 plot_format, width = 8, height = 5)
          files_created <- c(files_created, saved)
        }
      }
      if (verbose) cat("[ELEUTHIA] Enrichment plots saved to:", plot_dir, "/\n")
    }
  }

  # --------------------------------------------------------------------------
  # Save RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, ".rds"))
    saveRDS(gsea_result, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("[ELEUTHIA] RDS:", rds_file, "\n")
  }

  if (verbose) {
    cat("\n[ELEUTHIA] --- Export Summary ---\n")
    cat("    Total files created:", length(files_created), "\n")
    cat("    Gene sets tested:", nrow(results_df), "\n")
    cat("    Significant (FDR <", gsea_result$params$fdr_threshold, "):",
        nrow(gsea_result$significant), "\n")
  }

  invisible(files_created)
}
