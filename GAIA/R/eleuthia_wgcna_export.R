# GAIA/Eleuthia/wgcna_export.R
# Export and reporting functions for WGCNA results
#
# Functions to export module gene lists, save results, and generate reports.
# Note: ELEUTHIA_export_enrichment() is in export_shared.R (generic function)


#' Export WGCNA results to files
#'
#' Exports various WGCNA analysis results to CSV files, generates plots,
#' and optionally saves R objects for later use.
#'
#' @param output_dir Output directory path.
#' @param modules Optional. wgcna_modules object to export.
#' @param trait_cor Optional. wgcna_trait_cor object to export.
#' @param gene_sig Optional. wgcna_gene_sig object to export.
#' @param hubs Optional. wgcna_hubs object to export.
#' @param enrichment Optional. gost_enrichment object from APOLLO_enrich_gost().
#'   If provided, exports enrichment results to <prefix>_enrichment/ subdirectory.
#' @param power_result Optional. wgcna_power object for power selection plot.
#' @param cluster_result Optional. wgcna_cluster object for sample dendrogram.
#' @param wgcna_data Optional. wgcna_data object (needed for sample dendrogram traits).
#' @param save_plots Logical. Generate and save plots. Default: TRUE.
#' @param plot_format Character. Plot format: "png", "pdf", or "both". Default: "png".
#' @param prefix Prefix for output filenames. Default: "wgcna".
#' @param save_rds Logical. Save R objects as RDS files. Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible list of file paths created.
#'
#' @details
#' Creates the following files (depending on objects provided):
#' - <prefix>_gene_modules.csv: Gene-to-module assignments
#' - <prefix>_module_summary.csv: Module sizes
#' - <prefix>_module_eigengenes.csv: Module eigengene values
#' - <prefix>_trait_correlations.csv: Module-trait correlations
#' - <prefix>_trait_pvalues.csv: Correlation p-values
#' - <prefix>_significant_associations.csv: Significant module-trait pairs
#' - <prefix>_gene_significance.csv: Full gene info table
#' - <prefix>_hub_genes.csv: Hub gene list
#' - <prefix>_hub_summary.csv: Hub counts per module
#' - Module-specific gene lists in <prefix>_modules/ subdirectory
#' - Enrichment results in <prefix>_enrichment/ subdirectory (if enrichment provided)
#' - Plots in <prefix>_plots/ subdirectory (if save_plots = TRUE)
#'
#' @examples
#' \dontrun{
#' ELEUTHIA_export_wgcna_results(
#'   output_dir = "results/wgcna",
#'   modules = modules,
#'   trait_cor = trait_cor,
#'   gene_sig = gene_sig,
#'   hubs = hubs,
#'   enrichment = gost_results,
#'   power_result = power_result,
#'   cluster_result = cluster_result
#' )
#'
#' }
#' @export
ELEUTHIA_export_wgcna_results <- function(output_dir,
                                           modules = NULL,
                                           trait_cor = NULL,
                                           gene_sig = NULL,
                                           hubs = NULL,
                                           enrichment = NULL,
                                           power_result = NULL,
                                           cluster_result = NULL,
                                           wgcna_data = NULL,
                                           save_plots = TRUE,
                                           plot_format = "png",
                                           prefix = "wgcna",
                                           save_rds = TRUE,
                                           verbose = TRUE) {

  if (verbose) cat("=== Exporting WGCNA Results ===\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created directory:", output_dir, "\n")
  }

  files_created <- character(0)

  # --------------------------------------------------------------------------
  # Export modules
  # --------------------------------------------------------------------------
  if (!is.null(modules)) {
    if (verbose) cat("\nExporting module results...\n")

    # Gene-module assignments
    gene_modules_file <- file.path(output_dir, paste0(prefix, "_gene_modules.csv"))
    write.csv(modules$gene_module_df, gene_modules_file, row.names = FALSE)
    files_created <- c(files_created, gene_modules_file)
    if (verbose) cat("  Gene-module assignments:", gene_modules_file, "\n")

    # Module summary
    summary_file <- file.path(output_dir, paste0(prefix, "_module_summary.csv"))
    write.csv(modules$module_summary, summary_file, row.names = FALSE)
    files_created <- c(files_created, summary_file)
    if (verbose) cat("  Module summary:", summary_file, "\n")

    # Module eigengenes
    me_file <- file.path(output_dir, paste0(prefix, "_module_eigengenes.csv"))
    me_df <- as.data.frame(modules$module_eigengenes)
    me_df$sample <- rownames(me_df)
    me_df <- me_df[, c("sample", setdiff(names(me_df), "sample"))]
    write.csv(me_df, me_file, row.names = FALSE)
    files_created <- c(files_created, me_file)
    if (verbose) cat("  Module eigengenes:", me_file, "\n")

    # Individual module gene lists
    modules_dir <- file.path(output_dir, paste0(prefix, "_modules"))
    if (!dir.exists(modules_dir)) {
      dir.create(modules_dir)
    }

    unique_modules <- unique(modules$module_names)
    for (mod in unique_modules) {
      mod_genes <- names(modules$module_names)[modules$module_names == mod]
      mod_df <- data.frame(gene = mod_genes, module = mod, stringsAsFactors = FALSE)
      mod_file <- file.path(modules_dir, paste0(mod, ".csv"))
      write.csv(mod_df, mod_file, row.names = FALSE)
      files_created <- c(files_created, mod_file)
    }
    if (verbose) cat("  Module gene lists:", modules_dir, "/\n")

    # Export color map
    if (!is.null(modules$module_colors)) {
      color_map_file <- file.path(output_dir, paste0(prefix, "_color_map.csv"))
      color_map_df <- data.frame(
        module = names(modules$module_colors),
        display_color = unname(modules$module_colors),
        color_name = if (!is.null(modules$color_names))
          modules$color_names[names(modules$module_colors)] else NA,
        stringsAsFactors = FALSE
      )
      write.csv(color_map_df, color_map_file, row.names = FALSE)
      files_created <- c(files_created, color_map_file)
      if (verbose) cat("  Color map:", color_map_file, "\n")
    }

    # Save RDS
    if (save_rds) {
      rds_file <- file.path(output_dir, paste0(prefix, "_modules.rds"))
      saveRDS(modules, rds_file)
      files_created <- c(files_created, rds_file)
      if (verbose) cat("  RDS:", rds_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Export trait correlations
  # --------------------------------------------------------------------------
  if (!is.null(trait_cor)) {
    if (verbose) cat("\nExporting trait correlation results...\n")

    # Correlation matrix
    cor_file <- file.path(output_dir, paste0(prefix, "_trait_correlations.csv"))
    cor_df <- as.data.frame(trait_cor$cor_matrix)
    cor_df$module <- rownames(cor_df)
    cor_df <- cor_df[, c("module", setdiff(names(cor_df), "module"))]
    write.csv(cor_df, cor_file, row.names = FALSE)
    files_created <- c(files_created, cor_file)
    if (verbose) cat("  Correlations:", cor_file, "\n")

    # P-value matrix
    pval_file <- file.path(output_dir, paste0(prefix, "_trait_pvalues.csv"))
    pval_df <- as.data.frame(trait_cor$pvalue_matrix)
    pval_df$module <- rownames(pval_df)
    pval_df <- pval_df[, c("module", setdiff(names(pval_df), "module"))]
    write.csv(pval_df, pval_file, row.names = FALSE)
    files_created <- c(files_created, pval_file)
    if (verbose) cat("  P-values:", pval_file, "\n")

    # Significant associations
    if (nrow(trait_cor$significant_associations) > 0) {
      sig_file <- file.path(output_dir, paste0(prefix, "_significant_associations.csv"))
      write.csv(trait_cor$significant_associations, sig_file, row.names = FALSE)
      files_created <- c(files_created, sig_file)
      if (verbose) cat("  Significant associations:", sig_file, "\n")
    }

    if (save_rds) {
      rds_file <- file.path(output_dir, paste0(prefix, "_trait_cor.rds"))
      saveRDS(trait_cor, rds_file)
      files_created <- c(files_created, rds_file)
      if (verbose) cat("  RDS:", rds_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Export gene significance
  # --------------------------------------------------------------------------
  if (!is.null(gene_sig)) {
    if (verbose) cat("\nExporting gene significance results...\n")

    gs_file <- file.path(output_dir, paste0(prefix, "_gene_significance.csv"))
    write.csv(gene_sig$gene_info, gs_file, row.names = FALSE)
    files_created <- c(files_created, gs_file)
    if (verbose) cat("  Gene significance:", gs_file, "\n")

    if (save_rds) {
      rds_file <- file.path(output_dir, paste0(prefix, "_gene_sig.rds"))
      saveRDS(gene_sig, rds_file)
      files_created <- c(files_created, rds_file)
      if (verbose) cat("  RDS:", rds_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Export hub genes
  # --------------------------------------------------------------------------
  if (!is.null(hubs)) {
    if (verbose) cat("\nExporting hub gene results...\n")

    hub_file <- file.path(output_dir, paste0(prefix, "_hub_genes.csv"))
    write.csv(hubs$hub_genes, hub_file, row.names = FALSE)
    files_created <- c(files_created, hub_file)
    if (verbose) cat("  Hub genes:", hub_file, "\n")

    hub_summary_file <- file.path(output_dir, paste0(prefix, "_hub_summary.csv"))
    write.csv(hubs$hub_summary, hub_summary_file, row.names = FALSE)
    files_created <- c(files_created, hub_summary_file)
    if (verbose) cat("  Hub summary:", hub_summary_file, "\n")

    if (save_rds) {
      rds_file <- file.path(output_dir, paste0(prefix, "_hubs.rds"))
      saveRDS(hubs, rds_file)
      files_created <- c(files_created, rds_file)
      if (verbose) cat("  RDS:", rds_file, "\n")
    }
  }

  # --------------------------------------------------------------------------
  # Generate plots
  # --------------------------------------------------------------------------
  if (save_plots) {
    has_plots <- !is.null(power_result) || !is.null(trait_cor) || !is.null(cluster_result)

    if (has_plots) {
      plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
      }

      if (verbose) cat("\nGenerating plots...\n")

      # Helper function to save plots
      .save_plot <- function(p, name, width = 10, height = 8) {
        if (is.null(p)) return(NULL)
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

      # Power selection plot
      if (!is.null(power_result)) {
        p <- tryCatch({
          AETHER_plot_wgcna_power(power_result, return_plots = FALSE)
        }, error = function(e) {
          if (verbose) cat("  Warning: Could not create power plot -", e$message, "\n")
          NULL
        })
        if (!is.null(p)) {
          # AETHER_plot_wgcna_power returns NULL invisibly when return_plots = FALSE
          # but it creates the plot. We need return_plots = TRUE
          p <- tryCatch({
            AETHER_plot_wgcna_power(power_result, return_plots = TRUE)
          }, error = function(e) NULL)

          if (!is.null(p) && inherits(p, "list")) {
            # Combine the two plots
            combined <- tryCatch({
              gridExtra::grid.arrange(p$sft_plot, p$connectivity_plot, ncol = 2)
            }, error = function(e) NULL)

            if (!is.null(combined)) {
              if (plot_format %in% c("png", "both")) {
                png_file <- file.path(plot_dir, "power_selection.png")
                png(png_file, width = 12, height = 5, units = "in", res = 150)
                gridExtra::grid.arrange(p$sft_plot, p$connectivity_plot, ncol = 2)
                dev.off()
                files_created <- c(files_created, png_file)
              }
              if (plot_format %in% c("pdf", "both")) {
                pdf_file <- file.path(plot_dir, "power_selection.pdf")
                pdf(pdf_file, width = 12, height = 5)
                gridExtra::grid.arrange(p$sft_plot, p$connectivity_plot, ncol = 2)
                dev.off()
                files_created <- c(files_created, pdf_file)
              }
              if (verbose) cat("  Power selection plot\n")
            }
          }
        }
      }

      # Module-trait heatmap
      if (!is.null(trait_cor)) {
        p <- tryCatch({
          AETHER_plot_module_trait_heatmap_gg(trait_cor)
        }, error = function(e) {
          if (verbose) cat("  Warning: Could not create trait heatmap -", e$message, "\n")
          NULL
        })
        if (!is.null(p)) {
          # Calculate size based on number of modules and traits
          n_modules <- nrow(trait_cor$cor_matrix)
          n_traits <- ncol(trait_cor$cor_matrix)
          width <- max(6, n_traits * 1.2 + 3)
          height <- max(6, n_modules * 0.4 + 2)
          plot_files <- .save_plot(p, "module_trait_heatmap", width = width, height = height)
          files_created <- c(files_created, plot_files)
          if (verbose) cat("  Module-trait heatmap\n")
        }
      }
            if (!is.null(trait_cor)) {
        p <- tryCatch({
          AETHER_plot_module_trait_heatmap_gg(trait_cor,use_padj = TRUE)
        }, error = function(e) {
          if (verbose) cat("  Warning: Could not create trait heatmap -", e$message, "\n")
          NULL
        })
        if (!is.null(p)) {
          # Calculate size based on number of modules and traits
          n_modules <- nrow(trait_cor$cor_matrix)
          n_traits <- ncol(trait_cor$cor_matrix)
          width <- max(6, n_traits * 1.2 + 3)
          height <- max(6, n_modules * 0.4 + 2)
          plot_files <- .save_plot(p, "module_trait_heatmap_padj", width = width, height = height)
          files_created <- c(files_created, plot_files)
          if (verbose) cat("  Module-trait heatmap padjusted\n")
        }
      }

      # Sample dendrogram
      if (!is.null(cluster_result)) {
        # Get traits for color bar if available
        traits <- NULL
        if (!is.null(wgcna_data) && !is.null(wgcna_data$datTraits)) {
          traits <- wgcna_data$datTraits
        }

        # Sample dendrogram uses base R plotting
        if (plot_format %in% c("png", "both")) {
          png_file <- file.path(plot_dir, "sample_dendrogram.png")
          png(png_file, width = 12, height = 6, units = "in", res = 150)
          tryCatch({
            AETHER_plot_sample_dendrogram(cluster_result, traits = traits)
          }, error = function(e) {
            if (verbose) cat("  Warning: Could not create dendrogram -", e$message, "\n")
          })
          dev.off()
          files_created <- c(files_created, png_file)
        }
        if (plot_format %in% c("pdf", "both")) {
          pdf_file <- file.path(plot_dir, "sample_dendrogram.pdf")
          pdf(pdf_file, width = 12, height = 6)
          tryCatch({
            AETHER_plot_sample_dendrogram(cluster_result, traits = traits)
          }, error = function(e) NULL)
          dev.off()
          files_created <- c(files_created, pdf_file)
        }
        if (verbose) cat("  Sample dendrogram\n")
      }

      if (verbose) cat("  Plots saved to:", plot_dir, "/\n")
    }
  }

  # --------------------------------------------------------------------------
  # Export enrichment results
  # --------------------------------------------------------------------------
  if (!is.null(enrichment)) {
    if (verbose) cat("\nExporting enrichment results...\n")

    enrich_dir <- file.path(output_dir, paste0(prefix, "_enrichment"))
    enrich_files <- ELEUTHIA_export_enrichment(
      enrichment = enrichment,
      output_dir = enrich_dir,
      prefix = "enrichment",
      by_source = FALSE,
      by_module = TRUE,
      save_rds = save_rds,
      verbose = FALSE
    )
    files_created <- c(files_created, enrich_files)

    if (verbose) {
      cat("  Enrichment directory:", enrich_dir, "\n")
      cat("  Files created:", length(enrich_files), "\n")
    }
  }

  if (verbose) {
    cat("\n--- Export Summary ---\n")
    cat("Total files created:", length(files_created), "\n")
    cat("Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}


#' Generate WGCNA analysis report
#'
#' Creates a formatted text report summarizing the WGCNA analysis results.
#'
#' @param output_dir Output directory for report.
#' @param wgcna_data Optional. wgcna_data object (data summary).
#' @param power_result Optional. wgcna_power object (power selection).
#' @param modules Optional. wgcna_modules object (module detection).
#' @param trait_cor Optional. wgcna_trait_cor object (trait correlations).
#' @param gene_sig Optional. wgcna_gene_sig object (gene significance).
#' @param hubs Optional. wgcna_hubs object (hub genes).
#' @param prefix Filename prefix. Default: "wgcna".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible path to report file.
#'
#' @export
ELEUTHIA_wgcna_report <- function(output_dir,
                                   wgcna_data = NULL,
                                   power_result = NULL,
                                   modules = NULL,
                                   trait_cor = NULL,
                                   gene_sig = NULL,
                                   hubs = NULL,
                                   prefix = "wgcna",
                                   verbose = TRUE) {

  if (verbose) cat("=== Generating WGCNA Report ===\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  report_file <- file.path(output_dir, paste0(prefix, "_report.txt"))

  # Helper function
  add_line <- function(...) {
    cat(..., "\n", file = report_file, append = TRUE, sep = "")
  }

  # Initialize file
  cat("", file = report_file)

  add_line("================================================================================")
  add_line("WGCNA ANALYSIS REPORT")
  add_line("================================================================================")
  add_line("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add_line("")

  # --------------------------------------------------------------------------
  # Data summary
  # --------------------------------------------------------------------------
  if (!is.null(wgcna_data)) {
    add_line("================================================================================")
    add_line("DATA SUMMARY")
    add_line("================================================================================")
    add_line("Samples: ", wgcna_data$n_samples)
    add_line("Genes/Features: ", wgcna_data$n_genes)
    add_line("Traits: ", ncol(wgcna_data$datTraits))
    add_line("  ", paste(colnames(wgcna_data$datTraits), collapse = ", "))

    if (length(wgcna_data$removed_genes) > 0) {
      add_line("Genes removed (zero variance): ", length(wgcna_data$removed_genes))
    }
    if (length(wgcna_data$removed_samples) > 0) {
      add_line("Samples removed: ", length(wgcna_data$removed_samples))
    }
    add_line("")
  }

  # --------------------------------------------------------------------------
  # Power selection
  # --------------------------------------------------------------------------
  if (!is.null(power_result)) {
    add_line("================================================================================")
    add_line("POWER SELECTION")
    add_line("================================================================================")
    add_line("Selected power: ", power_result$power)
    add_line("Selection method: ", power_result$selection_method)
    add_line("R-squared at power: ", round(power_result$r2_at_power, 3))
    add_line("Mean connectivity: ", round(power_result$mean_k_at_power, 1))
    add_line("R-squared cutoff used: ", power_result$r2_cutoff)
    add_line("Network type: ", power_result$network_type)
    add_line("")
  }

  # --------------------------------------------------------------------------
  # Module detection
  # --------------------------------------------------------------------------
  if (!is.null(modules)) {
    add_line("================================================================================")
    add_line("MODULE DETECTION")
    add_line("================================================================================")
    add_line("Power used: ", modules$power)
    add_line("Network type: ", modules$network_type)
    add_line("Modules detected: ", modules$n_modules, " (excluding module_0/unassigned)")
    add_line("Total genes: ", modules$n_genes)
    add_line("")
    add_line("Module sizes:")
    for (i in seq_len(min(nrow(modules$module_summary), 15))) {
      row <- modules$module_summary[i, ]
      add_line("  ", row$module, ": ", row$n_genes, " genes")
    }
    if (nrow(modules$module_summary) > 15) {
      add_line("  ... (", nrow(modules$module_summary) - 15, " more modules)")
    }
    add_line("")
  }

  # --------------------------------------------------------------------------
  # Trait correlations
  # --------------------------------------------------------------------------
  if (!is.null(trait_cor)) {
    add_line("================================================================================")
    add_line("MODULE-TRAIT CORRELATIONS")
    add_line("================================================================================")
    add_line("Correlation method: ", trait_cor$cor_method)
    add_line("P-value adjustment: ", trait_cor$p_adjust)
    add_line("Significant associations: ", trait_cor$n_significant)
    add_line("")

    if (trait_cor$n_significant > 0) {
      add_line("Top significant associations:")
      top_assoc <- head(trait_cor$significant_associations, 10)
      for (i in seq_len(nrow(top_assoc))) {
        row <- top_assoc[i, ]
        add_line("  ", row$module, " ~ ", row$trait,
                 ": r = ", round(row$correlation, 3),
                 ", p = ", formatC(row$padj, format = "e", digits = 2))
      }
    }
    add_line("")
  }

  # --------------------------------------------------------------------------
  # Hub genes
  # --------------------------------------------------------------------------
  if (!is.null(hubs)) {
    add_line("================================================================================")
    add_line("HUB GENES")
    add_line("================================================================================")
    add_line("MM threshold: ", hubs$criteria$mm_threshold)
    if (!is.na(hubs$criteria$gs_threshold)) {
      add_line("GS threshold: ", hubs$criteria$gs_threshold,
               " (trait: ", hubs$criteria$trait_name, ")")
    }
    add_line("Total hub genes: ", nrow(hubs$hub_genes))
    add_line("")

    add_line("Hub genes per module:")
    for (i in seq_len(min(nrow(hubs$hub_summary), 10))) {
      row <- hubs$hub_summary[i, ]
      add_line("  ", row$module, ": ", row$n_hubs, " hubs (top: ", row$top_hub, ")")
    }
    add_line("")
  }

  # --------------------------------------------------------------------------
  # Footer
  # --------------------------------------------------------------------------
  add_line("================================================================================")
  add_line("END OF REPORT")
  add_line("================================================================================")

  if (verbose) {
    cat("Report saved to:", report_file, "\n")
  }

  invisible(report_file)
}
