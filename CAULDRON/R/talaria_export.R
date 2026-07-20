# ==============================================================================
# TALARIA - Export: First Steps Analysis
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Bundles the outputs of a completed first-steps analysis into a structured
# directory: the SCE object, analysis parameter summary, sweep result tables,
# QC plots, diagnostic sweep plots, embedding grid comparisons, and final
# embedding plots.
#
# All functions prefixed: TALARIA_
# ==============================================================================


#' Export a first-steps analysis
#'
#' Saves all outputs of a CAULDRON first-steps pipeline run into a structured
#' directory tree.  The function auto-detects what has been computed from the
#' SCE object and only exports what exists, so it works whether the user ran
#' the full pipeline or just loaded data and produced a single UMAP.  Sweep
#' objects (if any were run) are passed explicitly.
#'
#' \strong{Output directory structure:}
#' \preformatted{
#' output_dir/
#' ├── sce/
#' │   └── sce.rds
#' ├── plots/
#' │   ├── qc_metrics.png
#' │   ├── qc_metrics_by_sample.png  (if sample_col provided)
#' │   ├── pca_elbow.png
#' │   ├── sweep_hvg.png
#' │   ├── sweep_pc.png
#' │   ├── sweep_k.png
#' │   ├── sweep_resolution.png
#' │   ├── sweep_umap.png
#' │   ├── sweep_tsne.png
#' │   ├── grid_umap.png
#' │   ├── grid_tsne.png
#' │   ├── umap.png
#' │   └── tsne.png
#' └── tables/
#'     ├── analysis_parameters.csv
#'     ├── cell_metadata.csv
#'     ├── sweep_hvg_results.csv
#'     ├── sweep_pc_results.csv
#'     ├── sweep_k_results.csv
#'     ├── sweep_resolution_results.csv
#'     ├── sweep_umap_results.csv
#'     └── sweep_tsne_results.csv
#' }
#'
#' \strong{BPCells backend:} \code{saveRDS(sce, ...)} is used for both
#' in-memory and BPCells-backed objects — the BPCells matrix is already on
#' disk and the RDS stores only the SCE shell with a path reference.  The
#' BPCells directory (stored in \code{metadata(sce)$bpcells_dir}) must remain
#' accessible alongside the RDS for the object to be reloadable.
#'
#' \strong{Analysis parameters table:} a two-column CSV (\code{parameter},
#' \code{value}) summarising the key choices made during the analysis —
#' \code{n_hvgs}, \code{n_pcs_used}, \code{k}, \code{clustering_method},
#' \code{resolution}, \code{n_neighbors}, \code{min_dist}, and
#' \code{perplexity}.  Values are populated from sweep objects where provided;
#' parameters whose sweep was not run are recorded as \code{NA}.
#'
#' @param sce A \code{SingleCellExperiment} object — the main analysis object.
#'   The function inspects its state (assays, reducedDims, colData, metadata)
#'   to determine what has been computed and should be exported.
#' @param output_dir Character.  Path to the root output directory.  Created
#'   (recursively) if it does not exist.
#' @param hvg_sweep A \code{pyri_hvg_sweep} object from
#'   \code{\link{PYRI_tune_hvg}}, or \code{NULL} (default).
#' @param pc_sweep A \code{talos_pc_sweep} object from
#'   \code{\link{TALOS_tune_pc}}, or \code{NULL} (default).
#' @param k_sweep A \code{talos_k_sweep} object from
#'   \code{\link{TALOS_tune_k}}, or \code{NULL} (default).
#' @param resolution_sweep A \code{talos_resolution_sweep} object from
#'   \code{\link{TALOS_tune_resolution}}, or \code{NULL} (default).
#' @param umap_sweep A \code{talos_embedding_sweep} object (type \code{"umap"})
#'   from \code{\link{TALOS_tune_umap}}, or \code{NULL} (default).
#' @param tsne_sweep A \code{talos_embedding_sweep} object (type \code{"tsne"})
#'   from \code{\link{TALOS_tune_tsne}}, or \code{NULL} (default).
#' @param colour_by Character.  \code{colData} column used to colour embedding
#'   plots.  Default \code{"cluster"}.  If not found, falls back to the first
#'   factor or character column in \code{colData} with a console note.
#' @param sample_col Character or \code{NULL}.  If provided, an additional QC
#'   plot grouped by this \code{colData} column is saved alongside the
#'   ungrouped version.  Default \code{NULL}.
#' @param assay_name Character.  Assay used when \code{colour_by} is a gene
#'   name in the embedding plots.  Default \code{"logcounts"}.
#' @param format Character.  Output format for plots: \code{"png"} (default,
#'   600 DPI raster) or \code{"pdf"} (vector, no DPI).
#' @param verbose Logical.  Print a step-by-step progress log and a final
#'   file manifest summary.  Default \code{TRUE}.
#'
#' @return Invisibly returns a list with \code{output_dir} (character) and
#'   \code{files} (character vector of all paths successfully written).
#' @export
TALARIA_export_first_steps <- function(sce,
                                        output_dir,
                                        hvg_sweep        = NULL,
                                        pc_sweep         = NULL,
                                        k_sweep          = NULL,
                                        resolution_sweep = NULL,
                                        umap_sweep       = NULL,
                                        tsne_sweep       = NULL,
                                        colour_by        = "cluster",
                                        sample_col       = NULL,
                                        assay_name       = "logcounts",
                                        format           = c("png", "pdf"),
                                        verbose          = TRUE) {

  format <- match.arg(format)

  if (!inherits(sce, "SingleCellExperiment"))
    stop("'sce' must be a SingleCellExperiment object.", call. = FALSE)

  if (!is.null(sample_col) && !sample_col %in% names(colData(sce)))
    stop("sample_col '", sample_col, "' not found in colData.",
         call. = FALSE)

  # ── Create output directories ─────────────────────────────────────────────
  dir_sce    <- file.path(output_dir, "sce")
  dir_plots  <- file.path(output_dir, "plots")
  dir_tables <- file.path(output_dir, "tables")
  for (d in c(dir_sce, dir_plots, dir_tables))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # ── Auto-detect computed state ────────────────────────────────────────────
  is_bpcells <- any(vapply(assayNames(sce), function(a)
    inherits(assay(sce, a), "IterableMatrix"), logical(1L)))
  has_qc     <- all(c("sum", "detected") %in% names(colData(sce)))
  has_pca    <- "PCA"  %in% reducedDimNames(sce)
  has_umap   <- "UMAP" %in% reducedDimNames(sce)
  has_tsne   <- "tSNE" %in% reducedDimNames(sce)
  has_hvgs   <- !is.null(metadata(sce)$hvg)

  # Resolve colour_by — fall back gracefully if missing
  embed_col <- if (colour_by %in% names(colData(sce))) {
    colour_by
  } else {
    cd    <- colData(sce)
    fallb <- names(cd)[vapply(names(cd), function(x)
      is.factor(cd[[x]]) || is.character(cd[[x]]), logical(1L))]
    if (length(fallb)) {
      if (verbose)
        cat("  Note: '", colour_by, "' not in colData",
            " \u2014 colouring by '", fallb[1L], "'.\n", sep = "")
      fallb[1L]
    } else {
      if (verbose)
        cat("  Note: no discrete colData column found \u2014 embedding plots will be skipped.\n")
      NULL
    }
  }

  saved <- character(0L)

  # ── Internal save helper (ggplot only) ────────────────────────────────────
  .save_gg <- function(p, name, width, height) {
    path <- file.path(dir_plots, paste0(name, ".", format))
    tryCatch({
      if (format == "png") {
        ggsave(path, plot = p, width = width, height = height,
               dpi = 600, units = "in")
      } else {
        ggsave(path, plot = p, width = width, height = height,
               device = "pdf", units = "in")
      }
      saved <<- c(saved, path)
      if (verbose) cat("    ", basename(path), "\n", sep = "")
    }, error = function(e)
      cat("  Warning: could not save '", name, "': ",
          conditionMessage(e), "\n", sep = ""))
  }

  if (verbose)
    cat(sprintf(
      "\u2500\u2500 TALARIA: export first steps %s\n  Output     : %s\n  Backend    : %s\n  Cells      : %s  |  Genes: %s\n%s\n",
      strrep("\u2500", 24),
      output_dir,
      if (is_bpcells) "BPCells (on-disk)" else "in-memory",
      format(ncol(sce), big.mark = ","),
      format(nrow(sce), big.mark = ","),
      strrep("\u2500", 56)))

  # ── [1/5] SCE object ──────────────────────────────────────────────────────
  if (verbose) cat("  [1/5] Saving SCE object ...\n")
  tryCatch({
    sce_path <- file.path(dir_sce, "sce.rds")
    saveRDS(sce, sce_path)
    saved <- c(saved, sce_path)
    if (verbose) cat("    sce.rds\n")
    if (is_bpcells && verbose)
      cat("    Note: BPCells matrix lives at metadata(sce)$bpcells_dir \u2014 keep that directory alongside sce.rds.\n")
  }, error = function(e)
    cat("  Warning: could not save SCE: ", conditionMessage(e), "\n", sep = ""))

  # ── [2/5] Tables ──────────────────────────────────────────────────────────
  if (verbose) cat("  [2/5] Saving tables ...\n")

  # Cell metadata
  tryCatch({
    p <- file.path(dir_tables, "cell_metadata.csv")
    write.csv(as.data.frame(colData(sce)), p, row.names = TRUE)
    saved <- c(saved, p)
    if (verbose) cat("    cell_metadata.csv\n")
  }, error = function(e)
    cat("  Warning: cell_metadata: ", conditionMessage(e), "\n", sep = ""))

  # Analysis parameter summary
  tryCatch({
    params_df <- data.frame(
      parameter = c("n_hvgs", "n_pcs_used", "k",
                    "clustering_method", "resolution",
                    "n_neighbors", "min_dist", "perplexity"),
      value = c(
        if (has_hvgs)
          as.character(length(metadata(sce)$hvg))                      else NA_character_,
        if (!is.null(pc_sweep))
          as.character(pc_sweep$best_n_pcs)                            else NA_character_,
        if (!is.null(k_sweep))
          as.character(k_sweep$best_k)                                 else NA_character_,
        if (!is.null(resolution_sweep))
          as.character(resolution_sweep$best_params$method)            else NA_character_,
        if (!is.null(resolution_sweep))
          as.character(resolution_sweep$best_params$resolution)        else NA_character_,
        if (!is.null(umap_sweep))
          as.character(umap_sweep$best_params$n_neighbors)             else NA_character_,
        if (!is.null(umap_sweep))
          as.character(umap_sweep$best_params$min_dist)                else NA_character_,
        if (!is.null(tsne_sweep))
          as.character(tsne_sweep$best_params$perplexity)              else NA_character_
      ),
      stringsAsFactors = FALSE
    )
    p <- file.path(dir_tables, "analysis_parameters.csv")
    write.csv(params_df, p, row.names = FALSE)
    saved <- c(saved, p)
    if (verbose) cat("    analysis_parameters.csv\n")
  }, error = function(e)
    cat("  Warning: analysis_parameters: ", conditionMessage(e), "\n", sep = ""))

  # Sweep result tables
  sweep_tables <- list(
    sweep_hvg_results        = if (!is.null(hvg_sweep))        hvg_sweep$results        else NULL,
    sweep_pc_results         = if (!is.null(pc_sweep))         pc_sweep$results         else NULL,
    sweep_k_results          = if (!is.null(k_sweep))          k_sweep$results          else NULL,
    sweep_resolution_results = if (!is.null(resolution_sweep)) resolution_sweep$results else NULL,
    sweep_umap_results       = if (!is.null(umap_sweep))       umap_sweep$results       else NULL,
    sweep_tsne_results       = if (!is.null(tsne_sweep))       tsne_sweep$results       else NULL
  )
  for (nm in names(sweep_tables)) {
    if (is.null(sweep_tables[[nm]])) next
    tryCatch({
      p <- file.path(dir_tables, paste0(nm, ".csv"))
      write.csv(sweep_tables[[nm]], p, row.names = FALSE)
      saved <- c(saved, p)
      if (verbose) cat("    ", basename(p), "\n", sep = "")
    }, error = function(e)
      cat("  Warning: '", nm, "': ", conditionMessage(e), "\n", sep = ""))
  }

  # ── [3/5] QC plots ────────────────────────────────────────────────────────
  if (has_qc) {
    if (verbose) cat("  [3/5] Saving QC plots ...\n")
    qc_metrics <- intersect(
      c("sum", "detected", "subsets_mt_percent", "subsets_ribo_percent"),
      names(colData(sce)))
    n_rows <- ceiling(length(qc_metrics) / 3L)

    .save_gg(ASPIS_plot_qc(sce, metrics = qc_metrics),
             "qc_metrics", width = 10, height = 3.5 * n_rows)

    if (!is.null(sample_col)) {
      n_samples  <- nlevels(factor(colData(sce)[[sample_col]]))
      plot_width <- max(10, 1.5 * n_samples)
      .save_gg(ASPIS_plot_qc(sce, metrics = qc_metrics, group_by = sample_col),
               "qc_metrics_by_sample", width = plot_width, height = 3.5 * n_rows)
    }
  } else {
    if (verbose)
      cat("  [3/5] QC metrics not found in colData \u2014 skipping QC plots.\n")
  }

  # ── [4/5] Diagnostic & sweep plots ───────────────────────────────────────
  if (verbose) cat("  [4/5] Saving diagnostic and sweep plots ...\n")

  if (has_pca)
    .save_gg(ASPIS_plot_elbow(sce, suggest = "both"),
             "pca_elbow", width = 8, height = 5)

  if (!is.null(hvg_sweep))
    .save_gg(hvg_sweep$plot,        "sweep_hvg",        width = 6,  height = 8)
  if (!is.null(pc_sweep))
    .save_gg(pc_sweep$plot,         "sweep_pc",         width = 6,  height = 8)
  if (!is.null(k_sweep))
    .save_gg(k_sweep$plot,          "sweep_k",          width = 6,  height = 8)
  if (!is.null(resolution_sweep))
    .save_gg(resolution_sweep$plot, "sweep_resolution", width = 7,  height = 11)
  if (!is.null(umap_sweep))
    .save_gg(umap_sweep$plot,       "sweep_umap",       width = 8,  height = 6)
  if (!is.null(tsne_sweep))
    .save_gg(tsne_sweep$plot,       "sweep_tsne",       width = 6,  height = 8)

  # Embedding grids (require a resolvable colour column)
  if (!is.null(umap_sweep) && !is.null(embed_col))
    .save_gg(ASPIS_plot_embedding_grid(umap_sweep, sce,
                                        colour_by  = embed_col,
                                        assay_name = assay_name),
             "grid_umap", width = 14, height = 12)

  if (!is.null(tsne_sweep) && !is.null(embed_col))
    .save_gg(ASPIS_plot_embedding_grid(tsne_sweep, sce,
                                        colour_by  = embed_col,
                                        assay_name = assay_name),
             "grid_tsne", width = 14, height = 12)

  # ── [5/5] Final embedding plots ───────────────────────────────────────────
  if (verbose) cat("  [5/5] Saving final embedding plots ...\n")

  if (has_umap && !is.null(embed_col))
    .save_gg(ASPIS_plot_umap(sce,
                              colour_by      = embed_col,
                              label_clusters = TRUE,
                              assay_name     = assay_name),
             "umap", width = 7, height = 6)

  if (has_tsne && !is.null(embed_col))
    .save_gg(ASPIS_plot_tsne(sce,
                              colour_by      = embed_col,
                              label_clusters = TRUE,
                              assay_name     = assay_name),
             "tsne", width = 7, height = 6)

  # ── Summary ───────────────────────────────────────────────────────────────
  if (verbose) {
    n_plots  <- sum(endsWith(saved, paste0(".", format)))
    n_tables <- sum(endsWith(saved, ".csv"))
    n_rds    <- sum(endsWith(saved, ".rds"))
    cat(sprintf(
      "%s\n  Done  : %d plot(s)  |  %d table(s)  |  %d RDS\n  Saved to: %s\n%s\n",
      strrep("\u2500", 56),
      n_plots, n_tables, n_rds,
      output_dir,
      strrep("\u2500", 56)))
  }

  invisible(list(output_dir = output_dir, files = saved))
}


# ==============================================================================
# KERAUNOS result exports
# ==============================================================================
# Internal plot helpers and export functions for DE, GSEA, and ORA results
# produced by the KERAUNOS module.  Plot code is native to CAULDRON — no
# dependency on GAIA/AETHER; inspired by those implementations but standalone.
# ==============================================================================


# ── Internal plot helpers ──────────────────────────────────────────────────────

# Volcano plot for a single pseudobulk cluster result data.frame.
# Columns expected: gene, log2FoldChange, pvalue, <sig_col>
.talaria_plot_volcano <- function(df, sig_col, sig_thresh,
                                   l2fc_thresh = 1, top_n_label = 10,
                                   title = "Volcano Plot") {

  df <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), , drop = FALSE]

  sig     <- !is.na(df[[sig_col]]) & df[[sig_col]] < sig_thresh
  df$.cat <- "ns"
  df$.cat[sig & df$log2FoldChange >  l2fc_thresh] <- "up"
  df$.cat[sig & df$log2FoldChange < -l2fc_thresh] <- "down"
  df$.cat[sig & abs(df$log2FoldChange) <= l2fc_thresh] <- "low"

  n_up   <- sum(df$.cat == "up")
  n_down <- sum(df$.cat == "down")
  n_low  <- sum(df$.cat == "low")
  n_ns   <- sum(df$.cat == "ns")

  lvl_up   <- paste0("Up (n=",      n_up,   ")")
  lvl_down <- paste0("Down (n=",    n_down, ")")
  lvl_low  <- paste0("Low-reg (n=", n_low,  ")")
  lvl_ns   <- paste0("NS (n=",      n_ns,   ")")

  cat_map         <- c(up = lvl_up, down = lvl_down, low = lvl_low, ns = lvl_ns)
  df$Significance <- factor(cat_map[df$.cat],
                             levels = c(lvl_up, lvl_down, lvl_low, lvl_ns))
  df$.cat         <- NULL
  colors          <- c("#B31B21", "#1465AC", "#4DAF4A", "darkgray")
  names(colors)   <- c(lvl_up, lvl_down, lvl_low, lvl_ns)

  sig_rows   <- df[!is.na(df[[sig_col]]) & df[[sig_col]] < sig_thresh, ]
  sig_line_y <- if (nrow(sig_rows) > 0)
    -log10(max(sig_rows$pvalue, na.rm = TRUE)) else NA_real_

  df <- df[order(df$Significance, decreasing = TRUE), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange,
                                         y = -log10(pvalue),
                                         colour = Significance)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-l2fc_thresh, l2fc_thresh),
                        linetype = "dashed", colour = "black",
                        linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = colors, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 4, alpha = 1))) +
    ggplot2::xlab(expression("Log"[2] * "Fold Change")) +
    ggplot2::ylab(expression("-log"[10] * "(p-value)")) +
    ggplot2::ggtitle(title) +
    KHALKOS_theme_cauldron() +
    ggplot2::theme(legend.position  = "bottom",
                   legend.direction = "vertical")

  if (!is.na(sig_line_y))
    p <- p + ggplot2::geom_hline(yintercept = sig_line_y,
                                  linetype = "dashed", colour = "black",
                                  linewidth = 0.4)

  if (top_n_label > 0 && nrow(sig_rows) > 0) {
    lab_df <- head(sig_rows[order(sig_rows[[sig_col]]), ], top_n_label)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data    = lab_df,
        mapping = ggplot2::aes(label = gene),
        size    = 2.5, colour = "black", show.legend = FALSE,
        box.padding = 0.3, max.overlaps = 20
      )
    }
  }
  p
}


# MA plot for a single pseudobulk cluster result data.frame.
# Columns expected: gene, baseMean, log2FoldChange, pvalue, <sig_col>
.talaria_plot_ma <- function(df, sig_col, sig_thresh,
                               l2fc_thresh = 1, top_n_label = 10,
                               title = "MA Plot") {

  df <- df[!is.na(df$baseMean) & !is.na(df$log2FoldChange), , drop = FALSE]
  df$baseMean_log2 <- log2(df$baseMean + 1)

  sig     <- !is.na(df[[sig_col]]) & df[[sig_col]] < sig_thresh
  df$.cat <- "ns"
  df$.cat[sig & df$log2FoldChange >  l2fc_thresh] <- "up"
  df$.cat[sig & df$log2FoldChange < -l2fc_thresh] <- "down"
  df$.cat[sig & abs(df$log2FoldChange) <= l2fc_thresh] <- "low"

  n_up   <- sum(df$.cat == "up")
  n_down <- sum(df$.cat == "down")
  n_low  <- sum(df$.cat == "low")
  n_ns   <- sum(df$.cat == "ns")

  lvl_up   <- paste0("Up (n=",      n_up,   ")")
  lvl_down <- paste0("Down (n=",    n_down, ")")
  lvl_low  <- paste0("Low-reg (n=", n_low,  ")")
  lvl_ns   <- paste0("NS (n=",      n_ns,   ")")

  cat_map         <- c(up = lvl_up, down = lvl_down, low = lvl_low, ns = lvl_ns)
  df$Significance <- factor(cat_map[df$.cat],
                             levels = c(lvl_up, lvl_down, lvl_low, lvl_ns))
  df$.cat         <- NULL
  colors          <- c("#B31B21", "#1465AC", "#4DAF4A", "darkgray")
  names(colors)   <- c(lvl_up, lvl_down, lvl_low, lvl_ns)

  df       <- df[order(df$Significance, decreasing = TRUE), ]
  sig_rows <- df[!is.na(df[[sig_col]]) & df[[sig_col]] < sig_thresh, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = baseMean_log2,
                                         y = log2FoldChange,
                                         colour = Significance)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "firebrick", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-l2fc_thresh, l2fc_thresh),
                        linetype = "dashed", colour = "black",
                        linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = colors, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 4, alpha = 1))) +
    ggplot2::xlab(expression("log"[2] * "(baseMean + 1)")) +
    ggplot2::ylab(expression("log"[2] * "Fold Change")) +
    ggplot2::ggtitle(title) +
    KHALKOS_theme_cauldron() +
    ggplot2::theme(legend.position  = "bottom",
                   legend.direction = "vertical")

  if (top_n_label > 0 && nrow(sig_rows) > 0) {
    lab_df <- head(sig_rows[order(sig_rows[[sig_col]]), ], top_n_label)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data    = lab_df,
        mapping = ggplot2::aes(label = gene),
        size    = 2.5, colour = "black", show.legend = FALSE,
        box.padding = 0.3, max.overlaps = 20
      )
    }
  }
  p
}


# GSEA dotplot.  df: fgsea results (pathway, NES, padj, size).
# Top n_each enriched (NES > 0) and n_each depleted (NES < 0) by padj.
.talaria_gsea_dotplot <- function(df, top_n = 20, title = "GSEA") {

  df <- df[!is.na(df$padj) & !is.na(df$NES), , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)

  n_each  <- ceiling(top_n / 2)
  pos_df  <- df[df$NES > 0, , drop = FALSE]
  neg_df  <- df[df$NES < 0, , drop = FALSE]
  pos_top <- head(pos_df[order(pos_df$padj), ], n_each)
  neg_top <- head(neg_df[order(neg_df$padj), ], n_each)
  plot_df <- rbind(pos_top, neg_top)

  if (nrow(plot_df) == 0L) return(NULL)

  plot_df$label <- ifelse(
    nchar(plot_df$pathway) > 55,
    paste0(substr(plot_df$pathway, 1L, 52L), "\u2026"),
    plot_df$pathway
  )
  plot_df       <- plot_df[order(plot_df$NES), , drop = FALSE]
  plot_df$label <- factor(plot_df$label, levels = unique(plot_df$label))

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = NES, y = label,
                               size = size, colour = -log10(padj))) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.5) +
    ggplot2::scale_colour_gradient(low = "#deebf7", high = "#08519c",
                                    name = "-log10(padj)") +
    ggplot2::scale_size_continuous(name  = "Gene set\noverlap",
                                    range = c(2, 8)) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("Top ", nrow(pos_top), " enriched / ",
                        nrow(neg_top), " depleted (by padj)"),
      x = "NES", y = NULL
    ) +
    KHALKOS_theme_cauldron() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}


# ORA dotplot.  df: gprofiler2 results (term_name, p_value, intersection_size,
# term_size, query_size, effective_domain_size, source).
# Top top_n terms distributed across sources, ordered by fold enrichment.
.talaria_ora_dotplot <- function(df, top_n = 20, title = "ORA") {

  if (nrow(df) == 0L) return(NULL)

  df$fold_enrichment <- (df$intersection_size / df$query_size) /
                        (df$term_size / df$effective_domain_size)
  df$fold_enrichment[!is.finite(df$fold_enrichment)] <- NA_real_

  n_sources   <- length(unique(df$source))
  top_per_src <- max(2L, ceiling(top_n / n_sources))

  plot_df <- do.call(rbind, lapply(split(df, df$source), function(sub) {
    head(sub[order(sub$p_value), ], top_per_src)
  }))
  rownames(plot_df) <- NULL

  plot_df$label <- ifelse(
    nchar(plot_df$term_name) > 55,
    paste0(substr(plot_df$term_name, 1L, 52L), "\u2026"),
    plot_df$term_name
  )
  plot_df       <- plot_df[order(plot_df$source,
                                  plot_df$fold_enrichment), , drop = FALSE]
  plot_df$label <- factor(plot_df$label, levels = unique(plot_df$label))

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = fold_enrichment, y = label,
                               size = intersection_size,
                               colour = -log10(p_value))) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(source ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_gradient(low = "#fee5d9", high = "#a50f15",
                                    name = "-log10(p-value)") +
    ggplot2::scale_size_continuous(name  = "Intersection\nsize",
                                    range = c(2, 8)) +
    ggplot2::labs(title = title, x = "Fold Enrichment", y = NULL) +
    KHALKOS_theme_cauldron() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 8),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)
    )
}


# GSEA enrichment (running score) plot for a single pathway.
# Wraps fgsea::plotEnrichment() and applies CAULDRON theme.
.talaria_gsea_enrichment_plot <- function(pathway_genes, ranked_genes,
                                           title = NULL) {
  if (!requireNamespace("fgsea", quietly = TRUE)) return(NULL)
  p <- fgsea::plotEnrichment(pathway = pathway_genes, stats = ranked_genes)
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p + KHALKOS_theme_cauldron()
}


# ── Public export functions ────────────────────────────────────────────────────

#' Export DE results to files
#'
#' Writes CSV tables, an RDS object, and optional plots for a
#' \code{keraunos_markers}, \code{keraunos_contrast}, or
#' \code{keraunos_pseudobulk} object.
#'
#' Volcano and MA plots are generated for \code{keraunos_pseudobulk} only
#' (requires \code{log2FoldChange}, \code{pvalue}, and \code{baseMean} columns
#' from DESeq2).  Marker and contrast results receive only CSV + RDS.
#'
#' @param de_result A \code{keraunos_markers}, \code{keraunos_contrast}, or
#'   \code{keraunos_pseudobulk} object.
#' @param output_dir Character. Directory to write into (created if absent).
#' @param prefix Character. Filename prefix. Default \code{"de"}.
#' @param plots Logical. Generate volcano and MA plots (pseudobulk only).
#'   Default \code{TRUE}.
#' @param top_n Integer. Number of genes to label on each plot by significance.
#'   Requires \code{ggrepel} (Suggests). Default \code{10L}.
#' @param l2fc_thresh Numeric. \eqn{\log_2} fold-change threshold for
#'   up/down categories on plots. Default \code{1}.
#' @param plot_width Numeric. PNG width in inches. Default \code{8}.
#' @param plot_height Numeric. PNG height in inches. Default \code{6}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return Invisibly, a character vector of file paths written.
#' @export
TALARIA_export_de <- function(de_result,
                               output_dir,
                               prefix      = "de",
                               plots       = TRUE,
                               top_n       = 10L,
                               l2fc_thresh = 1,
                               plot_width  = 8,
                               plot_height = 6,
                               verbose     = TRUE) {

  dir_data  <- file.path(output_dir, "data")
  dir_plots <- file.path(output_dir, "plots")
  dir.create(dir_data,  recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(plots))
    dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

  saved <- character(0)
  clean <- function(x) gsub("[^A-Za-z0-9_]", "_", as.character(x))

  # Move named columns to the front of a data frame (silently skip missing ones)
  .front <- function(df, cols) {
    cols <- intersect(cols, names(df))
    if (length(cols)) df <- df[, c(cols, setdiff(names(df), cols)), drop = FALSE]
    df
  }

  write_csv <- function(df, fname) {
    path <- file.path(dir_data, fname)
    utils::write.csv(df, path, row.names = FALSE)
    saved <<- c(saved, path)
    invisible(path)
  }

  if (inherits(de_result, "keraunos_markers")) {
    for (cl in names(de_result$markers)) {
      df      <- as.data.frame(de_result$markers[[cl]])
      df$gene <- rownames(df)
      df      <- .front(df, c("gene", "cluster"))
      write_csv(df, paste0(prefix, "_markers_cluster_", clean(cl), ".csv"))
    }
    if (!is.null(de_result$top) && nrow(de_result$top) > 0)
      write_csv(.front(de_result$top, c("gene", "cluster")),
                paste0(prefix, "_markers_top.csv"))
    saveRDS(de_result, file.path(dir_data, paste0(prefix, "_markers.rds")))
    saved <- c(saved, file.path(dir_data, paste0(prefix, "_markers.rds")))

  } else if (inherits(de_result, "keraunos_contrast")) {
    if (is.null(de_result$params$cluster_col)) {
      for (g in names(de_result$results)) {
        df      <- as.data.frame(de_result$results[[g]])
        df$gene <- rownames(df)
        write_csv(.front(df, c("gene", "cluster")),
                  paste0(prefix, "_contrast_", clean(g), ".csv"))
      }
    } else {
      for (cl in names(de_result$results)) {
        for (g in names(de_result$results[[cl]])) {
          df      <- as.data.frame(de_result$results[[cl]][[g]])
          df$gene <- rownames(df)
          write_csv(.front(df, c("gene", "cluster")),
                    paste0(prefix, "_contrast_cluster_", clean(cl),
                           "_group_", clean(g), ".csv"))
        }
      }
    }
    if (!is.null(de_result$top) && nrow(de_result$top) > 0)
      write_csv(.front(de_result$top, c("gene", "cluster")),
                paste0(prefix, "_contrast_top.csv"))
    saveRDS(de_result, file.path(dir_data, paste0(prefix, "_contrast.rds")))
    saved <- c(saved, file.path(dir_data, paste0(prefix, "_contrast.rds")))

  } else if (inherits(de_result, "keraunos_pseudobulk")) {
    for (cl in names(de_result$results))
      write_csv(.front(de_result$results[[cl]], c("gene", "cluster")),
                paste0(prefix, "_pseudobulk_cluster_", clean(cl), ".csv"))
    if (!is.null(de_result$summary))
      write_csv(de_result$summary, paste0(prefix, "_pseudobulk_summary.csv"))
    saveRDS(de_result, file.path(dir_data, paste0(prefix, "_pseudobulk.rds")))
    saved <- c(saved, file.path(dir_data, paste0(prefix, "_pseudobulk.rds")))

    if (isTRUE(plots)) {
      sig_col    <- de_result$params$sig_col
      sig_thresh <- de_result$params$sig_thresh
      comparison <- de_result$params$comparison

      for (cl in names(de_result$results)) {
        df    <- de_result$results[[cl]]
        title <- paste0(comparison, " | Cluster ", cl)

        p_vol <- .talaria_plot_volcano(
          df, sig_col = sig_col, sig_thresh = sig_thresh,
          l2fc_thresh = l2fc_thresh, top_n_label = as.integer(top_n),
          title = paste0("Volcano \u2014 ", title)
        )
        p_ma  <- .talaria_plot_ma(
          df, sig_col = sig_col, sig_thresh = sig_thresh,
          l2fc_thresh = l2fc_thresh, top_n_label = as.integer(top_n),
          title = paste0("MA \u2014 ", title)
        )
        vol_path <- file.path(dir_plots,
                              paste0(prefix, "_pseudobulk_cluster_",
                                     clean(cl), "_volcano.png"))
        ma_path  <- file.path(dir_plots,
                              paste0(prefix, "_pseudobulk_cluster_",
                                     clean(cl), "_ma.png"))
        ggplot2::ggsave(vol_path, plot = p_vol, width = plot_width,
                        height = plot_height, dpi = 300, units = "in",
                        limitsize = FALSE)
        ggplot2::ggsave(ma_path,  plot = p_ma,  width = plot_width,
                        height = plot_height, dpi = 300, units = "in",
                        limitsize = FALSE)
        saved <- c(saved, vol_path, ma_path)
      }
    }

  } else {
    stop("'de_result' must be a keraunos_markers, keraunos_contrast, or ",
         "keraunos_pseudobulk object.", call. = FALSE)
  }

  if (isTRUE(verbose))
    message("Exported ", length(saved), " file(s) to: ", output_dir)
  invisible(saved)
}


#' Export GSEA results to files
#'
#' Writes CSV tables, an RDS object, NES dotplots (one per gene set collection),
#' and optional per-pathway enrichment score plots for a \code{keraunos_gsea}
#' or \code{keraunos_gsea_multi} object.  Dotplots are controlled by
#' \code{make_dotplots}; enrichment score plots are controlled by
#' \code{do_gsea_plots}. CSV/RDS export always happens regardless of either
#' setting.
#'
#' For \code{keraunos_gsea_multi} each cluster gets its own subdirectory.
#' Enrichment score plots require \code{$ranked_genes} and \code{$gene_sets}
#' stored in the object (present automatically when produced by
#' \code{\link[CAULDRON]{KERAUNOS_gsea}} or
#' \code{\link[CAULDRON]{KERAUNOS_gsea_pseudobulk}}).
#'
#' @param result A \code{keraunos_gsea} or \code{keraunos_gsea_multi} object.
#' @param output_dir Character. Directory to write into (created if absent).
#' @param prefix Character. Filename prefix. Default \code{"gsea"}.
#' @param make_dotplots Logical. Generate NES dotplots (one per gene set
#'   collection). Set \code{FALSE} when exporting many comparisons at once
#'   (e.g. all pairwise cluster combinations) to avoid writing one dotplot
#'   per comparison — CSV/RDS results are still written either way, so
#'   dotplots can be produced later for a chosen subset by calling this
#'   function again with \code{make_dotplots = TRUE} on just those results.
#'   Default \code{TRUE}.
#' @param do_gsea_plots Logical. Generate per-pathway enrichment score plots
#'   (running-score curves via \code{fgsea::plotEnrichment}).  Default
#'   \code{TRUE}.
#' @param top_n Integer. Pathways shown per dotplot (split equally between
#'   enriched and depleted) and number of enrichment score plots per cluster.
#'   Default \code{20L}.
#' @param plot_width Numeric. PNG width in inches. Default \code{8}.
#' @param plot_height Numeric. PNG height for dotplots in inches. Default \code{6}.
#' @param enrich_plot_height Numeric. PNG height for enrichment score plots.
#'   Default \code{4}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return Invisibly, a character vector of file paths written.
#' @export
TALARIA_export_gsea <- function(result,
                                 output_dir,
                                 prefix             = "gsea",
                                 make_dotplots      = TRUE,
                                 do_gsea_plots      = TRUE,
                                 top_n              = 20L,
                                 plot_width         = 8,
                                 plot_height        = 6,
                                 enrich_plot_height = 4,
                                 verbose            = TRUE) {

  if (!inherits(result, c("keraunos_gsea", "keraunos_gsea_multi")))
    stop("'result' must be a keraunos_gsea or keraunos_gsea_multi object.",
         call. = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saved <- character(0)
  clean <- function(x) gsub("[^A-Za-z0-9_]", "_", as.character(x))

  save_plot <- function(p, path, w = plot_width, h = plot_height) {
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi = 300, units = "in", limitsize = FALSE)
    saved <<- c(saved, path)
  }

  .write_cluster_gsea <- function(cl_result, cl_dir, cl_label) {
    utils::write.csv(
      cl_result$results,
      file.path(cl_dir, paste0(prefix, "_", cl_label, "_results.csv")),
      row.names = FALSE
    )
    saved <<- c(saved,
                file.path(cl_dir,
                          paste0(prefix, "_", cl_label, "_results.csv")))

    if (nrow(cl_result$results) == 0L) return(invisible(NULL))

    # ── Dotplots (gated by make_dotplots) — one per collection if present ────
    if (isTRUE(make_dotplots) &&
        "collection" %in% names(cl_result$results) &&
        !all(is.na(cl_result$results$collection))) {
      collections <- unique(cl_result$results$collection[
        !is.na(cl_result$results$collection)])
      for (col in collections) {
        col_df <- cl_result$results[
          !is.na(cl_result$results$collection) &
          cl_result$results$collection == col, , drop = FALSE]
        p_dot <- .talaria_gsea_dotplot(
          col_df, top_n = as.integer(top_n),
          title = paste0("GSEA \u2014 ", gsub("_", " ", cl_label),
                         " \u2014 ", col)
        )
        if (!is.null(p_dot))
          save_plot(p_dot,
                    file.path(cl_dir,
                              paste0(prefix, "_", cl_label, "_dotplot_",
                                     clean(col), ".png")))
      }
    } else if (isTRUE(make_dotplots)) {
      p_dot <- .talaria_gsea_dotplot(
        cl_result$results, top_n = as.integer(top_n),
        title = paste0("GSEA \u2014 ", gsub("_", " ", cl_label))
      )
      if (!is.null(p_dot))
        save_plot(p_dot,
                  file.path(cl_dir,
                            paste0(prefix, "_", cl_label, "_dotplot.png")))
    }

    # ── Enrichment score plots (optional) ─────────────────────────────────────
    if (isTRUE(do_gsea_plots) &&
        !is.null(cl_result$ranked_genes) &&
        !is.null(cl_result$gene_sets)) {
      top_paths <- head(
        cl_result$results[order(cl_result$results$padj, na.last = TRUE), ],
        as.integer(top_n)
      )
      for (i in seq_len(nrow(top_paths))) {
        pw   <- top_paths$pathway[i]
        gens <- cl_result$gene_sets[[pw]]
        if (is.null(gens)) next
        p_es <- .talaria_gsea_enrichment_plot(
          pathway_genes = gens,
          ranked_genes  = cl_result$ranked_genes,
          title         = paste0(pw, " \u2014 ", gsub("_", " ", cl_label))
        )
        if (!is.null(p_es))
          save_plot(p_es,
                    file.path(cl_dir,
                              paste0(prefix, "_", cl_label,
                                     "_enrichment_", clean(pw), ".png")),
                    h = enrich_plot_height)
      }
    }
  }

  if (inherits(result, "keraunos_gsea")) {
    saveRDS(result, file.path(output_dir, paste0(prefix, ".rds")))
    saved <- c(saved, file.path(output_dir, paste0(prefix, ".rds")))
    .write_cluster_gsea(result, output_dir, prefix)

  } else {
    if (!is.null(result$summary) && nrow(result$summary) > 0L) {
      utils::write.csv(result$summary,
                       file.path(output_dir, paste0(prefix, "_summary.csv")),
                       row.names = FALSE)
      saved <- c(saved,
                 file.path(output_dir, paste0(prefix, "_summary.csv")))
    }
    saveRDS(result, file.path(output_dir, paste0(prefix, "_multi.rds")))
    saved <- c(saved, file.path(output_dir, paste0(prefix, "_multi.rds")))

    for (cl in names(result$results)) {
      cl_dir <- file.path(output_dir, paste0("cluster_", clean(cl)))
      dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
      .write_cluster_gsea(result$results[[cl]], cl_dir,
                           paste0("cluster_", clean(cl)))
    }
  }

  if (isTRUE(verbose))
    message("Exported ", length(saved), " file(s) to: ", output_dir)
  invisible(saved)
}


#' Export ORA results to files
#'
#' Writes CSV tables, an RDS object, and optional fold-enrichment dotplots for
#' a \code{keraunos_ora} or \code{keraunos_ora_multi} object.
#'
#' For \code{keraunos_ora_multi} each cluster gets its own subdirectory.
#' Up- and down-regulated directions are exported separately.
#'
#' @param result A \code{keraunos_ora} or \code{keraunos_ora_multi} object.
#' @param output_dir Character. Directory to write into (created if absent).
#' @param prefix Character. Filename prefix. Default \code{"ora"}.
#' @param plots Logical. Generate dotplots. Default \code{TRUE}.
#' @param top_n Integer. Total terms shown per dotplot, distributed across
#'   sources. Default \code{20L}.
#' @param plot_width Numeric. PNG width in inches. Default \code{8}.
#' @param plot_height Numeric. PNG height in inches. Default \code{7}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return Invisibly, a character vector of file paths written.
#' @export
TALARIA_export_ora <- function(result,
                                output_dir,
                                prefix      = "ora",
                                plots       = TRUE,
                                top_n       = 20L,
                                plot_width  = 8,
                                plot_height = 7,
                                verbose     = TRUE) {

  if (!inherits(result, c("keraunos_ora", "keraunos_ora_multi")))
    stop("'result' must be a keraunos_ora or keraunos_ora_multi object.",
         call. = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saved <- character(0)
  clean <- function(x) gsub("[^A-Za-z0-9_]", "_", as.character(x))

  save_plot <- function(p, path) {
    ggplot2::ggsave(path, plot = p, width = plot_width, height = plot_height,
                    dpi = 300, units = "in", limitsize = FALSE)
    saved <<- c(saved, path)
  }

  .write_direction_ora <- function(ora_obj, out_dir, label) {
    if (is.null(ora_obj)) return(invisible(NULL))
    csv_path <- file.path(out_dir, paste0(prefix, "_", label, "_results.csv"))
    utils::write.csv(ora_obj$results, csv_path, row.names = FALSE)
    saved <<- c(saved, csv_path)

    if (isTRUE(plots) && nrow(ora_obj$significant) > 0L) {
      p_dot <- .talaria_ora_dotplot(
        ora_obj$significant, top_n = as.integer(top_n),
        title = paste0("ORA \u2014 ", gsub("_", " ", label))
      )
      if (!is.null(p_dot))
        save_plot(p_dot,
                  file.path(out_dir,
                            paste0(prefix, "_", label, "_dotplot.png")))
    }
  }

  if (inherits(result, "keraunos_ora")) {
    saveRDS(result, file.path(output_dir, paste0(prefix, ".rds")))
    saved <- c(saved, file.path(output_dir, paste0(prefix, ".rds")))
    .write_direction_ora(result, output_dir, prefix)

  } else {
    if (!is.null(result$summary) && nrow(result$summary) > 0L) {
      utils::write.csv(result$summary,
                       file.path(output_dir, paste0(prefix, "_summary.csv")),
                       row.names = FALSE)
      saved <- c(saved,
                 file.path(output_dir, paste0(prefix, "_summary.csv")))
    }
    saveRDS(result, file.path(output_dir, paste0(prefix, "_multi.rds")))
    saved <- c(saved, file.path(output_dir, paste0(prefix, "_multi.rds")))

    for (cl in names(result$results)) {
      cl_dir <- file.path(output_dir, paste0("cluster_", clean(cl)))
      dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
      for (direction in c("up", "down")) {
        .write_direction_ora(
          result$results[[cl]][[direction]],
          cl_dir,
          paste0("cluster_", clean(cl), "_", direction)
        )
      }
    }
  }

  if (isTRUE(verbose))
    message("Exported ", length(saved), " file(s) to: ", output_dir)
  invisible(saved)
}
