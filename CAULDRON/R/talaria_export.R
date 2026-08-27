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
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
}


# ORA dotplot.  df: gprofiler2 results (term_name, p_value, intersection_size,
# term_size, query_size, effective_domain_size, source).
# Top top_n terms distributed across sources, ordered by fold enrichment.
# low_colour/high_colour set the -log10(p-value) gradient -- used to give the
# "up" and "down" panels distinct colour families when combined side by side.
.talaria_ora_dotplot <- function(df, top_n = 20, title = "ORA",
                                  low_colour = "#fee5d9", high_colour = "#a50f15") {

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
    ggplot2::scale_colour_gradient(low = low_colour, high = high_colour,
                                    name = "-log10(p-value)") +
    ggplot2::scale_size_continuous(name  = "Intersection\nsize",
                                    range = c(2, 8)) +
    ggplot2::labs(title = title, x = "Fold Enrichment", y = NULL) +
    KHALKOS_theme_cauldron() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 10),
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
#' Up- and down-regulated directions are combined into a single results CSV
#' (distinguished by a \code{direction} column) and a single dotplot image
#' with the "up" and "down" panels shown side by side, each keeping its own
#' fold-enrichment axis and significance colour gradient (red for up, blue
#' for down) rather than a shared scale.
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
      cl_dir   <- file.path(output_dir, paste0("cluster_", clean(cl)))
      cl_label <- paste0("cluster_", clean(cl))
      dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

      up_obj   <- result$results[[cl]][["up"]]
      down_obj <- result$results[[cl]][["down"]]

      # Combined results table, tagged by direction
      combined_df <- do.call(rbind, c(
        if (!is.null(up_obj))   list(cbind(direction = "up",   up_obj$results))   else NULL,
        if (!is.null(down_obj)) list(cbind(direction = "down", down_obj$results)) else NULL
      ))
      if (!is.null(combined_df) && nrow(combined_df) > 0L) {
        csv_path <- file.path(cl_dir, paste0(prefix, "_", cl_label, "_results.csv"))
        utils::write.csv(combined_df, csv_path, row.names = FALSE)
        saved <- c(saved, csv_path)
      }

      # Combined dotplot: "up" (red) and "down" (blue) panels side by side,
      # each keeping its own fold-enrichment axis and colour scale.
      if (isTRUE(plots)) {
        p_up <- if (!is.null(up_obj) && nrow(up_obj$significant) > 0L)
          .talaria_ora_dotplot(
            up_obj$significant, top_n = as.integer(top_n),
            title = paste0("ORA \u2014 ", gsub("_", " ", cl_label), " \u2014 up"),
            low_colour = "#fee5d9", high_colour = "#a50f15"
          )

        p_down <- if (!is.null(down_obj) && nrow(down_obj$significant) > 0L)
          .talaria_ora_dotplot(
            down_obj$significant, top_n = as.integer(top_n),
            title = paste0("ORA \u2014 ", gsub("_", " ", cl_label), " \u2014 down"),
            low_colour = "#eff3ff", high_colour = "#08519c"
          )

        panels <- Filter(Negate(is.null), list(p_up, p_down))
        if (length(panels) > 0L) {
          if (!requireNamespace("gridExtra", quietly = TRUE))
            stop("Package 'gridExtra' is required for combined ORA dotplots. ",
                 "Install via: install.packages(\"gridExtra\")", call. = FALSE)

          combined_plot <- gridExtra::arrangeGrob(grobs = panels, ncol = length(panels))
          plot_path <- file.path(cl_dir, paste0(prefix, "_", cl_label, "_dotplot.png"))
          ggplot2::ggsave(plot_path, plot = combined_plot,
                          width = plot_width * length(panels), height = plot_height,
                          dpi = 300, units = "in", limitsize = FALSE)
          saved <- c(saved, plot_path)
        }
      }
    }
  }

  if (isTRUE(verbose))
    message("Exported ", length(saved), " file(s) to: ", output_dir)
  invisible(saved)
}


#' Export an interactive gene-expression-by-group widget for R Markdown
#'
#' Builds a self-contained, dependency-free HTML/JS widget that lets a report
#' viewer search any gene and see its mean expression across levels of one or
#' more grouping variables (e.g. cluster, genotype, cell label), rendered as
#' an SVG bar chart. Meant to be called directly inside an R Markdown chunk:
#' the return value knits as raw HTML automatically (see
#' \code{\link[knitr]{asis_output}}) — no \code{results='asis'} chunk option
#' or manual \code{cat()} needed.
#'
#' Expression is precomputed as a gene x group-level \strong{mean} matrix, one
#' per grouping variable, not exported per cell. For a full transcriptome
#' this keeps the embedded payload on the order of a few MB rather than the
#' multi-GB a per-cell export would require, at the cost of showing each
#' group's mean rather than per-cell resolution — plotting a per-cell UMAP
#' colored by a per-cluster average would imply resolution the data doesn't
#' have, so a bar chart is used instead.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param groupings Character vector of \code{colData(sce)} column names to
#'   expose as grouping options (e.g. \code{c("cluster", "genotype",
#'   "cell_label")}).
#' @param assay_name Character. Assay averaged per group. Default
#'   \code{"logcounts"}.
#' @param widget_id Character. DOM id prefix — must be unique if more than
#'   one widget is embedded in the same document. Default \code{"gene-expr"}.
#' @param digits Integer. Rounding applied to exported mean values. Default
#'   \code{4}.
#'
#' @return An object of class \code{knit_asis} (see
#'   \code{\link[knitr]{asis_output}}). Printing or auto-printing it inside
#'   an R Markdown chunk renders the widget; outside of knitr it just prints
#'   as plain HTML text.
#' @export
TALARIA_export_gene_expr_widget <- function(sce,
                                             groupings  = NULL,
                                             assay_name = "logcounts",
                                             widget_id  = "gene-expr",
                                             digits     = 4) {

  if (!inherits(sce, "SingleCellExperiment"))
    stop("'sce' must be a SingleCellExperiment.", call. = FALSE)
  if (is.null(groupings) || length(groupings) == 0L)
    stop("'groupings' must name at least one colData(sce) column.",
         call. = FALSE)
  missing_cols <- setdiff(groupings, names(colData(sce)))
  if (length(missing_cols) > 0L)
    stop("colData(sce) is missing: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  if (!assay_name %in% assayNames(sce))
    stop("'assay_name' (\"", assay_name, "\") not found in assayNames(sce).",
         call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for TALARIA_export_gene_expr_widget(). ",
         "Install it with install.packages(\"jsonlite\").", call. = FALSE)
  if (!requireNamespace("knitr", quietly = TRUE))
    stop("Package 'knitr' is required for TALARIA_export_gene_expr_widget(). ",
         "Install it with install.packages(\"knitr\").", call. = FALSE)

  expr_mat <- assay(sce, assay_name)

  grouping_export <- lapply(groupings, function(g) {
    grp <- as.character(colData(sce)[[g]])
    lv <- unique(grp)
    lv_num <- suppressWarnings(as.numeric(lv))
    lv <- if (!anyNA(lv_num)) as.character(sort(lv_num)) else sort(lv)
    mat <- sapply(lv, function(l) Matrix::rowMeans(expr_mat[, grp == l, drop = FALSE]))
    list(levels = lv, matrix = round(mat, digits))
  })
  names(grouping_export) <- groupings

  export_obj <- list(genes = rownames(sce), groupings = grouping_export)
  json_payload <- jsonlite::toJSON(export_obj, auto_unbox = TRUE)

  knitr::asis_output(.talaria_gene_expr_widget_html(widget_id, json_payload))
}


# Builds the widget's HTML/CSS/JS from a template, substituting the DOM id
# prefix (so multiple widgets can coexist on one page) and the precomputed
# JSON payload. Not exported -- called only from TALARIA_export_gene_expr_widget().
.talaria_gene_expr_widget_html <- function(widget_id, json_payload) {

  template <- r"-(
<div id="__ID__-widget">
<div class="filter-row">
<label for="__ID__-grouping-select">Group by</label>
<select id="__ID__-grouping-select"></select>
<label for="__ID__-gene-input">Gene</label>
<input id="__ID__-gene-input" placeholder="e.g. TH" autocomplete="off">
</div>
<div id="__ID__-status"></div>
<svg id="__ID__-chart" width="700" height="380"></svg>
</div>

<style>
#__ID__-widget .filter-row { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; }
#__ID__-widget label { font-size: 13px; }
#__ID__-status { font-size: 12px; color: #b00; min-height: 16px; }
#__ID__-chart text { fill: #333; }
</style>

<script id="__ID__-data" type="application/json">
__JSON__
</script>

<script>
(function () {
var payload = JSON.parse(document.getElementById('__ID__-data').textContent);
var geneIndex = new Map();
payload.genes.forEach(function (g, i) { geneIndex.set(g.toUpperCase(), i); });

var groupingSelect = document.getElementById('__ID__-grouping-select');
var geneInput = document.getElementById('__ID__-gene-input');
var statusEl = document.getElementById('__ID__-status');
var svg = document.getElementById('__ID__-chart');
var BAR_COLOR = '#2a78d6';
var GRID_COLOR = '#e1e0d9';
var AXIS_COLOR = '#c3c2b7';
var TEXT_COLOR = '#52514e';
var ns = 'http://www.w3.org/2000/svg';

function niceMax(value) {
if (value <= 0) return 1;
var exponent = Math.floor(Math.log10(value));
var fraction = value / Math.pow(10, exponent);
var niceFraction = fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
return niceFraction * Math.pow(10, exponent);
}

Object.keys(payload.groupings).forEach(function (g) {
var opt = document.createElement('option');
opt.value = g; opt.textContent = g;
groupingSelect.appendChild(opt);
});

function draw() {
while (svg.firstChild) svg.removeChild(svg.firstChild);
var gene = geneInput.value.trim();
if (!gene) { statusEl.textContent = ''; return; }
var idx = geneIndex.get(gene.toUpperCase());
if (idx === undefined) {
statusEl.textContent = 'Gene "' + gene + '" not found.';
return;
}
statusEl.textContent = '';

var grp = payload.groupings[groupingSelect.value];
var values = grp.matrix[idx];
var levels = grp.levels;

var ROT_DEG = 45;
var longestLabel = levels.reduce(function (a, b) { return String(b).length > a ? String(b).length : a; }, 0);
var estLabelPx = longestLabel * 6.2;
// Full label length, since text-anchor="start" below makes the whole label
// hang down-and-right from its pivot rather than splitting across it.
var padB = Math.max(50, estLabelPx * Math.sin(ROT_DEG * Math.PI / 180) + 24);

var plotH = 260;
var padL = 60, padT = 24;
var basePadR = 20;
var plotW = 700 - padL - basePadR; // bar layout is fixed, independent of margin

// The rightmost bar's label (text-anchor="start", rotated) hangs down-and-
// right from its pivot and can run past a fixed-width canvas -- widen just
// the right margin to fit it, without touching the bar layout above.
var labelOverhangX = estLabelPx * Math.cos(ROT_DEG * Math.PI / 180);
var padR = Math.max(basePadR, labelOverhangX + 15);

var width = padL + plotW + padR;
var height = padT + plotH + padB;
svg.setAttribute('width', width);
svg.setAttribute('height', height);

var rawMax = Math.max.apply(null, values.concat([0.001]));
var axisMax = niceMax(rawMax);
var barW = plotW / values.length;

var title = document.createElementNS(ns, 'text');
title.setAttribute('x', padL);
title.setAttribute('y', 14);
title.setAttribute('font-size', '12');
title.setAttribute('font-weight', '600');
title.textContent = gene + ' — mean expression by ' + groupingSelect.value;
svg.appendChild(title);

var nTicks = 4;
for (var t = 0; t <= nTicks; t++) {
var tickVal = (axisMax / nTicks) * t;
var ty = padT + plotH - (tickVal / axisMax) * plotH;

var grid = document.createElementNS(ns, 'line');
grid.setAttribute('x1', padL);
grid.setAttribute('x2', width - padR);
grid.setAttribute('y1', ty);
grid.setAttribute('y2', ty);
grid.setAttribute('stroke', GRID_COLOR);
grid.setAttribute('stroke-width', '1');
svg.appendChild(grid);

var tickLabel = document.createElementNS(ns, 'text');
tickLabel.setAttribute('x', padL - 8);
tickLabel.setAttribute('y', ty + 3);
tickLabel.setAttribute('text-anchor', 'end');
tickLabel.setAttribute('font-size', '10');
tickLabel.setAttribute('fill', TEXT_COLOR);
tickLabel.textContent = Number(tickVal.toFixed(2)).toString();
svg.appendChild(tickLabel);
}

var axisLine = document.createElementNS(ns, 'line');
axisLine.setAttribute('x1', padL);
axisLine.setAttribute('x2', padL);
axisLine.setAttribute('y1', padT);
axisLine.setAttribute('y2', padT + plotH);
axisLine.setAttribute('stroke', AXIS_COLOR);
axisLine.setAttribute('stroke-width', '1');
svg.appendChild(axisLine);

var axisTitle = document.createElementNS(ns, 'text');
axisTitle.setAttribute('x', -(padT + plotH / 2));
axisTitle.setAttribute('y', 14);
axisTitle.setAttribute('text-anchor', 'middle');
axisTitle.setAttribute('font-size', '10');
axisTitle.setAttribute('fill', TEXT_COLOR);
axisTitle.setAttribute('transform', 'rotate(-90)');
axisTitle.textContent = 'Mean logcounts expression';
svg.appendChild(axisTitle);

values.forEach(function (v, i) {
var barH = (v / axisMax) * plotH;
var x = padL + i * barW + barW * 0.15;
var y = padT + (plotH - barH);

var rect = document.createElementNS(ns, 'rect');
rect.setAttribute('x', x);
rect.setAttribute('y', y);
rect.setAttribute('width', barW * 0.7);
rect.setAttribute('height', barH);
rect.setAttribute('fill', BAR_COLOR);
svg.appendChild(rect);

var valLabel = document.createElementNS(ns, 'text');
valLabel.setAttribute('x', x + barW * 0.35);
valLabel.setAttribute('y', y - 4);
valLabel.setAttribute('text-anchor', 'middle');
valLabel.setAttribute('font-size', '10');
valLabel.textContent = v.toFixed(2);
svg.appendChild(valLabel);

// text-anchor="start" pins the pivot at the label's own top-left corner, so
// the whole string swings down-and-right from the axis when rotated -- with
// "middle", half the string would swing the other way and intrude upward
// into the bars instead of reading cleanly below them.
var labelX = x + barW * 0.35;
var labelY = height - padB + 8;
var label = document.createElementNS(ns, 'text');
label.setAttribute('x', labelX);
label.setAttribute('y', labelY);
label.setAttribute('text-anchor', 'start');
label.setAttribute('font-size', '11');
label.setAttribute('transform', 'rotate(' + ROT_DEG + ' ' + labelX + ' ' + labelY + ')');
label.textContent = levels[i];
svg.appendChild(label);
});
}

groupingSelect.addEventListener('change', draw);
geneInput.addEventListener('input', draw);
})();
</script>
)-"

  html <- gsub("__ID__", widget_id, template, fixed = TRUE)
  sub("__JSON__", json_payload, html, fixed = TRUE)
}


#' Export a searchable rowData(sce) table widget for R Markdown
#'
#' Builds a self-contained, dependency-free HTML/JS widget that lets a report
#' viewer search genes by name (matched against \code{rownames(sce)}) and see
#' their associated \code{rowData(sce)} columns in a live-filtered table.
#' Mirrors \code{\link{TALARIA_export_gene_expr_widget}}'s embedded-JSON +
#' vanilla-JS approach so the report stays free of external JS dependencies
#' (a package like DT or crosstalk would each pull in their own bundled
#' JS/CSS) and stays portable as a single, self-contained HTML file. Meant to
#' be called directly inside an R Markdown chunk: the return value knits as
#' raw HTML automatically (see \code{\link[knitr]{asis_output}}) -- no
#' \code{results='asis'} chunk option or manual \code{cat()} needed.
#'
#' Only the currently matching rows are ever rendered into the DOM (capped at
#' \code{max_rows}), rather than every gene at once -- with a full
#' transcriptome (tens of thousands of genes), building a table row for every
#' gene up front would bloat the page and make typing sluggish.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param columns Character vector of \code{rowData(sce)} column names to
#'   expose. \code{NULL} (default) exposes every column.
#' @param widget_id Character. DOM id prefix -- must be unique if more than
#'   one widget is embedded in the same document. Default \code{"gene-table"}.
#' @param max_rows Integer. Maximum number of matching rows rendered into the
#'   table at once. Default \code{200}.
#'
#' @return An object of class \code{knit_asis} (see
#'   \code{\link[knitr]{asis_output}}). Printing or auto-printing it inside
#'   an R Markdown chunk renders the widget; outside of knitr it just prints
#'   as plain HTML text.
#' @export
TALARIA_export_gene_table_widget <- function(sce,
                                              columns   = NULL,
                                              widget_id = "gene-table",
                                              max_rows  = 200) {

  if (!inherits(sce, "SingleCellExperiment"))
    stop("'sce' must be a SingleCellExperiment.", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for TALARIA_export_gene_table_widget(). ",
         "Install it with install.packages(\"jsonlite\").", call. = FALSE)
  if (!requireNamespace("knitr", quietly = TRUE))
    stop("Package 'knitr' is required for TALARIA_export_gene_table_widget(). ",
         "Install it with install.packages(\"knitr\").", call. = FALSE)

  rd <- as.data.frame(rowData(sce))
  if (!is.null(columns)) {
    missing_cols <- setdiff(columns, colnames(rd))
    if (length(missing_cols) > 0L)
      stop("rowData(sce) is missing: ", paste(missing_cols, collapse = ", "),
           call. = FALSE)
    rd <- rd[, columns, drop = FALSE]
  }
  rd <- cbind(gene = rownames(sce), rd)

  # Row-oriented array of objects (one JSON object per gene), not
  # column-oriented, so the JS side can filter/slice by row directly without
  # reshaping the payload.
  rows_json <- jsonlite::toJSON(rd, dataframe = "rows", na = "null")
  cols_json <- jsonlite::toJSON(colnames(rd))

  knitr::asis_output(
    .talaria_gene_table_widget_html(widget_id, rows_json, cols_json, max_rows)
  )
}


# Builds the searchable rowData table widget's HTML/CSS/JS from a template,
# substituting the DOM id prefix, the precomputed row/column JSON, and the
# render cap. Not exported -- called only from
# TALARIA_export_gene_table_widget().
.talaria_gene_table_widget_html <- function(widget_id, rows_json, cols_json, max_rows) {

  template <- r"-(
<div id="__ID__-widget">
<div class="filter-row">
<label for="__ID__-search">Search gene</label>
<input id="__ID__-search" type="text" placeholder="e.g. TH" autocomplete="off">
</div>
<div id="__ID__-status"></div>
<div class="table-scroll">
<table id="__ID__-table"><thead></thead><tbody></tbody></table>
</div>
</div>

<style>
#__ID__-widget .filter-row { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; }
#__ID__-widget label { font-size: 13px; }
#__ID__-status { font-size: 12px; color: #52514e; min-height: 16px; margin-bottom: 6px; }
#__ID__-widget .table-scroll { max-height: 480px; overflow: auto; border: 1px solid #e1e0d9; }
#__ID__-widget table { border-collapse: collapse; width: 100%; font-size: 12px; }
#__ID__-widget th, #__ID__-widget td { padding: 4px 8px; border-bottom: 1px solid #e1e0d9; text-align: left; white-space: nowrap; }
#__ID__-widget th { position: sticky; top: 0; background: #f5f4ef; }
</style>

<script id="__ID__-rows" type="application/json">
__ROWS_JSON__
</script>
<script id="__ID__-cols" type="application/json">
__COLS_JSON__
</script>

<script>
(function () {
var rows = JSON.parse(document.getElementById('__ID__-rows').textContent);
var cols = JSON.parse(document.getElementById('__ID__-cols').textContent);
var MAX_ROWS = __MAX_ROWS__;

var searchEl = document.getElementById('__ID__-search');
var statusEl = document.getElementById('__ID__-status');
var table = document.getElementById('__ID__-table');
var thead = table.querySelector('thead');
var tbody = table.querySelector('tbody');

var headRow = document.createElement('tr');
cols.forEach(function (c) {
var th = document.createElement('th');
th.textContent = c;
headRow.appendChild(th);
});
thead.appendChild(headRow);

function render(matches) {
while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
matches.slice(0, MAX_ROWS).forEach(function (row) {
var tr = document.createElement('tr');
cols.forEach(function (c) {
var td = document.createElement('td');
var val = row[c];
td.textContent = (val === null || val === undefined) ? '' : val;
tr.appendChild(td);
});
tbody.appendChild(tr);
});
if (matches.length > MAX_ROWS) {
statusEl.textContent = 'Showing ' + MAX_ROWS + ' of ' + matches.length + ' matches -- refine your search to see more.';
} else {
statusEl.textContent = matches.length + ' gene(s)';
}
}

function filter() {
var q = searchEl.value.trim().toUpperCase();
if (!q) { render(rows.slice(0, MAX_ROWS)); return; }
var matches = rows.filter(function (row) {
return String(row.gene).toUpperCase().indexOf(q) !== -1;
});
render(matches);
}

searchEl.addEventListener('input', filter);
filter();
})();
</script>
)-"

  html <- gsub("__ID__", widget_id, template, fixed = TRUE)
  html <- sub("__ROWS_JSON__", rows_json, html, fixed = TRUE)
  html <- sub("__COLS_JSON__", cols_json, html, fixed = TRUE)
  sub("__MAX_ROWS__", max_rows, html, fixed = TRUE)
}


# Shared "union of significant genes + per-comparison log2FC/padj matrices"
# logic, used by both TALARIA_export_de_comparison_widget() and
# TALARIA_build_de_comparison_table() so the two stay in sync when called
# with matching arguments. Not exported.
.talaria_de_comparison_core <- function(de_list, sig_col, sig_thresh, l2fc_thresh) {

  if (!is.list(de_list) || length(de_list) == 0L ||
      is.null(names(de_list)) || any(names(de_list) == ""))
    stop("'de_list' must be a non-empty named list of data.frames (one per comparison).",
         call. = FALSE)
  # 'padj' is always required (always reported regardless of which column
  # gates significance); sig_col is required too, and may be 'padj' itself
  # (the default) or a different column (e.g. 'pvalue', for pipelines that
  # call significance on the raw p-value instead of FDR).
  req_cols <- unique(c("gene", "log2FoldChange", "padj", sig_col))
  for (nm in names(de_list)) {
    missing_cols <- setdiff(req_cols, names(de_list[[nm]]))
    if (length(missing_cols) > 0L)
      stop("de_list[[\"", nm, "\"]] is missing column(s): ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  comparisons <- names(de_list)

  # Union of genes significant in at least one comparison -- every value in
  # a selected/listed gene's row is still shown for every comparison,
  # significant there or not. Significance is judged on 'sig_col' (padj by
  # default), but 'padj' itself is always reported regardless, since that's
  # the more standard value to report.
  sig_genes <- unique(unlist(lapply(de_list, function(df) {
    mask <- !is.na(df[[sig_col]]) & !is.na(df$log2FoldChange) &
      df[[sig_col]] < sig_thresh & abs(df$log2FoldChange) >= l2fc_thresh
    df$gene[mask]
  })))
  sig_genes <- sort(sig_genes)

  if (length(sig_genes) == 0L)
    stop("No gene was significant (", sig_col, " < ", sig_thresh, ", |log2FC| >= ",
         l2fc_thresh, ") in any comparison in 'de_list'.", call. = FALSE)

  # gene x comparison matrices, NA where a gene wasn't tested in that
  # comparison (e.g. filtered out by min_cells/min_samples upstream).
  log2fc_mat <- matrix(NA_real_, nrow = length(sig_genes), ncol = length(comparisons),
                        dimnames = list(sig_genes, comparisons))
  padj_mat   <- log2fc_mat

  for (nm in comparisons) {
    df  <- de_list[[nm]]
    idx <- match(sig_genes, df$gene)
    hit <- !is.na(idx)
    log2fc_mat[hit, nm] <- df$log2FoldChange[idx[hit]]
    padj_mat[hit, nm]   <- df$padj[idx[hit]]
  }

  list(genes = sig_genes, comparisons = comparisons,
       log2fc = log2fc_mat, padj = padj_mat)
}


#' Build a browsable table of DE results across several comparisons
#'
#' @description Returns a wide data.frame -- one row per gene significant in
#' at least one comparison, with a log2FoldChange and padj column pair per
#' comparison -- meant to be rendered with \code{DT::datatable()} directly
#' above \code{\link{TALARIA_export_de_comparison_widget}} in the same
#' report, so a viewer can browse (sort/search) the full list of searchable
#' genes before picking one to plot in the widget. Uses the same
#' significance/union logic as that function, so the two stay in sync when
#' called with matching \code{sig_col}/\code{sig_thresh}/\code{l2fc_thresh}.
#'
#' @param de_list Named list of data.frames, one per comparison (names
#'   become the comparison shown in each column pair's name, in the order
#'   given). Each data.frame must have \code{gene}, \code{log2FoldChange},
#'   and \code{padj} columns (plus \code{sig_col} itself, if different from
#'   \code{"padj"}) -- the same shape written by \code{TALARIA_export_de()}'s
#'   per-cluster pseudobulk DE CSVs.
#' @param sig_col Character. Column used to decide which genes are included
#'   (i.e. significant in at least one comparison). Default \code{"padj"}.
#'   Set to \code{"pvalue"} for pipelines that call significance on the raw
#'   p-value instead of FDR.
#' @param sig_thresh Numeric. Cutoff applied to \code{sig_col}. Default
#'   \code{0.05}.
#' @param l2fc_thresh Numeric. Absolute log2 fold change cutoff used to
#'   decide which genes are included. Default \code{0.5}.
#' @param digits Integer. Rounding applied to log2FoldChange columns (padj
#'   columns are kept to 3 significant figures regardless). Default \code{3}.
#'
#' @return A data.frame with one \code{Gene} column plus a
#'   \code{"<comparison> log2FC"} / \code{"<comparison> padj"} column pair
#'   per comparison in \code{de_list}, in the order given.
#'
#' @examples
#' \dontrun{
#' tbl <- TALARIA_build_de_comparison_table(de_list, sig_col = "pvalue")
#' DT::datatable(tbl, rownames = FALSE)
#' }
#' @export
TALARIA_build_de_comparison_table <- function(de_list,
                                               sig_col     = "padj",
                                               sig_thresh  = 0.05,
                                               l2fc_thresh = 0.5,
                                               digits      = 3) {

  core <- .talaria_de_comparison_core(de_list, sig_col, sig_thresh, l2fc_thresh)

  out <- data.frame(Gene = core$genes, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in core$comparisons) {
    out[[paste0(nm, " log2FC")]] <- round(core$log2fc[, nm], digits)
    out[[paste0(nm, " padj")]]   <- signif(core$padj[, nm], 3)
  }
  rownames(out) <- NULL
  out
}


#' Export an interactive DE cross-comparison bar-chart widget for R Markdown
#'
#' @description Builds a self-contained, dependency-free HTML/JS widget that
#' lets a report viewer search any gene that was called significant in at
#' least one of several differential expression comparisons, and see a bar
#' chart of that gene's log2 fold change across \emph{every} comparison
#' supplied -- including comparisons where the gene did not reach
#' significance -- so all comparisons are visible side by side for one gene.
#' Each bar is labelled with its log2 fold change and its adjusted p-value.
#' Mirrors \code{\link{TALARIA_export_gene_expr_widget}}'s embedded-JSON +
#' vanilla-JS approach so the report stays free of external JS dependencies
#' (a package like plotly or crosstalk would each pull in their own bundled
#' JS/CSS) and stays portable as a single, self-contained HTML file.
#'
#' @param de_list Named list of data.frames, one per comparison (names
#'   become the comparison labels shown on the x-axis, in the order given).
#'   Each data.frame must have \code{gene}, \code{log2FoldChange}, and
#'   \code{padj} columns (plus \code{sig_col} itself, if different from
#'   \code{"padj"}) -- the same shape written by \code{TALARIA_export_de()}'s
#'   per-cluster pseudobulk DE CSVs.
#' @param sig_col Character. Column used to decide which genes are
#'   searchable (i.e. significant in at least one comparison). Default
#'   \code{"padj"}. Set to \code{"pvalue"} for pipelines that call
#'   significance on the raw p-value instead of FDR (e.g.
#'   \code{KERAUNOS_de_pseudobulk(..., fdr_threshold = NULL)}) -- the bar
#'   label always shows \code{padj} regardless of this setting.
#' @param sig_thresh Numeric. Cutoff applied to \code{sig_col}. Default
#'   \code{0.05}.
#' @param l2fc_thresh Numeric. Absolute log2 fold change cutoff used only to
#'   decide which genes are searchable. Default \code{0.5}.
#' @param widget_id Character. DOM id prefix -- must be unique if more than
#'   one widget is embedded in the same document. Default
#'   \code{"de-comparison"}.
#' @param digits Integer. Rounding applied to exported log2FoldChange values
#'   (padj is kept to 3 significant figures regardless). Default \code{4}.
#'
#' @details
#' A gene not tested in a given comparison (e.g. filtered out upstream by
#' \code{min_cells}/\code{min_samples}) is shown as a "not tested" marker at
#' the zero line for that comparison rather than a bar, so its absence isn't
#' mistaken for a log2FoldChange of exactly zero.
#'
#' @return An object of class \code{knit_asis} (see
#'   \code{\link[knitr]{asis_output}}). Printing or auto-printing it inside
#'   an R Markdown chunk renders the widget; outside of knitr it just prints
#'   as plain HTML text.
#'
#' @examples
#' \dontrun{
#' TALARIA_export_de_comparison_widget(
#'   de_list = list(
#'     "SNc DA Neuron: AST23 vs CL21" = de_results[["SNc_DA_Neuron_CL21_vs_AST23"]]$results[["SNc_DA_Neuron"]],
#'     "SNc DA Neuron: AST23 vs CL18" = de_results[["SNc_DA_Neuron_CL18_vs_AST23"]]$results[["SNc_DA_Neuron"]],
#'     "DA Neuron: AST23 vs CL21"     = de_results[["DA_Neuron_CL21_vs_AST23"]]$results[["DA_Neuron"]],
#'     "DA Neuron: AST23 vs CL18"     = de_results[["DA_Neuron_CL18_vs_AST23"]]$results[["DA_Neuron"]]
#'   ),
#'   sig_col = "pvalue", sig_thresh = 0.05, l2fc_thresh = 0.5
#' )
#' }
#' @export
TALARIA_export_de_comparison_widget <- function(de_list,
                                                 sig_col     = "padj",
                                                 sig_thresh  = 0.05,
                                                 l2fc_thresh = 0.5,
                                                 widget_id   = "de-comparison",
                                                 digits      = 4) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for TALARIA_export_de_comparison_widget(). ",
         "Install it with install.packages(\"jsonlite\").", call. = FALSE)
  if (!requireNamespace("knitr", quietly = TRUE))
    stop("Package 'knitr' is required for TALARIA_export_de_comparison_widget(). ",
         "Install it with install.packages(\"knitr\").", call. = FALSE)

  core <- .talaria_de_comparison_core(de_list, sig_col, sig_thresh, l2fc_thresh)

  export_obj <- list(
    genes       = core$genes,
    comparisons = core$comparisons,
    log2fc      = round(unname(core$log2fc), digits),
    padj        = signif(unname(core$padj), 3)
  )
  json_payload <- jsonlite::toJSON(export_obj, auto_unbox = TRUE, na = "null")

  knitr::asis_output(.talaria_de_comparison_widget_html(widget_id, json_payload))
}


# Builds the DE cross-comparison widget's HTML/CSS/JS from a template,
# substituting the DOM id prefix (so multiple widgets can coexist on one
# page) and the precomputed JSON payload. Not exported -- called only from
# TALARIA_export_de_comparison_widget().
.talaria_de_comparison_widget_html <- function(widget_id, json_payload) {

  template <- r"-(
<div id="__ID__-widget">
<div class="filter-row">
<label for="__ID__-gene-input">Gene (significant in at least one comparison)</label>
<input id="__ID__-gene-input" list="__ID__-gene-list" placeholder="e.g. TH" autocomplete="off">
<datalist id="__ID__-gene-list"></datalist>
</div>
<div id="__ID__-status"></div>
<svg id="__ID__-chart" width="700" height="380"></svg>
</div>

<style>
#__ID__-widget .filter-row { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; }
#__ID__-widget label { font-size: 13px; }
#__ID__-status { font-size: 12px; color: #b00; min-height: 16px; }
#__ID__-chart text { fill: #333; }
</style>

<script id="__ID__-data" type="application/json">
__JSON__
</script>

<script>
(function () {
var payload = JSON.parse(document.getElementById('__ID__-data').textContent);
var geneIndex = new Map();
payload.genes.forEach(function (g, i) { geneIndex.set(g.toUpperCase(), i); });

var geneInput = document.getElementById('__ID__-gene-input');
var geneList  = document.getElementById('__ID__-gene-list');
var statusEl  = document.getElementById('__ID__-status');
var svg       = document.getElementById('__ID__-chart');
var UP_COLOR      = '#a50f15';
var DOWN_COLOR    = '#08519c';
var NOTTEST_COLOR = '#aaaaaa';
var GRID_COLOR    = '#e1e0d9';
var AXIS_COLOR    = '#c3c2b7';
var TEXT_COLOR    = '#52514e';
var ns = 'http://www.w3.org/2000/svg';

payload.genes.forEach(function (g) {
var opt = document.createElement('option');
opt.value = g;
geneList.appendChild(opt);
});

function niceMax(value) {
if (value <= 0) return 1;
var exponent = Math.floor(Math.log10(value));
var fraction = value / Math.pow(10, exponent);
var niceFraction = fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
return niceFraction * Math.pow(10, exponent);
}

function fmtPadj(p) {
if (p === null || p === undefined || isNaN(p)) return 'n/a';
if (p < 0.0001) return p.toExponential(1);
return p.toFixed(4);
}

function draw() {
while (svg.firstChild) svg.removeChild(svg.firstChild);
var gene = geneInput.value.trim();
if (!gene) { statusEl.textContent = ''; return; }
var idx = geneIndex.get(gene.toUpperCase());
if (idx === undefined) {
statusEl.textContent = 'Gene "' + gene + '" not found among genes significant in at least one comparison.';
return;
}
statusEl.textContent = '';

var comparisons = payload.comparisons;
var values = payload.log2fc[idx];
var padjs  = payload.padj[idx];

var ROT_DEG = 45;
var longestLabel = comparisons.reduce(function (a, b) { return String(b).length > a ? String(b).length : a; }, 0);
var estLabelPx = longestLabel * 6.2;
var padB = Math.max(50, estLabelPx * Math.sin(ROT_DEG * Math.PI / 180) + 24);

var plotH = 260;
var padL = 60, padT = 30;
var basePadR = 20;
var plotW = 700 - padL - basePadR;

var labelOverhangX = estLabelPx * Math.cos(ROT_DEG * Math.PI / 180);
var padR = Math.max(basePadR, labelOverhangX + 15);

var width = padL + plotW + padR;
var height = padT + plotH + padB;
svg.setAttribute('width', width);
svg.setAttribute('height', height);

var finiteVals = values.filter(function (v) { return v !== null && v !== undefined && !isNaN(v); });
var rawMax = finiteVals.length ? Math.max.apply(null, finiteVals.map(Math.abs)) : 0;
var axisMax = niceMax(Math.max(rawMax, 0.001));
var barW = plotW / values.length;
var zeroY = padT + plotH / 2;

var title = document.createElementNS(ns, 'text');
title.setAttribute('x', padL);
title.setAttribute('y', 14);
title.setAttribute('font-size', '12');
title.setAttribute('font-weight', '600');
title.textContent = gene + ' — log2 fold change by comparison';
svg.appendChild(title);

var nTicks = 4;
for (var t = -nTicks / 2; t <= nTicks / 2; t++) {
var tickVal = (axisMax / (nTicks / 2)) * t;
var ty = zeroY - (tickVal / axisMax) * (plotH / 2);

var grid = document.createElementNS(ns, 'line');
grid.setAttribute('x1', padL);
grid.setAttribute('x2', width - padR);
grid.setAttribute('y1', ty);
grid.setAttribute('y2', ty);
grid.setAttribute('stroke', t === 0 ? AXIS_COLOR : GRID_COLOR);
grid.setAttribute('stroke-width', t === 0 ? '1.2' : '1');
svg.appendChild(grid);

var tickLabel = document.createElementNS(ns, 'text');
tickLabel.setAttribute('x', padL - 8);
tickLabel.setAttribute('y', ty + 3);
tickLabel.setAttribute('text-anchor', 'end');
tickLabel.setAttribute('font-size', '10');
tickLabel.setAttribute('fill', TEXT_COLOR);
tickLabel.textContent = Number(tickVal.toFixed(2)).toString();
svg.appendChild(tickLabel);
}

var axisLine = document.createElementNS(ns, 'line');
axisLine.setAttribute('x1', padL);
axisLine.setAttribute('x2', padL);
axisLine.setAttribute('y1', padT);
axisLine.setAttribute('y2', padT + plotH);
axisLine.setAttribute('stroke', AXIS_COLOR);
axisLine.setAttribute('stroke-width', '1');
svg.appendChild(axisLine);

var axisTitle = document.createElementNS(ns, 'text');
axisTitle.setAttribute('x', -(padT + plotH / 2));
axisTitle.setAttribute('y', 14);
axisTitle.setAttribute('text-anchor', 'middle');
axisTitle.setAttribute('font-size', '10');
axisTitle.setAttribute('fill', TEXT_COLOR);
axisTitle.setAttribute('transform', 'rotate(-90)');
axisTitle.textContent = 'log2FoldChange';
svg.appendChild(axisTitle);

values.forEach(function (v, i) {
var x = padL + i * barW + barW * 0.15;
var barWidth = barW * 0.7;

var labelX = x + barWidth / 2;
var labelY = height - padB + 8;
var label = document.createElementNS(ns, 'text');
label.setAttribute('x', labelX);
label.setAttribute('y', labelY);
label.setAttribute('text-anchor', 'start');
label.setAttribute('font-size', '11');
label.setAttribute('transform', 'rotate(' + ROT_DEG + ' ' + labelX + ' ' + labelY + ')');
label.textContent = comparisons[i];
svg.appendChild(label);

if (v === null || v === undefined || isNaN(v)) {
var tick = document.createElementNS(ns, 'line');
tick.setAttribute('x1', x);
tick.setAttribute('x2', x + barWidth);
tick.setAttribute('y1', zeroY);
tick.setAttribute('y2', zeroY);
tick.setAttribute('stroke', NOTTEST_COLOR);
tick.setAttribute('stroke-width', '3');
svg.appendChild(tick);

var ntLabel = document.createElementNS(ns, 'text');
ntLabel.setAttribute('x', labelX);
ntLabel.setAttribute('y', zeroY - 6);
ntLabel.setAttribute('text-anchor', 'middle');
ntLabel.setAttribute('font-size', '9');
ntLabel.setAttribute('fill', NOTTEST_COLOR);
ntLabel.textContent = 'not tested';
svg.appendChild(ntLabel);
return;
}

var barH = Math.abs(v) / axisMax * (plotH / 2);
var y = v >= 0 ? zeroY - barH : zeroY;
var color = v >= 0 ? UP_COLOR : DOWN_COLOR;

var rect = document.createElementNS(ns, 'rect');
rect.setAttribute('x', x);
rect.setAttribute('y', y);
rect.setAttribute('width', barWidth);
rect.setAttribute('height', Math.max(barH, 0.5));
rect.setAttribute('fill', color);
svg.appendChild(rect);

var valLine = document.createElementNS(ns, 'text');
valLine.setAttribute('x', labelX);
valLine.setAttribute('text-anchor', 'middle');
valLine.setAttribute('font-size', '10');
valLine.setAttribute('font-weight', '600');
valLine.setAttribute('fill', color);
valLine.textContent = v.toFixed(2);

var padjLine = document.createElementNS(ns, 'text');
padjLine.setAttribute('x', labelX);
padjLine.setAttribute('text-anchor', 'middle');
padjLine.setAttribute('font-size', '9');
padjLine.setAttribute('fill', TEXT_COLOR);
padjLine.textContent = 'padj=' + fmtPadj(padjs[i]);

if (v >= 0) {
valLine.setAttribute('y', y - 16);
padjLine.setAttribute('y', y - 4);
} else {
valLine.setAttribute('y', y + barH + 12);
padjLine.setAttribute('y', y + barH + 24);
}
svg.appendChild(valLine);
svg.appendChild(padjLine);
});
}

geneInput.addEventListener('input', draw);
})();
</script>
)-"

  html <- gsub("__ID__", widget_id, template, fixed = TRUE)
  sub("__JSON__", json_payload, html, fixed = TRUE)
}
