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
  has_tsne   <- "TSNE" %in% reducedDimNames(sce)
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
