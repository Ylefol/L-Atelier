# ==============================================================================
# PANDORA - CellTypist annotation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Reference-based cell type annotation via CellTypist, running in an isolated
# Python environment managed by basilisk.  No manual conda/pip setup required.
#
# CellTypist models are downloaded on first use to ~/.celltypist/ (the
# celltypist default cache).  Models persist across sessions; use
# force_update = TRUE to refresh a specific model.
#
# All functions prefixed: PANDORA_
# ==============================================================================


#' List available CellTypist models
#'
#' Fetches the CellTypist model registry and prints a summary of all available
#' models.  Requires an internet connection (model index is hosted by the
#' Sanger Institute).
#'
#' @return Invisibly returns a \code{data.frame} of model information.
#' @export
PANDORA_list_celltypist_models <- function() {

  script <- system.file("python", "celltypist_list_models.py", package = "CAULDRON")

  models_df <- basilisk::basiliskRun(
    env  = .cauldron_celltypist_env,
    fun  = .pandora_celltypist_list_models,
    script = script
  )

  if (!is.data.frame(models_df))
    models_df <- as.data.frame(models_df)

  bar <- strrep("\u2500", 73)

  cat("\u2500\u2500 PANDORA: CellTypist models ", strrep("\u2500", 44), "\n\n", sep = "")

  # Print column names then each row
  cat(sprintf("  %-45s  %s\n", "Model", paste(names(models_df), collapse = "  |  ")))
  cat(bar, "\n", sep = "")
  for (i in seq_len(nrow(models_df))) {
    vals <- paste(as.character(models_df[i, ]), collapse = "  |  ")
    cat(sprintf("  %-45s  %s\n", rownames(models_df)[i], vals))
  }

  cat(bar, "\n\n", sep = "")
  cat("  Pass the model name (with .pkl) as model= in",
      "PANDORA_annotate_celltypist().\n")
  cat(bar, "\n", sep = "")

  invisible(models_df)
}


#' Cell type annotation via CellTypist
#'
#' Annotates cells using a pre-trained CellTypist model.  The model is
#' downloaded on first use to the celltypist cache (\code{~/.celltypist/}).
#' All Python execution is handled by basilisk in an isolated environment —
#' no manual Python setup is required.
#'
#' Raw counts are passed to Python and normalised to 10,000 counts per cell
#' followed by log1p transformation (scanpy's \code{normalize_total} +
#' \code{log1p}) inside the basilisk environment.  This is the exact
#' normalisation CellTypist requires, regardless of which R normalisation was
#' applied to \code{logcounts}.
#'
#' With \code{majority_voting = TRUE} (default), CellTypist over-clusters the
#' data internally and assigns a consensus label per cluster by majority vote,
#' which typically improves coherence.  With \code{FALSE}, each cell receives
#' an independent probabilistic prediction.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"counts"} assay.
#' @param model Character. Model filename, e.g.
#'   \code{"Mouse_Whole_Brain.pkl"}.  Use
#'   \code{\link{PANDORA_list_celltypist_models}} to browse available models.
#' @param majority_voting Logical. Use majority voting within over-clusters
#'   for more coherent labels.  Default \code{TRUE}.
#' @param force_update Logical. Re-download the model even if already cached.
#'   Default \code{FALSE}.
#' @param label_col Character. \code{colData} column for the result.
#'   Default \code{"celltypist_label"}.
#' @param conf_col Character or \code{NULL}.  If non-\code{NULL} and
#'   \code{majority_voting = FALSE}, stores the per-cell confidence score in
#'   this \code{colData} column.  Default \code{"celltypist_conf"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} populated.
#' @export
PANDORA_annotate_celltypist <- function(sce,
                                        model,
                                        majority_voting = TRUE,
                                        force_update    = FALSE,
                                        label_col       = "celltypist_label",
                                        conf_col        = "celltypist_conf",
                                        verbose         = TRUE) {

  if (!"counts" %in% assayNames(sce))
    stop("'counts' assay not found. PANDORA_annotate_celltypist() requires ",
         "raw counts — ensure the SCE was loaded via TALARIA.", call. = FALSE)

  # ── Extract and prepare raw counts ───────────────────────────────────────────
  mat <- counts(sce)
  if (inherits(mat, "IterableMatrix"))
    mat <- as(mat, "dgCMatrix")

  # Transpose to cells x genes (AnnData convention)
  mat_t <- Matrix::t(mat)

  genes <- rownames(sce)
  cells <- colnames(sce)

  if (isTRUE(verbose))
    cat(sprintf(
      "\u2500\u2500 PANDORA: CellTypist %s\n  Model      : %s\n  Cells      : %s  |  Genes: %s\n  Mode       : %s\n  Running Python...\n",
      strrep("\u2500", 34),
      model,
      format(length(cells), big.mark = ","),
      format(length(genes), big.mark = ","),
      if (majority_voting) "majority voting" else "per-cell"
    ))

  # ── Run CellTypist via basilisk ──────────────────────────────────────────────
  script <- system.file("python", "celltypist_run.py", package = "CAULDRON")

  result_df <- basilisk::basiliskRun(
    env             = .cauldron_celltypist_env,
    fun             = .pandora_celltypist_run,
    script          = script,
    mat_t           = mat_t,
    genes           = genes,
    cells           = cells,
    model_name      = model,
    majority_voting = majority_voting,
    force_update    = force_update
  )

  # ── Attach labels to colData ─────────────────────────────────────────────────
  # Align by cell name to guard against any row reordering
  cell_labels <- result_df[cells, "predicted_labels"]
  colData(sce)[[label_col]] <- cell_labels

  if (!majority_voting && !is.null(conf_col) && "conf_score" %in% names(result_df))
    colData(sce)[[conf_col]] <- as.numeric(result_df[cells, "conf_score"])

  # ── Summary ──────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    n_types <- length(unique(na.omit(cell_labels)))
    cat(sprintf(
      "  Cell types : %d unique types\n  Stored as  : colData(sce)[[\"%s\"]]\n%s\n",
      n_types, label_col, strrep("\u2500", 56)
    ))
  }

  sce
}


# ── Internal helpers ───────────────────────────────────────────────────────────
# Named package-level functions so that the CAULDRON namespace is in scope
# (ensures reticulate:: resolves correctly inside basiliskRun).
# Variables are injected into Python's __main__ via py_set_attr rather than
# reticulate::py$  to avoid S3 dispatch issues in the basilisk subprocess.

.pandora_celltypist_run <- function(script, mat_t, genes, cells,
                                    model_name, majority_voting, force_update) {
  main <- reticulate::import("__main__")
  reticulate::py_set_attr(main, "r_mat_t",           mat_t)
  reticulate::py_set_attr(main, "r_genes",            genes)
  reticulate::py_set_attr(main, "r_cells",            cells)
  reticulate::py_set_attr(main, "r_model_name",       model_name)
  reticulate::py_set_attr(main, "r_majority_voting",  majority_voting)
  reticulate::py_set_attr(main, "r_force_update",     force_update)
  reticulate::py_run_file(script)
  reticulate::py_to_r(reticulate::py_get_attr(main, "result_df"))
}

.pandora_celltypist_list_models <- function(script) {
  main <- reticulate::import("__main__")
  reticulate::py_run_file(script)
  reticulate::py_to_r(reticulate::py_get_attr(main, "result_df"))
}
