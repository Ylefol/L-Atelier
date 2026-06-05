# ==============================================================================
# TALOS - scVI Latent Embedding
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Variational autoencoder-based dimensionality reduction via scvi-tools.
# Produces a batch-corrected latent space suitable for use as a drop-in
# replacement for PCA in TALOS_run_umap() and TALOS_build_graph().
#
# Python execution is managed by basilisk (.cauldron_scvi_env).
#
# All functions prefixed: TALOS_
# ==============================================================================


#' scVI latent embedding
#'
#' Trains a scVI (single-cell Variational Inference) model on the raw count
#' matrix and extracts the per-cell latent representation.  The result is
#' stored in \code{reducedDims(sce)[["scVI"]]} and can be used as a drop-in
#' replacement for PCA in \code{\link{TALOS_run_umap}} and
#' \code{\link{TALOS_build_graph}} via their \code{use_rep} argument.
#'
#' scVI models counts with a negative binomial distribution and learns a
#' low-dimensional latent space that accounts for library size, technical
#' noise, and — when \code{batch_key} is provided — batch effects.  Unlike
#' PCA, no pre-normalisation is required: the function reads the raw
#' \code{"counts"} assay directly.
#'
#' Training is performed inside an isolated basilisk Python environment
#' (Python 3.11, scvi-tools 1.2.0, PyTorch 2.2).  The environment is
#' created automatically on first use (~1 GB download).
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"counts"} assay
#'   containing raw integer counts.
#' @param n_latent Integer. Dimensionality of the latent space. Default
#'   \code{10}. Higher values capture more variance but may include noise;
#'   10–20 is a typical range for scRNA-seq.
#' @param n_layers Integer. Number of encoder/decoder hidden layers. Default
#'   \code{2}.
#' @param n_hidden Integer. Number of nodes per hidden layer. Default
#'   \code{128}.
#' @param max_epochs Integer. Maximum training epochs. Early stopping is
#'   always enabled and will usually halt well before this limit. Default
#'   \code{400}.
#' @param batch_key Character or \code{NULL}. Column name in
#'   \code{colData(sce)} containing batch labels.  When provided, the scVI
#'   model conditions on batch during training, producing a batch-corrected
#'   latent space.  \code{NULL} (default) trains without batch correction.
#' @param seed Integer. Random seed for reproducibility. Default \code{42L}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with:
#'   \describe{
#'     \item{\code{reducedDims(sce)[["scVI"]]}}{Cells × \code{n_latent}
#'       matrix of latent coordinates.}
#'     \item{\code{metadata(sce)$scvi_params}}{List recording the training
#'       parameters for reproducibility.}
#'   }
#' @export
TALOS_run_scvi <- function(sce,
                            n_latent   = 10L,
                            n_layers   = 2L,
                            n_hidden   = 128L,
                            max_epochs = 400L,
                            batch_key  = NULL,
                            seed       = 42L,
                            verbose    = TRUE) {

  # ── Input validation ─────────────────────────────────────────────────────────
  if (!"counts" %in% assayNames(sce))
    stop("Assay 'counts' not found. TALOS_run_scvi() requires raw integer ",
         "counts — do not use a normalised assay.", call. = FALSE)

  if (!is.null(batch_key)) {
    if (!batch_key %in% names(colData(sce)))
      stop("batch_key '", batch_key, "' not found in colData(sce). ",
           "Available columns: ",
           paste(names(colData(sce)), collapse = ", "), call. = FALSE)
  }

  # ── Extract counts ───────────────────────────────────────────────────────────
  counts_mat <- assay(sce, "counts")
  if (!inherits(counts_mat, "dgCMatrix"))
    counts_mat <- methods::as(counts_mat, "dgCMatrix")

  batch_vec <- if (!is.null(batch_key)) {
    as.character(colData(sce)[[batch_key]])
  } else {
    NULL
  }

  # ── Progress header ──────────────────────────────────────────────────────────
  if (isTRUE(verbose))
    cat(sprintf(
      "── TALOS: scVI %s\n  Cells      : %s\n  Genes      : %s\n  n_latent   : %s  |  n_layers: %s  |  n_hidden: %s\n  max_epochs : %s  |  batch_key: %s\n%s\n",
      strrep("─", 43),
      format(ncol(sce), big.mark = ","),
      format(nrow(sce), big.mark = ","),
      n_latent, n_layers, n_hidden,
      max_epochs,
      if (!is.null(batch_key)) batch_key else "none (no batch correction)",
      strrep("─", 56)
    ))

  # ── Run via basilisk ─────────────────────────────────────────────────────────
  script <- system.file("python", "scvi_run.py", package = "CAULDRON")

  latent <- basilisk::basiliskRun(
    env        = .cauldron_scvi_env,
    fun        = .talos_scvi_run,
    script     = script,
    counts_mat = counts_mat,
    genes      = rownames(sce),
    cells      = colnames(sce),
    batch_vec  = batch_vec,
    n_latent   = as.integer(n_latent),
    n_layers   = as.integer(n_layers),
    n_hidden   = as.integer(n_hidden),
    max_epochs = as.integer(max_epochs),
    seed       = as.integer(seed)
  )

  # ── Store results ────────────────────────────────────────────────────────────
  rownames(latent) <- colnames(sce)
  colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
  reducedDim(sce, "scVI") <- latent

  metadata(sce)$scvi_params <- list(
    n_latent   = as.integer(n_latent),
    n_layers   = as.integer(n_layers),
    n_hidden   = as.integer(n_hidden),
    max_epochs = as.integer(max_epochs),
    batch_key  = batch_key,
    seed       = as.integer(seed)
  )

  if (isTRUE(verbose))
    cat(sprintf(
      "  Stored in  : reducedDims(sce)[[\"scVI\"]] (%s latent dims)\n%s\n",
      n_latent,
      strrep("─", 56)
    ))

  sce
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Package-level function (not anonymous lambda) so that the CAULDRON namespace
# is in scope inside basiliskRun — required for reticulate:: to resolve.
.talos_scvi_run <- function(script, counts_mat, genes, cells, batch_vec,
                              n_latent, n_layers, n_hidden, max_epochs, seed) {
  main <- reticulate::import("__main__")
  reticulate::py_set_attr(main, "r_counts",     counts_mat)
  reticulate::py_set_attr(main, "r_genes",      genes)
  reticulate::py_set_attr(main, "r_cells",      cells)
  reticulate::py_set_attr(main, "r_batch",      batch_vec)   # NULL → Python None
  reticulate::py_set_attr(main, "r_n_latent",   n_latent)
  reticulate::py_set_attr(main, "r_n_layers",   n_layers)
  reticulate::py_set_attr(main, "r_n_hidden",   n_hidden)
  reticulate::py_set_attr(main, "r_max_epochs", max_epochs)
  reticulate::py_set_attr(main, "r_seed",       seed)

  reticulate::py_run_file(script)

  reticulate::py_to_r(reticulate::py_get_attr(main, "result_latent"))
}
