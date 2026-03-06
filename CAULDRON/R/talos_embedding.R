# ==============================================================================
# TALOS - Embeddings
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Core dimensionality reduction: PCA, UMAP, and tSNE.
# Graph construction, clustering, and parameter tuning live in
# talos_graph.R and talos_tuning.R respectively.
#
# All functions prefixed: TALOS_
# ==============================================================================


#' Principal component analysis
#'
#' Runs PCA on the normalised \code{logcounts} assay, restricted to HVGs
#' selected by \code{\link{PYRI_select_hvg}} when available.  Results are
#' stored in \code{reducedDims(sce)[["PCA"]]}.
#'
#' PCA is computed via \code{irlba::irlba} directly for all matrix backends
#' (in-memory \code{dgCMatrix} and BPCells \code{IterableMatrix} alike).
#' This ensures identical numerical behaviour regardless of the backend chosen
#' at load time, making results comparable between the two modes.
#' \code{BiocGenerics::t()} is used for the transpose so that S4 dispatch
#' works correctly for all matrix classes.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"logcounts"} assay.
#' @param n_pcs Integer. Number of principal components to compute.
#'   Default \code{50}.
#' @param use_hvg Logical. If \code{TRUE} (default), restricts to genes in
#'   \code{metadata(sce)$hvg}.  Falls back to all genes with a warning if HVGs
#'   are not found.
#' @param assay_name Character. Assay to use. Default \code{"logcounts"}.
#' @param scale Logical. Scale genes to unit variance before PCA.
#'   Default \code{FALSE}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{reducedDims(sce)[["PCA"]]} populated.
#' @export
TALOS_run_pca <- function(sce,
                           n_pcs      = 50L,
                           use_hvg    = TRUE,
                           assay_name = "logcounts",
                           scale      = FALSE,
                           verbose    = TRUE) {

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  # ── HVG subset ──────────────────────────────────────────────────────────────
  hvg_genes <- NULL
  if (isTRUE(use_hvg)) {
    hvg_genes <- metadata(sce)$hvg
    if (is.null(hvg_genes)) {
      warning("use_hvg = TRUE but no HVGs found in metadata(sce)$hvg. ",
              "Running PCA on all genes. Run PYRI_select_hvg() first.",
              call. = FALSE)
    }
  }

  n_pcs <- min(as.integer(n_pcs), length(hvg_genes) %||% nrow(sce) - 1L)

  # Both BPCells IterableMatrix and in-memory dgCMatrix go through irlba
  # directly so that numerical results are identical regardless of backend.
  # (scater::runPCA wraps the matrix in a BiocSingular wrapper that applies
  # slightly different irlba defaults, producing different singular values.)
  mat <- assay(sce, assay_name)
  if (!is.null(hvg_genes))
    mat <- mat[hvg_genes, , drop = FALSE]
  reducedDim(sce, "PCA") <- .talos_irlba_pca(mat, n_pcs, scale)

  if (isTRUE(verbose)) {
    pct_var  <- attr(reducedDim(sce, "PCA"), "percentVar")
    n_input  <- if (!is.null(hvg_genes)) length(hvg_genes) else nrow(sce)
    top5     <- head(round(pct_var, 1), 5)
    top5_str <- paste0("PC", seq_along(top5), " ", top5, "%", collapse = ", ")
    cat(sprintf(
      "\u2500\u2500 TALOS: PCA %s\n  Input      : %s genes (%s)\n  Components : %s\n  Var. exp.  : %s ...\n%s\n",
      strrep("\u2500", 44),
      format(n_input, big.mark = ","),
      if (!is.null(hvg_genes)) "HVGs" else "all genes",
      n_pcs,
      top5_str,
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' UMAP embedding
#'
#' Computes a 2-dimensional UMAP embedding from the PCA reduced dimensions.
#' Requires \code{\link{TALOS_run_pca}} to have been run first.  Results are
#' stored in \code{reducedDims(sce)[["UMAP"]]}.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"PCA"} reduced
#'   dimension.
#' @param n_pcs Integer. Number of PCA components to use as input.
#'   Default \code{30}.
#' @param n_neighbors Integer. Number of nearest neighbours for the UMAP graph.
#'   Higher values produce a more global view; lower values emphasise local
#'   structure. Default \code{15}.
#' @param min_dist Numeric. Minimum distance between points in the embedding.
#'   Smaller values pack points more tightly. Default \code{0.1}.
#' @param seed Integer. Random seed for reproducibility. Default \code{42L}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{reducedDims(sce)[["UMAP"]]} populated.
#' @export
TALOS_run_umap <- function(sce,
                            n_pcs       = 30L,
                            n_neighbors = 15L,
                            min_dist    = 0.1,
                            seed        = 42L,
                            verbose     = TRUE) {

  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_pcs <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))

  set.seed(seed)
  sce <- scater::runUMAP(sce,
                          dimred      = "PCA",
                          n_dimred    = n_pcs,
                          n_neighbors = as.integer(n_neighbors),
                          min_dist    = min_dist,
                          name        = "UMAP")

  if (isTRUE(verbose))
    cat(sprintf(
      "\u2500\u2500 TALOS: UMAP %s\n  Input      : PCA (%s components)\n  n_neighbors: %s  |  min_dist: %s\n  Stored as  : reducedDims(sce)[[\"UMAP\"]]\n%s\n",
      strrep("\u2500", 43),
      n_pcs, n_neighbors, min_dist,
      strrep("\u2500", 56)
    ))

  sce
}


#' tSNE embedding
#'
#' Computes a 2-dimensional tSNE embedding from the PCA reduced dimensions.
#' Requires \code{\link{TALOS_run_pca}} to have been run first.  Results are
#' stored in \code{reducedDims(sce)[["TSNE"]]}.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"PCA"} reduced
#'   dimension.
#' @param n_pcs Integer. Number of PCA components to use. Default \code{30}.
#' @param perplexity Numeric. tSNE perplexity parameter — roughly the number
#'   of effective nearest neighbours.  Typical range 5–50. Default \code{30}.
#' @param max_iter Integer. Maximum number of tSNE iterations. Default
#'   \code{1000}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{reducedDims(sce)[["TSNE"]]} populated.
#' @export
TALOS_run_tsne <- function(sce,
                            n_pcs      = 30L,
                            perplexity = 30,
                            max_iter   = 1000L,
                            seed       = 42L,
                            verbose    = TRUE) {

  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  n_pcs <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))

  set.seed(seed)
  sce <- scater::runTSNE(sce,
                          dimred     = "PCA",
                          n_dimred   = n_pcs,
                          perplexity = perplexity,
                          max_iter   = as.integer(max_iter),
                          name       = "TSNE")

  if (isTRUE(verbose))
    cat(sprintf(
      "\u2500\u2500 TALOS: tSNE %s\n  Input      : PCA (%s components)\n  Perplexity : %s  |  Iterations: %s\n  Stored as  : reducedDims(sce)[[\"TSNE\"]]\n%s\n",
      strrep("\u2500", 43),
      n_pcs, perplexity, max_iter,
      strrep("\u2500", 56)
    ))

  sce
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Unified irlba PCA for all matrix backends.
#
# mat    : genes x cells — BPCells IterableMatrix or in-memory dgCMatrix/matrix.
#          Must already be subsetted to HVGs (if desired) before calling.
# n_pcs  : number of components to compute.
# scale  : logical; scale genes to unit variance before PCA.
#
# Returns a cells x n_pcs matrix with attributes "percentVar" and "rotation"
# matching the format produced by scater::runPCA, so the verbose block and all
# downstream functions (TALOS_run_umap, TALOS_build_graph) work unchanged.
#
# Both backends go through the same irlba call so results are numerically
# identical regardless of whether BPCells or in-memory storage is used.
.talos_irlba_pca <- function(mat, n_pcs, scale) {
  if (!requireNamespace("irlba", quietly = TRUE))
    stop("Package 'irlba' is required for PCA. ",
         "Install via: install.packages(\"irlba\")", call. = FALSE)

  gene_names <- rownames(mat)
  cell_names <- colnames(mat)
  n_cells    <- ncol(mat)
  n_genes    <- nrow(mat)
  n_pcs      <- min(n_pcs, n_genes - 1L, n_cells - 1L)

  # ── Row statistics — two paths for streaming vs in-memory ───────────────────
  if (inherits(mat, "IterableMatrix")) {
    # BPCells: two streaming passes from disk (one stat per call)
    gene_means <- BPCells::matrix_stats(mat, row_stats = "mean")$row_stats["mean", ]
    gene_vars  <- BPCells::matrix_stats(mat, row_stats = "variance")$row_stats["variance", ]
  } else {
    # In-memory: standard one-pass approach
    gene_means <- rowMeans(mat)
    gene_vars  <- (rowMeans(mat * mat) - gene_means^2) * n_cells / (n_cells - 1L)
  }
  names(gene_means) <- gene_names
  names(gene_vars)  <- gene_names

  scale_vec <- if (isTRUE(scale)) sqrt(pmax(gene_vars, .Machine$double.eps)) else FALSE

  # ── irlba SVD ───────────────────────────────────────────────────────────────
  # BiocGenerics::t() dispatches S4 methods: correct for IterableMatrix
  # (base::t.default() fails) and for dgCMatrix (returns correct dgCMatrix).
  result <- irlba::irlba(BiocGenerics::t(mat), nv = n_pcs,
                          center = gene_means, scale = scale_vec)

  # ── Cell embeddings: U %*% diag(d) — rows = cells, cols = PCs ──────────────
  pca_coords           <- sweep(result$u, 2, result$d, "*")
  rownames(pca_coords) <- cell_names
  colnames(pca_coords) <- paste0("PC", seq_len(n_pcs))

  # ── Percent variance: d_k^2 / (n-1) / total_var * 100 ─────────────────────
  total_var <- sum(gene_vars)
  pct_var   <- (result$d^2 / (n_cells - 1L)) / total_var * 100

  # ── Gene loadings (rotation) ─────────────────────────────────────────────────
  rot           <- result$v
  rownames(rot) <- gene_names
  colnames(rot) <- paste0("PC", seq_len(n_pcs))

  attr(pca_coords, "percentVar") <- pct_var
  attr(pca_coords, "rotation")   <- rot
  pca_coords
}


.talos_check_dimred <- function(sce, name, run_fn) {
  if (!name %in% reducedDimNames(sce))
    stop("'", name, "' not found in reducedDims. Run ", run_fn, " first.",
         call. = FALSE)
}

# NULL-coalescing operator (base R >= 4.4 has it natively; define for safety)
`%||%` <- function(x, y) if (!is.null(x)) x else y
