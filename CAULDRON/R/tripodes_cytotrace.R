# ==============================================================================
# TRIPODES - CytoTRACE v1 Stemness Scoring
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Pure R implementation of the CytoTRACE v1 algorithm (Gulati et al. 2020).
# Higher score = more stem-like / less differentiated.
#
# Algorithm:
#   1. Count expressed genes per cell (gene diversity)
#   2. Correlate each gene's expression with the gene diversity vector (Spearman)
#   3. Select top n_top_genes by absolute correlation
#   4. Raw score per cell = mean expression of those top genes
#   5. Smooth across neighbourhood graph
#   6. Normalise to [0, 1] and invert (high = stem-like)
#
# All public functions prefixed: TRIPODES_
# ==============================================================================


#' Score cells by differentiation state using CytoTRACE v1
#'
#' Implements the CytoTRACE v1 algorithm (Gulati et al. 2020) to assign a
#' stemness score to each cell.  A score of 1 indicates a highly undifferentiated
#' (stem-like) state; 0 indicates a fully differentiated state.
#'
#' The algorithm exploits the empirical observation that less differentiated cells
#' express a greater diversity of genes.  It correlates each gene's expression
#' with per-cell gene diversity, selects the most correlated genes, scores cells
#' by their mean expression of those genes, then smooths across the neighbourhood
#' graph.
#'
#' \strong{BPCells note:} If \code{assay_name} is backed by a BPCells
#' \code{IterableMatrix}, it is materialised to a dense matrix at the start of
#' the function.  CytoTRACE v1 requires full matrix operations (row-wise ranking
#' for Spearman correlation) that cannot be performed on a BPCells backend.
#'
#' \strong{Graph reuse:} By default the SNN graph built by
#' \code{\link{TALOS_build_graph}} and stored in
#' \code{metadata(sce)$snn_graph} is reused for neighbourhood smoothing.  This
#' ensures the CytoTRACE score is smoothed in the same neighbourhood structure
#' used for clustering.  Set \code{recompute_graph = TRUE} only if you have a
#' specific reason to use a different graph — a prominent warning will be issued.
#'
#' @param sce A \code{SingleCellExperiment}.  Should have been through the
#'   standard CAULDRON workflow (normalised, embedded, clustered).
#' @param assay_name Character. Assay to use for scoring. Default
#'   \code{"logcounts"}.
#' @param n_top_genes Integer. Number of top-correlated genes used to compute
#'   the raw score. Default \code{200}.
#' @param n_neighbours Integer. Number of nearest neighbours used for graph
#'   smoothing.  Only used when \code{recompute_graph = TRUE} or when no SNN
#'   graph is found. Default \code{10}.
#' @param use_dimred Character. Reduced dimension to use for kNN construction
#'   when \code{recompute_graph = TRUE} or no graph is available. Default
#'   \code{"PCA"}.
#' @param n_pcs Integer. Number of PCs to use from \code{use_dimred}. Default
#'   \code{30}.
#' @param recompute_graph Logical. If \code{TRUE}, ignore the existing SNN graph
#'   and build a fresh kNN graph for smoothing.  \strong{This is not recommended}
#'   — it breaks consistency with clustering. A prominent warning is always
#'   issued when this is set to \code{TRUE}. Default \code{FALSE}.
#' @param score_col Character. Name of the \code{colData} column to store the
#'   score in. Default \code{"cytotrace_v1"}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return The input \code{sce} with:
#'   \itemize{
#'     \item \code{colData(sce)[[score_col]]}: per-cell CytoTRACE score \eqn{[0, 1]}
#'     \item \code{metadata(sce)$cytotrace_v1_run = TRUE}
#'     \item \code{metadata(sce)$cytotrace_v1_col}: name of the score column
#'   }
#'
#' @references Gulati GS et al. (2020) Single-cell transcriptional diversity is
#'   a hallmark of developmental potential. \emph{Science} 367(6476):405–411.
#'
#' @export
TRIPODES_score_cytotrace_v1 <- function(sce,
                                          assay_name      = "logcounts",
                                          n_top_genes     = 200L,
                                          n_neighbours    = 10L,
                                          use_dimred      = "PCA",
                                          n_pcs           = 30L,
                                          recompute_graph = FALSE,
                                          score_col       = "cytotrace_v1",
                                          verbose         = TRUE) {

  # ── Input checks ─────────────────────────────────────────────────────────────
  if (!assay_name %in% assayNames(sce))
    stop("[TRIPODES] Assay '", assay_name, "' not found. Available: ",
         paste(assayNames(sce), collapse = ", "), ".", call. = FALSE)

  n_top_genes  <- as.integer(n_top_genes)
  n_neighbours <- as.integer(n_neighbours)
  n_pcs        <- as.integer(n_pcs)

  # ── Materialise expression matrix ────────────────────────────────────────────
  expr_raw <- assay(sce, assay_name)

  if (inherits(expr_raw, "IterableMatrix")) {
    message("[TRIPODES] Materialising BPCells-backed '", assay_name,
            "' assay — CytoTRACE v1 requires full matrix operations.")
    expr_mat <- as.matrix(expr_raw)
  } else {
    expr_mat <- as.matrix(expr_raw)
  }
  # expr_mat: genes × cells, dense

  n_genes <- nrow(expr_mat)
  n_cells <- ncol(expr_mat)

  # ── Step 1: gene diversity per cell ──────────────────────────────────────────
  gene_counts <- colSums(expr_mat > 0)

  if (all(gene_counts == 0))
    stop("[TRIPODES] All gene counts are zero — check that '", assay_name,
         "' contains non-zero expression values.", call. = FALSE)

  # ── Step 2: Spearman correlation of each gene vs gene_counts ─────────────────
  if (verbose)
    message("[TRIPODES] Computing gene-count correlations (",
            n_genes, " genes \u00d7 ", n_cells, " cells)...")

  gc_rank     <- rank(gene_counts, ties.method = "average")
  # Row-wise ranking via matrixStats (fast C implementation)
  expr_ranked <- matrixStats::rowRanks(expr_mat, ties.method = "average")

  gc_centered   <- gc_rank - mean(gc_rank)
  expr_centered <- expr_ranked - rowMeans(expr_ranked)

  gc_norm   <- sqrt(sum(gc_centered^2))
  cross     <- as.vector(expr_centered %*% gc_centered)     # n_genes vector
  row_norms <- sqrt(rowSums(expr_centered^2))

  # Guard against zero-variance genes or constant gene_counts
  gene_cors <- ifelse(row_norms > 0 & gc_norm > 0,
                      cross / (row_norms * gc_norm), 0)
  names(gene_cors) <- rownames(sce)

  # ── Step 3: select top n_top_genes by |correlation| ──────────────────────────
  n_top   <- min(n_top_genes, n_genes)
  top_idx <- order(abs(gene_cors), decreasing = TRUE)[seq_len(n_top)]

  if (verbose)
    message("[TRIPODES] Using top ", n_top, " genes (|Spearman r| range: [",
            round(min(abs(gene_cors[top_idx])), 3), ", ",
            round(max(abs(gene_cors[top_idx])), 3), "]).")

  # ── Step 4: raw CytoTRACE score = mean expression of top genes ───────────────
  raw_score <- colMeans(expr_mat[top_idx, , drop = FALSE])

  # ── Step 5: smooth over neighbourhood graph ───────────────────────────────────
  snn_graph   <- metadata(sce)$snn_graph
  has_snn     <- !is.null(snn_graph) &&
                  igraph::vcount(snn_graph) == n_cells

  if (!recompute_graph && has_snn) {
    if (verbose) message("[TRIPODES] Smoothing over existing SNN graph.")
    adj  <- igraph::as_adjacency_matrix(snn_graph, sparse = TRUE)
    deg  <- Matrix::rowSums(adj)
    deg[deg == 0] <- 1L                          # guard: isolated nodes
    norm_adj <- Matrix::Diagonal(x = 1 / deg) %*% adj
    smoothed <- as.vector(norm_adj %*% raw_score)

  } else {
    if (recompute_graph && has_snn) {
      message(
        "\n[TRIPODES] WARNING: recompute_graph = TRUE\n",
        "  A fresh kNN graph will be built instead of reusing the SNN graph\n",
        "  from TALOS_build_graph(). This breaks consistency with your\n",
        "  clustering neighbourhood structure. Use recompute_graph = FALSE\n",
        "  unless you have a specific reason to deviate.\n"
      )
    } else if (!has_snn) {
      if (!is.null(snn_graph))
        warning("[TRIPODES] SNN graph vertex count (", igraph::vcount(snn_graph),
                ") does not match ncol(sce) (", n_cells, "). ",
                "Building fresh kNN graph.", call. = FALSE)
      else
        message("[TRIPODES] No SNN graph in metadata. Building kNN from '",
                use_dimred, "' (k = ", n_neighbours, ").")
    }

    if (!use_dimred %in% reducedDimNames(sce))
      stop("[TRIPODES] Embedding '", use_dimred, "' not found. Available: ",
           paste(reducedDimNames(sce), collapse = ", "), ".", call. = FALSE)

    n_pcs_use <- min(n_pcs, ncol(reducedDim(sce, use_dimred)))
    pca_mat   <- reducedDim(sce, use_dimred)[, seq_len(n_pcs_use), drop = FALSE]
    knn_idx   <- BiocNeighbors::findKNN(pca_mat, k = n_neighbours)$index

    # Vectorised smoothing: each row of knn_idx is the k neighbours of cell i
    smoothed <- rowMeans(matrix(raw_score[knn_idx], nrow = nrow(knn_idx)))
  }

  # ── Step 6: normalise [0, 1] and invert ──────────────────────────────────────
  s_min <- min(smoothed)
  s_max <- max(smoothed)
  ct_score <- if (s_max > s_min) {
    1 - (smoothed - s_min) / (s_max - s_min)
  } else {
    rep(0.5, n_cells)
  }

  # ── Store results ─────────────────────────────────────────────────────────────
  colData(sce)[[score_col]]       <- ct_score
  metadata(sce)$cytotrace_v1_run  <- TRUE
  metadata(sce)$cytotrace_v1_col  <- score_col

  if (verbose) {
    message("[TRIPODES] CytoTRACE v1 complete.")
    message("  Score range: [0, 1]  (1 = most stem-like / undifferentiated)")
    message("  Median score: ", round(median(ct_score), 3))
    message("  Stored in colData(sce)$", score_col)
  }

  sce
}
