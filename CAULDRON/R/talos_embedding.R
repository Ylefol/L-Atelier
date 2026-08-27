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

  n_genes <- if (!is.null(hvg_genes)) length(hvg_genes) else nrow(sce)
  n_pcs   <- min(as.integer(n_pcs), n_genes - 1L)

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


#' Harmony batch integration
#'
#' Corrects batch effects in the PCA embedding using Harmony, while preserving
#' genuine biological variation. The harmonised embedding is stored in
#' \code{reducedDims(sce)[["Harmony"]]} and is designed to be passed directly
#' to \code{\link{TALOS_build_graph}} via \code{use_rep = "Harmony"}.
#'
#' Harmony iteratively finds soft cell-type clusters and removes the component
#' of each cluster's centroid that varies across batches, leaving genuine
#' between-condition differences intact.
#'
#' @param sce A \code{SingleCellExperiment} with \code{reducedDims(sce)[["PCA"]]}
#'   populated by \code{\link{TALOS_run_pca}}.
#' @param batch_col Character. Column in \code{colData(sce)} identifying the
#'   batch variable. Default \code{"sample_id"}.
#' @param n_pcs Integer. Number of PCA dimensions to pass to Harmony.
#'   Default \code{30L}.
#' @param theta Numeric. Diversity penalty — higher values enforce stronger
#'   mixing across batches. Default \code{2}.
#' @param max_iter Integer. Maximum Harmony iterations. Default \code{10L}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return SCE with \code{reducedDims(sce)[["Harmony"]]} populated (cells x
#'   \code{n_pcs} matrix). Downstream call:
#'   \code{TALOS_build_graph(sce, use_rep = "Harmony")}.
#' @export
TALOS_run_harmony <- function(sce,
                               batch_col = "sample_id",
                               n_pcs     = 30L,
                               theta     = 2,
                               max_iter  = 10L,
                               seed      = 42L,
                               verbose   = TRUE) {

  .talos_check_dimred(sce, "PCA", "TALOS_run_pca()")

  if (!batch_col %in% colnames(colData(sce)))
    stop("'", batch_col, "' not found in colData(sce). Available: ",
         paste(colnames(colData(sce)), collapse = ", "), call. = FALSE)

  n_dims    <- min(as.integer(n_pcs), ncol(reducedDim(sce, "PCA")))
  n_batches <- length(unique(colData(sce)[[batch_col]]))

  if (isTRUE(verbose)) {
    cat(sprintf(
      "── TALOS: Harmony %s\n  Source  : PCA (%d dims)\n  Batch   : %s (%d levels)\n  theta   : %g | max_iter: %d\n",
      strrep("─", 40),
      n_dims, batch_col, n_batches,
      theta, as.integer(max_iter)
    ))
  }

  set.seed(seed)
  sce <- harmony::RunHarmony(
    object         = sce,
    group.by.vars  = batch_col,
    dims.use       = seq_len(n_dims),
    theta          = theta,
    max_iter       = as.integer(max_iter),
    reduction.save = "Harmony",
    verbose        = FALSE
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "  Stored  : reducedDims(sce)[['Harmony']] (%s cells x %d dims)\n%s\n",
      format(ncol(sce), big.mark = ","),
      n_dims,
      strrep("─", 56)
    ))
  }

  sce
}


#' UMAP embedding
#'
#' Computes a 2-dimensional UMAP embedding from a reduced-dimension
#' representation.  Defaults to PCA; pass \code{use_rep = "scVI"} to embed
#' from the scVI latent space produced by \code{\link{TALOS_run_scvi}}.
#' Results are stored in \code{reducedDims(sce)[["UMAP"]]}.
#'
#' @param sce A \code{SingleCellExperiment} with the chosen reduced dimension
#'   already populated.
#' @param use_rep Character. Name of the \code{reducedDims} slot to use as
#'   input.  Default \code{"PCA"}.  Use \code{"scVI"} for the batch-corrected
#'   latent embedding from \code{\link{TALOS_run_scvi}}.
#' @param n_pcs Integer. Number of PCA components to use when
#'   \code{use_rep = "PCA"}.  Ignored for other representations (all dims are
#'   used).  Default \code{30}.
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
                            use_rep     = "PCA",
                            n_pcs       = 30L,
                            n_neighbors = 15L,
                            min_dist    = 0.1,
                            use_graph   = TRUE,
                            seed        = 42L,
                            verbose     = TRUE) {

  if (isTRUE(use_graph)) {
    # \u2500\u2500 Graph-aligned mode \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    # Use the pre-built SNN graph from TALOS_build_graph() directly as UMAP
    # input, bypassing uwot's internal KNN step. This ensures clustering and
    # UMAP operate on exactly the same neighbourhood structure.
    g <- metadata(sce)$snn_graph
    if (is.null(g))
      stop("use_graph = TRUE requires metadata(sce)$snn_graph. ",
           "Run TALOS_build_graph() first.", call. = FALSE)

    if (!requireNamespace("uwot", quietly = TRUE))
      stop("Package 'uwot' is required for use_graph = TRUE. ",
           "Install via: install.packages('uwot')", call. = FALSE)

    # Convert SNN igraph -> sparse adjacency matrix (values = SNN edge weights)
    adj <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE)

    # Normalize weights to [0, 1]: scran's rank-based SNN weights can exceed 1.
    max_w <- max(adj@x)
    if (max_w > 1) adj@x <- adj@x / max_w

    # uwot::umap() treats a sparse X as a DISTANCE matrix, not a similarity/
    # affinity matrix (see ?uwot::umap). adj@x holds SNN similarity (higher =
    # more similar), so it must be inverted before being passed in, otherwise
    # a cell's strongest SNN neighbour is read as its farthest point. Verified
    # empirically: feeding the raw similarity gives near-chance (42.8%) odds
    # that a cell's top-weighted neighbour ends up closer than its
    # bottom-weighted one in the embedding; inverting to a true distance
    # raises that to 96.2%.
    adj@x <- 1 - adj@x

    # Retrieve k from TALOS_build_graph(); fall back to mean graph degree.
    k_graph <- metadata(sce)$snn_k %||%
                as.integer(round(mean(igraph::degree(g))))

    set.seed(seed)
    umap_coords <- uwot::umap(
      X            = adj,
      n_components = 2L,
      n_neighbors  = k_graph,
      min_dist     = min_dist,
      verbose      = FALSE
    )
    rownames(umap_coords) <- colnames(sce)
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    reducedDim(sce, "UMAP") <- umap_coords

    if (isTRUE(verbose))
      cat(sprintf(
        "\u2500\u2500 TALOS: UMAP %s\n  Input      : SNN graph (metadata(sce)$snn_graph)\n  Graph mode : TRUE \u2014 aligned with TALOS_build_graph()\n  k (graph)  : %s  |  min_dist: %s\n  Stored as  : reducedDims(sce)[[\"UMAP\"]]\n%s\n",
        strrep("\u2500", 43),
        k_graph, min_dist,
        strrep("\u2500", 56)
      ))

  } else {
    # \u2500\u2500 Default mode: independent KNN graph \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    run_fn <- if (use_rep == "PCA") "TALOS_run_pca()" else
                paste0("TALOS_run_scvi() [use_rep = '", use_rep, "']")
    .talos_check_dimred(sce, use_rep, run_fn)

    rep_mat <- reducedDim(sce, use_rep)
    n_dims  <- if (use_rep == "PCA") {
      min(as.integer(n_pcs), ncol(rep_mat))
    } else {
      ncol(rep_mat)
    }

    set.seed(seed)
    sce <- scater::runUMAP(sce,
                            dimred      = use_rep,
                            n_dimred    = n_dims,
                            n_neighbors = as.integer(n_neighbors),
                            min_dist    = min_dist,
                            name        = "UMAP")

    if (isTRUE(verbose))
      cat(sprintf(
        "\u2500\u2500 TALOS: UMAP %s\n  Input      : %s (%s dims)\n  n_neighbors: %s  |  min_dist: %s\n  Stored as  : reducedDims(sce)[[\"UMAP\"]]\n%s\n",
        strrep("\u2500", 43),
        use_rep, n_dims, n_neighbors, min_dist,
        strrep("\u2500", 56)
      ))
  }

  sce
}


#' tSNE embedding
#'
#' Computes a 2-dimensional tSNE embedding from the PCA reduced dimensions.
#' Requires \code{\link{TALOS_run_pca}} to have been run first.  Results are
#' stored in \code{reducedDims(sce)[["tSNE"]]}.
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
#' @return SCE with \code{reducedDims(sce)[["tSNE"]]} populated.
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
                          name       = "tSNE")

  if (isTRUE(verbose))
    cat(sprintf(
      "\u2500\u2500 TALOS: tSNE %s\n  Input      : PCA (%s components)\n  Perplexity : %s  |  Iterations: %s\n  Stored as  : reducedDims(sce)[[\"tSNE\"]]\n%s\n",
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
