# ==============================================================================
# TRIPODES - Monocle3 Pseudotime & Trajectory Inference
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Wraps the Monocle3 trajectory inference workflow, re-using pre-computed
# embeddings from TALOS to ensure consistency with the rest of the CAULDRON
# pipeline.  Three root selection modes are supported:
#   - root_cells:        explicit cell barcodes
#   - root_cluster:      all cells from a named TALOS cluster
#   - root_by_cytotrace: top N% cells by CytoTRACE v1 stemness score
#
# Principal graph node coordinates and the graph topology are stored in SCE
# metadata for downstream use by ASPIS_plot_trajectory().
#
# All public functions prefixed: TRIPODES_
# ==============================================================================


#' Pseudotime and trajectory inference via Monocle3
#'
#' Runs the Monocle3 trajectory inference pipeline on a \code{SingleCellExperiment},
#' re-using an existing 2D embedding from TALOS rather than re-computing it.
#' Pseudotime and branch assignments are written back into \code{colData}.
#' The principal graph is stored in \code{metadata} for downstream plotting.
#'
#' \strong{Embedding reuse:} The pre-computed \code{"UMAP"} embedding from
#' TALOS is injected directly into the Monocle3 \code{cell_data_set},
#' bypassing Monocle3's internal normalisation and PCA steps.  This ensures
#' the trajectory is learned on exactly the same layout used for clustering.
#'
#' \strong{UMAP only:} Monocle3's \code{learn_graph()} and
#' \code{cluster_cells()} functions are hardcoded to require a \code{"UMAP"}
#' reducedDim regardless of any user-supplied \code{reduction_method} argument.
#' This is a known, unfixed limitation of the monocle3 package.
#' \code{TRIPODES_run_monocle()} therefore requires a UMAP to be present (run
#' \code{\link{TALOS_run_umap}} first) and will error immediately if it is
#' absent rather than silently failing inside Monocle3.
#'
#' \strong{Graph partitions:} Monocle3 clusters the cells internally on the
#' injected embedding (via \code{cluster_cells()}) to determine graph
#' partitions.  These internal partitions are used only by \code{learn_graph()}
#' and do not replace or override the TALOS cluster labels.  Root selection by
#' \code{root_cluster} refers to your TALOS cluster labels, not Monocle3's
#' internal partitions.
#'
#' \strong{Inf pseudotime:} Cells in partitions that do not contain a root cell
#' will receive \code{Inf} pseudotime from Monocle3.  These are stored as
#' \code{NA} in \code{colData}.
#'
#' @param sce A \code{SingleCellExperiment} that has been through the standard
#'   CAULDRON workflow (normalised, embedded, clustered via TALOS).
#' @param root_cells Character vector of cell barcodes to use as root.
#'   Exactly one of \code{root_cells}, \code{root_cluster}, or
#'   \code{root_by_cytotrace = TRUE} must be provided.
#' @param root_cluster Character. A single cluster label from
#'   \code{colData(sce)[[cluster_col]]}.  All cells in that cluster are used as
#'   the root.
#' @param root_by_cytotrace Logical. Use the top \code{cytotrace_top_pct}\%
#'   cells by CytoTRACE v1 score as root.  Requires
#'   \code{\link{TRIPODES_score_cytotrace_v1}} to have been run. Default
#'   \code{FALSE}.
#' @param cytotrace_top_pct Numeric. Percentage of top-scoring cells to use as
#'   root when \code{root_by_cytotrace = TRUE}. Default \code{5}.
#' @param cytotrace_col Character. \code{colData} column containing CytoTRACE
#'   scores. Default \code{"cytotrace_v1"}.
#' @param cluster_col Character. \code{colData} column containing TALOS cluster
#'   labels.  Used when \code{root_cluster} is specified. Default
#'   \code{"cluster"}.
#' @param use_partition Logical. If \code{TRUE}, Monocle3 learns separate
#'   graphs per partition (recommended for data with disconnected trajectories).
#'   If \code{FALSE}, a single graph is learned across all cells. Default
#'   \code{TRUE}.
#' @param close_loop Logical. Whether Monocle3 should attempt to close loops in
#'   the principal graph.  Set \code{TRUE} only if your biology has cyclic
#'   dynamics. Default \code{FALSE}.
#' @param pseudotime_col Character. Name of the \code{colData} column for
#'   pseudotime values. Default \code{"monocle_pseudotime"}.
#' @param branch_col Character. Name of the \code{colData} column for
#'   Monocle3 partition assignments. Default \code{"monocle_branch"}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return The input \code{sce} with:
#'   \itemize{
#'     \item \code{colData(sce)[[pseudotime_col]]}: per-cell pseudotime;
#'           \code{NA} for cells in partitions with no root.
#'     \item \code{colData(sce)[[branch_col]]}: Monocle3 partition assignment.
#'     \item \code{metadata(sce)$monocle_run = TRUE}
#'     \item \code{metadata(sce)$monocle_graph}: principal graph as an
#'           \code{igraph} object (topology only).
#'     \item \code{metadata(sce)$monocle_graph_nodes}: matrix of principal
#'           graph node coordinates in UMAP space (nodes \eqn{\times} 2).
#'     \item \code{metadata(sce)$monocle_dimred}: always \code{"UMAP"}.
#'   }
#'
#' @export
TRIPODES_run_monocle <- function(sce,
                                  root_cells          = NULL,
                                  root_cluster        = NULL,
                                  root_by_cytotrace   = FALSE,
                                  cytotrace_top_pct   = 5,
                                  cytotrace_col       = "cytotrace_v1",
                                  cluster_col         = "cluster",
                                  use_partition       = TRUE,
                                  close_loop          = FALSE,
                                  pseudotime_col      = "monocle_pseudotime",
                                  branch_col          = "monocle_branch",
                                  verbose             = TRUE) {

  # ── Check monocle3 available ─────────────────────────────────────────────────
  if (!requireNamespace("monocle3", quietly = TRUE))
    stop("[TRIPODES] monocle3 is not installed.\n",
         "  Install with: BiocManager::install('cole-trapnell-lab/monocle3')",
         call. = FALSE)

  # ── Validate root mode — exactly one ─────────────────────────────────────────
  n_root_modes <- sum(!is.null(root_cells), !is.null(root_cluster),
                       isTRUE(root_by_cytotrace))
  if (n_root_modes == 0)
    stop("[TRIPODES] Specify exactly one root mode: root_cells, root_cluster, ",
         "or root_by_cytotrace = TRUE.", call. = FALSE)
  if (n_root_modes > 1)
    stop("[TRIPODES] Only one root mode may be active at a time.", call. = FALSE)

  # ── Validate UMAP presence ───────────────────────────────────────────────────
  # Monocle3's learn_graph() and cluster_cells() are hardcoded to require a
  # "UMAP" reducedDim — this is a known, unfixed limitation of the package.
  # TRIPODES_run_monocle() therefore requires a pre-computed UMAP from TALOS.
  if (!"UMAP" %in% reducedDimNames(sce))
    stop("[TRIPODES] Monocle3 requires a 'UMAP' embedding.\n",
         "  Run TALOS_run_umap() before calling TRIPODES_run_monocle().\n",
         "  Note: Monocle3 does not support non-UMAP embeddings (e.g. tSNE) ",
         "due to hardcoded internal checks in learn_graph() — this is a ",
         "known limitation of the monocle3 package.",
         call. = FALSE)

  # ── Resolve root cell IDs ────────────────────────────────────────────────────
  if (!is.null(root_cells)) {
    missing <- setdiff(root_cells, colnames(sce))
    if (length(missing) > 0)
      stop("[TRIPODES] ", length(missing),
           " root_cells not found in colnames(sce).", call. = FALSE)
    root_ids <- root_cells

  } else if (!is.null(root_cluster)) {
    if (!cluster_col %in% names(colData(sce)))
      stop("[TRIPODES] cluster_col '", cluster_col,
           "' not found in colData(sce).", call. = FALSE)
    root_ids <- colnames(sce)[
      as.character(colData(sce)[[cluster_col]]) == as.character(root_cluster)
    ]
    if (length(root_ids) == 0)
      stop("[TRIPODES] No cells found in cluster '", root_cluster, "'. ",
           "Check cluster labels in colData(sce)$", cluster_col, ".",
           call. = FALSE)
    if (verbose)
      message("[TRIPODES] Root: ", length(root_ids), " cells from cluster '",
              root_cluster, "'.")

  } else {
    # root_by_cytotrace
    if (!isTRUE(metadata(sce)$cytotrace_v1_run))
      stop("[TRIPODES] root_by_cytotrace = TRUE requires ",
           "TRIPODES_score_cytotrace_v1() to be run first.", call. = FALSE)
    if (!cytotrace_col %in% names(colData(sce)))
      stop("[TRIPODES] cytotrace_col '", cytotrace_col,
           "' not found in colData(sce).", call. = FALSE)
    ct     <- colData(sce)[[cytotrace_col]]
    thresh <- stats::quantile(ct, 1 - cytotrace_top_pct / 100,
                               na.rm = TRUE)
    root_ids <- colnames(sce)[!is.na(ct) & ct >= thresh]
    if (length(root_ids) == 0)
      stop("[TRIPODES] No root cells selected at top ", cytotrace_top_pct,
           "% CytoTRACE threshold.", call. = FALSE)
    if (verbose)
      message("[TRIPODES] Root: top ", cytotrace_top_pct,
              "% by CytoTRACE (", length(root_ids), " cells, score >= ",
              round(thresh, 3), ").")
  }

  # ── Build Monocle3 CDS ───────────────────────────────────────────────────────
  if (verbose) message("[TRIPODES] Building Monocle3 cell_data_set...")

  counts_mat <- assay(sce, "counts")
  if (inherits(counts_mat, "IterableMatrix")) {
    message("[TRIPODES] Materialising BPCells-backed counts assay for Monocle3.")
    counts_mat <- as(counts_mat, "dgCMatrix")
  } else if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }

  gene_meta <- as.data.frame(rowData(sce))
  if (!"gene_short_name" %in% colnames(gene_meta))
    gene_meta$gene_short_name <- rownames(sce)

  cds <- monocle3::new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata   = as.data.frame(colData(sce)),
    gene_metadata   = gene_meta
  )

  # Inject pre-computed UMAP — bypasses monocle3's internal preprocess/PCA
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <-
    reducedDim(sce, "UMAP")

  # ── Cluster cells (monocle3-internal, for graph partitions only) ──────────────
  if (verbose)
    message("[TRIPODES] Clustering on 'UMAP' for graph partition assignment...")
  cds <- suppressMessages(
    monocle3::cluster_cells(cds, reduction_method = "UMAP")
  )

  # ── Learn principal graph ─────────────────────────────────────────────────────
  if (verbose) message("[TRIPODES] Learning principal graph...")
  cds <- monocle3::learn_graph(cds,
                                use_partition = use_partition,
                                close_loop    = close_loop,
                                verbose       = FALSE)

  # ── Order cells from root ─────────────────────────────────────────────────────
  if (verbose) message("[TRIPODES] Ordering cells from root...")
  cds <- monocle3::order_cells(cds, root_cells = root_ids)

  # ── Extract results ───────────────────────────────────────────────────────────
  pt   <- monocle3::pseudotime(cds)
  part <- as.character(monocle3::partitions(cds))

  # Inf → NA (cells in partitions with no root)
  pt[is.infinite(pt)] <- NA_real_

  colData(sce)[[pseudotime_col]] <- pt[colnames(sce)]
  colData(sce)[[branch_col]]     <- part[colnames(sce)]

  # Store principal graph info for ASPIS_plot_trajectory.
  pg_graph       <- monocle3::principal_graph(cds)[["UMAP"]]
  pg_node_coords <- t(monocle3::principal_graph_aux(cds)[["UMAP"]]$dp_mst)

  metadata(sce)$monocle_run          <- TRUE
  metadata(sce)$monocle_graph        <- pg_graph
  metadata(sce)$monocle_graph_nodes  <- pg_node_coords
  metadata(sce)$monocle_dimred       <- "UMAP"

  if (verbose) {
    n_ordered <- sum(!is.na(pt))
    message("[TRIPODES] Monocle3 complete.")
    message("  ", n_ordered, " / ", ncol(sce), " cells ordered ",
            "(", ncol(sce) - n_ordered, " in rootless partitions → NA).")
    message("  Pseudotime stored in colData(sce)$", pseudotime_col)
    message("  Branch stored in colData(sce)$", branch_col)
  }

  sce
}
