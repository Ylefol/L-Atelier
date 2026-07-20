# ==============================================================================
# KERAUNOS - Differential Expression & Gene Set Enrichment
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The thunderbolts of Zeus, forged by Hephaestus and the Cyclopes in his
# forge. The thunderbolt is the definitive, decisive statement — it strikes
# and the matter is settled. KERAUNOS delivers that same finality through
# statistical testing.
#
# Core responsibilities:
#   - Cluster marker identification
#   - Condition-level group contrasts (cell-level, quick path)
#   - Pseudobulk differential expression (DESeq2)
#   - Gene set enrichment analysis (fgsea / gprofiler2)
#   - Result summarisation and export
#
# All functions prefixed: KERAUNOS_
# ==============================================================================


# ==============================================================================
# Internal helpers
# ==============================================================================

# Extract a tidy top-N table from a scran findMarkers list.
# Detects AUC or logFC columns and computes mean_effect across pairwise comps.
.keraunos_extract_top <- function(markers_list, top_n, fdr_threshold) {
  rows <- lapply(names(markers_list), function(cl) {
    df         <- as.data.frame(markers_list[[cl]])
    df$gene    <- rownames(df)
    df$cluster <- cl

    auc_cols <- grep("^AUC\\.",   names(df), value = TRUE)
    lfc_cols <- grep("^logFC\\.", names(df), value = TRUE)

    if (length(auc_cols) > 0)
      df$mean_effect <- rowMeans(df[, auc_cols, drop = FALSE], na.rm = TRUE)
    else if (length(lfc_cols) > 0)
      df$mean_effect <- rowMeans(df[, lfc_cols, drop = FALSE], na.rm = TRUE)
    else
      df$mean_effect <- NA_real_

    sig <- df[!is.na(df$FDR) & df$FDR < fdr_threshold, , drop = FALSE]
    sig <- sig[order(sig$Top, sig$FDR), , drop = FALSE]
    if (nrow(sig) > top_n) sig <- sig[seq_len(top_n), , drop = FALSE]

    keep <- intersect(c("cluster", "gene", "Top", "p.value", "FDR",
                        "mean_effect"), names(sig))
    sig[, keep, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  if (!is.null(out) && nrow(out) > 0) rownames(out) <- NULL
  out
}


# ==============================================================================
# Cluster marker identification
# ==============================================================================

#' Identify cluster marker genes
#'
#' Finds differentially expressed marker genes for each cluster using pairwise
#' statistical tests via \code{scran::findMarkers()}.  By default runs pairwise
#' Wilcoxon rank-sum tests and returns genes that are significantly upregulated
#' in at least one pairwise comparison (\code{pval_type = "any"}).
#'
#' The result object contains both the full ranked gene list per cluster
#' (\code{$markers}) and a tidy top-N summary table (\code{$top}).
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or other)
#'   assay.
#' @param cluster_col Character. \code{colData} column with cluster labels.
#'   Default \code{"cluster"}.
#' @param assay_name Character. Assay to test. Default \code{"logcounts"}.
#' @param test_type Character. Test to use: \code{"wilcox"} (default), \code{"t"},
#'   or \code{"binom"}.
#' @param pval_type Character. How to combine pairwise p-values:
#'   \code{"any"} (default — gene significant in any comparison),
#'   \code{"some"} (at least \code{min_prop} comparisons), or
#'   \code{"all"} (all comparisons).
#' @param direction Character. Direction of effect: \code{"up"} (default),
#'   \code{"down"}, or \code{"any"}.
#' @param min_prop Numeric. Minimum proportion of comparisons required when
#'   \code{pval_type = "some"}. Default \code{0.5}.
#' @param block_col Character or \code{NULL}. \code{colData} column to use as a
#'   blocking factor (e.g. sample ID).  When provided, each pairwise comparison
#'   is stratified by block, controlling for batch or sample effects.
#'   Default \code{NULL} (no blocking).
#' @param top_n Integer. Maximum number of top markers per cluster in the tidy
#'   summary table. Default \code{10L}.
#' @param fdr_threshold Numeric. FDR cut-off for the tidy summary table.
#'   Default \code{0.05}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{keraunos_markers} object (list) with:
#'   \itemize{
#'     \item \code{$markers} — named list of DataFrames from
#'       \code{scran::findMarkers()}, one per cluster.
#'     \item \code{$top} — tidy data.frame of top markers (columns:
#'       cluster, gene, Top, p.value, FDR, mean_effect).
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_find_markers <- function(sce,
                                   cluster_col   = "cluster",
                                   assay_name    = "logcounts",
                                   test_type     = c("wilcox", "t", "binom"),
                                   pval_type     = c("any", "some", "all"),
                                   direction     = c("up", "down", "any"),
                                   min_prop      = 0.5,
                                   block_col     = NULL,
                                   top_n         = 10L,
                                   fdr_threshold = 0.05,
                                   verbose       = TRUE) {

  test_type <- match.arg(test_type)
  pval_type <- match.arg(pval_type)
  direction <- match.arg(direction)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.  Run PYRI_normalize() first.",
         call. = FALSE)
  if (!cluster_col %in% names(colData(sce)))
    stop("cluster_col '", cluster_col, "' not found in colData(sce).",
         call. = FALSE)
  if (!is.null(block_col) && !block_col %in% names(colData(sce)))
    stop("block_col '", block_col, "' not found in colData(sce).", call. = FALSE)

  clusters   <- colData(sce)[[cluster_col]]
  n_clusters <- length(unique(clusters))
  block      <- if (!is.null(block_col)) colData(sce)[[block_col]] else NULL

  if (isTRUE(verbose))
    cat("Finding markers across ", n_clusters, " clusters (",
            test_type, " test, pval_type=", pval_type,
            if (!is.null(block_col)) paste0(", blocked by ", block_col) else "",
            ")...")

  markers_list <- scran::findMarkers(
    sce,
    groups     = clusters,
    assay.type = assay_name,
    test.type  = test_type,
    pval.type  = pval_type,
    direction  = direction,
    min.prop   = min_prop,
    block      = block
  )

  top_df <- .keraunos_extract_top(markers_list, as.integer(top_n),
                                   fdr_threshold)

  if (isTRUE(verbose)) {
    cl_levels   <- names(markers_list)
    n_sig_per   <- vapply(cl_levels, function(cl) {
      df <- as.data.frame(markers_list[[cl]])
      sum(!is.na(df$FDR) & df$FDR < fdr_threshold)
    }, integer(1))

    cat(sprintf(
      "\u2500\u2500 KERAUNOS: find_markers %s\n  Clusters   : %d\n  FDR < %.2f : %s significant markers (range %d\u2013%d per cluster)\n  Stored as  : keraunos_markers object\n%s\n",
      strrep("\u2500", 33),
      n_clusters,
      fdr_threshold,
      format(sum(n_sig_per), big.mark = ","),
      min(n_sig_per), max(n_sig_per),
      strrep("\u2500", 56)
    ))
  }

  structure(
    list(
      markers = markers_list,
      top     = top_df,
      params  = list(
        cluster_col = cluster_col, assay_name = assay_name,
        test_type = test_type, pval_type = pval_type,
        direction = direction, min_prop = min_prop,
        block_col = block_col,
        top_n = top_n, fdr_threshold = fdr_threshold,
        method = "findMarkers"
      )
    ),
    class = "keraunos_markers"
  )
}

#' @export
print.keraunos_markers <- function(x, ...) {
  cat("keraunos_markers object\n")
  cat("  Method    :", x$params$method %||% "findMarkers", "\n")
  cat("  Clusters  :", length(x$markers), "\n")
  if (identical(x$params$method, "scoreMarkers")) {
    cat("  Top genes :", nrow(x$top), "rows (top", x$params$top_n,
        "per cluster, median AUC >=", x$params$min_auc, ")\n")
  } else {
    cat("  Top genes :", nrow(x$top), "rows (top", x$params$top_n,
        "per cluster at FDR <", x$params$fdr_threshold, ")\n")
  }
  cat("  Access    : $markers (list), $top (data.frame), $params\n")
  invisible(x)
}


# ==============================================================================
# Annotation-focused marker scoring
# ==============================================================================

# Internal helper: extract top-N genes per cluster from scran::scoreMarkers() output.
# restrict_to: optional character vector — only genes in this set are considered.
.keraunos_extract_top_score <- function(score_list, top_n, min_auc,
                                         restrict_to = NULL) {
  rows <- lapply(names(score_list), function(cl) {
    df         <- as.data.frame(score_list[[cl]])
    df$gene    <- rownames(df)
    df$cluster <- cl

    if (!is.null(restrict_to))
      df <- df[df$gene %in% restrict_to, , drop = FALSE]

    if (!is.null(min_auc) && "median.AUC" %in% names(df))
      df <- df[!is.na(df$median.AUC) & df$median.AUC >= min_auc, , drop = FALSE]

    if (nrow(df) == 0L) return(NULL)

    if ("rank.AUC" %in% names(df) && "median.AUC" %in% names(df))
      df <- df[order(df$rank.AUC, -df$median.AUC), , drop = FALSE]

    if (nrow(df) > top_n) df <- df[seq_len(top_n), , drop = FALSE]

    data.frame(
      cluster     = df$cluster,
      gene        = df$gene,
      Top         = if ("rank.AUC"   %in% names(df)) df$rank.AUC   else NA_integer_,
      median_auc  = if ("median.AUC" %in% names(df)) df$median.AUC else NA_real_,
      mean_effect = if ("median.AUC" %in% names(df)) df$median.AUC else NA_real_,
      FDR         = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (!is.null(out) && nrow(out) > 0L) rownames(out) <- NULL
  out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Annotation-focused cluster marker scoring
#'
#' Scores cluster markers using pairwise effect sizes via
#' \code{scran::scoreMarkers()}.  Unlike \code{\link{KERAUNOS_find_markers}},
#' this function does not compute p-values; instead it returns AUC and Cohen's d
#' summaries across all pairwise cluster comparisons.  This is the recommended
#' approach when the goal is cell type annotation rather than hypothesis testing,
#' as p-values from single-cell tests are heavily inflated due to non-independence
#' of cells.
#'
#' Genes are ranked by \code{rank.AUC} (the minimum rank across all pairwise
#' AUC comparisons — lower is better) and filtered by \code{min_auc}
#' (\code{median.AUC >= min_auc} across pairwise comparisons, where 0.5 = random
#' and 1.0 = perfect classifier).
#'
#' The returned object is a \code{keraunos_markers} and is fully compatible with
#' \code{ASPIS_plot_marker_dotplot()} and \code{ASPIS_plot_marker_heatmap()}.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or other)
#'   assay.
#' @param cluster_col Character. \code{colData} column with cluster labels.
#'   Default \code{"cluster"}.
#' @param assay_name Character. Assay to score. Default \code{"logcounts"}.
#' @param block_col Character or \code{NULL}. \code{colData} column to use as a
#'   blocking factor (e.g. sample ID).  Stratifies each pairwise comparison by
#'   block, controlling for batch or sample effects.  Default \code{NULL}.
#' @param restrict_to Character vector or \code{NULL}.  When provided, only genes
#'   present in this set are considered when building the tidy \code{$top} table.
#'   The full \code{$markers} list is unaffected (all genes are scored).  Use
#'   \code{\link{KERAUNOS_fetch_marker_genes}} to obtain a validated whitelist from
#'   MSigDB C8 or PanglaoDB, which filters out uninformative features such as
#'   lncRNAs, ribosomal genes, and ubiquitous housekeeping genes.
#'   Default \code{NULL} (no restriction).
#' @param min_auc Numeric. Minimum \code{median.AUC} for a gene to appear in
#'   the tidy summary table.  AUC of 0.5 = no better than random; 1.0 = perfect
#'   classifier.  Default \code{0.5}.
#' @param top_n Integer. Maximum number of top markers per cluster in the tidy
#'   summary table, ranked by \code{rank.AUC}.  Default \code{10L}.
#' @param full_stats Logical.  Passed through to \code{scran::scoreMarkers()}'s
#'   \code{full.stats} argument.  When \code{TRUE}, each per-cluster DataFrame
#'   in \code{$markers} additionally carries a \code{full.AUC} (and
#'   \code{full.logFC.cohen}, \code{full.logFC.detected}) nested column giving
#'   the per-gene effect size against every *individual* other cluster,
#'   not just the \code{median}/\code{mean}/\code{min}/\code{max} summary
#'   across all of them.  Required by
#'   \code{\link{KERAUNOS_rank_pairwise_markers}} to build a ranked gene list
#'   for one specific pair of clusters (e.g. for GSEA between two clusters
#'   without needing pseudobulk replicate samples).  Increases result object
#'   size roughly in proportion to the number of clusters.  Default
#'   \code{FALSE} (unchanged from previous behaviour).
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{keraunos_markers} object (list) with:
#'   \itemize{
#'     \item \code{$markers} — named list of DataFrames from
#'       \code{scran::scoreMarkers()}, one per cluster.  Each DataFrame contains
#'       per-gene effect size summaries (\code{median.AUC}, \code{rank.AUC},
#'       \code{median.logFC.cohen}, etc.), plus \code{full.*} pairwise columns
#'       when \code{full_stats = TRUE}.
#'     \item \code{$top} — tidy data.frame of top markers (columns:
#'       cluster, gene, Top, median_auc, mean_effect, FDR).
#'       \code{Top = rank.AUC}; \code{FDR = NA} (no p-values).
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_score_markers <- function(sce,
                                    cluster_col  = "cluster",
                                    assay_name   = "logcounts",
                                    block_col    = NULL,
                                    restrict_to  = NULL,
                                    min_auc      = 0.5,
                                    top_n        = 10L,
                                    full_stats   = FALSE,
                                    verbose      = TRUE) {

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.  Run PYRI_normalize() first.",
         call. = FALSE)
  if (!cluster_col %in% names(colData(sce)))
    stop("cluster_col '", cluster_col, "' not found in colData(sce).",
         call. = FALSE)
  if (!is.null(block_col) && !block_col %in% names(colData(sce)))
    stop("block_col '", block_col, "' not found in colData(sce).", call. = FALSE)

  clusters   <- colData(sce)[[cluster_col]]
  n_clusters <- length(unique(clusters))
  block      <- if (!is.null(block_col)) colData(sce)[[block_col]] else NULL

  if (isTRUE(verbose)) {
    cat("[KERAUNOS] Scoring markers across ", n_clusters, " clusters",
        if (!is.null(block_col)) paste0(" (blocked by ", block_col, ")") else "",
        "...\n", sep = "")
    if (!is.null(restrict_to)) {
      n_in_sce <- sum(restrict_to %in% rownames(sce))
      cat("    Restricting $top to known marker genes: ",
          format(length(restrict_to), big.mark = ","), " provided, ",
          format(n_in_sce, big.mark = ","), " present in SCE.\n", sep = "")
    }
  }

  score_list <- scran::scoreMarkers(
    sce,
    groups     = clusters,
    assay.type = assay_name,
    block      = block,
    full.stats = full_stats
  )

  top_df <- .keraunos_extract_top_score(score_list, as.integer(top_n),
                                         min_auc, restrict_to)

  if (isTRUE(verbose)) {
    n_per <- vapply(names(score_list), function(cl) {
      sum(!is.na(top_df$cluster) & top_df$cluster == cl)
    }, integer(1L))

    restrict_line <- if (!is.null(restrict_to))
      paste0("\n  Restricted : known marker genes whitelist (",
             format(length(restrict_to), big.mark = ","), " genes)")
    else ""

    cat(sprintf(
      "── KERAUNOS: score_markers %s\n  Clusters   : %d\n  min AUC    : %.2f%s\n  Top genes  : %s total (range %d–%d per cluster)\n  Stored as  : keraunos_markers object\n%s\n",
      strrep("─", 28),
      n_clusters,
      if (is.null(min_auc)) 0 else min_auc,
      restrict_line,
      format(nrow(top_df), big.mark = ","),
      min(n_per), max(n_per),
      strrep("─", 56)
    ))
  }

  structure(
    list(
      markers = score_list,
      top     = top_df,
      params  = list(
        cluster_col  = cluster_col,
        assay_name   = assay_name,
        block_col    = block_col,
        restrict_to  = restrict_to,
        min_auc      = min_auc,
        top_n        = as.integer(top_n),
        full_stats   = isTRUE(full_stats),
        method       = "scoreMarkers"
      )
    ),
    class = "keraunos_markers"
  )
}


#' Pairwise ranked gene list from marker scores (for pairwise GSEA)
#'
#' Extracts a per-gene ranking statistic for one specific pair of clusters
#' from a \code{\link{KERAUNOS_score_markers}} result run with
#' \code{full_stats = TRUE}.  Unlike \code{$top}/\code{$markers}' default
#' summary columns (\code{median.AUC}, etc., aggregated across *every* other
#' cluster), this pulls the effect size against exactly \code{cluster2},
#' giving a proper pairwise ranking for two clusters without requiring
#' pseudobulk replicate samples — effect sizes are computed at single-cell
#' resolution by \code{scran::scoreMarkers()}. Feed the result directly into
#' \code{\link{KERAUNOS_gsea}}'s \code{ranked_genes} argument.
#'
#' @param markers_result A \code{keraunos_markers} object from
#'   \code{KERAUNOS_score_markers(..., full_stats = TRUE)}.
#' @param cluster1 Character. The focal cluster (positive ranking values mean
#'   higher in this cluster).
#' @param cluster2 Character. The comparison cluster (negative ranking values,
#'   for \code{stat = "AUC"}, mean higher in this cluster).
#' @param stat Character. Which pairwise effect size to use:
#'   \code{"AUC"} (default; scaled cell-detection-rank statistic, shifted by
#'   \code{-0.5} here so 0 = no difference, matching \code{scran}'s convention
#'   that 0.5 = random/no difference — the shift makes the vector suitable
#'   as a signed GSEA ranking statistic), \code{"logFC.cohen"} (Cohen's d,
#'   already signed and centred at 0), or \code{"logFC.detected"} (log-fold
#'   change in detection rate, already signed).
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A named numeric vector (names = gene symbols), NA-removed and
#'   sorted descending — ready to pass to \code{KERAUNOS_gsea(ranked_genes = )}.
#'
#' @export
KERAUNOS_rank_pairwise_markers <- function(markers_result,
                                            cluster1,
                                            cluster2,
                                            stat    = c("AUC", "logFC.cohen", "logFC.detected"),
                                            verbose = TRUE) {

  stat <- match.arg(stat)

  if (!inherits(markers_result, "keraunos_markers"))
    stop("'markers_result' must be a keraunos_markers object from ",
         "KERAUNOS_score_markers().", call. = FALSE)

  if (!isTRUE(markers_result$params$full_stats))
    stop("'markers_result' was scored without full_stats = TRUE, so no ",
         "per-cluster-pair effect sizes are available (only aggregated ",
         "median/mean/min/max across all other clusters).  Re-run ",
         "KERAUNOS_score_markers(sce, ..., full_stats = TRUE).", call. = FALSE)

  cl_names <- names(markers_result$markers)
  for (cl in c(cluster1, cluster2)) {
    if (!cl %in% cl_names)
      stop("Cluster '", cl, "' not found. Available: ",
           paste(cl_names, collapse = ", "), call. = FALSE)
  }

  df       <- markers_result$markers[[cluster1]]
  full_col <- paste0("full.", stat)
  if (!full_col %in% names(df))
    stop("Column '", full_col, "' not found in markers_result$markers[['",
         cluster1, "']].", call. = FALSE)

  pairwise <- df[[full_col]]
  if (!cluster2 %in% names(pairwise))
    stop("Cluster '", cluster2, "' not found among the pairwise comparisons ",
         "for '", cluster1, "'. Available: ",
         paste(names(pairwise), collapse = ", "), call. = FALSE)

  values <- pairwise[[cluster2]]
  names(values) <- rownames(df)

  if (identical(stat, "AUC")) values <- values - 0.5

  values <- values[!is.na(values)]
  values <- sort(values, decreasing = TRUE)

  if (isTRUE(verbose))
    cat("[KERAUNOS] Pairwise ranked genes: ", cluster1, " vs ", cluster2,
        " (", stat, "): ", length(values), " genes\n", sep = "")

  values
}


# ==============================================================================
# Known cell type marker gene databases
# ==============================================================================

#' Fetch known cell type marker genes from a curated database
#'
#' Retrieves the union of all gene symbols annotated as cell type markers across
#' a curated database.  The resulting character vector is intended to be passed
#' to the \code{restrict_to} parameter of \code{\link{KERAUNOS_score_markers}},
#' limiting reported markers to biologically validated genes and filtering out
#' uninformative features such as long non-coding RNAs, ribosomal genes, and
#' ubiquitously expressed housekeeping genes.
#'
#' @section Sources:
#' \describe{
#'   \item{\code{"C8"} (default)}{MSigDB collection C8 — cell type signature
#'     gene sets aggregated from PanglaoDB, CellMarker, DICE, Blueprint, Monaco,
#'     HPCA, and others (~15,000 unique genes in human).  Requires the
#'     \code{msigdbr} package, which is already used by
#'     \code{\link{KERAUNOS_gsea}}.}
#'   \item{\code{"PanglaoDB"}}{PanglaoDB marker database (~8,000 marker
#'     associations across >1,100 cell types in human and mouse).  Reads from a
#'     local TSV file (download from \url{https://panglaodb.se/markers.html}).
#'     Default path: \code{data/PanglaoDB_markers.tsv} relative to the current
#'     working directory.  Override with the \code{path} argument.  Supports
#'     tissue filtering via the \code{tissue} argument.}
#' }
#'
#' @param source Character. Database to query: \code{"C8"} (default) or
#'   \code{"PanglaoDB"}.
#' @param species Character.  Species for gene symbol resolution.
#'   For \code{source = "C8"}: full species name passed to \code{msigdbr}
#'   (e.g. \code{"Homo sapiens"}, \code{"Mus musculus"}).
#'   For \code{source = "PanglaoDB"}: \code{"Hs"} (human, default),
#'   \code{"Mm"} (mouse), or \code{"both"} (no species filter).
#'   Default \code{"Homo sapiens"}.
#' @param tissue Character or \code{NULL}.  For \code{source = "PanglaoDB"}
#'   only: restrict to markers annotated for a specific tissue/organ
#'   (e.g. \code{"Brain"}).  Case-insensitive partial match against the
#'   \code{organ} column.  \code{NULL} (default) returns all tissues.
#' @param path Character or \code{NULL}.  For \code{source = "PanglaoDB"}
#'   only: path to the local PanglaoDB markers TSV file.  When \code{NULL}
#'   (default), looks for \code{data/PanglaoDB_markers.tsv} relative to the
#'   current working directory.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A sorted character vector of unique gene symbols known to be cell
#'   type markers in the chosen database.
#'
#' @examples
#' \dontrun{
#' known_markers <- KERAUNOS_fetch_marker_genes(source = "C8",
#'                                               species = "Homo sapiens")
#' markers <- KERAUNOS_score_markers(sce, restrict_to = known_markers)
#'
#' # PanglaoDB — brain markers only
#' brain_markers <- KERAUNOS_fetch_marker_genes(source = "PanglaoDB",
#'                                               tissue = "Brain")
#' }
#' @export
KERAUNOS_fetch_marker_genes <- function(source  = c("C8", "PanglaoDB"),
                                         species = "Homo sapiens",
                                         tissue  = NULL,
                                         path    = NULL,
                                         verbose = TRUE) {
  source <- match.arg(source)

  genes <- if (source == "C8") {

    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Package 'msigdbr' is required for source = 'C8'.\n",
           "  Install via: install.packages(\"msigdbr\")", call. = FALSE)

    if (isTRUE(verbose))
      cat("[KERAUNOS] Fetching C8 cell type marker genes from MSigDB (",
          species, ")...\n", sep = "")

    gs_df <- tryCatch(
      as.data.frame(msigdbr::msigdbr(species = species, collection = "C8")),
      error = function(e)
        stop("msigdbr failed for C8: ", conditionMessage(e), call. = FALSE)
    )

    if (!"gene_symbol" %in% names(gs_df))
      stop("Unexpected msigdbr output — 'gene_symbol' column not found.",
           call. = FALSE)

    unique(gs_df$gene_symbol)

  } else {

    # Resolve file path
    tsv_path <- if (!is.null(path)) path else
      file.path(getwd(), "data", "PanglaoDB_markers.tsv")

    if (!file.exists(tsv_path))
      stop("PanglaoDB TSV file not found at: ", tsv_path, "\n",
           "  Download from https://panglaodb.se/markers.html and place at\n",
           "  data/PanglaoDB_markers.tsv (or pass path= explicitly).",
           call. = FALSE)

    if (isTRUE(verbose))
      cat("[KERAUNOS] Reading PanglaoDB markers from: ", tsv_path,
          if (!is.null(tissue)) paste0(" (tissue = '", tissue, "')") else "",
          "\n", sep = "")

    db <- tryCatch(
      read.delim(tsv_path, stringsAsFactors = FALSE),
      error = function(e)
        stop("Failed to read PanglaoDB TSV: ", conditionMessage(e), call. = FALSE)
    )

    # Column names after read.delim: spaces become dots.
    # Expected: official.gene.symbol, species, organ
    if (!"official.gene.symbol" %in% names(db))
      stop("Expected column 'official gene symbol' not found in PanglaoDB TSV.\n",
           "  Available columns: ", paste(names(db), collapse = ", "),
           call. = FALSE)

    # Species filter — PanglaoDB encodes as "Hs", "Mm", or "Mm Hs" / "Hs Mm"
    sp_code <- switch(
      tolower(trimws(species)),
      "homo sapiens" = "Hs",
      "mus musculus" = "Mm",
      species  # pass through "Hs" / "Mm" / "both" directly
    )

    if (!identical(sp_code, "both") && "species" %in% names(db))
      db <- db[grepl(sp_code, db$species, fixed = TRUE), , drop = FALSE]

    if (!is.null(tissue) && "organ" %in% names(db))
      db <- db[grepl(tissue, db$organ, ignore.case = TRUE), , drop = FALSE]

    if (nrow(db) == 0L)
      stop("No markers found in PanglaoDB for the given species/tissue filter.",
           call. = FALSE)

    unique(db$official.gene.symbol)
  }

  genes <- sort(unique(genes[!is.na(genes) & nzchar(genes)]))

  if (isTRUE(verbose)) {
    cat(sprintf(
      "── KERAUNOS: fetch_marker_genes %s\n  Source     : %s\n  Species    : %s%s\n  Result     : %s unique gene symbols\n%s\n",
      strrep("─", 22),
      source,
      species,
      if (!is.null(tissue)) paste0("\n  Tissue     : ", tissue) else "",
      format(length(genes), big.mark = ","),
      strrep("─", 56)
    ))
  }

  genes
}


# ==============================================================================
# Cell-level condition contrast
# ==============================================================================

#' Cell-level contrast between groups (quick path)
#'
#' Tests for differentially expressed genes between two (or more) groups at
#' the single-cell level using \code{scran::findMarkers()}.  Can run globally
#' or within each cluster independently.
#'
#' \strong{Important:} This is a cell-level test.  With multiple biological
#' samples per group, individual cells from the same sample are not independent
#' — this is the pseudoreplication problem.  For a statistically rigorous
#' multi-sample comparison, use \code{\link{KERAUNOS_de_pseudobulk}} instead.
#' Use this function for exploratory analysis or single-sample datasets.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or other)
#'   assay.
#' @param group_col Character. \code{colData} column defining the groups to
#'   contrast (must have \eqn{\geq 2} levels).
#' @param cluster_col Character or \code{NULL}.  If \code{NULL} (default), a
#'   single global contrast is performed.  If a \code{colData} column name is
#'   provided, the contrast is run independently within each cluster.
#' @param assay_name Character. Assay to test. Default \code{"logcounts"}.
#' @param test_type Character. Test: \code{"wilcox"} (default), \code{"t"},
#'   or \code{"binom"}.
#' @param direction Character. Direction of effect: \code{"any"} (default),
#'   \code{"up"}, or \code{"down"}.
#' @param fdr_threshold Numeric. FDR cut-off for the tidy summary.
#'   Default \code{0.05}.
#' @param top_n Integer. Maximum top hits per group/cluster in the tidy table.
#'   Default \code{10L}.
#' @param verbose Logical. Print a summary and pseudoreplication warning.
#'   Default \code{TRUE}.
#'
#' @return A \code{keraunos_contrast} object (list) with:
#'   \itemize{
#'     \item \code{$results} — named list of findMarkers DataFrames (one per
#'       group if global; named list \code{[[cluster]][[group]]} if clustered).
#'     \item \code{$top} — tidy data.frame of top hits.
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_contrast_groups <- function(sce,
                                      group_col,
                                      cluster_col   = NULL,
                                      assay_name    = "logcounts",
                                      test_type     = c("wilcox", "t", "binom"),
                                      direction     = c("any", "up", "down"),
                                      fdr_threshold = 0.05,
                                      top_n         = 10L,
                                      verbose       = TRUE) {

  test_type <- match.arg(test_type)
  direction <- match.arg(direction)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.  Run PYRI_normalize() first.",
         call. = FALSE)
  if (!group_col %in% names(colData(sce)))
    stop("group_col '", group_col, "' not found in colData(sce).",
         call. = FALSE)

  groups       <- as.character(colData(sce)[[group_col]])
  group_levels <- sort(unique(stats::na.omit(groups)))

  if (length(group_levels) < 2L)
    stop("group_col '", group_col, "' has fewer than 2 levels: ",
         paste(group_levels, collapse = ", "), call. = FALSE)

  if (isTRUE(verbose))
    cat(
      "KERAUNOS_contrast_groups: cell-level test between groups [",
      paste(group_levels, collapse = ", "), "].\n",
      "  Note: no sample aggregation — pseudoreplication risk for ",
      "multi-sample data.\n",
      "  For rigorous multi-sample DE use KERAUNOS_de_pseudobulk()."
    )

  # Helper: run findMarkers for a group contrast on a sub-SCE
  .run_fm <- function(sub_sce, sub_groups) {
    scran::findMarkers(
      sub_sce,
      groups     = sub_groups,
      assay.type = assay_name,
      test.type  = test_type,
      pval.type  = "any",
      direction  = direction
    )
  }

  # Helper: build tidy top table from a findMarkers list, optionally tagging cl
  .tidy_top <- function(res_list, top_n, fdr_threshold, cluster_tag = NULL) {
    rows <- lapply(names(res_list), function(g) {
      df         <- as.data.frame(res_list[[g]])
      df$gene    <- rownames(df)
      df$group   <- g
      if (!is.null(cluster_tag)) df$cluster <- cluster_tag

      auc_cols <- grep("^AUC\\.",   names(df), value = TRUE)
      lfc_cols <- grep("^logFC\\.", names(df), value = TRUE)
      if (length(auc_cols) > 0)
        df$mean_effect <- rowMeans(df[, auc_cols, drop = FALSE], na.rm = TRUE)
      else if (length(lfc_cols) > 0)
        df$mean_effect <- rowMeans(df[, lfc_cols, drop = FALSE], na.rm = TRUE)
      else
        df$mean_effect <- NA_real_

      sig <- df[!is.na(df$FDR) & df$FDR < fdr_threshold, , drop = FALSE]
      sig <- sig[order(sig$Top, sig$FDR), , drop = FALSE]
      if (nrow(sig) > top_n) sig <- sig[seq_len(top_n), , drop = FALSE]

      keep <- intersect(
        c("cluster", "group", "gene", "Top", "p.value", "FDR", "mean_effect"),
        names(sig)
      )
      sig[, keep, drop = FALSE]
    })
    out <- do.call(rbind, rows)
    if (!is.null(out) && nrow(out) > 0) rownames(out) <- NULL
    out
  }

  if (is.null(cluster_col)) {
    # ── Global contrast ───────────────────────────────────────────────────────
    results <- .run_fm(sce, groups)
    top_df  <- .tidy_top(results, as.integer(top_n), fdr_threshold)

    if (isTRUE(verbose)) {
      cat(sprintf(
        "\u2500\u2500 KERAUNOS: contrast_groups %s\n  Groups     : %s\n  Mode       : global\n  FDR < %.2f : %s significant genes\n%s\n",
        strrep("\u2500", 29),
        paste(group_levels, collapse = " vs "),
        fdr_threshold,
        format(nrow(top_df), big.mark = ","),
        strrep("\u2500", 56)
      ))
    }

  } else {
    # ── Per-cluster contrast ──────────────────────────────────────────────────
    if (!cluster_col %in% names(colData(sce)))
      stop("cluster_col '", cluster_col, "' not found in colData(sce).",
           call. = FALSE)

    clusters       <- as.character(colData(sce)[[cluster_col]])
    cluster_levels <- sort(unique(clusters))
    n_skipped      <- 0L

    results <- stats::setNames(
      lapply(cluster_levels, function(cl) {
        keep     <- clusters == cl
        sub_sce  <- sce[, keep]
        sub_grp  <- groups[keep]

        tbl <- table(sub_grp)
        if (any(tbl < 2L)) {
          if (isTRUE(verbose))
            cat("  Cluster '", cl, "': skipping — a group has < 2 cells.")
          n_skipped <<- n_skipped + 1L
          return(NULL)
        }
        .run_fm(sub_sce, sub_grp)
      }),
      cluster_levels
    )
    results <- Filter(Negate(is.null), results)

    top_list <- lapply(names(results), function(cl)
      .tidy_top(results[[cl]], as.integer(top_n), fdr_threshold,
                cluster_tag = cl))
    top_df <- do.call(rbind, top_list)
    if (is.null(top_df)) top_df <- data.frame()
    if (nrow(top_df) > 0) rownames(top_df) <- NULL

    if (isTRUE(verbose)) {
      cat(sprintf(
        "\u2500\u2500 KERAUNOS: contrast_groups %s\n  Groups     : %s\n  Clusters   : %d tested  |  %d skipped\n  FDR < %.2f : %s significant genes\n%s\n",
        strrep("\u2500", 29),
        paste(group_levels, collapse = " vs "),
        length(results), n_skipped,
        fdr_threshold,
        format(nrow(top_df), big.mark = ","),
        strrep("\u2500", 56)
      ))
    }
  }

  structure(
    list(
      results = results,
      top     = top_df,
      params  = list(
        group_col = group_col, cluster_col = cluster_col,
        assay_name = assay_name, test_type = test_type,
        direction = direction, fdr_threshold = fdr_threshold,
        top_n = top_n
      )
    ),
    class = "keraunos_contrast"
  )
}

#' @export
print.keraunos_contrast <- function(x, ...) {
  cat("keraunos_contrast object\n")
  cat("  Group col  :", x$params$group_col, "\n")
  cat("  Cluster col:",
      ifelse(is.null(x$params$cluster_col), "none", x$params$cluster_col), "\n")
  cat("  Top genes  :", nrow(x$top), "rows\n")
  cat("  Access     : $results (list), $top (data.frame), $params\n")
  invisible(x)
}


# ==============================================================================
# Pseudobulk differential expression
# ==============================================================================

#' Pseudobulk differential expression via DESeq2
#'
#' Aggregates raw counts per sample per cluster (pseudobulk), then runs
#' \code{DESeq2} on each cluster's pseudobulk matrix to identify genes
#' differentially expressed between two conditions.
#'
#' This is the statistically rigorous multi-sample approach: each
#' biological replicate contributes one pseudobulk column, eliminating
#' pseudoreplication.  Requires at least \code{min_samples} replicates per
#' condition per cluster after cell-count filtering.
#'
#' @param sce A \code{SingleCellExperiment} with a raw \code{counts} assay.
#' @param condition_col Character. \code{colData} column with condition labels
#'   (must have exactly 2 levels).
#' @param sample_col Character. \code{colData} column with sample/replicate IDs
#'   used for pseudobulk aggregation (one pseudobulk column per sample).
#' @param cluster_col Character. \code{colData} column with cluster labels.
#'   Default \code{"cluster"}.
#' @param assay_name Character. Assay to aggregate.  Must be raw integer counts.
#'   Default \code{"counts"}.
#' @param reference_level Character or \code{NULL}.  Which condition level is
#'   the reference (denominator) in the DESeq2 contrast.  When \code{NULL}
#'   (default), the alphabetically first level is used.
#' @param min_cells Integer. Minimum cells a sample must contribute to a
#'   cluster to be included as a pseudobulk replicate.  Samples below this
#'   threshold are dropped. Default \code{10L}.
#' @param min_samples Integer. Minimum number of replicates required per
#'   condition per cluster after \code{min_cells} filtering.  Clusters not
#'   meeting this requirement are skipped. Default \code{2L}.
#' @param fdr_threshold Numeric or \code{NULL}. FDR (adjusted p-value)
#'   threshold used to define significant genes in the summary table and for
#'   DESeq2 independent filtering.  Default \code{0.05}.  Set to \code{NULL}
#'   when using \code{pval_threshold} instead.
#' @param pval_threshold Numeric or \code{NULL}. Raw p-value threshold.
#'   Only used when \code{fdr_threshold = NULL}.  Intended for niche cases
#'   where FDR adjustment is not appropriate (e.g. very few tests).  Providing
#'   both \code{fdr_threshold} and \code{pval_threshold} raises an error.
#'   Default \code{NULL}.
#' @param lfc_threshold Numeric. Absolute log2 fold-change threshold passed to
#'   \code{DESeq2::results()}.  Default \code{0} (no LFC filter).
#' @param verbose Logical. Print per-cluster progress and a summary.
#'   Default \code{TRUE}.
#'
#' @return A \code{keraunos_pseudobulk} object (list) with:
#'   \itemize{
#'     \item \code{$results} — named list of data.frames, one per cluster
#'       (columns: gene, cluster, baseMean, log2FoldChange, lfcSE, stat,
#'       pvalue, padj).
#'     \item \code{$dds} — named list of fitted \code{DESeqDataSet} objects
#'       (for downstream use, e.g. \code{DESeq2::lfcShrink()}).
#'     \item \code{$summary} — data.frame: cluster | n_tested | n_sig |
#'       n_up | n_down.
#'     \item \code{$skipped} — character vector of cluster names that were
#'       skipped.
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_de_pseudobulk <- function(sce,
                                    condition_col,
                                    sample_col,
                                    cluster_col     = "cluster",
                                    assay_name      = "counts",
                                    reference_level = NULL,
                                    min_cells       = 10L,
                                    min_samples     = 2L,
                                    fdr_threshold   = 0.05,
                                    pval_threshold  = NULL,
                                    lfc_threshold   = 0,
                                    verbose         = TRUE) {

  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Package 'DESeq2' is required.  ",
         "Install via: BiocManager::install(\"DESeq2\")", call. = FALSE)

  if (!is.null(fdr_threshold) && !is.null(pval_threshold))
    stop("Provide either fdr_threshold or pval_threshold, not both.",
         call. = FALSE)
  if (is.null(fdr_threshold) && is.null(pval_threshold))
    stop("One of fdr_threshold or pval_threshold must be non-NULL.",
         call. = FALSE)

  # Resolve which column and threshold to use for filtering / summary
  use_fdr      <- !is.null(fdr_threshold)
  sig_col      <- if (use_fdr) "padj"   else "pvalue"
  sig_thresh   <- if (use_fdr) fdr_threshold else pval_threshold
  thresh_label <- if (use_fdr) paste0("FDR < ", sig_thresh) else
                               paste0("p < ",   sig_thresh)
  # alpha for DESeq2::results() independent filtering (FDR-based internally);
  # when the user has chosen raw p-value filtering we use the DESeq2 default.
  deseq_alpha  <- if (use_fdr) fdr_threshold else 0.1

  for (col in c(condition_col, sample_col, cluster_col)) {
    if (!col %in% names(colData(sce)))
      stop("Column '", col, "' not found in colData(sce).", call. = FALSE)
  }
  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.  Pseudobulk requires raw counts.",
         call. = FALSE)

  conditions <- as.character(colData(sce)[[condition_col]])
  samples    <- as.character(colData(sce)[[sample_col]])
  clusters   <- as.character(colData(sce)[[cluster_col]])

  cond_levels <- sort(unique(stats::na.omit(conditions)))
  if (length(cond_levels) != 2L)
    stop("condition_col '", condition_col,
         "' must have exactly 2 levels.  Found: ",
         paste(cond_levels, collapse = ", "), call. = FALSE)

  if (!is.null(reference_level)) {
    if (!reference_level %in% cond_levels)
      stop("reference_level '", reference_level,
           "' not found in condition levels: ",
           paste(cond_levels, collapse = ", "), call. = FALSE)
    cond_levels <- c(reference_level, setdiff(cond_levels, reference_level))
  }

  cluster_levels <- sort(unique(clusters))
  is_bpcells     <- inherits(assay(sce, assay_name), "IterableMatrix")

  if (isTRUE(verbose))
    cat("Running pseudobulk DE across ", length(cluster_levels),
            " clusters [", paste(cond_levels, collapse = " vs "), "]...")

  results_list <- list()
  dds_list     <- list()
  skipped      <- character(0)

  for (cl in cluster_levels) {

    cl_idx        <- which(clusters == cl)
    cl_mat        <- assay(sce, assay_name)[, cl_idx, drop = FALSE]
    if (is_bpcells) cl_mat <- as(cl_mat, "dgCMatrix")

    cl_samples    <- samples[cl_idx]
    cl_conditions <- conditions[cl_idx]

    smp_unique  <- sort(unique(cl_samples))
    cell_counts <- table(cl_samples)
    low_samples <- names(cell_counts)[cell_counts < min_cells]
    smp_use     <- setdiff(smp_unique, low_samples)

    if (length(low_samples) > 0 && isTRUE(verbose))
      cat("  Cluster '", cl, "': dropped ", length(low_samples),
              " sample(s) with < ", min_cells, " cells.")

    if (length(smp_use) == 0) {
      if (isTRUE(verbose))
        cat("  Cluster '", cl, "': skipped — no samples remain after ",
                "min_cells filter.")
      skipped <- c(skipped, cl)
      next
    }

    # Condition per remaining sample
    con_per_smp <- vapply(smp_use, function(s)
      cl_conditions[cl_samples == s][1L], character(1))

    cond_counts <- table(con_per_smp)
    has_both <- all(cond_levels %in% names(cond_counts)) &&
      all(cond_counts[cond_levels] >= min_samples)

    if (!has_both) {
      if (isTRUE(verbose))
        cat("  Cluster '", cl, "': skipped — need >= ", min_samples,
                " replicates per condition after filtering.")
      skipped <- c(skipped, cl)
      next
    }

    # Build pseudobulk matrix: genes x samples
    pb_mat <- do.call(cbind, lapply(smp_use, function(s) {
      Matrix::rowSums(cl_mat[, cl_samples == s, drop = FALSE])
    }))
    colnames(pb_mat) <- smp_use
    rownames(pb_mat) <- rownames(sce)

    sample_info <- data.frame(
      condition = factor(con_per_smp, levels = cond_levels),
      row.names = smp_use
    )

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = as.matrix(pb_mat),
      colData   = sample_info,
      design    = ~condition
    )
    # Low-count gene filter
    dds <- dds[rowSums(DESeq2::counts(dds)) >= 10L, ]

    dds <- tryCatch(
      DESeq2::DESeq(dds, quiet = !verbose),
      error = function(e) {
        if (isTRUE(verbose))
          cat("  Cluster '", cl, "': DESeq2 failed — ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(dds)) {
      skipped <- c(skipped, cl)
      next
    }

    res <- DESeq2::results(
      dds,
      contrast     = c("condition", cond_levels[2], cond_levels[1]),
      alpha        = deseq_alpha,
      lfcThreshold = lfc_threshold
    )

    res_df         <- as.data.frame(res)
    res_df$gene    <- rownames(res_df)
    res_df$cluster <- cl
    rownames(res_df) <- NULL

    results_list[[cl]] <- res_df
    dds_list[[cl]]     <- dds
  }

  # Build summary table
  summary_df <- if (length(results_list) > 0) {
    do.call(rbind, lapply(names(results_list), function(cl) {
      df  <- results_list[[cl]]
      sig <- df[!is.na(df[[sig_col]]) & df[[sig_col]] < sig_thresh, ,
                drop = FALSE]
      data.frame(
        cluster  = cl,
        n_tested = nrow(df),
        n_sig    = nrow(sig),
        n_up     = sum(sig$log2FoldChange > 0L, na.rm = TRUE),
        n_down   = sum(sig$log2FoldChange < 0L, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  } else NULL

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 KERAUNOS: de_pseudobulk %s\n  Contrast   : %s vs %s  (%s = reference)\n  Clusters   : %d tested  |  %d skipped\n  Threshold  : %s  |  see $summary\n%s\n",
      strrep("\u2500", 31),
      cond_levels[2], cond_levels[1], cond_levels[1],
      length(results_list), length(skipped),
      thresh_label,
      strrep("\u2500", 56)
    ))
    if (!is.null(summary_df) && nrow(summary_df) > 0)
      print(summary_df, row.names = FALSE)
  }

  structure(
    list(
      results  = results_list,
      dds      = dds_list,
      summary  = summary_df,
      skipped  = skipped,
      params   = list(
        condition_col   = condition_col,
        sample_col      = sample_col,
        cluster_col     = cluster_col,
        assay_name      = assay_name,
        reference_level = cond_levels[1],
        comparison      = paste0(cond_levels[2], " vs ", cond_levels[1]),
        min_cells       = min_cells,
        min_samples     = min_samples,
        fdr_threshold   = fdr_threshold,
        pval_threshold  = pval_threshold,
        sig_col         = sig_col,
        sig_thresh      = sig_thresh,
        lfc_threshold   = lfc_threshold
      )
    ),
    class = "keraunos_pseudobulk"
  )
}

#' @export
print.keraunos_pseudobulk <- function(x, ...) {
  cat("keraunos_pseudobulk object\n")
  cat("  Contrast   :", x$params$comparison, "\n")
  cat("  Clusters   :", length(x$results), "tested,",
      length(x$skipped), "skipped\n")
  cat("  Access     : $results (list), $dds (list), $summary, $params\n")
  invisible(x)
}


# ==============================================================================
# Gene set helpers
# ==============================================================================

#' List available MSigDB gene set collections
#'
#' Wraps \code{msigdbr::msigdbr_collections()} and prints a formatted table
#' grouped by database (human vs mouse).  Use this to discover which collection
#' codes to pass to \code{\link{KERAUNOS_gsea}}.
#'
#' Human collections use plain codes (\code{H}, \code{C2}, \code{C5}, …).
#' Mouse collections (\code{db_species = "MM"}) use \code{M}-prefixed codes
#' (\code{MH}, \code{M2}, \code{M5}, …).  \code{\link{KERAUNOS_gsea}} adjusts
#' these automatically when \code{db_species = "MM"} is set, so you can still
#' pass \code{collection = c("H", "C2")} and the correct codes will be used.
#'
#' @param db_species Character or \code{NULL}.  Filter to a specific database:
#'   \code{"HS"} (human), \code{"MM"} (mouse), or \code{NULL} (default — show
#'   both).
#' @param verbose Logical. Print the formatted table. Default \code{TRUE}.
#'
#' @return Invisibly returns the full collections \code{data.frame} from
#'   \code{msigdbr::msigdbr_collections()}, with an added \code{database}
#'   column (\code{"HS"} or \code{"MM"}).
#' @export
KERAUNOS_list_gene_sets <- function(db_species = NULL, verbose = TRUE) {

  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required.  ",
         "Install via: install.packages(\"msigdbr\")", call. = FALSE)

  db_upper <- if (!is.null(db_species)) toupper(trimws(db_species)) else NULL

  # In newer msigdbr, msigdbr_collections() accepts db_species to return the
  # correct species-specific collection codes.  Try that first; fall back to
  # calling without args (older msigdbr returns all collections in one call).
  cols <- tryCatch(
    as.data.frame(
      if (!is.null(db_upper))
        msigdbr::msigdbr_collections(db_species = db_upper)
      else
        msigdbr::msigdbr_collections()
    ),
    error = function(e) as.data.frame(msigdbr::msigdbr_collections())
  )

  # Normalise column names: msigdbr >= 7.5 uses gs_collection/gs_subcollection;
  # older versions use gs_cat/gs_subcat.
  if (!"gs_collection" %in% names(cols) && "gs_cat" %in% names(cols)) {
    cols$gs_collection    <- cols$gs_cat
    cols$gs_subcollection <- cols$gs_subcat
  }
  if (!"gs_subcollection" %in% names(cols))
    cols$gs_subcollection <- NA_character_

  # When db_species was explicitly requested, tag all rows with that database;
  # otherwise infer from code prefix (mouse codes start with "M").
  if (!is.null(db_upper)) {
    cols$database <- db_upper
  } else {
    cols$database <- ifelse(grepl("^M", cols$gs_collection), "MM", "HS")
  }

  if (nrow(cols) == 0L) {
    cat(sprintf("No collections found for db_species = '%s'.  Valid values: 'HS', 'MM'.\n",
                db_upper))
  }

  if (isTRUE(verbose) && nrow(cols) > 0L) {
    bar <- strrep("\u2500", 60)
    for (db in intersect(c("HS", "MM"), unique(cols$database))) {
      label <- if (db == "MM")
        "Mouse  (db_species = \"MM\")"
      else
        "Human  (db_species = NULL or \"HS\")"
      cat(sprintf("\n\u2500\u2500 MSigDB: %s %s\n", label,
                  strrep("\u2500", max(0L, 44L - nchar(label)))))
      cat(sprintf("  %-12s  %10s  %s\n", "Code", "Gene sets", "Subcategory"))
      cat(bar, "\n")

      sub <- cols[cols$database == db, , drop = FALSE]
      sub <- sub[order(sub$gs_collection, sub$gs_subcollection), , drop = FALSE]
      for (i in seq_len(nrow(sub))) {
        subc <- if (is.na(sub$gs_subcollection[i]) ||
                    !nzchar(sub$gs_subcollection[i]))
          "" else sub$gs_subcollection[i]
        cat(sprintf("  %-12s  %10s  %s\n",
                    sub$gs_collection[i],
                    format(sub$num_genesets[i], big.mark = ","),
                    subc))
      }
    }
    cat(sprintf(
      "\n  Pass the Code to collection= in KERAUNOS_gsea().\n%s\n", bar
    ))
  }

  invisible(cols)
}


#' List species available in MSigDB
#'
#' Wraps \code{msigdbr::msigdbr_species()} and prints a formatted table of
#' available species names and their two-letter \code{db_species} codes.  Use
#' this to find the correct \code{species} and \code{db_species} values for
#' \code{\link{KERAUNOS_gsea}} and \code{\link{KERAUNOS_list_gene_sets}}.
#'
#' @param verbose Logical. Print the formatted table. Default \code{TRUE}.
#'
#' @return Invisibly returns the species \code{data.frame} from
#'   \code{msigdbr::msigdbr_species()}.
#' @export
KERAUNOS_list_species <- function(verbose = TRUE) {

  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required.  ",
         "Install via: install.packages(\"msigdbr\")", call. = FALSE)

  sp <- as.data.frame(msigdbr::msigdbr_species())

  # Detect species name and db_species code columns — case-insensitive search
  # across all names msigdbr has used across versions.
  col_lower  <- tolower(names(sp))
  name_col   <- names(sp)[col_lower %in% c("species_name", "organism_name",
                                            "name", "species")][1L]
  code_col   <- names(sp)[col_lower %in% c("db_species", "species_db",
                                            "code", "database")][1L]

  # If the code column is absent, annotate with known HS/MM codes.
  # Vertebrate MSigDB supports HS (Homo sapiens) and MM (Mus musculus);
  # the two-letter code is derived from the genus + species initials.
  if (is.na(code_col) && !is.na(name_col)) {
    sp$db_species_code <- .keraunos_species_code(sp[[name_col]])
    code_col <- "db_species_code"
  }

  if (isTRUE(verbose)) {
    bar <- strrep("\u2500", 60)
    cat(sprintf("\n\u2500\u2500 MSigDB: available species %s\n",
                strrep("\u2500", 33)))

    if (!is.na(name_col)) {
      cat(sprintf("  %-6s  %s\n", "Code", "Species"))
      cat(bar, "\n")
      sp_ord <- sp[order(sp[[name_col]]), , drop = FALSE]
      for (i in seq_len(nrow(sp_ord))) {
        code <- if (!is.na(code_col) && !is.na(sp_ord[[code_col]][i]))
          sp_ord[[code_col]][i] else ""
        cat(sprintf("  %-6s  %s\n", code, sp_ord[[name_col]][i]))
      }
    } else {
      # Last resort: print the raw data frame so the user can see what's there
      print(sp)
    }

    cat(sprintf(
      "\n  Pass Code to db_species= and the full name to species= in KERAUNOS_gsea().\n%s\n",
      bar
    ))
  }

  invisible(sp)
}

# Internal: derive two-letter MSigDB db_species code from species name.
# Only covers species known to be in MSigDB; returns NA for unknowns.
.keraunos_species_code <- function(species_names) {
  known <- c(
    "Homo sapiens"              = "HS",
    "Mus musculus"              = "MM",
    "Rattus norvegicus"         = "RN",
    "Danio rerio"               = "DR",
    "Drosophila melanogaster"   = "DM",
    "Caenorhabditis elegans"    = "CE",
    "Saccharomyces cerevisiae"  = "SC",
    "Sus scrofa"                = "SS",
    "Bos taurus"                = "BT",
    "Canis lupus familiaris"    = "CF",
    "Macaca mulatta"            = "MQ",
    "Pan troglodytes"           = "PT",
    "Gallus gallus"             = "GG",
    "Xenopus tropicalis"        = "XT"
  )
  unname(known[species_names])
}


# ==============================================================================
# Gene set enrichment analysis (GSEA)
# ==============================================================================

#' Gene set enrichment analysis via fgsea
#'
#' Runs \code{fgsea::fgseaMultilevel()} on a pre-ranked gene list against a
#' collection of gene sets.  Gene sets can be supplied as a named list or
#' fetched automatically from MSigDB via \code{msigdbr}.
#'
#' @param ranked_genes Named numeric vector.  Names are gene symbols (or IDs
#'   matching the gene sets); values are the ranking statistic (e.g. signed
#'   \eqn{-\log_{10}(p)} \eqn{\times} sign(LFC), or the DESeq2 \code{stat}
#'   column).  Must have no missing names.
#' @param gene_sets Named list of character vectors (gene set name → member
#'   genes), or \code{NULL} to fetch from MSigDB via \code{msigdbr}.
#' @param species Character. Species name passed to \code{msigdbr::msigdbr()}
#'   when \code{gene_sets = NULL}.  Default \code{"Homo sapiens"}.
#' @param db_species Character or \code{NULL}.  MSigDB species database to
#'   query.  \code{NULL} (default) lets msigdbr use its default (human-centric,
#'   with orthologs for non-human species).  Use \code{"MM"} to query the
#'   native mouse database (gene sets defined in mouse studies with native
#'   Mus musculus gene symbols, not human orthologs).  Only relevant when
#'   \code{gene_sets = NULL}.
#' @param collection Character vector.  MSigDB collection codes to fetch when
#'   \code{gene_sets = NULL}.  Default \code{c("H", "C2", "C5")}.  Subcategory
#'   can be appended with a colon (e.g. \code{"C2:CP:REACTOME"}).
#' @param min_size Integer. Minimum gene set size (after overlap with ranked
#'   genes).  Default \code{15L}.
#' @param max_size Integer. Maximum gene set size.  Default \code{500L}.
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant}.
#'   Default \code{0.05}.
#' @param verbose Logical. Print progress and summary. Default \code{TRUE}.
#'
#' @return A \code{keraunos_gsea} object (list) with:
#'   \itemize{
#'     \item \code{$results} — full fgsea result data.frame (ordered by padj).
#'     \item \code{$significant} — filtered to \code{fdr_threshold}.
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_gsea <- function(ranked_genes,
                           gene_sets     = NULL,
                           species       = "Homo sapiens",
                           db_species    = NULL,
                           collection    = c("H", "C2", "C5"),
                           min_size      = 15L,
                           max_size      = 500L,
                           fdr_threshold = 0.05,
                           verbose       = TRUE) {

  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Package 'fgsea' is required.  ",
         "Install via: BiocManager::install(\"fgsea\")", call. = FALSE)

  if (is.null(names(ranked_genes)) || length(ranked_genes) == 0L)
    stop("'ranked_genes' must be a named numeric vector (names = gene IDs).",
         call. = FALSE)

  # Remove NA, deduplicate names (keep highest absolute value), sort descending
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  if (anyDuplicated(names(ranked_genes))) {
    ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
    if (isTRUE(verbose))
      cat("  Duplicate gene names removed (kept first occurrence).")
  }
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  # Resolve gene sets
  pathway_to_col <- NULL   # populated below when gene_sets come from msigdbr
  if (is.null(gene_sets)) {
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Package 'msigdbr' is required when gene_sets = NULL.  ",
           "Install via: install.packages(\"msigdbr\")", call. = FALSE)

    # Warn if db_species="MM" is set but collection codes look like human codes.
    # Mouse MSigDB uses "M"-prefixed codes (MH, M2, M5…); human uses plain
    # codes (H, C2, C5…). The user must supply the correct codes explicitly.
    if (!is.null(db_species) && toupper(db_species) == "MM") {
      wrong <- collection[!grepl("^M", collection)]
      if (length(wrong) > 0)
        warning(
          "db_species='MM' (mouse) but collection code(s) do not use the ",
          "mouse prefix: ", paste(wrong, collapse = ", "), ".\n",
          "  Mouse MSigDB collections use uppercase M-prefixed codes ",
          "(e.g. 'MH', 'M2', 'M5').\n",
          "  Run KERAUNOS_list_gene_sets(db_species='MM') to see all available codes.",
          call. = FALSE
        )
    }

    if (isTRUE(verbose)) {
      db_str <- if (!is.null(db_species)) paste0(" [db_species=", db_species, "]") else ""
      cat("Fetching MSigDB gene sets (",
              paste(collection, collapse = " + "),
              ") for ", species, db_str, "...")
    }

    # msigdbr::msigdbr()'s db_species argument requires a character vector --
    # unlike msigdbr_collections(), it errors on NULL rather than falling
    # back to its own default ("HS"), so NULL must be resolved here.
    msigdbr_db_species <- if (is.null(db_species)) "HS" else db_species

    gs_dfs <- lapply(collection, function(col) {
      # subcollection codes can themselves contain a colon (e.g. "GO:BP",
      # "CP:REACTOME"), so everything after the first ":" must be rejoined
      # rather than truncated to the second token only -- taking just
      # cat_sub[2] silently fetched the wrong (broader) gene set for
      # "C2:CP:REACTOME" (all of C2:CP, not just REACTOME) and an invalid
      # subcollection for "C5:GO:BP" ("GO" instead of "GO:BP").
      cat_sub <- strsplit(col, ":", fixed = TRUE)[[1]]
      cat  <- cat_sub[1]
      subc <- if (length(cat_sub) > 1) paste(cat_sub[-1], collapse = ":") else NULL

      tryCatch(
        msigdbr::msigdbr(species = species, db_species = msigdbr_db_species,
                         collection = cat, subcollection = subc),
        error = function(e) {
          warning("msigdbr failed for collection '", col, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    })
    gs_dfs <- Filter(Negate(is.null), gs_dfs)

    if (length(gs_dfs) == 0L)
      stop("No gene sets retrieved from MSigDB.  ",
           "Check species name and collection codes.", call. = FALSE)

    gs_df     <- do.call(rbind, gs_dfs)
    gene_sets <- split(gs_df$gene_symbol, gs_df$gs_name)

    # Build pathway → collection map for annotating fgsea results later
    unique_gs      <- gs_df[!duplicated(gs_df$gs_name), ]
    pathway_to_col <- setNames(unique_gs$gs_collection, unique_gs$gs_name)
  }

  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("'gene_sets' must be a named list of character vectors.", call. = FALSE)

  if (isTRUE(verbose))
    cat("Running fgsea on ", length(gene_sets), " gene sets, ",
            format(length(ranked_genes), big.mark = ","), " ranked genes...")

  res <- fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats    = ranked_genes,
    minSize  = as.integer(min_size),
    maxSize  = as.integer(max_size),
    eps      = 0
  )

  # Convert to data.frame; collapse leadingEdge list column
  res_df <- as.data.frame(res)
  res_df$leadingEdge <- vapply(
    res$leadingEdge,
    function(x) paste(x, collapse = ", "),
    character(1)
  )
  res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
  rownames(res_df) <- NULL

  # Annotate with gene set collection when fetched from msigdbr
  if (!is.null(pathway_to_col))
    res_df$collection <- pathway_to_col[res_df$pathway]

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < fdr_threshold, ,
                drop = FALSE]

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 KERAUNOS: gsea %s\n  Gene sets  : %d tested\n  FDR < %.2f : %d significant  (%d enriched, %d depleted)\n%s\n",
      strrep("\u2500", 39),
      nrow(res_df),
      fdr_threshold,
      nrow(sig),
      sum(sig$NES > 0, na.rm = TRUE),
      sum(sig$NES < 0, na.rm = TRUE),
      strrep("\u2500", 56)
    ))
  }

  # Keep only tested gene sets (those that passed min/max size and appear in
  # results) — needed for enrichment plots without re-fetching from MSigDB.
  tested_sets <- gene_sets[names(gene_sets) %in% res_df$pathway]

  structure(
    list(
      results      = res_df,
      significant  = sig,
      ranked_genes = ranked_genes,
      gene_sets    = tested_sets,
      params       = list(
        species = species, db_species = db_species,
        collection = collection,
        min_size = min_size, max_size = max_size,
        n_genes = length(ranked_genes),
        n_sets  = nrow(res_df),
        fdr_threshold = fdr_threshold
      )
    ),
    class = "keraunos_gsea"
  )
}

#' @export
print.keraunos_gsea <- function(x, ...) {
  cat("keraunos_gsea object\n")
  cat("  Gene sets :", x$params$n_sets, "tested\n")
  cat("  Significant:", nrow(x$significant),
      paste0("(FDR < ", x$params$fdr_threshold, ")\n"))
  cat("  Access    : $results, $significant, $params\n")
  invisible(x)
}


# ==============================================================================
# Over-representation analysis (ORA)
# ==============================================================================

#' Over-representation analysis via gprofiler2
#'
#' Tests whether a set of genes is enriched for annotated biological terms
#' using \code{gprofiler2::gost()}.  Supports GO terms, KEGG, Reactome,
#' and other databases available through g:Profiler.
#'
#' @param genes Character vector of gene symbols (or Ensembl IDs, Entrez IDs,
#'   etc. — whatever gprofiler2 accepts for the chosen organism).
#' @param background Character vector of background genes, or \code{NULL}
#'   (default) to use the gprofiler2 reference genome as background.
#' @param organism Character. gprofiler2 organism code.
#'   Default \code{"hsapiens"}.  Use \code{"mmusculus"} for mouse, etc.
#'   See \code{gprofiler2::gost()} for full list.
#' @param sources Character vector. Databases to query.
#'   Default \code{c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")}.
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant}.
#'   Default \code{0.05}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{keraunos_ora} object (list) with:
#'   \itemize{
#'     \item \code{$results} — full gprofiler2 result data.frame (ordered by
#'       p-value).
#'     \item \code{$significant} — filtered to \code{fdr_threshold}.
#'     \item \code{$params} — list of parameters used.
#'   }
#' @export
KERAUNOS_ora <- function(genes,
                          background    = NULL,
                          organism      = "hsapiens",
                          sources       = c("GO:BP", "GO:MF", "GO:CC",
                                            "KEGG", "REAC"),
                          fdr_threshold = 0.05,
                          verbose       = TRUE) {

  if (!requireNamespace("gprofiler2", quietly = TRUE))
    stop("Package 'gprofiler2' is required.  ",
         "Install via: install.packages(\"gprofiler2\")", call. = FALSE)

  if (length(genes) == 0L)
    stop("'genes' is empty.", call. = FALSE)

  if (isTRUE(verbose))
    cat("Running ORA via gprofiler2 for ",
            format(length(genes), big.mark = ","), " genes...")

  gost_res <- gprofiler2::gost(
    query             = genes,
    organism          = organism,
    custom_bg         = background,
    sources           = sources,
    correction_method = "fdr",
    user_threshold    = fdr_threshold,
    significant       = FALSE   # return all; filter ourselves
  )

  if (is.null(gost_res) || is.null(gost_res$result) ||
      nrow(gost_res$result) == 0L) {
    if (isTRUE(verbose))
      cat("  No results returned from gprofiler2.")
    res_df <- data.frame()
    sig    <- data.frame()

  } else {
    res_df <- gost_res$result

    # Collapse list columns (parents field)
    list_cols <- vapply(res_df, is.list, logical(1))
    res_df[list_cols] <- lapply(res_df[list_cols], function(col)
      vapply(col, function(x) paste(x, collapse = ","), character(1)))

    res_df <- res_df[order(res_df$p_value), , drop = FALSE]
    rownames(res_df) <- NULL
    sig <- res_df[!is.na(res_df$p_value) & res_df$p_value < fdr_threshold, ,
                  drop = FALSE]
  }

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 KERAUNOS: ora %s\n  Query      : %s genes\n  Organism   : %s\n  Sources    : %s\n  FDR < %.2f : %d significant terms\n%s\n",
      strrep("\u2500", 42),
      format(length(genes), big.mark = ","),
      organism,
      paste(sources, collapse = ", "),
      fdr_threshold,
      nrow(sig),
      strrep("\u2500", 56)
    ))
  }

  structure(
    list(
      results     = res_df,
      significant = sig,
      params      = list(
        organism       = organism,
        sources        = sources,
        n_genes        = length(genes),
        has_background = !is.null(background),
        fdr_threshold  = fdr_threshold
      )
    ),
    class = "keraunos_ora"
  )
}

#' @export
print.keraunos_ora <- function(x, ...) {
  cat("keraunos_ora object\n")
  cat("  Organism  :", x$params$organism, "\n")
  cat("  Genes     :", x$params$n_genes, "\n")
  cat("  Significant:", nrow(x$significant),
      paste0("(FDR < ", x$params$fdr_threshold, ")\n"))
  cat("  Access    : $results, $significant, $params\n")
  invisible(x)
}


# ==============================================================================
# Result summarisation & export
# ==============================================================================

#' Summarise differential expression results
#'
#' Produces a concise table (n_tested, n_sig, n_up, n_down per cluster/group)
#' from any KERAUNOS DE result object.
#'
#' @param de_result A \code{keraunos_markers}, \code{keraunos_contrast}, or
#'   \code{keraunos_pseudobulk} object.
#' @param fdr_threshold Numeric. FDR cut-off applied to the results.
#'   Default \code{0.05}.
#' @param lfc_threshold Numeric. Minimum absolute log2 fold-change required
#'   (only applied to pseudobulk results, which have an LFC column).
#'   Default \code{0}.
#'
#' @return A data.frame with one row per cluster/group.
#' @export
KERAUNOS_summarise_de <- function(de_result,
                                   fdr_threshold = 0.05,
                                   lfc_threshold = 0) {

  if (inherits(de_result, "keraunos_markers")) {
    do.call(rbind, lapply(names(de_result$markers), function(cl) {
      df  <- as.data.frame(de_result$markers[[cl]])
      sig <- df[!is.na(df$FDR) & df$FDR < fdr_threshold, , drop = FALSE]
      data.frame(cluster = cl, n_tested = nrow(df), n_sig = nrow(sig),
                 stringsAsFactors = FALSE)
    }))

  } else if (inherits(de_result, "keraunos_contrast")) {

    if (is.null(de_result$params$cluster_col)) {
      # Global: one row per group in the findMarkers list
      do.call(rbind, lapply(names(de_result$results), function(g) {
        df  <- as.data.frame(de_result$results[[g]])
        sig <- df[!is.na(df$FDR) & df$FDR < fdr_threshold, , drop = FALSE]
        data.frame(group = g, n_tested = nrow(df), n_sig = nrow(sig),
                   stringsAsFactors = FALSE)
      }))
    } else {
      # Per-cluster: iterate clusters then groups
      do.call(rbind, lapply(names(de_result$results), function(cl) {
        do.call(rbind, lapply(names(de_result$results[[cl]]), function(g) {
          df  <- as.data.frame(de_result$results[[cl]][[g]])
          sig <- df[!is.na(df$FDR) & df$FDR < fdr_threshold, , drop = FALSE]
          data.frame(cluster = cl, group = g,
                     n_tested = nrow(df), n_sig = nrow(sig),
                     stringsAsFactors = FALSE)
        }))
      }))
    }

  } else if (inherits(de_result, "keraunos_pseudobulk")) {
    # Respect whichever threshold was used when the object was created
    pb_sig_col    <- de_result$params$sig_col
    pb_sig_thresh <- de_result$params$sig_thresh
    do.call(rbind, lapply(names(de_result$results), function(cl) {
      df  <- de_result$results[[cl]]
      sig <- df[!is.na(df[[pb_sig_col]]) &
                  df[[pb_sig_col]] < pb_sig_thresh &
                  abs(df$log2FoldChange) > lfc_threshold, , drop = FALSE]
      data.frame(
        cluster  = cl,
        n_tested = nrow(df),
        n_sig    = nrow(sig),
        n_up     = sum(sig$log2FoldChange >  lfc_threshold, na.rm = TRUE),
        n_down   = sum(sig$log2FoldChange < -lfc_threshold, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))

  } else {
    stop("'de_result' must be a keraunos_markers, keraunos_contrast, or ",
         "keraunos_pseudobulk object.", call. = FALSE)
  }
}


# NOTE: Export functions for KERAUNOS results live in TALARIA (talaria_export.R)
#       in line with the CAULDRON architecture where TALARIA owns all I/O.
#       See: TALARIA_export_de(), TALARIA_export_gsea(), TALARIA_export_ora().


# ==============================================================================
# Multi-cluster GSEA on pseudobulk DE results
# ==============================================================================

#' Run GSEA across all clusters from pseudobulk DE results
#'
#' A convenience wrapper around \code{\link{KERAUNOS_gsea}} that accepts a
#' \code{keraunos_pseudobulk} object and runs GSEA for every cluster
#' automatically.  Gene sets are fetched \strong{once} from MSigDB (not per
#' cluster), making multi-cluster runs efficient.
#'
#' The ranking statistic is derived from the DESeq2 result columns:
#' \describe{
#'   \item{\code{"stat"}}{DESeq2 Wald test statistic (default).  Most
#'     informative ranking; captures both magnitude and significance.}
#'   \item{\code{"lfc_pval"}}{\eqn{\mathrm{sign}(\log_2 FC) \times
#'     -\log_{10}(p\text{-value})}.  Useful when \code{stat} has many
#'     \code{NA}s (e.g. genes filtered by DESeq2 independent filtering).}
#' }
#' Genes with \code{NA} ranking values are removed before the GSEA call.
#' Clusters with fewer than \code{min_ranked_genes} valid ranking values are
#' skipped with a message.
#'
#' @param de_result A \code{keraunos_pseudobulk} object (output of
#'   \code{\link{KERAUNOS_de_pseudobulk}}).
#' @param gene_sets Named list of character vectors, or \code{NULL} (default)
#'   to fetch from MSigDB via \code{msigdbr}.  Fetched once and reused across
#'   all clusters.
#' @param species Character. Species for MSigDB gene set retrieval.
#'   Default \code{"Homo sapiens"}.
#' @param db_species Character or \code{NULL}.  MSigDB database to query.
#'   Use \code{"MM"} for native mouse gene sets (requires M-prefixed
#'   collection codes, e.g. \code{c("MH", "M2", "M5")}).  See
#'   \code{\link{KERAUNOS_list_gene_sets}} to browse available codes.
#' @param collection Character vector. MSigDB collection codes.
#'   Default \code{c("H", "C2", "C5")}.  Subcategory can be appended with a
#'   colon, e.g. \code{"C2:CP:REACTOME"}.
#' @param rank_by Character. Ranking statistic to derive from DESeq2 results.
#'   \code{"stat"} (default) or \code{"lfc_pval"}.
#' @param min_ranked_genes Integer. Minimum number of non-\code{NA} ranked
#'   genes required to run GSEA for a cluster.  Clusters below this threshold
#'   are skipped.  Default \code{10L}.
#' @param min_size Integer. Minimum gene set size after overlap with ranked
#'   genes.  Default \code{15L}.
#' @param max_size Integer. Maximum gene set size.  Default \code{500L}.
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant}.
#'   Default \code{0.05}.
#' @param verbose Logical. Print per-cluster progress and a summary table.
#'   Default \code{TRUE}.
#'
#' @return A \code{keraunos_gsea_multi} object (list) with:
#'   \itemize{
#'     \item \code{$results}     — named list of \code{keraunos_gsea} objects,
#'       one per cluster that was run.
#'     \item \code{$significant} — named list of significant-only
#'       data.frames, one per cluster.
#'     \item \code{$summary}     — data.frame:
#'       cluster / n_ranked / n_sets / n_sig / n_enriched / n_depleted.
#'     \item \code{$skipped}     — character vector of cluster names skipped.
#'     \item \code{$params}      — parameters used.
#'   }
#' @export
KERAUNOS_gsea_pseudobulk <- function(de_result,
                                      gene_sets        = NULL,
                                      species          = "Homo sapiens",
                                      db_species       = NULL,
                                      collection       = c("H", "C2", "C5"),
                                      rank_by          = c("stat", "lfc_pval"),
                                      min_ranked_genes = 10L,
                                      min_size         = 15L,
                                      max_size         = 500L,
                                      fdr_threshold    = 0.05,
                                      verbose          = TRUE) {

  if (!inherits(de_result, "keraunos_pseudobulk"))
    stop("'de_result' must be a keraunos_pseudobulk object (output of ",
         "KERAUNOS_de_pseudobulk()).", call. = FALSE)

  if (length(de_result$results) == 0L)
    stop("No cluster results found in 'de_result'.", call. = FALSE)

  rank_by <- match.arg(rank_by)

  # ── Fetch gene sets once ────────────────────────────────────────────────────
  if (is.null(gene_sets)) {
    if (!requireNamespace("fgsea", quietly = TRUE))
      stop("Package 'fgsea' is required.  ",
           "Install via: BiocManager::install(\"fgsea\")", call. = FALSE)
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Package 'msigdbr' is required when gene_sets = NULL.  ",
           "Install via: install.packages(\"msigdbr\")", call. = FALSE)

    # Warn if db_species="MM" but collection codes look like human codes
    if (!is.null(db_species) && toupper(db_species) == "MM") {
      wrong <- collection[!grepl("^M", collection)]
      if (length(wrong) > 0)
        warning(
          "db_species='MM' (mouse) but collection code(s) do not use the ",
          "mouse prefix: ", paste(wrong, collapse = ", "), ".\n",
          "  Mouse MSigDB collections use uppercase M-prefixed codes ",
          "(e.g. 'MH', 'M2', 'M5').\n",
          "  Run KERAUNOS_list_gene_sets(db_species='MM') to see all available codes.",
          call. = FALSE
        )
    }

    if (isTRUE(verbose)) {
      db_str <- if (!is.null(db_species))
        paste0(" [db_species=", db_species, "]") else ""
      cat("Fetching MSigDB gene sets (", paste(collection, collapse = " + "),
              ") for ", species, db_str, "...")
    }

    # msigdbr::msigdbr()'s db_species argument requires a character vector --
    # unlike msigdbr_collections(), it errors on NULL rather than falling
    # back to its own default ("HS"), so NULL must be resolved here.
    msigdbr_db_species <- if (is.null(db_species)) "HS" else db_species

    gs_dfs <- lapply(collection, function(col) {
      # subcollection codes can themselves contain a colon (e.g. "GO:BP",
      # "CP:REACTOME") -- see the matching comment in KERAUNOS_gsea() above
      # for why cat_sub[2] alone silently truncates multi-colon codes.
      cat_sub <- strsplit(col, ":", fixed = TRUE)[[1]]
      cat  <- cat_sub[1]
      subc <- if (length(cat_sub) > 1) paste(cat_sub[-1], collapse = ":") else NULL
      tryCatch(
        msigdbr::msigdbr(species = species, db_species = msigdbr_db_species,
                         collection = cat, subcollection = subc),
        error = function(e) {
          warning("msigdbr failed for collection '", col, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    })
    gs_dfs <- Filter(Negate(is.null), gs_dfs)
    if (length(gs_dfs) == 0L)
      stop("No gene sets retrieved from MSigDB.  ",
           "Check species name and collection codes.", call. = FALSE)

    gs_df     <- do.call(rbind, gs_dfs)
    gene_sets <- split(gs_df$gene_symbol, gs_df$gs_name)
    if (isTRUE(verbose))
      cat("  Retrieved ", format(length(gene_sets), big.mark = ","),
              " gene sets.")
  }

  # ── Per-cluster GSEA ────────────────────────────────────────────────────────
  cluster_names <- names(de_result$results)
  results_list  <- list()
  sig_list      <- list()
  skipped       <- character(0)

  summary_rows <- lapply(cluster_names, function(cl) {

    df <- de_result$results[[cl]]

    # Build ranked vector from DESeq2 columns
    ranked <- if (rank_by == "stat") {
      stats::setNames(df$stat, df$gene)
    } else {
      stats::setNames(
        sign(df$log2FoldChange) * (-log10(df$pvalue)),
        df$gene
      )
    }
    ranked <- ranked[!is.na(ranked)]

    if (length(ranked) < as.integer(min_ranked_genes)) {
      if (isTRUE(verbose))
        cat("  Cluster ", cl, ": skipped (",
                length(ranked), " valid ranked genes < min_ranked_genes = ",
                as.integer(min_ranked_genes), ").")
      skipped <<- c(skipped, cl)
      return(data.frame(
        cluster    = cl,
        n_ranked   = length(ranked),
        n_sets     = NA_integer_,
        n_sig      = NA_integer_,
        n_enriched = NA_integer_,
        n_depleted = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    if (isTRUE(verbose))
      cat("  Cluster ", cl, ": ",
              format(length(ranked), big.mark = ","), " ranked genes...")

    gsea_cl <- KERAUNOS_gsea(
      ranked_genes  = ranked,
      gene_sets     = gene_sets,
      species       = species,
      db_species    = db_species,
      collection    = collection,
      min_size      = min_size,
      max_size      = max_size,
      fdr_threshold = fdr_threshold,
      verbose       = FALSE   # suppress per-call header; summary printed below
    )

    results_list[[cl]] <<- gsea_cl
    sig_list[[cl]]     <<- gsea_cl$significant

    data.frame(
      cluster    = cl,
      n_ranked   = length(ranked),
      n_sets     = nrow(gsea_cl$results),
      n_sig      = nrow(gsea_cl$significant),
      n_enriched = sum(gsea_cl$significant$NES > 0, na.rm = TRUE),
      n_depleted = sum(gsea_cl$significant$NES < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 KERAUNOS: gsea_pseudobulk %s\n  Clusters   : %d run, %d skipped\n  FDR < %.2f : %s sig. pathways total\n%s\n",
      strrep("\u2500", 26),
      length(results_list),
      length(skipped),
      fdr_threshold,
      format(sum(summary_df$n_sig, na.rm = TRUE), big.mark = ","),
      strrep("\u2500", 56)
    ))
    if (nrow(summary_df) > 0)
      print(summary_df, row.names = FALSE)
  }

  structure(
    list(
      results     = results_list,
      significant = sig_list,
      summary     = summary_df,
      skipped     = skipped,
      params      = list(
        species          = species,
        db_species       = db_species,
        collection       = collection,
        rank_by          = rank_by,
        min_ranked_genes = as.integer(min_ranked_genes),
        min_size         = min_size,
        max_size         = max_size,
        fdr_threshold    = fdr_threshold,
        n_clusters_run   = length(results_list),
        n_gene_sets      = length(gene_sets),
        comparison       = de_result$params$comparison
      )
    ),
    class = "keraunos_gsea_multi"
  )
}

#' @export
print.keraunos_gsea_multi <- function(x, ...) {
  cat("keraunos_gsea_multi object\n")
  cat("  Contrast   :", x$params$comparison, "\n")
  cat("  Clusters   :", x$params$n_clusters_run, "run,",
      length(x$skipped), "skipped\n")
  cat("  Gene sets  :", x$params$n_gene_sets, "tested per cluster\n")
  cat("  Significant:",
      format(sum(x$summary$n_sig, na.rm = TRUE), big.mark = ","),
      paste0("total (FDR < ", x$params$fdr_threshold, ")\n"))
  cat("  Access     : $results (list), $significant (list), $summary, $params\n")
  invisible(x)
}


# ==============================================================================
# Pairwise cluster-vs-cluster GSEA from marker scores (no replicates needed)
# ==============================================================================

#' Run GSEA between pairs of clusters using marker effect sizes
#'
#' A marker-score-based alternative to \code{\link{KERAUNOS_gsea_pseudobulk}}
#' for comparing specific clusters to each other rather than testing a
#' condition within each cluster. Ranking statistics come from
#' \code{scran::scoreMarkers()} pairwise effect sizes (via
#' \code{\link{KERAUNOS_rank_pairwise_markers}}), computed at single-cell
#' resolution — unlike DESeq2 pseudobulk DE, this does not require \eqn{\ge2}
#' biological replicate samples per cluster, making it suitable for comparing
#' two specific clusters directly (e.g. two developmentally-related or
#' visually-adjacent clusters) even when sample composition is skewed or a
#' cluster is dominated by one sample/clone.
#'
#' @param markers_result A \code{keraunos_markers} object from
#'   \code{KERAUNOS_score_markers(sce, ..., full_stats = TRUE)}.
#' @param cluster_pairs \code{NULL} (default) or a list of length-2 character
#'   vectors, e.g. \code{list(c("1","2"), c("8","13"))}, restricting the run
#'   to specific pairs. \code{NULL} runs every pairwise combination of
#'   clusters present in \code{markers_result}.
#' @param stat Character. Which pairwise effect size to rank on — see
#'   \code{\link{KERAUNOS_rank_pairwise_markers}}. Default \code{"AUC"}.
#' @param gene_sets Named list of character vectors, or \code{NULL} to fetch
#'   from MSigDB via msigdbr (fetched once and reused across all pairs).
#' @param species Character. Species for MSigDB gene set retrieval.
#'   Default \code{"Homo sapiens"}.
#' @param db_species Character or \code{NULL}. MSigDB database to query.
#'   \code{NULL} (default): human gene sets; \code{"MM"}: native mouse (use
#'   \code{M}-prefixed collection codes, e.g. \code{c("MH", "M2", "M5")}).
#'   See \code{\link{KERAUNOS_list_gene_sets}} to browse available codes.
#' @param collection Character vector. MSigDB collection codes.
#'   Default \code{c("H", "C2", "C5")}. Subcollection can be appended with a
#'   colon, e.g. \code{"C2:CP:REACTOME"}.
#' @param min_ranked_genes Integer. Minimum number of non-\code{NA} ranked
#'   genes required to run GSEA for a pair. Pairs below this threshold are
#'   skipped. Default \code{10L}.
#' @param min_size Integer. Minimum gene set size after overlap with ranked
#'   genes. Default \code{15L}.
#' @param max_size Integer. Maximum gene set size. Default \code{500L}.
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant}.
#'   Default \code{0.05}.
#' @param verbose Logical. Print per-pair progress and a summary table.
#'   Default \code{TRUE}.
#'
#' @return A \code{keraunos_gsea_pairwise} object (also inherits
#'   \code{keraunos_gsea_multi}, so it is directly compatible with
#'   \code{\link{TALARIA_export_gsea}}) with:
#'   \itemize{
#'     \item \code{$results}     — named list of \code{keraunos_gsea} objects,
#'       one per pair run (names are \code{"<cluster1>_vs_<cluster2>"}).
#'     \item \code{$significant} — named list of significant-only
#'       data.frames, one per pair.
#'     \item \code{$summary}     — data.frame: comparison / cluster1 /
#'       cluster2 / n_ranked / n_sets / n_sig / n_enriched / n_depleted.
#'     \item \code{$skipped}     — character vector of pair labels skipped.
#'     \item \code{$params}      — parameters used.
#'   }
#' @export
KERAUNOS_gsea_markers_pairwise <- function(markers_result,
                                            cluster_pairs    = NULL,
                                            stat             = c("AUC", "logFC.cohen", "logFC.detected"),
                                            gene_sets        = NULL,
                                            species          = "Homo sapiens",
                                            db_species       = NULL,
                                            collection       = c("H", "C2", "C5"),
                                            min_ranked_genes = 10L,
                                            min_size         = 15L,
                                            max_size         = 500L,
                                            fdr_threshold    = 0.05,
                                            verbose          = TRUE) {

  if (!inherits(markers_result, "keraunos_markers"))
    stop("'markers_result' must be a keraunos_markers object from ",
         "KERAUNOS_score_markers().", call. = FALSE)

  if (!isTRUE(markers_result$params$full_stats))
    stop("'markers_result' was scored without full_stats = TRUE, so no ",
         "per-cluster-pair effect sizes are available. Re-run ",
         "KERAUNOS_score_markers(sce, ..., full_stats = TRUE).", call. = FALSE)

  stat <- match.arg(stat)

  cl_names <- names(markers_result$markers)
  if (is.null(cluster_pairs)) {
    if (length(cl_names) < 2L)
      stop("'markers_result' has fewer than 2 clusters — nothing to compare.",
           call. = FALSE)
    pair_mat     <- utils::combn(cl_names, 2L)
    cluster_pairs <- lapply(seq_len(ncol(pair_mat)), function(i) pair_mat[, i])
  } else {
    bad <- Filter(function(p) length(p) != 2L, cluster_pairs)
    if (length(bad) > 0L)
      stop("Every element of 'cluster_pairs' must be a length-2 vector.",
           call. = FALSE)
    unknown <- setdiff(unlist(cluster_pairs), cl_names)
    if (length(unknown) > 0L)
      stop("Cluster(s) not found in 'markers_result': ",
           paste(unknown, collapse = ", "), call. = FALSE)
  }

  # ── Fetch gene sets once ────────────────────────────────────────────────────
  if (is.null(gene_sets)) {
    if (!requireNamespace("fgsea", quietly = TRUE))
      stop("Package 'fgsea' is required.  ",
           "Install via: BiocManager::install(\"fgsea\")", call. = FALSE)
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Package 'msigdbr' is required when gene_sets = NULL.  ",
           "Install via: install.packages(\"msigdbr\")", call. = FALSE)

    if (!is.null(db_species) && toupper(db_species) == "MM") {
      wrong <- collection[!grepl("^M", collection)]
      if (length(wrong) > 0)
        warning(
          "db_species='MM' (mouse) but collection code(s) do not use the ",
          "mouse prefix: ", paste(wrong, collapse = ", "), ".\n",
          "  Mouse MSigDB collections use uppercase M-prefixed codes ",
          "(e.g. 'MH', 'M2', 'M5').\n",
          "  Run KERAUNOS_list_gene_sets(db_species='MM') to see all available codes.",
          call. = FALSE
        )
    }

    if (isTRUE(verbose)) {
      db_str <- if (!is.null(db_species))
        paste0(" [db_species=", db_species, "]") else ""
      cat("Fetching MSigDB gene sets (", paste(collection, collapse = " + "),
              ") for ", species, db_str, "...")
    }

    # msigdbr::msigdbr()'s db_species argument requires a character vector --
    # unlike msigdbr_collections(), it errors on NULL rather than falling
    # back to its own default ("HS"), so NULL must be resolved here.
    msigdbr_db_species <- if (is.null(db_species)) "HS" else db_species

    gs_dfs <- lapply(collection, function(col) {
      cat_sub <- strsplit(col, ":", fixed = TRUE)[[1]]
      cat  <- cat_sub[1]
      subc <- if (length(cat_sub) > 1) paste(cat_sub[-1], collapse = ":") else NULL
      tryCatch(
        msigdbr::msigdbr(species = species, db_species = msigdbr_db_species,
                         collection = cat, subcollection = subc),
        error = function(e) {
          warning("msigdbr failed for collection '", col, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    })
    gs_dfs <- Filter(Negate(is.null), gs_dfs)
    if (length(gs_dfs) == 0L)
      stop("No gene sets retrieved from MSigDB.  ",
           "Check species name and collection codes.", call. = FALSE)

    gs_df     <- do.call(rbind, gs_dfs)
    gene_sets <- split(gs_df$gene_symbol, gs_df$gs_name)
    if (isTRUE(verbose))
      cat("  Retrieved ", format(length(gene_sets), big.mark = ","),
              " gene sets.")
  }

  # ── Per-pair GSEA ────────────────────────────────────────────────────────────
  results_list <- list()
  sig_list     <- list()
  skipped      <- character(0)

  summary_rows <- lapply(cluster_pairs, function(p) {

    cl1   <- p[1L]; cl2 <- p[2L]
    label <- paste0(cl1, "_vs_", cl2)

    ranked <- KERAUNOS_rank_pairwise_markers(markers_result, cl1, cl2,
                                              stat = stat, verbose = FALSE)

    if (length(ranked) < as.integer(min_ranked_genes)) {
      if (isTRUE(verbose))
        cat("  ", label, ": skipped (", length(ranked),
                " valid ranked genes < min_ranked_genes = ",
                as.integer(min_ranked_genes), ").")
      skipped <<- c(skipped, label)
      return(data.frame(
        comparison = label, cluster1 = cl1, cluster2 = cl2,
        n_ranked = length(ranked), n_sets = NA_integer_,
        n_sig = NA_integer_, n_enriched = NA_integer_, n_depleted = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    if (isTRUE(verbose))
      cat("  ", label, ": ", format(length(ranked), big.mark = ","),
              " ranked genes...")

    gsea_p <- KERAUNOS_gsea(
      ranked_genes  = ranked,
      gene_sets     = gene_sets,
      species       = species,
      db_species    = db_species,
      collection    = collection,
      min_size      = min_size,
      max_size      = max_size,
      fdr_threshold = fdr_threshold,
      verbose       = FALSE
    )

    results_list[[label]] <<- gsea_p
    sig_list[[label]]     <<- gsea_p$significant

    data.frame(
      comparison = label, cluster1 = cl1, cluster2 = cl2,
      n_ranked   = length(ranked),
      n_sets     = nrow(gsea_p$results),
      n_sig      = nrow(gsea_p$significant),
      n_enriched = sum(gsea_p$significant$NES > 0, na.rm = TRUE),
      n_depleted = sum(gsea_p$significant$NES < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)

  if (isTRUE(verbose)) {
    cat(sprintf(
      "── KERAUNOS: gsea_markers_pairwise %s\n  Pairs      : %d run, %d skipped\n  FDR < %.2f : %s sig. pathways total\n%s\n",
      strrep("─", 15),
      length(results_list),
      length(skipped),
      fdr_threshold,
      format(sum(summary_df$n_sig, na.rm = TRUE), big.mark = ","),
      strrep("─", 56)
    ))
    if (nrow(summary_df) > 0)
      print(summary_df, row.names = FALSE)
  }

  structure(
    list(
      results     = results_list,
      significant = sig_list,
      summary     = summary_df,
      skipped     = skipped,
      params      = list(
        stat             = stat,
        species          = species,
        db_species       = db_species,
        collection       = collection,
        min_ranked_genes = as.integer(min_ranked_genes),
        min_size         = min_size,
        max_size         = max_size,
        fdr_threshold    = fdr_threshold,
        n_pairs_run      = length(results_list),
        n_gene_sets      = length(gene_sets),
        comparison       = paste0("pairwise marker-based GSEA (", stat, ")")
      )
    ),
    class = c("keraunos_gsea_pairwise", "keraunos_gsea_multi")
  )
}

#' @export
print.keraunos_gsea_pairwise <- function(x, ...) {
  cat("keraunos_gsea_pairwise object\n")
  cat("  Ranking    :", x$params$comparison, "\n")
  cat("  Pairs      :", x$params$n_pairs_run, "run,",
      length(x$skipped), "skipped\n")
  cat("  Gene sets  :", x$params$n_gene_sets, "tested per pair\n")
  cat("  Significant:",
      format(sum(x$summary$n_sig, na.rm = TRUE), big.mark = ","),
      paste0("total (FDR < ", x$params$fdr_threshold, ")\n"))
  cat("  Access     : $results (list), $significant (list), $summary, $params\n")
  invisible(x)
}


# ==============================================================================
# Multi-cluster ORA on pseudobulk DE results (up/down split)
# ==============================================================================

#' Run ORA across all clusters from pseudobulk DE results
#'
#' A convenience wrapper around \code{\link{KERAUNOS_ora}} that accepts a
#' \code{keraunos_pseudobulk} object and runs ORA separately on upregulated
#' and downregulated significant genes for every cluster.
#'
#' Significance is determined using whichever threshold was applied when
#' \code{\link{KERAUNOS_de_pseudobulk}} was run (FDR or raw p-value, stored in
#' \code{de_result$params}).  Direction is determined by \code{log2FoldChange}:
#' genes with \eqn{LFC > lfc\_threshold} are "up"; genes with
#' \eqn{LFC < -lfc\_threshold} are "down".
#'
#' Clusters, or directions within a cluster, with fewer than
#' \code{min_sig_genes} significant genes are skipped with a message.
#'
#' @param de_result A \code{keraunos_pseudobulk} object (output of
#'   \code{\link{KERAUNOS_de_pseudobulk}}).
#' @param organism Character. gprofiler2 organism code.  Default
#'   \code{"hsapiens"}.  Use \code{"mmusculus"} for mouse.
#' @param sources Character vector. Databases to query.  Default
#'   \code{c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")}.
#' @param lfc_threshold Numeric. Minimum absolute \eqn{\log_2} fold-change for
#'   a gene to be considered directionally significant.  Genes with
#'   \eqn{|LFC| \leq lfc\_threshold} are excluded from both up and down sets.
#'   Default \code{0} (all significant genes included).
#' @param background Character. Background gene set to pass to gprofiler2.
#'   \code{"tested"} (default) — all genes tested in the cluster (statistically
#'   correct for pseudobulk); \code{"genome"} — gprofiler2 reference genome
#'   (no custom background).
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant} ORA
#'   terms.  Default \code{0.05}.
#' @param min_sig_genes Integer. Minimum number of significant DE genes
#'   required to run ORA for a given cluster/direction.  Default \code{5L}.
#' @param verbose Logical. Print per-cluster progress and a summary table.
#'   Default \code{TRUE}.
#'
#' @return A \code{keraunos_ora_multi} object (list) with:
#'   \itemize{
#'     \item \code{$results}     — named list (cluster → list with \code{$up}
#'       and \code{$down}, each a \code{keraunos_ora} object or \code{NULL}
#'       if skipped).
#'     \item \code{$significant} — named list (cluster → list with \code{$up}
#'       and \code{$down} significant-only data.frames).
#'     \item \code{$summary}     — data.frame:
#'       cluster / n_sig_genes / n_up_genes / n_down_genes /
#'       n_sig_terms_up / n_sig_terms_down.
#'     \item \code{$skipped}     — character vector of cluster names where
#'       neither direction had enough genes.
#'     \item \code{$params}      — parameters used.
#'   }
#' @export
KERAUNOS_ora_pseudobulk <- function(de_result,
                                     organism      = "hsapiens",
                                     sources       = c("GO:BP", "GO:MF",
                                                       "GO:CC", "KEGG",
                                                       "REAC"),
                                     lfc_threshold = 0,
                                     background    = c("tested", "genome"),
                                     fdr_threshold = 0.05,
                                     min_sig_genes = 5L,
                                     verbose       = TRUE) {

  if (!inherits(de_result, "keraunos_pseudobulk"))
    stop("'de_result' must be a keraunos_pseudobulk object (output of ",
         "KERAUNOS_de_pseudobulk()).", call. = FALSE)

  if (length(de_result$results) == 0L)
    stop("No cluster results found in 'de_result'.", call. = FALSE)

  background <- match.arg(background)

  # Resolve the significance filter used when pseudobulk was run
  pb_sig_col    <- de_result$params$sig_col
  pb_sig_thresh <- de_result$params$sig_thresh

  # ── Per-cluster ORA ─────────────────────────────────────────────────────────
  cluster_names <- names(de_result$results)
  results_list  <- list()
  sig_list      <- list()
  skipped       <- character(0)

  summary_rows <- lapply(cluster_names, function(cl) {

    df <- de_result$results[[cl]]

    # Background: all genes tested in this cluster
    bg <- if (background == "tested") df$gene else NULL

    # Identify significant DE genes respecting original threshold + LFC filter
    is_sig <- !is.na(df[[pb_sig_col]]) & df[[pb_sig_col]] < pb_sig_thresh &
              !is.na(df$log2FoldChange)

    up_genes   <- df$gene[is_sig & df$log2FoldChange >  lfc_threshold]
    down_genes <- df$gene[is_sig & df$log2FoldChange < -lfc_threshold]

    n_up   <- length(up_genes)
    n_down <- length(down_genes)
    n_sig  <- n_up + n_down

    if (n_up < as.integer(min_sig_genes) && n_down < as.integer(min_sig_genes)) {
      if (isTRUE(verbose))
        cat("  Cluster ", cl, ": skipped (", n_up, " up / ", n_down,
                " down genes — both below min_sig_genes = ",
                as.integer(min_sig_genes), ").")
      skipped <<- c(skipped, cl)
      return(data.frame(
        cluster         = cl,
        n_sig_genes     = n_sig,
        n_up_genes      = n_up,
        n_down_genes    = n_down,
        n_sig_terms_up  = NA_integer_,
        n_sig_terms_down = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    if (isTRUE(verbose))
      cat("  Cluster ", cl, ": ", n_up, " up / ", n_down, " down genes.")

    # Run ORA for upregulated genes
    ora_up <- NULL
    if (n_up >= as.integer(min_sig_genes)) {
      ora_up <- KERAUNOS_ora(
        genes         = up_genes,
        background    = bg,
        organism      = organism,
        sources       = sources,
        fdr_threshold = fdr_threshold,
        verbose       = FALSE
      )
    } else if (isTRUE(verbose)) {
      cat("    Up: skipped (", n_up, " genes < min_sig_genes = ",
              as.integer(min_sig_genes), ").")
    }

    # Run ORA for downregulated genes
    ora_down <- NULL
    if (n_down >= as.integer(min_sig_genes)) {
      ora_down <- KERAUNOS_ora(
        genes         = down_genes,
        background    = bg,
        organism      = organism,
        sources       = sources,
        fdr_threshold = fdr_threshold,
        verbose       = FALSE
      )
    } else if (isTRUE(verbose)) {
      cat("    Down: skipped (", n_down, " genes < min_sig_genes = ",
              as.integer(min_sig_genes), ").")
    }

    results_list[[cl]] <<- list(up = ora_up, down = ora_down)
    sig_list[[cl]] <<- list(
      up   = if (!is.null(ora_up))   ora_up$significant   else data.frame(),
      down = if (!is.null(ora_down)) ora_down$significant else data.frame()
    )

    data.frame(
      cluster          = cl,
      n_sig_genes      = n_sig,
      n_up_genes       = n_up,
      n_down_genes     = n_down,
      n_sig_terms_up   = if (!is.null(ora_up))   nrow(ora_up$significant)   else NA_integer_,
      n_sig_terms_down = if (!is.null(ora_down)) nrow(ora_down$significant) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)

  if (isTRUE(verbose)) {
    total_up   <- sum(summary_df$n_sig_terms_up,   na.rm = TRUE)
    total_down <- sum(summary_df$n_sig_terms_down, na.rm = TRUE)
    cat(sprintf(
      "\u2500\u2500 KERAUNOS: ora_pseudobulk %s\n  Clusters   : %d run, %d skipped\n  Organism   : %s\n  FDR < %.2f : %s sig. terms up | %s sig. terms down\n%s\n",
      strrep("\u2500", 27),
      length(results_list),
      length(skipped),
      organism,
      fdr_threshold,
      format(total_up,   big.mark = ","),
      format(total_down, big.mark = ","),
      strrep("\u2500", 56)
    ))
    if (nrow(summary_df) > 0)
      print(summary_df, row.names = FALSE)
  }

  structure(
    list(
      results     = results_list,
      significant = sig_list,
      summary     = summary_df,
      skipped     = skipped,
      params      = list(
        organism       = organism,
        sources        = sources,
        lfc_threshold  = lfc_threshold,
        background     = background,
        fdr_threshold  = fdr_threshold,
        min_sig_genes  = as.integer(min_sig_genes),
        pb_sig_col     = pb_sig_col,
        pb_sig_thresh  = pb_sig_thresh,
        n_clusters_run = length(results_list),
        comparison     = de_result$params$comparison
      )
    ),
    class = "keraunos_ora_multi"
  )
}

#' @export
print.keraunos_ora_multi <- function(x, ...) {
  cat("keraunos_ora_multi object\n")
  cat("  Contrast   :", x$params$comparison, "\n")
  cat("  Clusters   :", x$params$n_clusters_run, "run,",
      length(x$skipped), "skipped\n")
  cat("  Organism   :", x$params$organism, "\n")
  total_up   <- sum(x$summary$n_sig_terms_up,   na.rm = TRUE)
  total_down <- sum(x$summary$n_sig_terms_down, na.rm = TRUE)
  cat("  Sig. terms :", format(total_up, big.mark = ","), "up |",
      format(total_down, big.mark = ","),
      paste0("down (FDR < ", x$params$fdr_threshold, ")\n"))
  cat("  Access     : $results (list), $significant (list), $summary, $params\n")
  invisible(x)
}
