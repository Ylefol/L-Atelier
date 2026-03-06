# ==============================================================================
# PANDORA - Cell Type Annotation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Fashioned from clay by Hephaestus at Zeus's command; every god bestowed
# a gift upon her (pan = all, dora = gifts). PANDORA assigns identities
# drawn from all available reference knowledge — each cell receives its
# gifts of cell type labels and markers.
#
# Core responsibilities:
#   - Reference-based annotation (SingleR / celldex)
#   - Marker-based annotation (scType, bundled ScTypeDB)
#   - Manual label assignment (planned: PANDORA_assign_labels)
#
# All functions prefixed: PANDORA_
# ==============================================================================


# ── celldex shorthand lookup table ────────────────────────────────────────────
.pandora_celldex_refs <- data.frame(
  shorthand  = c("HumanPrimaryCellAtlas", "BlueprintEncode", "DICE",
                 "Monaco", "Novershtern", "ImmGen", "MouseRNAseq"),
  fn         = c("HumanPrimaryCellAtlasData", "BlueprintEncodeData",
                 "DatabaseImmuneCellExpressionData", "MonacoImmuneData",
                 "NovershternHematopoieticData", "ImmGenData",
                 "MouseRNAseqData"),
  organism   = c("Human", "Human", "Human", "Human", "Human",
                 "Mouse", "Mouse"),
  coverage   = c("Broad", "Immune + stromal", "Immune", "Immune",
                 "Hematopoietic", "Immune", "Broad"),
  label_main = c(37L, 24L, 5L, 11L, 17L, 20L, 18L),
  label_fine = c(157L, 43L, 17L, 29L, 38L, 253L, NA_integer_),
  stringsAsFactors = FALSE
)


#' List available annotation references
#'
#' Prints a formatted summary of supported \code{celldex} reference shorthands
#' for \code{\link{PANDORA_annotate_singler}} and available tissue types from
#' the bundled ScTypeDB for \code{\link{PANDORA_annotate_sctype}}.  Use this
#' to quickly confirm whether either approach covers your dataset before
#' committing to an annotation strategy.
#'
#' @param source Character. Which sources to list: \code{"all"} (default),
#'   \code{"singler"}, or \code{"sctype"}.
#' @param db \code{NULL} (use bundled \code{ScTypeDB_full.xlsx}), a file path
#'   to a custom Excel or CSV, or a \code{data.frame} already loaded in R.
#'   Only relevant when \code{source} includes \code{"sctype"}.
#'
#' @return Invisibly returns a named list with \code{$singler} (reference info
#'   \code{data.frame}) and \code{$sctype} (character vector of tissue types).
#' @export
PANDORA_list_references <- function(source = c("all", "singler", "sctype"),
                                    db     = NULL) {
  source <- match.arg(source)
  bar    <- strrep("\u2500", 73)

  singler_out <- NULL
  sctype_out  <- NULL

  # ── SingleR / celldex ────────────────────────────────────────────────────────
  if (source %in% c("all", "singler")) {
    refs     <- .pandora_celldex_refs
    fine_str <- ifelse(is.na(refs$label_fine), "\u2014",
                       as.character(refs$label_fine))
    type_str <- paste0(refs$label_main, " main / ", fine_str, " fine")

    cat("\u2500\u2500 PANDORA: Available annotation references ",
        strrep("\u2500", 27), "\n\n", sep = "")
    cat("SingleR / celldex  (reference= in PANDORA_annotate_singler)\n")
    cat(bar, "\n", sep = "")
    cat(sprintf("  %-30s %-8s %-18s %s\n",
                "Shorthand", "Organism", "Coverage", "Cell types"))
    cat(bar, "\n", sep = "")
    for (i in seq_len(nrow(refs))) {
      cat(sprintf("  %-30s %-8s %-18s %s\n",
                  refs$shorthand[i], refs$organism[i],
                  refs$coverage[i],  type_str[i]))
    }
    cat("\n  Pass labels_col = \"label.main\" (default) or \"label.fine\"",
        "for finer resolution.\n")
    cat("  Or pass a pre-loaded SummarizedExperiment as reference=.\n")
    cat(bar, "\n\n", sep = "")
    singler_out <- refs
  }

  # ── scType tissues ───────────────────────────────────────────────────────────
  if (source %in% c("all", "sctype")) {
    db_df   <- .pandora_load_sctype_db(db)
    tissues <- sort(unique(db_df$tissueType))

    cat("scType tissues  (tissue= in PANDORA_annotate_sctype)\n")
    cat(bar, "\n", sep = "")
    tissue_str <- paste(tissues, collapse = ", ")
    cat(paste(strwrap(tissue_str, width = 71, indent = 2, exdent = 2),
              collapse = "\n"), "\n\n", sep = "")

    ct_counts <- table(db_df$tissueType)
    cat("  Cell types per tissue:\n")
    for (tis in tissues)
      cat(sprintf("    %-24s %d types\n", tis, ct_counts[[tis]]))
    cat(bar, "\n", sep = "")

    sctype_out <- tissues
  }

  invisible(list(singler = singler_out, sctype = sctype_out))
}


#' Reference-based cell type annotation via SingleR
#'
#' Annotates cells (or clusters) by correlating expression profiles against a
#' reference dataset using \code{SingleR::SingleR()}.  A pre-loaded
#' \code{SummarizedExperiment} reference or a shorthand string for a
#' \code{celldex} dataset can be supplied.  Run
#' \code{\link{PANDORA_list_references}} to see available shorthand values.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"logcounts"} (or
#'   other) assay.
#' @param reference A \code{SummarizedExperiment} reference, or a character
#'   shorthand: one of \code{"HumanPrimaryCellAtlas"},
#'   \code{"BlueprintEncode"}, \code{"DICE"}, \code{"Monaco"},
#'   \code{"Novershtern"} (human), or \code{"ImmGen"},
#'   \code{"MouseRNAseq"} (mouse).
#' @param labels_col Character. Column in the reference \code{colData} to use
#'   as labels.  Default \code{"label.main"}.  Use \code{"label.fine"} for
#'   finer resolution.
#' @param cluster_col Character or \code{NULL}.  If \code{NULL} (default),
#'   annotation is at the cell level.  If a \code{colData} column name is
#'   provided, expression is aggregated per cluster before annotation —
#'   one label per cluster is broadcast back to all cells in that cluster.
#' @param assay_name Character. Assay to use. Default \code{"logcounts"}.
#' @param label_col Character. \code{colData} column for the result.
#'   Default \code{"singler_label"}.
#' @param prune Logical. Replace low-confidence labels with \code{NA} using
#'   SingleR's delta-score pruning.  Default \code{TRUE}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} populated.
#' @export
PANDORA_annotate_singler <- function(sce,
                                     reference,
                                     labels_col  = "label.main",
                                     cluster_col = NULL,
                                     assay_name  = "logcounts",
                                     label_col   = "singler_label",
                                     prune       = TRUE,
                                     verbose     = TRUE) {

  if (!requireNamespace("SingleR", quietly = TRUE))
    stop("Package 'SingleR' is required. ",
         "Install via: BiocManager::install(\"SingleR\")", call. = FALSE)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  # ── Resolve reference ────────────────────────────────────────────────────────
  ref <- .pandora_resolve_celldex_ref(reference)

  if (!labels_col %in% colnames(colData(ref)))
    stop("'", labels_col, "' not found in reference colData. ",
         "Available: ", paste(colnames(colData(ref)), collapse = ", "),
         call. = FALSE)

  # ── Cluster vector (NULL = cell-level) ───────────────────────────────────────
  clusters <- NULL
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(colData(sce)))
      stop("cluster_col '", cluster_col, "' not found in colData(sce).",
           call. = FALSE)
    clusters <- colData(sce)[[cluster_col]]
  }

  # ── Run SingleR ──────────────────────────────────────────────────────────────
  result <- SingleR::SingleR(
    test            = sce,
    ref             = ref,
    labels          = ref[[labels_col]],
    clusters        = clusters,
    assay.type.test = assay_name
  )

  # ── Extract labels — broadcast cluster results back to cells if needed ───────
  label_vec <- if (isTRUE(prune)) result$pruned.labels else result$labels

  if (!is.null(clusters)) {
    cell_labels         <- label_vec[as.character(clusters)]
    names(cell_labels)  <- colnames(sce)
  } else {
    cell_labels <- label_vec
  }

  colData(sce)[[label_col]] <- cell_labels

  # ── Summary ──────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    ref_name    <- if (is.character(reference)) reference else "custom"
    n_annotated <- sum(!is.na(cell_labels))
    n_pruned    <- sum(is.na(cell_labels))
    n_types     <- length(unique(na.omit(cell_labels)))
    mode_str    <- if (is.null(clusters)) "cell-level" else
      paste0("cluster-level (via '", cluster_col, "')")
    cat(sprintf(
      "\u2500\u2500 PANDORA: SingleR annotation %s\n  Reference  : %s (%s)\n  Mode       : %s\n  Annotated  : %s cells  |  Pruned (low-confidence): %s\n  Cell types : %d unique types\n  Stored as  : colData(sce)[[\"%s\"]]\n%s\n",
      strrep("\u2500", 26),
      ref_name, labels_col,
      mode_str,
      format(n_annotated, big.mark = ","),
      format(n_pruned,    big.mark = ","),
      n_types,
      label_col,
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' Marker-based cell type annotation via scType
#'
#' Scores clusters against a curated marker gene database to assign cell type
#' labels.  Annotation is cluster-level: expression is aggregated per cluster,
#' scored against each cell type's positive and negative markers (from the
#' bundled \code{ScTypeDB_full.xlsx}), and the highest-scoring type is assigned
#' to all cells in each cluster.
#'
#' Use \code{\link{PANDORA_list_references}(source = "sctype")} to see
#' available tissues and cell type counts before running.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{"logcounts"} (or
#'   other) assay.
#' @param tissue Character. Tissue type string matching an entry in the
#'   scTypeDB (case-sensitive).
#' @param cluster_col Character. \code{colData} column with cluster
#'   assignments.  Default \code{"cluster"}.
#' @param db \code{NULL} (bundled \code{ScTypeDB_full.xlsx}), a file path to a
#'   custom Excel or CSV, or a \code{data.frame} with columns
#'   \code{tissueType}, \code{cellName}, \code{geneSymbolmore1} (positive
#'   markers, comma-separated), \code{geneSymbolmore2} (negative markers,
#'   comma-separated; \code{NA} = none).
#' @param assay_name Character. Assay to score. Default \code{"logcounts"}.
#' @param label_col Character. \code{colData} column for the result.
#'   Default \code{"sctype_label"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} populated.
#' @export
PANDORA_annotate_sctype <- function(sce,
                                    tissue,
                                    cluster_col = "cluster",
                                    db          = NULL,
                                    assay_name  = "logcounts",
                                    label_col   = "sctype_label",
                                    verbose     = TRUE) {

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  if (!cluster_col %in% names(colData(sce)))
    stop("cluster_col '", cluster_col, "' not found in colData(sce).",
         call. = FALSE)

  # ── Load and filter DB ───────────────────────────────────────────────────────
  db_df <- .pandora_load_sctype_db(db)
  db_df <- db_df[db_df$tissueType == tissue, , drop = FALSE]

  if (nrow(db_df) == 0) {
    avail <- paste(sort(unique(.pandora_load_sctype_db(db)$tissueType)),
                   collapse = ", ")
    stop("Tissue '", tissue, "' not found in scTypeDB.\n",
         "  Available: ", avail, call. = FALSE)
  }

  # ── Parse marker lists ───────────────────────────────────────────────────────
  pos_markers <- lapply(db_df$geneSymbolmore1,
                        function(x) unlist(strsplit(x, ",")))
  neg_markers <- lapply(db_df$geneSymbolmore2,
                        function(x) {
                          if (is.na(x) || nchar(trimws(x)) == 0)
                            return(character(0))
                          unlist(strsplit(x, ","))
                        })
  cell_types <- db_df$cellName

  # ── Subset to marker genes (safe for BPCells — small subset) ────────────────
  all_markers <- unique(c(unlist(pos_markers), unlist(neg_markers)))
  all_markers <- intersect(all_markers, rownames(sce))

  if (length(all_markers) == 0)
    stop("None of the DB marker genes are present in the SCE. ",
         "Check that rownames(sce) are HGNC gene symbols.", call. = FALSE)

  mat <- assay(sce, assay_name)[all_markers, , drop = FALSE]
  if (inherits(mat, "IterableMatrix"))
    mat <- as(mat, "dgCMatrix")
  mat <- as.matrix(mat)

  # ── Scale per gene (mean 0, sd 1) capped at ±10 ─────────────────────────────
  gene_means <- rowMeans(mat)
  gene_sds   <- apply(mat, 1, stats::sd)
  gene_sds   <- pmax(gene_sds, .Machine$double.eps)
  scaled_mat <- (mat - gene_means) / gene_sds
  scaled_mat <- pmin(pmax(scaled_mat, -10), 10)

  # ── Score cells against each cell type ───────────────────────────────────────
  scores <- .pandora_sctype_score(scaled_mat, cell_types,
                                  pos_markers, neg_markers, all_markers)

  # ── Aggregate by cluster and assign best-scoring label ───────────────────────
  clusters        <- as.character(colData(sce)[[cluster_col]])
  unique_clusters <- sort(unique(clusters))

  cluster_labels <- vapply(unique_clusters, function(cl) {
    cells     <- which(clusters == cl)
    cl_scores <- rowSums(scores[, cells, drop = FALSE])
    cell_types[which.max(cl_scores)]
  }, character(1))
  names(cluster_labels) <- unique_clusters

  colData(sce)[[label_col]] <- cluster_labels[clusters]

  # ── Summary ──────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    n_queried <- length(unique(c(unlist(pos_markers), unlist(neg_markers))))
    n_found   <- length(all_markers)
    n_types   <- length(unique(cluster_labels))
    cat(sprintf(
      "\u2500\u2500 PANDORA: scType annotation %s\n  Tissue     : %s\n  Clusters   : %d  \u2192  %d unique cell types\n  Markers    : %d found in SCE (of %d queried)\n  Stored as  : colData(sce)[[\"%s\"]]\n%s\n",
      strrep("\u2500", 27),
      tissue,
      length(unique_clusters), n_types,
      n_found, n_queried,
      label_col,
      strrep("\u2500", 56)
    ))
  }

  sce
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Load scTypeDB from bundled file, user-supplied path, or data.frame.
.pandora_load_sctype_db <- function(db) {
  if (is.data.frame(db)) {
    required <- c("tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2")
    missing  <- setdiff(required, names(db))
    if (length(missing) > 0)
      stop("Custom db missing columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    return(db)
  }

  if (is.null(db))
    db <- system.file("extdata", "ScTypeDB_full.xlsx", package = "CAULDRON")

  if (!file.exists(db))
    stop("scTypeDB file not found: ", db, call. = FALSE)

  ext <- tolower(tools::file_ext(db))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Package 'readxl' is required to read .xlsx files. ",
           "Install via: install.packages(\"readxl\")", call. = FALSE)
    df <- readxl::read_excel(db)
  } else if (ext == "csv") {
    df <- utils::read.csv(db, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported scTypeDB format '.", ext, "'. Use .xlsx or .csv.",
         call. = FALSE)
  }

  as.data.frame(df, stringsAsFactors = FALSE)
}


# Resolve a celldex shorthand string to a SummarizedExperiment.
# Returns reference unchanged if already a SummarizedExperiment.
.pandora_resolve_celldex_ref <- function(reference) {
  if (inherits(reference, "SummarizedExperiment")) return(reference)

  if (!is.character(reference) || length(reference) != 1L)
    stop("'reference' must be a SummarizedExperiment or a single character ",
         "shorthand.  Run PANDORA_list_references() to see options.",
         call. = FALSE)

  if (!requireNamespace("celldex", quietly = TRUE))
    stop("Package 'celldex' is required to use reference shorthands. ",
         "Install via: BiocManager::install(\"celldex\")", call. = FALSE)

  idx <- match(reference, .pandora_celldex_refs$shorthand)
  if (is.na(idx))
    stop("Unknown reference shorthand '", reference, "'. ",
         "Run PANDORA_list_references() to see valid options.", call. = FALSE)

  fn <- getExportedValue("celldex", .pandora_celldex_refs$fn[idx])
  fn()
}


# Score cells against cell type marker sets.
# Returns a cell_types x cells matrix of scType scores.
.pandora_sctype_score <- function(scaled_mat, cell_types,
                                  pos_markers, neg_markers, all_markers) {
  n_types <- length(cell_types)
  n_cells <- ncol(scaled_mat)
  scores  <- matrix(0, nrow = n_types, ncol = n_cells,
                    dimnames = list(cell_types, colnames(scaled_mat)))

  for (i in seq_len(n_types)) {
    pos       <- intersect(pos_markers[[i]], all_markers)
    neg       <- intersect(neg_markers[[i]], all_markers)
    n_markers <- length(pos) + length(neg)
    if (n_markers == 0L) next

    pos_score <- if (length(pos) > 0)
      colSums(scaled_mat[pos, , drop = FALSE]) else 0
    neg_score <- if (length(neg) > 0)
      colSums(scaled_mat[neg, , drop = FALSE]) else 0

    scores[i, ] <- (pos_score - neg_score) / sqrt(n_markers)
  }

  scores
}
