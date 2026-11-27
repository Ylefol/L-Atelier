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
#   - Anchor-based label transfer (Seurat FindTransferAnchors/TransferData)
#   - Manual label assignment (PANDORA_assign_labels)
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
    names(label_vec)    <- rownames(result)
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
#' The bundled ScTypeDB uses human HGNC gene symbols.  For non-human data,
#' set \code{species} to the target organism and orthologous gene symbols will
#' be resolved via \code{orthogene} before scoring.  The conversion report
#' (database used, mapping rate) is printed to allow you to judge the
#' reliability of the annotation before interpreting results.
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
#' @param species Character. Target organism.  \code{"human"} (default) uses
#'   the DB markers directly.  Any other value (e.g. \code{"mouse"},
#'   \code{"Mus musculus"}) triggers ortholog conversion via \code{orthogene}.
#'   Accepts any species string recognised by \code{orthogene}.
#' @param ortholog_method Character. Ortholog database passed to
#'   \code{orthogene::convert_orthologs()}.  Default \code{"homologene"}
#'   (offline, fast, covers common model organisms).  Use \code{"gprofiler"}
#'   for broader species coverage (requires internet).
#' @param assay_name Character. Assay to score. Default \code{"logcounts"}.
#' @param label_col Character. \code{colData} column for the result.
#'   Default \code{"sctype_label"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} populated.
#' @export
PANDORA_annotate_sctype <- function(sce,
                                    tissue,
                                    cluster_col      = "cluster",
                                    db               = NULL,
                                    species          = "human",
                                    ortholog_method  = "homologene",
                                    assay_name       = "logcounts",
                                    label_col        = "sctype_label",
                                    verbose          = TRUE) {

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

  # ── Ortholog conversion (non-human species) ──────────────────────────────────
  ortho_info <- NULL
  if (!identical(tolower(trimws(species)), "human")) {
    ortho_info  <- .pandora_convert_sctype_markers(
      pos_markers, neg_markers, species, ortholog_method, verbose
    )
    pos_markers <- ortho_info$pos_markers
    neg_markers <- ortho_info$neg_markers
  }

  # ── Subset to marker genes (safe for BPCells — small subset) ────────────────
  n_db_markers <- length(unique(c(unlist(pos_markers), unlist(neg_markers))))
  all_markers  <- unique(c(unlist(pos_markers), unlist(neg_markers)))
  all_markers  <- intersect(all_markers, rownames(sce))

  if (length(all_markers) == 0) {
    species_hint <- if (!is.null(ortho_info))
      "Check species name and ortholog_method." else
      "Check that rownames(sce) are HGNC gene symbols."
    stop("None of the marker genes are present in the SCE. ", species_hint,
         call. = FALSE)
  }

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
    n_types  <- length(unique(cluster_labels))
    cat(sprintf(
      "\u2500\u2500 PANDORA: scType annotation %s\n  Tissue     : %s\n  Species    : %s\n",
      strrep("\u2500", 27), tissue, species
    ))
    if (!is.null(ortho_info)) {
      cat(sprintf(
        "  Orthologs  : %d / %d human markers mapped (%s%%)  |  DB: %s\n",
        ortho_info$n_mapped, ortho_info$n_total,
        ortho_info$pct, ortholog_method
      ))
      if (ortho_info$n_mapped < ortho_info$n_total)
        cat(sprintf("  Unmapped   : %d markers excluded\n",
                    ortho_info$n_total - ortho_info$n_mapped))
      if (ortho_info$pct_num < 50)
        cat("  WARNING    : Conversion rate <50% - interpret results with caution.\n")
    }
    cat(sprintf(
      "  Clusters   : %d  \u2192  %d unique cell types\n  Markers    : %d found in SCE (of %d post-conversion)\n  Stored as  : colData(sce)[[\"%s\"]]\n%s\n",
      length(unique_clusters), n_types,
      length(all_markers), n_db_markers,
      label_col,
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' Anchor-based label transfer via Seurat
#'
#' Annotates cells by projecting them into a pre-reduced reference's embedding
#' and transferring labels using \code{Seurat::FindTransferAnchors()} /
#' \code{Seurat::TransferData()}.  Unlike \code{\link{PANDORA_annotate_singler}}
#' (correlation-based) and \code{\link{PANDORA_annotate_sctype}}
#' (marker-score-based), this performs anchor-based transfer — matching the
#' Seurat label-transfer method used alongside SingleR-style correlation and
#' marker-score approaches in multi-method annotation designs (e.g. Mathys
#' et al. 2024).
#'
#' This function validates but never prepares the reference: \code{reference}
#' must already carry the dimensionality reduction named by
#' \code{reference_reduction} (e.g. from your own \code{Seurat::RunPCA()} or
#' \code{Seurat::SCTransform()} recipe) before calling this function. Baking a
#' generic prep recipe into the toolkit risks silently applying defaults that
#' don't suit a given reference's composition — that decision is left to you,
#' in the calling script. The query (\code{sce}) is likewise never
#' normalized here.
#'
#' @param sce A \code{SingleCellExperiment} (the query).
#' @param reference A \code{Seurat} object (already normalized/reduced), or a
#'   \code{SingleCellExperiment}/\code{SummarizedExperiment} that will be
#'   converted via \code{\link{TALARIA_to_seurat}} (its \code{reducedDims}
#'   carry over in that conversion, so the equivalent prep recipe can be run
#'   on the SCE beforehand instead).
#' @param reference_reduction Character. Name of the dimensionality reduction
#'   on \code{reference} to project the query into.  Default \code{"pca"}.
#'   Must already exist — not computed by this function.
#' @param dims Integer vector or \code{NULL}.  Which components of
#'   \code{reference_reduction} to use.  If \code{NULL} (default), uses every
#'   component already present (\code{1:ncol(Embeddings(reference,
#'   reference_reduction))}) rather than guessing a count.
#' @param labels_col Character. Column in the reference's metadata holding
#'   the labels to transfer.  Default \code{"label.main"}.
#' @param cluster_col Character or \code{NULL}.  If \code{NULL} (default),
#'   annotation is at the cell level.  If a \code{colData} column name is
#'   provided, the per-cell transferred labels are aggregated per cluster by
#'   majority vote (mean confidence score) and broadcast back to all cells in
#'   that cluster — same semantics as \code{cluster_col} in
#'   \code{\link{PANDORA_annotate_singler}}.
#' @param assay_name Character. Normalized-data assay on \code{sce} (and, if
#'   \code{reference} needs converting, on \code{reference}) used for anchor
#'   finding.  Default \code{"logcounts"}.
#' @param counts_assay Character. Raw-counts assay used only when
#'   \code{sce}/\code{reference} need converting via
#'   \code{\link{TALARIA_to_seurat}}.  Default \code{"counts"}.
#' @param min_score Numeric or \code{NULL}.  If set, (cluster- or cell-level)
#'   confidence scores below this threshold get \code{NA} in \code{label_col}
#'   instead of a forced low-confidence label — \code{TransferData()} has no
#'   built-in delta-score pruning like SingleR's, so thresholding is the only
#'   option here.  Default \code{NULL} (no thresholding).
#' @param label_col Character. \code{colData} column for the result.
#'   Default \code{"transfer_label"}.  The per-cell confidence score is always
#'   stored alongside it in \code{paste0(label_col, "_score")}, regardless of
#'   whether \code{min_score} is set — so it's available for QC/diagnostic use
#'   even on an unthresholded run.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} and
#'   \code{colData(sce)[[paste0(label_col, "_score")]]} populated.
#'
#' @seealso \code{\link{PANDORA_annotate_singler}},
#'   \code{\link{PANDORA_annotate_sctype}}, \code{\link{TALARIA_to_seurat}}
#'
#' @export
PANDORA_annotate_transfer <- function(sce,
                                       reference,
                                       reference_reduction = "pca",
                                       dims                 = NULL,
                                       labels_col           = "label.main",
                                       cluster_col          = NULL,
                                       assay_name           = "logcounts",
                                       counts_assay         = "counts",
                                       min_score            = NULL,
                                       label_col            = "transfer_label",
                                       verbose              = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required. ",
         "Install via: install.packages(\"Seurat\")", call. = FALSE)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Run PYRI_normalize() first.",
         call. = FALSE)

  # ── Resolve reference and query to Seurat objects (reuse TALARIA converters,
  # rather than duplicating SCE<->Seurat logic) ────────────────────────────────
  ref_sobj   <- .pandora_resolve_transfer_reference(reference, counts_assay,
                                                     assay_name)
  query_sobj <- TALARIA_to_seurat(sce, counts = counts_assay, data = assay_name)

  # ── Validate the reference is already prepared -- never auto-prepared ───────
  if (!reference_reduction %in% Seurat::Reductions(ref_sobj))
    stop(
      "Reduction '", reference_reduction, "' not found on reference ",
      "(available: ",
      paste(Seurat::Reductions(ref_sobj), collapse = ", "), "). ",
      "PANDORA_annotate_transfer() does not prepare the reference for you -- ",
      "run your own normalization / feature-selection / dimensionality-",
      "reduction recipe on it first, e.g.\n",
      "    Seurat::NormalizeData() -> Seurat::FindVariableFeatures() -> ",
      "Seurat::ScaleData() -> Seurat::RunPCA()\n",
      "  or\n",
      "    Seurat::SCTransform() -> Seurat::RunPCA()\n",
      "then pass the prepared object as reference=.",
      call. = FALSE
    )

  if (!labels_col %in% colnames(ref_sobj[[]]))
    stop("'", labels_col, "' not found in reference metadata. Available: ",
         paste(colnames(ref_sobj[[]]), collapse = ", "), call. = FALSE)

  # ── Resolve dims -- default to however many components the reference
  # already has, rather than guessing a count ─────────────────────────────────
  if (is.null(dims))
    dims <- seq_len(ncol(Seurat::Embeddings(ref_sobj, reference_reduction)))

  # ── Cluster vector (NULL = cell-level) ───────────────────────────────────────
  clusters <- NULL
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(colData(sce)))
      stop("cluster_col '", cluster_col, "' not found in colData(sce).",
           call. = FALSE)
    clusters <- as.character(colData(sce)[[cluster_col]])
  }

  # ── Run anchor-based label transfer ──────────────────────────────────────────
  anchors <- Seurat::FindTransferAnchors(
    reference           = ref_sobj,
    query               = query_sobj,
    dims                = dims,
    reference.reduction = reference_reduction
  )

  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata   = ref_sobj[[labels_col, drop = TRUE]],
    dims      = dims
  )

  predicted_id     <- stats::setNames(predictions$predicted.id, rownames(predictions))
  prediction_score <- stats::setNames(predictions$prediction.score.max, rownames(predictions))
  predicted_id     <- predicted_id[colnames(sce)]
  prediction_score <- prediction_score[colnames(sce)]

  # ── Cluster-level majority vote (broadcast back to cells) -- same semantics
  # as cluster_col handling in PANDORA_annotate_singler() ──────────────────────
  if (!is.null(clusters)) {
    unique_clusters <- sort(unique(clusters))

    cluster_label <- vapply(unique_clusters, function(cl) {
      votes <- table(predicted_id[clusters == cl])
      names(votes)[which.max(votes)]
    }, character(1))

    cluster_score <- vapply(unique_clusters, function(cl) {
      mean(prediction_score[clusters == cl], na.rm = TRUE)
    }, numeric(1))

    names(cluster_label) <- names(cluster_score) <- unique_clusters
    predicted_id     <- stats::setNames(cluster_label[clusters], colnames(sce))
    prediction_score <- stats::setNames(cluster_score[clusters], colnames(sce))
  }

  # ── min_score thresholding -- NA below threshold; no forced low-confidence
  # label, since TransferData() has no delta-score pruning like SingleR's ─────
  if (!is.null(min_score))
    predicted_id[prediction_score < min_score] <- NA

  score_col <- paste0(label_col, "_score")
  colData(sce)[[label_col]] <- predicted_id
  colData(sce)[[score_col]] <- prediction_score

  # ── Summary ──────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    ref_name    <- if (is.character(reference)) reference else class(reference)[1]
    n_annotated <- sum(!is.na(predicted_id))
    n_na        <- sum(is.na(predicted_id))
    n_types     <- length(unique(stats::na.omit(predicted_id)))
    mode_str    <- if (is.null(clusters)) "cell-level" else
      paste0("cluster-level (via '", cluster_col, "')")
    cat(sprintf(
      "── PANDORA: label transfer annotation %s\n  Reference  : %s (%s)\n  Dims used  : %d (1:%d)\n  Mode       : %s\n  Annotated  : %s cells  |  Below min_score (NA): %s\n  Cell types : %d unique types\n  Stored as  : colData(sce)[[\"%s\"]], colData(sce)[[\"%s\"]]\n%s\n",
      strrep("─", 19),
      ref_name, labels_col,
      length(dims), max(dims),
      mode_str,
      format(n_annotated, big.mark = ","),
      format(n_na, big.mark = ","),
      n_types,
      label_col, score_col,
      strrep("─", 56)
    ))
  }

  sce
}


# ==============================================================================
# Manual label assignment
# ==============================================================================

#' Manually assign cell type labels to clusters
#'
#' Maps cluster IDs to user-defined cell type labels using a named character
#' vector, storing the result as a new \code{colData} column.  Clusters
#' without a matching entry in \code{labels} receive \code{NA}.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param labels Named character vector.  Names are cluster IDs (must match
#'   levels in \code{cluster_col}); values are cell type labels.
#'   Example: \code{c("0" = "T cell", "1" = "B cell", "2" = "Myeloid")}.
#' @param cluster_col Character. \code{colData} column containing cluster IDs.
#'   Default \code{"cluster"}.
#' @param label_col Character. Name of the new \code{colData} column to create.
#'   Must be supplied explicitly; no default to prevent accidental overwrites.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return SCE with \code{colData(sce)[[label_col]]} populated.
#' @export
PANDORA_assign_labels <- function(sce,
                                   labels,
                                   cluster_col = "cluster",
                                   label_col,
                                   verbose     = TRUE) {

  if (!is.character(labels) || is.null(names(labels)))
    stop("'labels' must be a named character vector.  ",
         "Example: c('0' = 'T cell', '1' = 'B cell').", call. = FALSE)

  if (!cluster_col %in% names(colData(sce)))
    stop("cluster_col '", cluster_col, "' not found in colData(sce).",
         call. = FALSE)

  if (missing(label_col) || !nzchar(label_col))
    stop("'label_col' must be specified - name for the new colData column.",
         call. = FALSE)

  clusters        <- as.character(colData(sce)[[cluster_col]])
  unique_clusters <- sort(unique(clusters))

  unmatched <- setdiff(unique_clusters, names(labels))
  extra     <- setdiff(names(labels), unique_clusters)

  if (length(unmatched) > 0)
    warning("Clusters without labels (assigned NA): ",
            paste(unmatched, collapse = ", "), call. = FALSE)

  if (length(extra) > 0 && isTRUE(verbose))
    message("  Note: labels provided for clusters not present in SCE: ",
            paste(extra, collapse = ", "))

  cell_labels         <- labels[clusters]
  names(cell_labels)  <- colnames(sce)

  colData(sce)[[label_col]] <- cell_labels

  if (isTRUE(verbose)) {
    n_assigned <- sum(!is.na(cell_labels))
    n_types    <- length(unique(stats::na.omit(cell_labels)))
    cat(sprintf(
      "\u2500\u2500 PANDORA: assign_labels %s\n  Source    : %s  \u2192  %s\n  Assigned  : %s cells  |  NA: %s\n  Cell types: %d unique\n  Stored as : colData(sce)[[\"%s\"]]\n%s\n",
      strrep("\u2500", 32),
      cluster_col, label_col,
      format(n_assigned, big.mark = ","),
      format(sum(is.na(cell_labels)), big.mark = ","),
      n_types,
      label_col,
      strrep("\u2500", 56)
    ))
  }

  sce
}


#' Summarise cell type label composition
#'
#' Returns overall label frequencies and optionally a breakdown by cluster
#' and/or sample.  Prints a formatted composition table to the console.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param label_col Character. \code{colData} column with cell type labels
#'   (output of any PANDORA annotation function or
#'   \code{\link{PANDORA_assign_labels}}).
#' @param cluster_col Character or \code{NULL}.  If provided, also tabulates
#'   label composition within each cluster. Default \code{NULL}.
#' @param sample_col Character or \code{NULL}.  If provided, also tabulates
#'   label composition within each sample. Default \code{NULL}.
#' @param verbose Logical. Print the composition table. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#'   \itemize{
#'     \item \code{$overall} — data.frame of overall label frequencies.
#'     \item \code{$by_cluster} — data.frame or \code{NULL}.
#'     \item \code{$by_sample} — data.frame or \code{NULL}.
#'   }
#' @export
PANDORA_summarise_labels <- function(sce,
                                      label_col,
                                      cluster_col = NULL,
                                      sample_col  = NULL,
                                      verbose     = TRUE) {

  if (!label_col %in% names(colData(sce)))
    stop("label_col '", label_col, "' not found in colData(sce).",
         call. = FALSE)

  labels <- as.character(colData(sce)[[label_col]])

  tbl    <- table(labels)
  out_df <- data.frame(
    label      = names(tbl),
    n_cells    = as.integer(tbl),
    proportion = as.numeric(tbl) / sum(tbl),
    stringsAsFactors = FALSE
  )
  out_df <- out_df[order(-out_df$n_cells), , drop = FALSE]
  rownames(out_df) <- NULL

  # Optional cluster breakdown
  cluster_df <- NULL
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(colData(sce)))
      stop("cluster_col '", cluster_col, "' not found in colData(sce).",
           call. = FALSE)
    clusters   <- as.character(colData(sce)[[cluster_col]])
    ct         <- table(cluster = clusters, label = labels)
    cluster_df <- as.data.frame(ct, stringsAsFactors = FALSE)
    cluster_df$proportion <- cluster_df$Freq /
      ave(cluster_df$Freq, cluster_df$cluster, FUN = sum)
    cluster_df <- cluster_df[cluster_df$Freq > 0L, , drop = FALSE]
    cluster_df <- cluster_df[order(cluster_df$cluster, -cluster_df$Freq), ,
                              drop = FALSE]
    rownames(cluster_df) <- NULL
  }

  # Optional sample breakdown
  sample_df <- NULL
  if (!is.null(sample_col)) {
    if (!sample_col %in% names(colData(sce)))
      stop("sample_col '", sample_col, "' not found in colData(sce).",
           call. = FALSE)
    samples   <- as.character(colData(sce)[[sample_col]])
    st        <- table(sample = samples, label = labels)
    sample_df <- as.data.frame(st, stringsAsFactors = FALSE)
    sample_df$proportion <- sample_df$Freq /
      ave(sample_df$Freq, sample_df$sample, FUN = sum)
    sample_df <- sample_df[sample_df$Freq > 0L, , drop = FALSE]
    sample_df <- sample_df[order(sample_df$sample, -sample_df$Freq), ,
                            drop = FALSE]
    rownames(sample_df) <- NULL
  }

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 PANDORA: summarise_labels %s\n  Label col : %s\n  Cell types: %d unique (%s cells total)\n\n",
      strrep("\u2500", 29),
      label_col,
      nrow(out_df),
      format(ncol(sce), big.mark = ",")
    ))
    cat("  Overall composition:\n")
    show_n <- min(nrow(out_df), 15L)
    for (i in seq_len(show_n)) {
      cat(sprintf("    %-32s %6s cells  (%4.1f%%)\n",
                  out_df$label[i],
                  format(out_df$n_cells[i], big.mark = ","),
                  out_df$proportion[i] * 100))
    }
    if (nrow(out_df) > 15L)
      cat("    ... and", nrow(out_df) - 15L, "more.\n")
    cat(strrep("\u2500", 56), "\n")
  }

  invisible(list(overall = out_df, by_cluster = cluster_df,
                 by_sample = sample_df))
}


# ==============================================================================
# PANDORA_convert_orthologs
# ==============================================================================

#' Convert gene symbols to orthologs
#'
#' Converts a character vector of gene symbols from one species to another
#' using \pkg{orthogene}.  By default the returned vector is the same length
#' as the input, with \code{NA} for genes that could not be mapped.  Set
#' \code{drop_na = TRUE} to silently remove unmapped entries (used internally
#' by \code{\link{PANDORA_annotate_sctype}}).
#'
#' @param genes Character vector of gene symbols to convert.
#' @param from Character. Source species.  Default \code{"human"}.
#' @param to Character. Target species (e.g. \code{"mouse"},
#'   \code{"Mus musculus"}, \code{"zebrafish"}).  Any species string accepted
#'   by \code{orthogene::convert_orthologs()} is valid.
#' @param method Character. Ortholog database passed to
#'   \code{orthogene::convert_orthologs()}.  Default \code{"homologene"}.
#'   Use \code{"gprofiler"} for broader species coverage.
#' @param drop_na Logical. If \code{FALSE} (default), unmapped genes are
#'   returned as \code{NA} (output length equals input length).  If
#'   \code{TRUE}, unmapped genes are dropped from the result.
#' @param verbose Logical. Print mapping statistics.  Default \code{TRUE}.
#'
#' @return A named character vector.  Names are the input gene symbols; values
#'   are the corresponding ortholog symbols (or \code{NA} when unmapped and
#'   \code{drop_na = FALSE}).
#'
#' @export
PANDORA_convert_orthologs <- function(genes,
                                       from     = "human",
                                       to,
                                       method   = "homologene",
                                       drop_na  = FALSE,
                                       verbose  = TRUE) {

  if (!requireNamespace("orthogene", quietly = TRUE))
    stop("Package 'orthogene' is required for ortholog conversion.\n",
         "  Install via: BiocManager::install(\"orthogene\")", call. = FALSE)

  genes        <- as.character(genes)
  genes_unique <- unique(genes)
  n_input      <- length(genes_unique)

  ortho_df <- suppressMessages(
    orthogene::convert_orthologs(
      gene_df        = genes_unique,
      input_species  = from,
      output_species = to,
      method         = method,
      drop_nonorths  = TRUE,
      verbose        = FALSE
    )
  )

  if (!is.data.frame(ortho_df) || nrow(ortho_df) == 0) {
    warning("orthogene returned no orthologs for '", to,
            "' using method '", method, "'. ",
            "Try method = \"gprofiler\" for broader species coverage.",
            call. = FALSE)
    result <- stats::setNames(rep(NA_character_, length(genes)), genes)
    if (isTRUE(drop_na)) result <- result[!is.na(result)]
    return(result)
  }

  if (!"input_gene" %in% colnames(ortho_df))
    stop("Unexpected orthogene output format: 'input_gene' column not found. ",
         "This may indicate an orthogene version incompatibility.", call. = FALSE)

  # Named map: input symbol -> ortholog symbol (rownames = target species)
  gene_map <- stats::setNames(rownames(ortho_df), ortho_df[["input_gene"]])
  n_mapped <- length(gene_map)

  # Build output: same length as input, NA for unmapped genes
  result        <- gene_map[genes]
  names(result) <- genes

  if (isTRUE(verbose)) {
    pct <- round(100 * n_mapped / n_input, 1)
    cat("[PANDORA] Ortholog conversion:", from, "->", to, "\n")
    cat("    Method   :", method, "\n")
    cat("    Mapped   :", n_mapped, "/", n_input,
        paste0("(", format(pct, nsmall = 1), "%)\n"))
    if (isTRUE(drop_na) && n_mapped < n_input)
      cat("    Unmapped :", n_input - n_mapped, "genes dropped\n")
    else if (n_mapped < n_input)
      cat("    Unmapped :", n_input - n_mapped, "genes set to NA\n")
    if (pct < 50)
      cat("    WARNING  : Conversion rate <50% - interpret results with caution.\n")
  }

  if (isTRUE(drop_na))
    result <- result[!is.na(result)]

  result
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


# Resolve an annotation reference (Seurat, SingleCellExperiment, or generic
# SummarizedExperiment) to a Seurat object for PANDORA_annotate_transfer().
# Reuses TALARIA_to_seurat() rather than duplicating SCE<->Seurat conversion
# logic. A Seurat object is passed through unchanged -- trusted as already
# prepared (normalized + reduced) by the caller.
.pandora_resolve_transfer_reference <- function(reference, counts_assay, data_assay) {
  if (inherits(reference, "Seurat")) return(reference)

  if (inherits(reference, "SummarizedExperiment") &&
      !inherits(reference, "SingleCellExperiment"))
    reference <- as(reference, "SingleCellExperiment")

  if (!inherits(reference, "SingleCellExperiment"))
    stop("'reference' must be a Seurat, SingleCellExperiment, or ",
         "SummarizedExperiment object.", call. = FALSE)

  if (!data_assay %in% assayNames(reference))
    stop("Assay '", data_assay, "' not found in reference. Available: ",
         paste(assayNames(reference), collapse = ", "), call. = FALSE)

  TALARIA_to_seurat(reference, counts = counts_assay, data = data_assay)
}


# Remap scType human marker lists to target species orthologs.
# Delegates conversion to PANDORA_convert_orthologs(); reporting is handled
# by the caller (PANDORA_annotate_sctype).
#
# Returns a list:
#   $pos_markers  — remapped positive marker lists (target species symbols)
#   $neg_markers  — remapped negative marker lists (target species symbols)
#   $n_total      — total unique human markers queried
#   $n_mapped     — number successfully mapped to target species
#   $pct          — mapping rate as formatted string ("87.1")
#   $pct_num      — mapping rate as numeric (for threshold checks)
.pandora_convert_sctype_markers <- function(pos_markers, neg_markers,
                                            target_species, method, verbose) {
  all_human <- unique(c(unlist(pos_markers), unlist(neg_markers)))
  n_total   <- length(all_human)

  gene_map <- PANDORA_convert_orthologs(
    genes    = all_human,
    from     = "human",
    to       = target_species,
    method   = method,
    drop_na  = TRUE,
    verbose  = FALSE   # reporting handled by PANDORA_annotate_sctype
  )

  n_mapped <- length(gene_map)
  pct_num  <- round(100 * n_mapped / n_total, 1)

  # Remap marker lists — genes without an ortholog are silently dropped
  remap <- function(markers) {
    lapply(markers, function(gs) unname(stats::na.omit(gene_map[gs])))
  }

  list(
    pos_markers = remap(pos_markers),
    neg_markers = remap(neg_markers),
    n_total     = n_total,
    n_mapped    = n_mapped,
    pct         = format(pct_num, nsmall = 1),
    pct_num     = pct_num
  )
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
