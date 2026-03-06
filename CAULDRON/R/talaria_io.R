# ==============================================================================
# TALARIA - Import, Export & Format Conversion
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The winged sandals forged by Hephaestus for Hermes. TALARIA carries data
# in and out of the toolkit — pure transport, bridging external formats and
# the internal SingleCellExperiment object model.
#
# Core responsibilities:
#   - Load 10x Genomics MEX directories (barcodes/features/matrix)
#   - Load H5 and H5AD (AnnData) files
#   - Load Loom files
#   - Convert between SingleCellExperiment and Seurat objects
#   - Export results (count matrices, metadata, reduced dims)
#   - Report generation
#
# All functions prefixed: TALARIA_
# ==============================================================================


#' Load an AnnData H5AD file into a SingleCellExperiment
#'
#' Reads a \code{.h5ad} file produced by Python's AnnData / Scanpy and returns
#' a \code{SingleCellExperiment} object.  Conversion is handled by
#' \pkg{zellkonverter}, which maps:
#'
#' \itemize{
#'   \item \code{X} (and \code{raw.X} if present) → assay(s) of the SCE
#'   \item \code{layers} → additional named assays
#'   \item \code{obs} columns → \code{colData}
#'   \item \code{var} columns → \code{rowData}
#'   \item \code{obsm} entries (e.g. \code{X_pca}, \code{X_umap}) →
#'         \code{reducedDims}
#' }
#'
#' @param path Character. Path to the \code{.h5ad} file.
#' @param use_hdf5 Logical. If \code{TRUE}, count matrices are kept on-disk as
#'   \code{HDF5Array} objects rather than loaded into RAM.  Useful for very
#'   large datasets. Default \code{FALSE}.
#' @param reader Character. \code{"R"} (default, pure R via \pkg{rhdf5}) or
#'   \code{"python"} (requires a working Python / basilisk environment).
#' @param backend Character. \code{"memory"} (default) loads the counts assay
#'   into RAM as a standard sparse matrix.  \code{"bpcells"} writes the counts
#'   assay to an on-disk BPCells store and replaces the in-memory matrix with a
#'   lazy \code{BPCells::IterableMatrix} that streams from disk.  Requires
#'   \pkg{BPCells}.
#' @param bpcells_dir Character. Path for the BPCells on-disk store directory.
#'   Required (and created) when \code{backend = "bpcells"}; ignored otherwise.
#'   Choose a permanent location — the path is embedded in the returned SCE and
#'   must remain accessible after saving/reloading the object.
#' @param verbose Logical. Print a summary of the loaded object. Default
#'   \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment} object.
#' @export
TALARIA_load_h5ad <- function(path,
                               use_hdf5    = FALSE,
                               reader      = c("R", "python"),
                               backend     = c("memory", "bpcells"),
                               bpcells_dir = NULL,
                               verbose     = TRUE) {

  reader  <- match.arg(reader)
  backend <- match.arg(backend)

  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)

  if (!requireNamespace("zellkonverter", quietly = TRUE))
    stop("Package 'zellkonverter' is required. ",
         "Install via: BiocManager::install(\"zellkonverter\")", call. = FALSE)

  sce <- zellkonverter::readH5AD(path,
                                  use_hdf5 = use_hdf5,
                                  reader   = reader,
                                  verbose  = FALSE)

  sce <- .talaria_normalise_assay_names(sce)

  if (backend == "bpcells")
    sce <- .talaria_to_bpcells_backend(sce, bpcells_dir)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(sce, path)

  sce
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Ensure the raw count assay is named "counts" — the Bioconductor convention
# expected by scDblFinder, scran, scater, etc.
#
# zellkonverter behaviour:
#   h5ad with raw.X  →  raw.X becomes "counts", X becomes "logcounts"  (ok)
#   h5ad without raw →  X stays as "X"  (needs normalisation)
#
# If "X" is present and appears to contain raw integer counts, rename it.
# If "X" is present but appears to be log-normalised (continuous values),
# stop with a descriptive error — raw counts are not available.
.talaria_normalise_assay_names <- function(sce) {
  an <- assayNames(sce)

  if ("counts" %in% an)
    return(sce)   # already correct

  if (!"X" %in% an)
    stop("No 'counts' or 'X' assay found in the loaded object. ",
         "Available assays: ", paste(an, collapse = ", "), ".",
         call. = FALSE)

  if (.talaria_is_integer_matrix(assay(sce, "X"))) {
    assayNames(sce)[assayNames(sce) == "X"] <- "counts"
    cat("Note: assay 'X' contained raw integer counts and has been renamed to 'counts' (Bioconductor convention).\n")
  } else {
    stop(
      "Assay 'X' does not appear to contain raw integer counts — values are ",
      "continuous, suggesting the matrix has already been log-normalised.\n",
      "Raw counts are required for downstream QC and doublet detection.\n",
      "Re-export the h5ad with raw counts stored in 'adata.raw' so they are ",
      "read as the 'counts' assay, or ensure 'adata.X' holds the unnormalised ",
      "count matrix before exporting.",
      call. = FALSE
    )
  }

  sce
}


# Sample a small block of non-zero values and test whether they are all
# non-negative integers.  Using a block slice (rather than @x) works for
# both in-memory sparse matrices and on-disk HDF5Array objects.
.talaria_is_integer_matrix <- function(mat, n_rows = 200, n_cols = 50) {
  nr   <- min(nrow(mat), n_rows)
  nc   <- min(ncol(mat), n_cols)
  vals <- as.vector(as.matrix(mat[seq_len(nr), seq_len(nc)]))
  vals <- vals[vals != 0]

  if (length(vals) == 0)
    return(TRUE)   # empty block — assume integers (e.g. very sparse data)

  all(vals >= 0) && all(abs(vals - round(vals)) < 1e-6)
}


# ==============================================================================
# Loading functions
# ==============================================================================


#' Load a 10x Genomics MEX directory into a SingleCellExperiment
#'
#' Reads the three-file MEX format produced by Cell Ranger
#' (\code{barcodes.tsv.gz}, \code{features.tsv.gz}, \code{matrix.mtx.gz})
#' using \pkg{DropletUtils}.
#'
#' @param path Character. Path to the MEX directory (the folder containing the
#'   three Cell Ranger output files).
#' @param sample_names Character. Sample name(s) assigned to the
#'   \code{colData} \code{Sample} column.  Defaults to the directory name.
#' @param backend Character. \code{"memory"} (default) or \code{"bpcells"}.
#'   See \code{\link{TALARIA_load_h5ad}} for details.
#' @param bpcells_dir Character. Path for the BPCells on-disk store directory.
#'   Required when \code{backend = "bpcells"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment} with a \code{"counts"} assay.
#' @export
TALARIA_load_10x <- function(path,
                              sample_names = NULL,
                              backend      = c("memory", "bpcells"),
                              bpcells_dir  = NULL,
                              verbose      = TRUE) {

  backend <- match.arg(backend)

  if (!dir.exists(path))
    stop("Directory not found: ", path, call. = FALSE)

  if (!requireNamespace("DropletUtils", quietly = TRUE))
    stop("Package 'DropletUtils' is required. ",
         "Install via: BiocManager::install(\"DropletUtils\")", call. = FALSE)

  samples <- if (!is.null(sample_names)) sample_names else basename(path)

  sce <- DropletUtils::read10xCounts(path, sample.names = samples)
  sce <- .talaria_normalise_assay_names(sce)

  if (backend == "bpcells")
    sce <- .talaria_to_bpcells_backend(sce, bpcells_dir)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(sce, path)

  sce
}


#' Load a 10x Genomics HDF5 file into a SingleCellExperiment
#'
#' Reads a \code{.h5} file produced by Cell Ranger using \pkg{DropletUtils}.
#'
#' @param path Character. Path to the \code{.h5} file.
#' @param sample_names Character. Sample name assigned to the \code{colData}
#'   \code{Sample} column.  Defaults to the file name without extension.
#' @param backend Character. \code{"memory"} (default) or \code{"bpcells"}.
#'   See \code{\link{TALARIA_load_h5ad}} for details.
#' @param bpcells_dir Character. Path for the BPCells on-disk store directory.
#'   Required when \code{backend = "bpcells"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment} with a \code{"counts"} assay.
#' @export
TALARIA_load_h5 <- function(path,
                             sample_names = NULL,
                             backend      = c("memory", "bpcells"),
                             bpcells_dir  = NULL,
                             verbose      = TRUE) {

  backend <- match.arg(backend)

  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)

  if (!requireNamespace("DropletUtils", quietly = TRUE))
    stop("Package 'DropletUtils' is required. ",
         "Install via: BiocManager::install(\"DropletUtils\")", call. = FALSE)

  samples <- if (!is.null(sample_names)) sample_names else
    tools::file_path_sans_ext(basename(path))

  sce <- DropletUtils::read10xCounts(path, sample.names = samples)
  sce <- .talaria_normalise_assay_names(sce)

  if (backend == "bpcells")
    sce <- .talaria_to_bpcells_backend(sce, bpcells_dir)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(sce, path)

  sce
}


#' Load a Loom file into a SingleCellExperiment
#'
#' Reads a \code{.loom} file using \pkg{LoomExperiment} and returns a
#' \code{SingleCellExperiment}.  Loom is a legacy format; prefer \code{.h5ad}
#' for new projects.
#'
#' @param path Character. Path to the \code{.loom} file.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment}.
#' @export
TALARIA_load_loom <- function(path, verbose = TRUE) {

  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)

  if (!requireNamespace("LoomExperiment", quietly = TRUE))
    stop("Package 'LoomExperiment' is required. ",
         "Install via: BiocManager::install(\"LoomExperiment\")", call. = FALSE)

  # import() with type="SingleCellLoomExperiment" returns an SCE subclass
  sce <- LoomExperiment::import(path, type = "SingleCellLoomExperiment")
  sce <- .talaria_normalise_assay_names(sce)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(sce, path)

  sce
}


# ==============================================================================
# BPCells backend
# ==============================================================================


#' Convert a SingleCellExperiment counts assay to an on-disk BPCells matrix
#'
#' Writes the specified assay to a \pkg{BPCells} bitpacked directory store and
#' replaces the in-memory matrix with a lazy \code{IterableMatrix} that streams
#' from disk.  Most downstream CAULDRON operations (QC metrics, log-normalisation,
#' HVG selection, PCA via irlba) are natively compatible with BPCells and will
#' not load the full matrix into RAM.
#'
#' \strong{Note on \code{AEGIS_detect_doublets}:} \pkg{scDblFinder} may
#' internally coerce the BPCells matrix to a standard sparse matrix, loading
#' the full counts into memory for doublet detection.  This is a limitation of
#' the scDblFinder implementation.  For very large datasets consider running
#' doublet detection on a per-sample subset.
#'
#' \strong{Persistence:} The \code{bpcells_dir} path is stored in
#' \code{metadata(sce)$bpcells_dir}.  When saving the SCE with
#' \code{\link{TALARIA_save_sce}} and reloading with
#' \code{\link{TALARIA_load_sce}}, the lazy matrix is automatically
#' re-connected — provided \code{bpcells_dir} still exists at the same path.
#' Use a permanent, session-independent path.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param bpcells_dir Character. Directory path for the BPCells on-disk store.
#'   Created (recursively) if it does not exist.  Must be a permanent path —
#'   temporary directories will break the object after session restart.
#' @param assay_name Character. Assay to convert.  Default \code{"counts"}.
#' @param verbose Logical. Print a summary message.  Default \code{TRUE}.
#'
#' @return The input \code{SingleCellExperiment} with the specified assay
#'   replaced by a \pkg{BPCells} \code{IterableMatrix}.
#'   \code{metadata(sce)$bpcells_dir} is set to \code{bpcells_dir}.
#' @export
TALARIA_to_bpcells <- function(sce,
                                bpcells_dir,
                                assay_name = "counts",
                                verbose    = TRUE) {

  if (!requireNamespace("BPCells", quietly = TRUE))
    stop("Package 'BPCells' is required. ",
         "Install via: devtools::install_github(\"bnprks/BPCells/r\")",
         call. = FALSE)

  if (!inherits(sce, "SingleCellExperiment"))
    stop("sce must be a SingleCellExperiment object.", call. = FALSE)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(assayNames(sce), collapse = ", "), call. = FALSE)

  dir.create(dirname(bpcells_dir), recursive = TRUE, showWarnings = FALSE)

  mat    <- assay(sce, assay_name)
  BPCells::write_matrix_dir(mat, bpcells_dir)
  bp_mat <- BPCells::open_matrix_dir(bpcells_dir)

  assay(sce, assay_name)       <- bp_mat
  metadata(sce)$bpcells_dir    <- bpcells_dir

  if (isTRUE(verbose))
    cat("BPCells: '", assay_name, "' assay written to: ", bpcells_dir,
        "\n  Class: ", class(bp_mat),
        "\n  Dims : ", nrow(sce), " genes x ", ncol(sce), " cells\n", sep = "")

  sce
}


# Internal dispatcher called by the three loaders when backend = "bpcells".
# Validates bpcells_dir and delegates to TALARIA_to_bpcells().
.talaria_to_bpcells_backend <- function(sce, bpcells_dir, assay_name = "counts") {
  if (is.null(bpcells_dir))
    stop("bpcells_dir must be provided when backend = \"bpcells\".",
         call. = FALSE)
  TALARIA_to_bpcells(sce, bpcells_dir, assay_name = assay_name, verbose = TRUE)
}


# ==============================================================================
# Conversion functions
# ==============================================================================


#' Convert a SingleCellExperiment to a Seurat object
#'
#' Wraps \code{Seurat::as.Seurat()} to convert a \code{SingleCellExperiment}
#' to a \code{Seurat} object, transferring counts, normalised data,
#' dimensionality reductions, and cell metadata.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param counts Character. Assay to use as the raw counts slot in Seurat.
#'   Default \code{"counts"}.
#' @param data Character. Assay to use as the normalised data slot.
#'   Default \code{"logcounts"} if present, otherwise \code{NULL}.
#'
#' @return A \code{Seurat} object.
#' @export
TALARIA_to_seurat <- function(sce, counts = "counts", data = NULL) {

  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required. ",
         "Install via: install.packages(\"Seurat\")", call. = FALSE)

  if (!counts %in% assayNames(sce))
    stop("Assay '", counts, "' not found. Available: ",
         paste(assayNames(sce), collapse = ", "), call. = FALSE)

  # Default data slot to logcounts if present
  if (is.null(data))
    data <- if ("logcounts" %in% assayNames(sce)) "logcounts" else NULL

  sobj <- Seurat::as.Seurat(sce, counts = counts, data = data)

  cat("Converted SCE (", ncol(sce), " cells, ", nrow(sce),
      " genes) to Seurat object.\n", sep = "")
  sobj
}


#' Convert a Seurat object to a SingleCellExperiment
#'
#' Wraps \code{Seurat::as.SingleCellExperiment()} to convert a \code{Seurat}
#' object to a \code{SingleCellExperiment}, preserving assays, embeddings,
#' and cell metadata.
#'
#' @param sobj A \code{Seurat} object.
#' @param assay Character. Which Seurat assay to convert.  Default
#'   \code{"RNA"}.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment}.
#' @export
TALARIA_from_seurat <- function(sobj, assay = "RNA", verbose = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required. ",
         "Install via: install.packages(\"Seurat\")", call. = FALSE)

  if (!inherits(sobj, "Seurat"))
    stop("sobj must be a Seurat object.", call. = FALSE)

  sce <- Seurat::as.SingleCellExperiment(sobj, assay = assay)
  sce <- .talaria_normalise_assay_names(sce)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(sce, path = paste0("<Seurat:", assay, ">"))

  sce
}


# ==============================================================================
# Export / persistence functions
# ==============================================================================


#' Export a count matrix to disk
#'
#' Writes the specified assay to disk in CSV or 10x MEX format.
#'
#' \strong{Note:} CSV export calls \code{as.matrix()} on the assay, which
#' loads the full matrix into RAM.  For large datasets prefer \code{"mex"}
#' (requires \pkg{DropletUtils}).
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param output_dir Character. Directory to write output into.  Created if it
#'   does not exist.
#' @param format Character. \code{"csv"} (default) or \code{"mex"}.
#' @param assay_name Character. Assay to export.  Default \code{"counts"}.
#' @param verbose Logical. Print confirmation. Default \code{TRUE}.
#'
#' @return Output path(s), invisibly.
#' @export
TALARIA_export_matrix <- function(sce,
                                   output_dir,
                                   format     = c("csv", "mex"),
                                   assay_name = "counts",
                                   verbose    = TRUE) {

  format <- match.arg(format)

  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(assayNames(sce), collapse = ", "), call. = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (format == "csv") {
    n_cells <- ncol(sce)
    n_genes <- nrow(sce)
    if (n_cells * n_genes > 5e7)
      warning("Matrix has ", format(n_cells * n_genes, big.mark = ","),
              " elements — CSV export will be large and slow. ",
              "Consider format = \"mex\".", call. = FALSE)

    out_path <- file.path(output_dir, paste0(assay_name, "_matrix.csv"))
    mat      <- as.matrix(assay(sce, assay_name))
    write.csv(mat, out_path)

    if (isTRUE(verbose))
      cat("Matrix written to: ", out_path, "\n", sep = "")
    invisible(out_path)

  } else {
    if (!requireNamespace("DropletUtils", quietly = TRUE))
      stop("Package 'DropletUtils' is required for MEX export. ",
           "Install via: BiocManager::install(\"DropletUtils\")", call. = FALSE)

    DropletUtils::write10xCounts(output_dir, assay(sce, assay_name),
                                  gene.id   = rownames(sce),
                                  gene.symbol = rownames(sce),
                                  barcodes  = colnames(sce),
                                  overwrite = TRUE)

    if (isTRUE(verbose))
      cat("MEX files written to: ", output_dir, "\n", sep = "")
    invisible(output_dir)
  }
}


#' Export cell metadata to CSV
#'
#' Writes \code{colData(sce)} to a CSV file, with cell barcodes as the first
#' column.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param output_path Character. Full path for the output CSV file.
#' @param verbose Logical. Print confirmation. Default \code{TRUE}.
#'
#' @return \code{output_path}, invisibly.
#' @export
TALARIA_export_metadata <- function(sce, output_path, verbose = TRUE) {

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  meta           <- as.data.frame(colData(sce))
  meta           <- cbind(cell_barcode = rownames(meta), meta)
  rownames(meta) <- NULL

  write.csv(meta, output_path, row.names = FALSE)

  if (isTRUE(verbose))
    cat("Metadata (", nrow(meta), " cells, ", ncol(meta) - 1L,
        " columns) written to: ", output_path, "\n", sep = "")
  invisible(output_path)
}


#' Save a SingleCellExperiment object to disk
#'
#' Serialises the SCE to an RDS file.  Use \code{\link{TALARIA_load_sce}} to
#' reload.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param output_path Character. Full path for the \code{.rds} file.
#' @param verbose Logical. Print confirmation. Default \code{TRUE}.
#'
#' @return \code{output_path}, invisibly.
#' @export
TALARIA_save_sce <- function(sce, output_path, verbose = TRUE) {

  if (!inherits(sce, "SingleCellExperiment"))
    stop("sce must be a SingleCellExperiment object.", call. = FALSE)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sce, output_path)

  if (isTRUE(verbose))
    cat("SCE (", ncol(sce), " cells, ", nrow(sce),
        " genes) saved to: ", output_path, "\n", sep = "")
  invisible(output_path)
}


#' Load a SingleCellExperiment object from disk
#'
#' Reads an RDS file previously saved with \code{\link{TALARIA_save_sce}} and
#' validates that it contains a \code{SingleCellExperiment}.
#'
#' @param path Character. Path to the \code{.rds} file.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A \code{SingleCellExperiment} object.
#' @export
TALARIA_load_sce <- function(path, verbose = TRUE) {

  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)

  obj <- readRDS(path)

  if (!inherits(obj, "SingleCellExperiment"))
    stop("Object loaded from '", path, "' is not a SingleCellExperiment ",
         "(class: ", paste(class(obj), collapse = ", "), ").", call. = FALSE)

  if (isTRUE(verbose))
    .talaria_print_sce_summary(obj, path)

  obj
}


# ==============================================================================
# Internal helpers
# ==============================================================================


.talaria_print_sce_summary <- function(sce, path) {
  bar <- strrep("\u2500", 56)

  n_cells    <- ncol(sce)
  n_genes    <- nrow(sce)
  assays     <- assayNames(sce)
  red_dims   <- reducedDimNames(sce)
  col_cols   <- names(colData(sce))
  n_col_cols <- length(col_cols)

  # Truncate long colData column lists
  col_preview <- if (n_col_cols <= 8) {
    paste(col_cols, collapse = ", ")
  } else {
    paste0(paste(head(col_cols, 8), collapse = ", "),
           " ... (", n_col_cols, " total)")
  }

  cat("\u2500\u2500 TALARIA: H5AD loaded ", bar, "\n", sep = "")
  cat(sprintf("  File       : %s\n", basename(path)))
  cat(sprintf("  Cells      : %s\n", format(n_cells, big.mark = ",")))
  cat(sprintf("  Genes      : %s\n", format(n_genes, big.mark = ",")))
  cat(sprintf("  Assays     : %s\n",
                  if (length(assays)) paste(assays, collapse = ", ") else "(none)"))
  cat(sprintf("  Reductions : %s\n",
                  if (length(red_dims)) paste(red_dims, collapse = ", ") else "(none)"))
  cat(sprintf("  Metadata   : %s\n", col_preview))
  cat(bar, "\n", sep = "")
}
