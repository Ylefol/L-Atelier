#' Build a Bowtie2 genome index
#'
#' Wraps \code{Rbowtie2::bowtie2_build()} to create a Bowtie2 index from one
#' or more reference FASTA files.  The index only needs to be built once per
#' reference genome; the resulting \code{index_name} path is passed directly
#' to \code{\link{HORIZON_run_bowtie2}}.
#'
#' @param reference Character vector. Path(s) to reference FASTA file(s).
#'   Multiple FASFTAs are concatenated by bowtie2-build.
#' @param index_dir Character. Directory in which to write the index files.
#'   Created automatically if it does not exist.
#' @param index_name Character or \code{NULL}. Basename for the index files
#'   (without extension).  Defaults to the filename of
#'   \code{reference[1]} without its extension.
#' @param threads Integer. Number of threads passed to \code{--threads}.
#'   Default 4.
#' @param force Logical. If \code{TRUE}, rebuild the index even when it
#'   already exists.  Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{Rbowtie2::bowtie2_build}.
#'
#' @return Character. The full index basename path (i.e.
#'   \code{<index_dir>/<index_name>}), invisibly.  Pass this to
#'   \code{index} in \code{\link{HORIZON_run_bowtie2}}.
#' @export
HORIZON_build_bowtie2_index <- function(reference,
                                         index_dir,
                                         index_name = NULL,
                                         threads    = 4L,
                                         force      = FALSE,
                                         ...) {

  if (!requireNamespace("Rbowtie2", quietly = TRUE))
    stop("Package 'Rbowtie2' is required. Install via BiocManager::install('Rbowtie2').",
         call. = FALSE)

  if (length(reference) == 0 || !all(file.exists(reference)))
    stop("One or more reference FASTA files not found.", call. = FALSE)

  # Default index_name from first reference
  if (is.null(index_name))
    index_name <- tools::file_path_sans_ext(basename(reference[1]))

  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  full_index <- file.path(index_dir, index_name)

  # Check for existing index
  if (!isTRUE(force) && file.exists(paste0(full_index, ".1.bt2"))) {
    cat("Bowtie2 index already exists at: ", full_index,
            "\n  Use force=TRUE to rebuild.")
    return(invisible(full_index))
  }

  cat("Building Bowtie2 index: ", full_index,
          "\n  Index type: small (32-bit offsets, supports genomes up to ~4 GB)",
          "\n  For genomes exceeding 4 GB pass \"--large-index\" via ...")
  Rbowtie2::bowtie2_build(
    references = reference,
    bt2Index   = full_index,
    paste0("--threads ", as.integer(threads)),
    ...
  )

  cat("Bowtie2 index complete: ", full_index)
  invisible(full_index)
}
