#' Build a Rsubread genome index
#'
#' Wrapper around \code{\link[Rsubread]{buildindex}} for constructing the genome
#' index required for alignment. This is a one-time operation per reference
#' genome; the index can be stored on any accessible drive (including external
#' storage).
#'
#' RAM note: index building is the most RAM-intensive step. The \code{memory}
#' parameter directly controls how much RAM Rsubread uses during construction.
#' With \code{index_split = TRUE} (default), the resulting index is stored in
#' chunks, which significantly reduces per-alignment RAM usage at the cost of
#' slightly slower alignment.
#'
#' @param basename Character. Path prefix for the index files
#'   (e.g., \code{"/data/refs/hg38/hg38_subread"}). The directory must already
#'   exist.
#' @param reference Character. Path to the reference genome FASTA file.
#' @param memory Integer. RAM to use for index building in MB. Default 8000
#'   (~8 GB). Reduce on memory-constrained machines.
#' @param index_split Logical. Split the index into chunks to reduce RAM usage
#'   during alignment. Recommended for large genomes. Default TRUE.
#' @param gapped Logical. Build a gapped (splice-aware) index suitable for
#'   RNA-seq alignment. Default TRUE.
#' @param force Logical. If TRUE, rebuild the index even if one already exists
#'   at \code{basename}. Default FALSE.
#' @param ... Additional arguments passed to \code{\link[Rsubread]{buildindex}}.
#'
#' @return Character. The index basename path, invisibly.
#' @export
HORIZON_build_index <- function(basename,
                                reference,
                                memory      = 8000,
                                index_split = TRUE,
                                gapped      = TRUE,
                                force       = FALSE,
                                ...) {
  if (!file.exists(reference)) stop("Reference FASTA not found: ", reference)

  index_dir <- dirname(basename)
  if (!dir.exists(index_dir)) {
    stop("Index output directory does not exist: ", index_dir,
         "\nCreate it first or check the path.")
  }

  # Check whether an index already exists at this basename
  existing <- Sys.glob(paste0(basename, ".*"))
  if (length(existing) > 0) {
    if (!force) {
      cat("Index already exists at: ", basename)
      cat("  (", length(existing), " file(s) found — use force = TRUE to rebuild)")
      return(invisible(basename))
    }
    cat("Overwriting existing index at: ", basename, " (force = TRUE)")
  }

  cat("Building Rsubread genome index — this may take a while for large genomes.")
  cat("  Reference  : ", reference)
  cat("  Index path : ", basename)
  cat("  Memory     : ", memory, " MB")
  cat("  Split index: ", index_split)

  Rsubread::buildindex(
    basename    = basename,
    reference   = reference,
    memory      = memory,
    indexSplit  = index_split,
    gappedIndex = gapped,
    ...
  )

  cat("Index build complete.")
  invisible(basename)
}
