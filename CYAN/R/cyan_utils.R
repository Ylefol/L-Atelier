#' Validate a BAM file
#'
#' Checks that a BAM file exists and has an accompanying index (\code{.bai}).
#' Stops with an informative message if either is missing.
#'
#' @param bam Character. Path to the BAM file.
#'
#' @return Invisibly returns \code{bam} (normalised absolute path).
#'
#' @export
CYAN_check_bam <- function(bam) {
  bam <- normalizePath(bam, mustWork = FALSE)
  if (!file.exists(bam))
    stop("[CYAN] BAM file not found: ", bam, call. = FALSE)
  bai <- paste0(bam, ".bai")
  bai_alt <- sub("\\.bam$", ".bai", bam)
  if (!file.exists(bai) && !file.exists(bai_alt))
    stop("[CYAN] BAM index not found. Expected: ", bai,
         "\n  Run 'samtools index ", bam, "' to create it.", call. = FALSE)
  invisible(bam)
}

#' Validate reference files
#'
#' Checks that a reference FASTA exists and has an accompanying FASTA index
#' (\code{.fai}). Optionally validates additional VCF or BED paths.
#'
#' @param reference Character. Path to the reference FASTA.
#' @param ...        Additional file paths to check for existence (e.g. VCFs,
#'                   BED files). Checked for existence only, not for indexing.
#'
#' @return Invisibly returns \code{reference} (normalised absolute path).
#'
#' @export
CYAN_check_references <- function(reference, ...) {
  reference <- normalizePath(reference, mustWork = FALSE)
  if (!file.exists(reference))
    stop("[CYAN] Reference FASTA not found: ", reference, call. = FALSE)
  fai <- paste0(reference, ".fai")
  if (!file.exists(fai))
    stop("[CYAN] Reference FASTA index not found: ", fai,
         "\n  Run 'samtools faidx ", reference, "' to create it.", call. = FALSE)

  extras <- list(...)
  for (path in extras) {
    if (!is.null(path) && !file.exists(normalizePath(path, mustWork = FALSE)))
      stop("[CYAN] Required file not found: ", path, call. = FALSE)
  }
  invisible(reference)
}

#' Resolve the GATK executable path
#'
#' Resolves the GATK executable in the following order:
#' \enumerate{
#'   \item Explicit \code{gatk_path} argument
#'   \item \code{GATK_PATH} environment variable
#'   \item \code{gatk} on the system \code{PATH} (via \code{Sys.which})
#' }
#' Stops with an informative error if GATK cannot be found.
#'
#' @param gatk_path Character or \code{NULL}. Explicit path to the GATK
#'   executable (e.g. \code{"/opt/gatk-4.4/gatk"}). If \code{NULL} the
#'   function falls back to the environment variable and then \code{PATH}.
#'
#' @return Character. Resolved path to the GATK executable.
#'
#' @export
CYAN_resolve_gatk <- function(gatk_path = NULL) {
  if (!is.null(gatk_path)) {
    if (!file.exists(gatk_path))
      stop("[CYAN] Supplied gatk_path does not exist: ", gatk_path, call. = FALSE)
    return(gatk_path)
  }

  env_path <- Sys.getenv("GATK_PATH", unset = NA)
  if (!is.na(env_path) && nzchar(env_path)) {
    if (!file.exists(env_path))
      stop("[CYAN] GATK_PATH env var points to a non-existent path: ", env_path,
           call. = FALSE)
    return(env_path)
  }

  which_path <- Sys.which("gatk")
  if (nzchar(which_path))
    return(unname(which_path))

  stop(
    "[CYAN] Could not locate GATK.\n",
    "  Supply it via the 'gatk_path' argument, the GATK_PATH environment variable,\n",
    "  or ensure 'gatk' is on your PATH.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Ensure a directory exists, creating it with a message if absent
#' @keywords internal
.cyan_ensure_dir <- function(path, label = "directory") {
  if (!dir.exists(path)) {
    message("[CYAN] Creating ", label, ": ", path)
    dir.create(path, recursive = TRUE)
  }
  normalizePath(path)
}

#' Derive a sample name from a BAM file path
#' Strip directory and extension, use as a safe file-system label.
#' @keywords internal
.cyan_sample_name <- function(bam, sample_name = NULL) {
  if (!is.null(sample_name)) return(sample_name)
  tools::file_path_sans_ext(basename(bam))
}
