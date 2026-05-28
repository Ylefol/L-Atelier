# ==============================================================================
# HORIZON - Chromosome Size Utilities
# ==============================================================================


#' Get chromosome sizes
#'
#' Returns a two-column data.frame of chromosome names and sizes.  Accepts
#' three input types so callers can use whichever file is most convenient:
#'
#' \describe{
#'   \item{FASTA (\code{.fa}, \code{.fasta}, \code{.fa.gz},
#'     \code{.fasta.gz})}{Runs \code{samtools faidx} via the registered conda
#'     environment to generate a \code{.fai} index next to the FASTA, then
#'     reads the first two columns.  Requires \code{\link{HORIZON_set_conda_env}}
#'     to have been called first.}
#'   \item{FAI (\code{.fai})}{Reads the samtools fasta index directly.
#'     Columns 1 and 2 are the chromosome name and length.}
#'   \item{Two-column text file}{Any tab- or whitespace-separated file whose
#'     first column is the chromosome name and second column is the size.
#'     Typical extensions: \code{.sizes}, \code{.chrom.sizes}, \code{.tsv},
#'     \code{.txt}.}
#' }
#'
#' The return value uses the same \code{data.frame(chr, size)} format as
#' \code{APOLLO_get_chromosome_sizes()} in the GAIA package so that the two
#' outputs are interchangeable.
#'
#' @param genome Character.  Path to a FASTA file, \code{.fai} index, or
#'   pre-made two-column chromosome sizes file.
#' @param chromosomes Character vector or \code{NULL}.  If supplied, only rows
#'   whose \code{chr} value is in this vector are returned.
#' @param cache_dir Character or \code{NULL}.  Directory for the cached RDS
#'   file.  Defaults to a \code{data/} subdirectory relative to the current
#'   working directory.
#' @param cache_name Character or \code{NULL}.  Stem of the cache filename
#'   (without extension).  Defaults to the input filename stem with
#'   \code{_chrom_sizes} appended.
#' @param force Logical.  If \code{TRUE}, regenerate the cache even when it
#'   already exists.  Default \code{FALSE}.
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{chr} (character) and
#'   \code{size} (integer), one row per chromosome.
#' @export
HORIZON_get_chrom_sizes <- function(genome,
                                     chromosomes = NULL,
                                     cache_dir   = NULL,
                                     cache_name  = NULL,
                                     force       = FALSE,
                                     verbose     = TRUE) {

  if (!file.exists(genome))
    stop("File not found: ", genome, call. = FALSE)

  # ── Cache setup ──────────────────────────────────────────────────────────────
  if (is.null(cache_dir))
    cache_dir <- file.path(getwd(), "data")
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(cache_name)) {
    stem       <- basename(genome)
    stem       <- sub("\\.gz$",    "", stem, ignore.case = TRUE)
    stem       <- sub("\\.[^.]+$", "", stem)
    cache_name <- paste0(stem, "_chrom_sizes")
  }

  cache_path <- file.path(cache_dir, paste0(cache_name, ".rds"))

  if (file.exists(cache_path) && !isTRUE(force)) {
    if (isTRUE(verbose))
      cat("[HORIZON] Loading cached chromosome sizes from:",
          basename(cache_path), "\n")
    sizes <- readRDS(cache_path)
    if (!is.null(chromosomes))
      sizes <- sizes[sizes$chr %in% chromosomes, , drop = FALSE]
    return(sizes)
  }

  # ── Determine input type and load ────────────────────────────────────────────
  is_fasta <- grepl("\\.(fa|fasta)(\\.gz)?$", genome, ignore.case = TRUE)
  is_fai   <- grepl("\\.fai$",                genome, ignore.case = TRUE)

  if (is_fasta) {
    fai_path <- paste0(genome, ".fai")
    if (!file.exists(fai_path) || isTRUE(force)) {
      if (isTRUE(verbose))
        cat("[HORIZON] Generating .fai index with samtools faidx...\n")
      exit_code <- .horizon_run_cli("samtools", c("faidx", genome))
      if (exit_code != 0)
        stop("samtools faidx failed for: ", genome, call. = FALSE)
    }
    genome <- fai_path
    is_fai <- TRUE
  }

  if (isTRUE(verbose))
    cat("[HORIZON] Reading chromosome sizes from:", basename(genome), "\n")

  raw <- utils::read.table(genome, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE)

  sizes <- data.frame(
    chr  = as.character(raw[[1]]),
    size = as.integer(raw[[2]]),
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose))
    cat("    Chromosomes read:", nrow(sizes), "\n")

  # ── Cache ────────────────────────────────────────────────────────────────────
  saveRDS(sizes, cache_path)
  if (isTRUE(verbose))
    cat("    Cached to:", basename(cache_path), "\n")

  if (!is.null(chromosomes))
    sizes <- sizes[sizes$chr %in% chromosomes, , drop = FALSE]

  sizes
}
