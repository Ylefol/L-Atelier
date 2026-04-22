#' Create a sample sheet template
#'
#' Generates a CSV template with all required columns for a HORIZON pipeline
#' run. Fill in the paths and metadata, then pass to
#' \code{\link{HORIZON_validate_sample_sheet}} before running any pipeline steps.
#'
#' Required columns:
#' \describe{
#'   \item{sample_id}{Unique sample identifier. Used to name output files and folders.}
#'   \item{fastq_r1}{Absolute path to R1 (or single-end) FASTQ file.}
#'   \item{fastq_r2}{Absolute path to R2 FASTQ file. Leave empty for single-end.}
#'   \item{paired_end}{Logical. TRUE for paired-end libraries.}
#'   \item{output_dir}{Absolute path to the root output directory for this sample.
#'     A sub-folder named \code{sample_id} will be created inside it.}
#'   \item{strandedness}{Library strandedness: \code{"unstranded"}, \code{"forward"},
#'     or \code{"reverse"}.}
#'   \item{condition}{Experimental condition or group label. Passed through to the
#'     aggregated count matrix as metadata.}
#' }
#'
#' Additional columns can be added freely; they are carried through as metadata.
#'
#' @param output_path Character or NULL. Path to save the template CSV.
#'   If NULL, the data.frame is returned without writing to disk. Default NULL.
#' @param n_samples Integer. Number of placeholder rows to include. Default 1.
#'
#' @return A data.frame with the sample sheet template, invisibly.
#' @export
HORIZON_create_sample_sheet <- function(output_path = NULL, n_samples = 1) {
  template <- data.frame(
    sample_id    = paste0("sample_", seq_len(n_samples)),
    fastq_r1     = rep("", n_samples),
    fastq_r2     = rep("", n_samples),
    paired_end   = rep(TRUE, n_samples),
    output_dir   = rep("", n_samples),
    strandedness = rep("unstranded", n_samples),
    condition    = rep("", n_samples),
    stringsAsFactors = FALSE
  )

  if (!is.null(output_path)) {
    write.csv(template, output_path, row.names = FALSE)
    message("Sample sheet template written to: ", output_path)
  }

  invisible(template)
}


#' Validate and load a sample sheet
#'
#' Reads a HORIZON sample sheet CSV (or accepts a data.frame) and validates it
#' before use. Checks for required columns, unique sample IDs, valid strandedness
#' values, and optionally whether all FASTQ paths exist on disk.
#'
#' @param path Character or data.frame. Path to a sample sheet CSV, or an
#'   already-loaded data.frame.
#' @param check_files Logical. Whether to verify that FASTQ file paths exist
#'   on disk. Default TRUE. Set to FALSE if files are on a drive not currently
#'   mounted.
#'
#' @return The validated sample sheet as a data.frame.
#' @export
HORIZON_validate_sample_sheet <- function(path, check_files = TRUE) {
  if (is.character(path)) {
    if (!file.exists(path)) stop("Sample sheet file not found: ", path)
    ss <- read.csv(path, stringsAsFactors = FALSE)
  } else if (is.data.frame(path)) {
    ss <- path
  } else {
    stop("'path' must be a file path string or a data.frame.")
  }

  # Required columns
  required_cols <- c("sample_id", "fastq_r1", "paired_end", "output_dir", "strandedness")
  missing_cols  <- setdiff(required_cols, colnames(ss))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Unique sample IDs
  if (anyDuplicated(ss$sample_id)) {
    dups <- unique(ss$sample_id[duplicated(ss$sample_id)])
    stop("Duplicate sample_id values: ", paste(dups, collapse = ", "))
  }

  # Valid strandedness
  valid_strand <- c("unstranded", "forward", "reverse")
  bad_strand   <- setdiff(unique(ss$strandedness), valid_strand)
  if (length(bad_strand) > 0) {
    stop("Invalid strandedness value(s): ", paste(bad_strand, collapse = ", "),
         ". Must be one of: ", paste(valid_strand, collapse = ", "))
  }

  # PARSE columns: validate chemistry and kit if present
  if (all(c("chemistry", "kit") %in% colnames(ss)))
    .validate_chemistry_kit(ss)

  # Add fastq_r2 column if missing (single-end only sheet)
  if (!"fastq_r2" %in% colnames(ss)) {
    ss$fastq_r2 <- NA_character_
  }

  if (check_files) {
    # R1 paths
    missing_r1 <- ss$fastq_r1[!file.exists(ss$fastq_r1)]
    if (length(missing_r1) > 0) {
      stop("R1 FASTQ file(s) not found:\n", paste(" ", missing_r1, collapse = "\n"))
    }

    # R2 paths for paired-end samples
    pe_rows <- ss[as.logical(ss$paired_end), ]
    if (nrow(pe_rows) > 0) {
      r2_paths   <- pe_rows$fastq_r2
      r2_missing <- r2_paths[is.na(r2_paths) | r2_paths == "" | !file.exists(r2_paths)]
      if (length(r2_missing) > 0) {
        stop("R2 FASTQ file(s) missing or not found for paired-end samples:\n",
             paste(" ", r2_missing, collapse = "\n"))
      }
    }
  }

  message("Sample sheet OK — ", nrow(ss), " sample(s) validated.")
  ss
}
