#' Concatenate multi-lane FASTQ files into one file per read direction
#'
#' Some aligners/tools accept multiple FASTQ files per sample natively (e.g.
#' \code{\link{HORIZON_run_salmon}}'s \code{fastq_r1}/\code{fastq_r2} vectors),
#' but most of the pipeline (\code{\link{HORIZON_run_qc_trim}},
#' \code{\link{HORIZON_run_align}}) expects exactly one FASTQ per sample per
#' read direction. When a sample was sequenced across multiple
#' flowcells/lanes, this function merges those lane-level FASTQs into a single
#' R1 (and R2, for paired-end data) file per sample, so the rest of the
#' pipeline can treat it like any single-lane sample.
#'
#' Concatenation is done at the byte level (raw \code{.fastq.gz} bytes
#' appended in sequence), never decompressing/recompressing. A gzip file is a
#' series of independent "members" (RFC 1952); concatenating gzip files this
#' way produces another valid gzip file that decompresses to the full,
#' correctly-ordered read set — this is the standard way multi-lane FASTQs are
#' merged (equivalent to \code{cat lane1.fastq.gz lane2.fastq.gz > merged.fastq.gz}),
#' and every downstream reader here (Rfastp, Rsubread, zcat) handles
#' multi-member gzip transparently.
#'
#' @param fastq_r1 Character vector. Paths to R1 FASTQ files across lanes, for
#'   one sample.
#' @param fastq_r2 Character vector or \code{NULL}. Paths to R2 FASTQ files
#'   across lanes, positionally matched to \code{fastq_r1} (same length, same
#'   lane order — \code{fastq_r1[i]} and \code{fastq_r2[i]} must be the same
#'   lane). \code{NULL} for single-end data.
#' @param output_dir Character. Directory to write the concatenated FASTQ(s)
#'   into. Created if it doesn't exist.
#' @param sample_id Character. Sample ID; output files are named
#'   \code{<sample_id>_R1.fastq.gz} / \code{<sample_id>_R2.fastq.gz}.
#' @param force Logical. If \code{FALSE} (default), skip concatenation when
#'   the output file(s) already exist and return their paths directly. Set to
#'   \code{TRUE} to re-run and overwrite.
#'
#' @return Named list with elements \code{fastq_r1} and \code{fastq_r2} (the
#'   latter \code{NULL} for single-end data) giving the paths to the
#'   concatenated FASTQ(s), invisibly.
#' @export
HORIZON_concat_lanes <- function(fastq_r1,
                                  fastq_r2   = NULL,
                                  output_dir,
                                  sample_id,
                                  force      = FALSE) {

  if (length(fastq_r1) == 0) stop("'fastq_r1' must contain at least one file.")
  if (!all(file.exists(fastq_r1)))
    stop("R1 FASTQ(s) not found: ",
         paste(fastq_r1[!file.exists(fastq_r1)], collapse = ", "))

  paired <- !is.null(fastq_r2)
  if (paired) {
    if (length(fastq_r2) != length(fastq_r1))
      stop("'fastq_r1' and 'fastq_r2' must have the same length (",
           length(fastq_r1), " vs ", length(fastq_r2), ").")
    if (!all(file.exists(fastq_r2)))
      stop("R2 FASTQ(s) not found: ",
           paste(fastq_r2[!file.exists(fastq_r2)], collapse = ", "))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  r1_out <- file.path(output_dir, paste0(sample_id, "_R1.fastq.gz"))
  r2_out <- if (paired) file.path(output_dir, paste0(sample_id, "_R2.fastq.gz")) else NULL

  outputs_exist <- file.exists(r1_out) && (!paired || file.exists(r2_out))
  if (!isTRUE(force) && outputs_exist) {
    cat("Concatenated FASTQ(s) already exist for: ", sample_id,
            " — skipping (use force=TRUE to overwrite)")
    return(invisible(list(fastq_r1 = r1_out, fastq_r2 = r2_out)))
  }

  cat("Concatenating ", length(fastq_r1), " lane(s) for: ", sample_id)
  .horizon_concat_bytes(fastq_r1, r1_out)
  if (paired) .horizon_concat_bytes(fastq_r2, r2_out)

  cat("Concatenation complete for: ", sample_id)
  invisible(list(fastq_r1 = r1_out, fastq_r2 = r2_out))
}


#' Concatenate files at the byte level [internal]
#'
#' @param input_files Character vector of input file paths, written to
#'   \code{output_file} in order.
#' @param output_file Character. Destination file path.
#' @keywords internal
.horizon_concat_bytes <- function(input_files, output_file) {
  con_out <- file(output_file, "wb")
  on.exit(close(con_out), add = TRUE)
  for (f in input_files) {
    con_in <- file(f, "rb")
    repeat {
      chunk <- readBin(con_in, "raw", n = 1e8)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }
    close(con_in)
  }
}
