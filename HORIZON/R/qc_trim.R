#' Run QC and adapter trimming with fastp
#'
#' Wraps \code{\link[Rfastp]{rfastp}} to perform quality control and adapter
#' trimming in a single pass. fastp auto-detects adapters by default:
#' for paired-end data it uses read-pair overlap analysis; for single-end data
#' it falls back to a built-in list of common Illumina adapters. Explicit
#' adapter sequences can be provided to override auto-detection.
#'
#' Outputs are written to \code{<output_dir>/<sample_id>/qc/}. Rfastp uses
#' \code{outputFastq} as a path prefix; trimmed FASTQs are named
#' \code{<prefix>_R1.fastq.gz} / \code{<prefix>_R2.fastq.gz} and QC reports
#' are \code{<prefix>.html} / \code{<prefix>.json}.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process. Must match a value in the
#'   \code{sample_id} column of \code{sample_sheet}.
#' @param threads Integer. Number of threads for fastp. Default 4.
#' @param adapter_r1 Character or NULL. Explicit R1 adapter sequence. Default
#'   NULL (auto-detect via \code{adapterSequenceRead1 = "auto"}).
#' @param adapter_r2 Character or NULL. Explicit R2 adapter sequence. Default
#'   NULL (auto-detect). Ignored for single-end samples.
#' @param force Logical. If \code{FALSE} (default), skip trimming when the
#'   output FASTQs already exist and return their paths directly.  Set to
#'   \code{TRUE} to re-run and overwrite existing files.
#' @param ... Additional arguments passed directly to
#'   \code{\link[Rfastp]{rfastp}} (e.g., \code{qualityFiltering},
#'   \code{qualityFilterPhred}, \code{lengthFiltering}, \code{minReadLength}).
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{trimmed_r1}{Path to trimmed R1 FASTQ.}
#'     \item{trimmed_r2}{Path to trimmed R2 FASTQ, or NULL for single-end.}
#'     \item{report_html}{Path to HTML QC report.}
#'     \item{report_json}{Path to JSON QC report.}
#'     \item{json_data}{The JSON report object returned by rfastp, or NULL
#'       when skipped.}
#'   }
#' @export
HORIZON_run_qc_trim <- function(sample_sheet,
                                sample_id,
                                threads    = 4,
                                adapter_r1 = NULL,
                                adapter_r2 = NULL,
                                force      = FALSE,
                                ...) {
  row    <- .get_sample_row(sample_sheet, sample_id)
  paired <- as.logical(row$paired_end)

  out_dir <- file.path(row$output_dir, sample_id, "qc")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # rfastp takes a path prefix; it appends _R1.fastq.gz / _R2.fastq.gz
  # and generates <prefix>.html / <prefix>.json automatically
  out_prefix <- file.path(out_dir, paste0(sample_id, "_trimmed"))
  r1_out <- paste0(out_prefix, "_R1.fastq.gz")
  r2_out <- paste0(out_prefix, "_R2.fastq.gz")

  # Skip if outputs already exist and force=FALSE
  outputs_exist <- file.exists(r1_out) && (!paired || file.exists(r2_out))
  if (!isTRUE(force) && outputs_exist) {
    cat("Trimmed FASTQs already exist for: ", sample_id,
            " — skipping (use force=TRUE to overwrite)")
    return(invisible(list(
      trimmed_r1  = r1_out,
      trimmed_r2  = if (paired) r2_out else NULL,
      report_html = paste0(out_prefix, ".html"),
      report_json = paste0(out_prefix, ".json"),
      json_data   = NULL
    )))
  }

  cat("Running fastp QC + trimming for: ", sample_id)

  args <- list(
    read1       = row$fastq_r1,
    outputFastq = out_prefix,
    thread      = threads,
    ...
  )

  if (paired) args$read2 <- row$fastq_r2

  if (!is.null(adapter_r1)) args$adapterSequenceRead1 <- adapter_r1
  if (!is.null(adapter_r2) && paired) args$adapterSequenceRead2 <- adapter_r2

  json_data <- do.call(Rfastp::rfastp, args)

  cat("QC + trimming complete for: ", sample_id)

  list(
    trimmed_r1  = paste0(out_prefix, "_R1.fastq.gz"),
    trimmed_r2  = if (paired) paste0(out_prefix, "_R2.fastq.gz") else NULL,
    report_html = paste0(out_prefix, ".html"),
    report_json = paste0(out_prefix, ".json"),
    json_data   = json_data
  )
}
