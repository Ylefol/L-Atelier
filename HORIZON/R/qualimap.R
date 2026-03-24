# HORIZON_run_qualimap removed — qualimap has heavy Java/conda dependencies
# that conflict with other packages in the basilisk environment.
# BAM QC is covered by flagstat output from HORIZON_process_bam and
# insert size distribution from HORIZON_run_fragment_size.


#' Compute insert size distribution with bamPEFragmentSize
#'
#' Wraps \code{deeptools bamPEFragmentSize} to produce an insert-size
#' histogram.  For ATAC-seq data, the nucleosomal periodicity (peaks near
#' ~200 bp, ~400 bp, ~600 bp) is visible in this plot.
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param threads Integer. Number of processors (\code{--numberOfProcessors}).
#'   Default 4.
#' @param max_length Integer. Maximum fragment length shown on the x-axis
#'   (\code{--maxFragmentLength}).  Default 1000.
#' @param force Logical. If \code{FALSE} (default), skip when the output PNG
#'   already exists.  Set \code{TRUE} to rerun.
#'
#' @return Character. Path to the output PNG, invisibly.
#' @export
HORIZON_run_fragment_size <- function(sample_sheet,
                                       sample_id,
                                       threads    = 4L,
                                       max_length = 1000L,
                                       force      = FALSE) {

  bam_path <- .latest_bam(sample_sheet, sample_id)
  row      <- .get_sample_row(sample_sheet, sample_id)
  out_dir  <- file.path(row$output_dir, sample_id, "aligned")
  out_png  <- file.path(out_dir, paste0(sample_id, "_fragment_size.png"))

  if (!isTRUE(force) && file.exists(out_png)) {
    message("Fragment size plot already exists for: ", sample_id,
            " — skipping (use force=TRUE to rerun)")
    return(invisible(out_png))
  }

  args <- c(
    "-b",                    bam_path,
    "-hist",                 out_png,
    "--numberOfProcessors",  as.integer(threads),
    "--maxFragmentLength",   as.integer(max_length),
    "--plotTitle",           sample_id
  )

  message("[", sample_id, "] Running bamPEFragmentSize")
  exit <- .horizon_run_cli("bamPEFragmentSize", args)
  if (exit != 0)
    stop("bamPEFragmentSize failed (exit code ", exit, ").", call. = FALSE)

  message("Fragment size plot complete: ", out_png)
  invisible(out_png)
}
