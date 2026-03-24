#' Call peaks with MACS3
#'
#' Wraps the \code{macs3 callpeak} command via the user-managed conda
#' environment registered with \code{\link{HORIZON_set_conda_env}}.
#'
#' \strong{Recommended settings by assay:}
#' \itemize{
#'   \item ATAC-seq — \code{nomodel=TRUE, shift=-100, extsize=200}
#'     (shift model is not appropriate for open-chromatin data)
#'   \item ChIP-seq (transcription factor) — defaults
#'     (MACS3 builds a shift model from the data; provide \code{control_bed})
#'   \item ChIP-seq (histone mark) — \code{broad=TRUE}
#'   \item CUT&TAG / CUT&RUN — \code{nomodel=TRUE},
#'     optionally \code{extra_flags="--keep-dup all"} if duplicates were not
#'     removed upstream (CUT&TAG can generate high local duplicate rates from
#'     genuine signal rather than PCR artefact)
#' }
#'
#' @param treatment_bed Character. Path to the treatment fragment BED file
#'   (from \code{\link{HORIZON_bam_to_bed}}).
#' @param control_bed Character or \code{NULL}. Path to a control (input/IgG)
#'   BED file.  Recommended for ChIP-seq; optional for CUT&TAG/CUT&RUN; not
#'   used for ATAC-seq.
#' @param output_dir Character. Directory for MACS3 output files.
#' @param sample_name Character. Prefix applied to all output filenames.
#' @param genome_size Character or numeric. Effective genome size for MACS3
#'   (\code{-g}).  Shortcuts: \code{"hs"}, \code{"mm"}, \code{"ce"},
#'   \code{"dm"}.  Pass a numeric value for other genomes.
#' @param q_value Numeric. FDR threshold (\code{-q}). Default 0.05.
#' @param broad Logical. Call broad peaks (\code{--broad}).  Recommended for
#'   histone mark ChIP-seq.  Default \code{FALSE}.
#' @param nomodel Logical. Skip the MACS3 shift-model building step
#'   (\code{--nomodel}).  Required when \code{shift} and \code{extsize} are
#'   set manually.  Default \code{FALSE}.
#' @param shift Integer or \code{NULL}. Read shift in bp (\code{--shift}).
#'   Only used when \code{nomodel=TRUE}.  Use \code{-100} for ATAC-seq.
#'   Default \code{NULL} (not passed).
#' @param extsize Integer or \code{NULL}. Fragment extension size in bp
#'   (\code{--extsize}).  Only used when \code{nomodel=TRUE}.  Use \code{200}
#'   for ATAC-seq.  Default \code{NULL} (not passed).
#' @param extra_flags Character vector of additional MACS3 flags to append
#'   verbatim (e.g. \code{c("--keep-dup", "all")} for CUT&TAG/CUT&RUN).
#'   Default \code{character(0)}.
#' @param force Logical. If \code{FALSE} (default), skip peak calling when the
#'   output peak file already exists.  Set \code{TRUE} to rerun.
#'
#' @return Character. Path to the narrowPeak (or broadPeak) file, invisibly.
#' @export
HORIZON_call_peaks <- function(treatment_bed,
                                control_bed  = NULL,
                                output_dir,
                                sample_name,
                                genome_size  = "hs",
                                q_value      = 0.05,
                                broad        = FALSE,
                                nomodel      = FALSE,
                                shift        = NULL,
                                extsize      = NULL,
                                extra_flags  = character(0),
                                force        = FALSE) {

  peak_ext       <- if (isTRUE(broad)) "_peaks.broadPeak" else "_peaks.narrowPeak"
  peak_file_check <- file.path(output_dir, paste0(sample_name, peak_ext))
  if (!isTRUE(force) && file.exists(peak_file_check)) {
    message("Peak file already exists for: ", sample_name,
            " — skipping (use force=TRUE to rerun)")
    return(invisible(peak_file_check))
  }

  if (!file.exists(treatment_bed))
    stop("Treatment BED not found: ", treatment_bed, call. = FALSE)
  if (!is.null(control_bed) && !file.exists(control_bed))
    stop("Control BED not found: ", control_bed, call. = FALSE)
  if (!is.null(shift) && !isTRUE(nomodel))
    warning("'shift' is set but nomodel=FALSE — MACS3 ignores --shift without --nomodel.",
            call. = FALSE)
  if (!is.null(extsize) && !isTRUE(nomodel))
    warning("'extsize' is set but nomodel=FALSE — MACS3 ignores --extsize without --nomodel.",
            call. = FALSE)

  # Assemble model flags
  model_flags <- character(0)
  if (isTRUE(nomodel)) {
    model_flags <- c(model_flags, "--nomodel")
    if (!is.null(shift))   model_flags <- c(model_flags, "--shift",   as.character(shift))
    if (!is.null(extsize)) model_flags <- c(model_flags, "--extsize", as.character(extsize))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "callpeak",
    "-t", treatment_bed,
    "-f", "BEDPE",
    "-n", sample_name,
    "--outdir", output_dir,
    "-g", as.character(genome_size),
    "-q", as.character(q_value),
    model_flags,
    extra_flags
  )
  if (!is.null(control_bed))
    args <- c(args, "-c", control_bed)
  if (isTRUE(broad))
    args <- c(args, "--broad")

  message("[", sample_name, "] Calling peaks with MACS3")
  exit_code <- .horizon_run_cli("macs3", args)
  if (exit_code != 0)
    stop("MACS3 failed for '", sample_name, "' (exit code ", exit_code, ").",
         call. = FALSE)

  peak_file <- file.path(output_dir, paste0(sample_name, peak_ext))

  if (!file.exists(peak_file))
    stop("MACS3 ran but output peak file not found: ", peak_file, call. = FALSE)

  message("Peak calling complete for: ", sample_name,
          "\n  Peaks: ", peak_file)
  invisible(peak_file)
}
