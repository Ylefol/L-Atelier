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
    cat("Peak file already exists for: ", sample_name,
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

  cat("[", sample_name, "] Calling peaks with MACS3")
  exit_code <- .horizon_run_cli("macs3", args)
  if (exit_code != 0)
    stop("MACS3 failed for '", sample_name, "' (exit code ", exit_code, ").",
         call. = FALSE)

  peak_file <- file.path(output_dir, paste0(sample_name, peak_ext))

  if (!file.exists(peak_file))
    stop("MACS3 ran but output peak file not found: ", peak_file, call. = FALSE)

  cat("Peak calling complete for: ", sample_name,
          "\n  Peaks: ", peak_file)
  invisible(peak_file)
}


# ==============================================================================
# SEACR — pure-R reimplementation
# ==============================================================================

# Convert a fragment BED data.frame to a per-base RleList coverage.
# chrom_sizes_vec: named integer vector (names = chr, values = length).
.horizon_bed_to_coverage <- function(bed_df, chrom_sizes_vec) {
  # Clip coordinates to chromosome boundaries before building coverage so that
  # out-of-bounds fragments (e.g. near centromeres) do not cause errors.
  chr_col   <- as.character(bed_df[[1]])
  start_col <- as.integer(bed_df[[2]])
  end_col   <- as.integer(bed_df[[3]])

  keep <- chr_col %in% names(chrom_sizes_vec)
  if (!all(keep)) {
    n_drop <- sum(!keep)
    warning(n_drop, " fragment(s) on unrecognised chromosomes dropped.",
            call. = FALSE)
    chr_col   <- chr_col[keep]
    start_col <- start_col[keep]
    end_col   <- end_col[keep]
  }

  chr_lens  <- chrom_sizes_vec[chr_col]
  start_col <- pmax(1L,         start_col)
  end_col   <- pmin(chr_lens,   end_col)

  gr <- GenomicRanges::GRanges(
    seqnames = chr_col,
    ranges   = IRanges::IRanges(start = start_col, end = end_col),
    seqlengths = chrom_sizes_vec
  )

  GenomicRanges::coverage(gr)
}


# Find contiguous non-zero signal blocks in an RleList and compute AUC and
# max signal for each block.
# Returns data.frame(chr, start, end, total_signal, max_signal).
.horizon_seacr_blocks <- function(cov_rle, chrom_sizes_vec) {
  chrs <- intersect(names(chrom_sizes_vec), names(cov_rle))

  block_list <- lapply(chrs, function(chr) {
    rle <- cov_rle[[chr]]
    # IRanges::slice finds runs of signal > 0
    slices <- IRanges::slice(rle, lower = 1L, rangesOnly = FALSE)
    if (length(slices) == 0L)
      return(NULL)

    data.frame(
      chr          = chr,
      start        = IRanges::start(slices),
      end          = IRanges::end(slices),
      total_signal = IRanges::viewSums(slices),
      max_signal   = IRanges::viewMaxs(slices),
      stringsAsFactors = FALSE
    )
  })

  blocks <- do.call(rbind, Filter(Negate(is.null), block_list))
  if (is.null(blocks) || nrow(blocks) == 0L)
    return(data.frame(chr = character(), start = integer(), end = integer(),
                      total_signal = numeric(), max_signal = numeric(),
                      stringsAsFactors = FALSE))
  rownames(blocks) <- NULL
  blocks
}


#' Call peaks with SEACR (pure-R reimplementation)
#'
#' A native R reimplementation of the SEACR (Sparse Enrichment Analysis for
#' CUT\&RUN) peak-calling method published by Meers, Tenenbaum and Bhatt
#' (2019, \doi{10.1186/s13072-019-0287-4}).  \strong{This function
#' reimplements their published algorithm; it is not an independent method.}
#' Results should be validated against the original SEACR tool
#' (\url{https://github.com/FredHutch/SEACR}) before use in publications.
#'
#' \strong{Algorithm summary (Meers et al. 2019):}
#' \enumerate{
#'   \item Build per-base fragment coverage from the treatment BED file.
#'   \item Identify contiguous blocks of non-zero signal.
#'   \item Compute the total signal (AUC) and maximum signal height per block.
#'   \item If a control BED is provided: subtract the per-block control AUC
#'     from the treatment AUC; retain blocks with positive net signal.
#'   \item If no control is provided: rank blocks by AUC and retain the top
#'     \code{threshold} fraction.
#'   \item Apply stringency: \code{"stringent"} requires a block to exceed the
#'     threshold on both total AUC and max signal; \code{"relaxed"} requires
#'     either criterion.
#' }
#'
#' The interface mirrors \code{\link{HORIZON_call_peaks}} (MACS3) so the two
#' functions can be called with the same core arguments in a comparison
#' workflow.
#'
#' @param treatment_bed Character.  Path to the treatment fragment BED file
#'   (from \code{\link{HORIZON_bam_to_bed}}).
#' @param control_bed Character or \code{NULL}.  Path to an IgG / input
#'   fragment BED file.  When \code{NULL} (default), threshold-based calling
#'   is used (\code{threshold}).
#' @param output_dir Character.  Directory for the output peak file.
#' @param sample_name Character.  Prefix applied to the output filename.
#' @param chrom_sizes A \code{data.frame} with columns \code{chr} and
#'   \code{size} (from \code{\link{HORIZON_get_chrom_sizes}}).
#' @param threshold Numeric between 0 and 1.  Fraction of blocks to retain
#'   by total AUC when no control is provided.  Default \code{0.05} (top 5\%).
#'   Ignored when \code{control_bed} is supplied.
#' @param normalize Character.  \code{"non"} (default): use raw fragment
#'   counts.  \code{"norm"}: scale coverage by total fragments (per million)
#'   before computing AUC, matching SEACR's \code{norm} mode.
#' @param stringency Character.  \code{"stringent"} (default): a block must
#'   exceed the threshold on both total AUC and max signal.
#'   \code{"relaxed"}: a block passes if either criterion is met.
#' @param force Logical.  Re-run even if the output file already exists.
#'   Default \code{FALSE}.
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return Character.  Path to the output BED file
#'   (\code{<sample_name>_seacr_peaks.bed}), invisibly.  Columns:
#'   \code{peak_id}, \code{chr}, \code{start}, \code{end},
#'   \code{total_signal}, \code{max_signal}.
#' @export
HORIZON_call_peaks_seacr <- function(treatment_bed,
                                      control_bed  = NULL,
                                      output_dir,
                                      sample_name,
                                      chrom_sizes,
                                      threshold    = 0.05,
                                      normalize    = c("non", "norm"),
                                      stringency   = c("stringent", "relaxed"),
                                      force        = FALSE,
                                      verbose      = TRUE) {

  normalize  <- match.arg(normalize)
  stringency <- match.arg(stringency)

  if (!file.exists(treatment_bed))
    stop("Treatment BED not found: ", treatment_bed, call. = FALSE)
  if (!is.null(control_bed) && !file.exists(control_bed))
    stop("Control BED not found: ", control_bed, call. = FALSE)
  if (!is.data.frame(chrom_sizes) ||
      !all(c("chr", "size") %in% colnames(chrom_sizes)))
    stop("chrom_sizes must be a data.frame with 'chr' and 'size' columns ",
         "(from HORIZON_get_chrom_sizes()).", call. = FALSE)
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1)
    stop("threshold must be a number strictly between 0 and 1.", call. = FALSE)

  out_file <- file.path(output_dir, paste0(sample_name, "_seacr_peaks.bed"))
  if (!isTRUE(force) && file.exists(out_file)) {
    if (isTRUE(verbose))
      cat("Peak file already exists for: ", sample_name,
              " — skipping (use force=TRUE to rerun)")
    return(invisible(out_file))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  chrom_sizes_vec <- setNames(as.integer(chrom_sizes$size),
                               as.character(chrom_sizes$chr))

  # ── Build treatment coverage ─────────────────────────────────────────────────
  if (isTRUE(verbose)) cat("[SEACR] Reading treatment BED:", basename(treatment_bed), "\n")
  t_bed <- utils::read.table(treatment_bed, header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE)
  n_t <- nrow(t_bed)
  if (isTRUE(verbose)) cat("    Fragments:", format(n_t, big.mark = ","), "\n")

  t_cov <- .horizon_bed_to_coverage(t_bed, chrom_sizes_vec)
  if (normalize == "norm") t_cov <- t_cov / (n_t / 1e6)

  # ── Find treatment blocks ────────────────────────────────────────────────────
  if (isTRUE(verbose)) cat("[SEACR] Identifying signal blocks...\n")
  t_blocks <- .horizon_seacr_blocks(t_cov, chrom_sizes_vec)

  if (nrow(t_blocks) == 0L)
    stop("No signal blocks found in treatment BED. Check input data.",
         call. = FALSE)

  if (isTRUE(verbose))
    cat("    Total blocks:", format(nrow(t_blocks), big.mark = ","), "\n")

  # ── Threshold derivation ─────────────────────────────────────────────────────
  if (!is.null(control_bed)) {

    if (isTRUE(verbose)) cat("[SEACR] Reading control BED:", basename(control_bed), "\n")
    c_bed <- utils::read.table(control_bed, header = FALSE, sep = "\t",
                                stringsAsFactors = FALSE)
    n_c <- nrow(c_bed)
    if (isTRUE(verbose)) cat("    Fragments:", format(n_c, big.mark = ","), "\n")

    c_cov <- .horizon_bed_to_coverage(c_bed, chrom_sizes_vec)
    if (normalize == "norm") c_cov <- c_cov / (n_c / 1e6)

    # For each treatment block, query total control signal over the same region
    t_gr <- GenomicRanges::GRanges(
      seqnames   = t_blocks$chr,
      ranges     = IRanges::IRanges(start = t_blocks$start, end = t_blocks$end),
      seqlengths = chrom_sizes_vec
    )
    c_views <- GenomicRanges::binnedAverage(
      bins     = t_gr,
      numvar   = c_cov,
      varname  = "ctrl_avg"
    )
    ctrl_auc <- c_views$ctrl_avg * IRanges::width(t_gr)
    net_auc  <- t_blocks$total_signal - ctrl_auc

    auc_thresh     <- 0
    max_thresh     <- 0
    t_blocks$net_auc <- net_auc

    if (stringency == "stringent") {
      pass <- net_auc > auc_thresh &
              t_blocks$max_signal > max_thresh
    } else {
      pass <- net_auc > auc_thresh |
              t_blocks$max_signal > max_thresh
    }

  } else {

    # Threshold-based: top `threshold` fraction by total AUC
    auc_cutoff <- stats::quantile(t_blocks$total_signal,
                                   probs = 1 - threshold,
                                   names = FALSE)
    max_cutoff <- stats::quantile(t_blocks$max_signal,
                                   probs = 1 - threshold,
                                   names = FALSE)

    if (isTRUE(verbose))
      cat("    AUC cutoff (top", scales::percent(threshold), "):",
          round(auc_cutoff, 2), "\n")

    if (stringency == "stringent") {
      pass <- t_blocks$total_signal >= auc_cutoff &
              t_blocks$max_signal   >= max_cutoff
    } else {
      pass <- t_blocks$total_signal >= auc_cutoff |
              t_blocks$max_signal   >= max_cutoff
    }
  }

  peaks <- t_blocks[pass, , drop = FALSE]

  if (isTRUE(verbose))
    cat("[SEACR] Peaks called:", format(nrow(peaks), big.mark = ","),
        "(", stringency, "mode)\n")

  if (nrow(peaks) == 0L)
    warning("No peaks passed the threshold. Consider raising 'threshold' or ",
            "switching to stringency = 'relaxed'.", call. = FALSE)

  # ── Write output ─────────────────────────────────────────────────────────────
  peaks <- peaks[order(peaks$chr, peaks$start), ]
  out <- data.frame(
    peak_id      = paste0(sample_name, "_peak_", seq_len(nrow(peaks))),
    chr          = peaks$chr,
    start        = peaks$start,
    end          = peaks$end,
    total_signal = round(peaks$total_signal, 4),
    max_signal   = round(peaks$max_signal,   4),
    stringsAsFactors = FALSE
  )

  utils::write.table(out, out_file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  if (isTRUE(verbose))
    cat("[SEACR] Done. Output: ", out_file)

  invisible(out_file)
}
