#' Process a BAM file: collate, fixmate, sort, mark/remove duplicates
#'
#' Runs the standard samtools duplicate-removal pipeline on the MAPQ-filtered
#' BAM produced by \code{\link{HORIZON_run_bowtie2}}.  Steps performed in
#' order:
#' \enumerate{
#'   \item \strong{collate} — name-sort for fixmate
#'   \item \strong{fixmate} — fill mate coordinates (required by markdup)
#'   \item \strong{sort}    — coordinate sort
#'   \item \strong{markdup} — mark and remove PCR/optical duplicates
#'   \item \strong{chrM removal} (optional) — remove mitochondrial reads
#'   \item \strong{index}   — index the final BAM
#'   \item \strong{flagstat} — write alignment statistics
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param threads Integer. Threads for all samtools commands. Default 4.
#' @param memory_per_thread Integer. Memory in MB per thread for
#'   \code{samtools sort} (\code{-m}). Default 2000.
#' @param remove_chrM Logical. Remove mitochondrial reads from the final BAM.
#'   Default \code{TRUE}.
#' @param remove_tmp Logical. Remove intermediate BAM files after processing.
#'   Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip processing when
#'   \code{_processed.bam} already exists.  Set \code{TRUE} to reprocess.
#'
#' @return Character. Path to the processed BAM file, invisibly.
#'   Pass to \code{\link{HORIZON_filter_blacklist}} or
#'   \code{\link{HORIZON_bam_to_bigwig}}.
#' @export
HORIZON_process_bam <- function(sample_sheet,
                                 sample_id,
                                 threads           = 4L,
                                 memory_per_thread = 2000L,
                                 remove_chrM       = TRUE,
                                 remove_tmp        = TRUE,
                                 force             = FALSE) {

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")

  input_bam   <- file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))
  collate_bam <- file.path(out_dir, paste0(sample_id, "_collated.bam"))
  fixmate_bam <- file.path(out_dir, paste0(sample_id, "_fixmate.bam"))
  sorted_bam  <- file.path(out_dir, paste0(sample_id, "_sorted.bam"))
  markdup_bam <- file.path(out_dir, paste0(sample_id, "_markdup.bam"))
  markdup_stat <- file.path(out_dir, paste0(sample_id, "_markdup_stats.txt"))
  final_bam   <- file.path(out_dir, paste0(sample_id, "_processed.bam"))
  flagstat_f  <- file.path(out_dir, paste0(sample_id, "_flagstat.txt"))

  if (!isTRUE(force) && file.exists(final_bam)) {
    message("Processed BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to reprocess)")
    return(invisible(final_bam))
  }

  if (!file.exists(input_bam))
    stop("MAPQ-filtered BAM not found: ", input_bam,
         ". Run HORIZON_run_bowtie2() first.", call. = FALSE)

  thr <- as.integer(threads)
  mem <- paste0(as.integer(memory_per_thread), "M")

  # Helper to run a samtools command via the Herper CLI env
  .run_samtools <- function(...) {
    exit <- .horizon_run_cli("samtools", c(...))
    if (exit != 0)
      stop("samtools exited with code ", exit, call. = FALSE)
    invisible(exit)
  }

  # Step 1: collate (name-sort) ------------------------------------------------
  message("[", sample_id, "] samtools collate")
  .run_samtools("collate", "-O", "-u",
                "--threads", thr,
                input_bam,
                file.path(out_dir, paste0(sample_id, "_collate_tmp")),
                "-o", collate_bam)

  # Step 2: fixmate ------------------------------------------------------------
  message("[", sample_id, "] samtools fixmate")
  .run_samtools("fixmate", "-m", "-u",
                "--threads", thr,
                collate_bam, fixmate_bam)
  if (isTRUE(remove_tmp)) file.remove(collate_bam)

  # Step 3: sort (coordinate) --------------------------------------------------
  message("[", sample_id, "] samtools sort")
  .run_samtools("sort",
                "--threads", thr,
                "-m", mem,
                "-T", file.path(out_dir, paste0(sample_id, "_sort_tmp")),
                fixmate_bam, "-o", sorted_bam)
  if (isTRUE(remove_tmp)) file.remove(fixmate_bam)

  # Step 4: markdup + stats ---------------------------------------------------
  message("[", sample_id, "] samtools markdup (removing duplicates)")
  .run_samtools("markdup",
                "--threads", thr,
                "-r",             # remove duplicates
                "-s",             # write stats to stderr (captured below)
                "-f", markdup_stat,
                sorted_bam, markdup_bam)
  if (isTRUE(remove_tmp)) file.remove(sorted_bam)

  # Step 5: chrM removal -------------------------------------------------------
  if (isTRUE(remove_chrM)) {
    message("[", sample_id, "] Removing chrM reads")
    # Index markdup BAM first (required for idxstats)
    .run_samtools("index", "--threads", thr, markdup_bam)

    # Get all chromosome names except chrM and MT
    idxstats_out <- .horizon_run_cli(
      "samtools",
      c("idxstats", markdup_bam),
      stdout = TRUE
    )
    keep_chroms <- grep("^(chrM|MT)\\b", idxstats_out, invert = TRUE,
                        value = TRUE, perl = TRUE)
    keep_chroms <- sub("\t.*", "", keep_chroms)   # first column only
    keep_chroms <- keep_chroms[nzchar(keep_chroms)]

    .run_samtools("view", "-b",
                  "--threads", thr,
                  "-o", final_bam,
                  markdup_bam,
                  keep_chroms)

    if (isTRUE(remove_tmp)) {
      file.remove(markdup_bam)
      bai <- paste0(markdup_bam, ".bai")
      if (file.exists(bai)) file.remove(bai)
    }
  } else {
    file.rename(markdup_bam, final_bam)
  }

  # Step 6: index final BAM ----------------------------------------------------
  message("[", sample_id, "] Indexing processed BAM")
  .run_samtools("index", "--threads", thr, final_bam)

  # Step 7: flagstat -----------------------------------------------------------
  message("[", sample_id, "] Running flagstat")
  .horizon_run_cli(
    "samtools",
    c("flagstat", "--threads", thr, final_bam),
    stdout = flagstat_f
  )

  message("BAM processing complete for: ", sample_id,
          "\n  Output: ", final_bam)
  invisible(final_bam)
}


#' Remove reads overlapping a genomic blacklist
#'
#' Filters a processed BAM file using \code{bedtools intersect -v} to exclude
#' reads that overlap a blacklist BED file (e.g. ENCODE blacklist regions).
#' For assemblies without an ENCODE blacklist (such as T2T-CHM13), pass
#' \code{blacklist_bed = NULL} to skip this step.
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param blacklist_bed Character. Path to blacklist BED file.
#'   Pass \code{NULL} to skip filtering (a warning is issued).
#' @param threads Integer. Threads for re-indexing the filtered BAM. Default 4.
#' @param force Logical. If \code{FALSE} (default), skip filtering when
#'   \code{_blacklist_filtered.bam} already exists.  Set \code{TRUE} to
#'   re-filter and overwrite.
#'
#' @return Character. Path to the blacklist-filtered BAM file (or the
#'   unfiltered processed BAM if \code{blacklist_bed = NULL}), invisibly.
#' @export
HORIZON_filter_blacklist <- function(sample_sheet,
                                      sample_id,
                                      blacklist_bed,
                                      threads = 4L,
                                      force   = FALSE) {

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")
  input_bam <- .latest_bam(sample_sheet, sample_id)
  out_bam   <- file.path(out_dir, paste0(sample_id, "_blacklist_filtered.bam"))

  if (!isTRUE(force) && !is.null(blacklist_bed) && file.exists(out_bam)) {
    message("Blacklist-filtered BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to re-filter)")
    return(invisible(out_bam))
  }

  if (is.null(blacklist_bed)) {
    warning("blacklist_bed is NULL — skipping blacklist filtering for: ",
            sample_id, "\n",
            "  For T2T-CHM13 or other assemblies without an ENCODE blacklist,\n",
            "  this is expected. Provide a custom BED to suppress this warning.",
            call. = FALSE)
    return(invisible(input_bam))
  }

  if (!file.exists(blacklist_bed))
    stop("Blacklist BED file not found: ", blacklist_bed, call. = FALSE)

  message("[", sample_id, "] Filtering blacklist regions")
  exit <- .horizon_run_cli(
    "bedtools",
    c("intersect", "-v", "-abam", input_bam, "-b", blacklist_bed),
    stdout = out_bam
  )
  if (exit != 0)
    stop("bedtools intersect failed (exit code ", exit, ").", call. = FALSE)

  # Index filtered BAM
  .horizon_run_cli(
    "samtools",
    c("index", "--threads", as.integer(threads), out_bam)
  )

  message("Blacklist filtering complete for: ", sample_id,
          "\n  Output: ", out_bam)
  invisible(out_bam)
}


#' Downsample a BAM file by a spike-in derived scale factor
#'
#' Randomly subsamples the most downstream BAM for a sample using
#' \code{samtools view -s}, producing a \code{<sample_id>_downsampled.bam}
#' that is picked up automatically by all downstream pipeline functions via
#' \code{.latest_bam()}.
#'
#' The \code{scale_factor} should come from the \code{scale_factor_bw} column
#' of \code{\link{HORIZON_compute_spike_in_factors}}: it is the ratio
#' \eqn{\min(N) / N_i} (\eqn{\leq 1}), where the sample with the fewest
#' spike-in reads receives a factor of 1 (no downsampling) and all others are
#' scaled down proportionally.
#'
#' When \code{scale_factor = 1}, no reads are discarded.  A copy of the input
#' BAM is written under the \code{_downsampled.bam} name so that all samples
#' in a batch have a consistent output file and \code{.latest_bam()} behaves
#' identically for every sample.
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param scale_factor Numeric in \code{(0, 1]}.  Downsampling fraction.
#'   Supply the \code{scale_factor_bw} value from
#'   \code{\link{HORIZON_compute_spike_in_factors}}.
#' @param seed Integer. Random seed for \code{samtools view -s}.  Using the
#'   same seed across samples ensures reproducibility.  Default \code{42L}.
#' @param threads Integer. Threads for samtools. Default \code{4L}.
#' @param remove_input Logical. Remove the input BAM (and its index) after
#'   downsampling.  Default \code{FALSE}.
#' @param force Logical. If \code{FALSE} (default), skip downsampling when
#'   \code{_downsampled.bam} already exists.  Set \code{TRUE} to redo.
#'
#' @return Character. Path to the downsampled BAM file, invisibly.
#' @export
HORIZON_downsample_bam <- function(sample_sheet,
                                    sample_id,
                                    scale_factor,
                                    seed         = 42L,
                                    threads      = 4L,
                                    remove_input = FALSE,
                                    force        = FALSE) {

  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      scale_factor <= 0 || scale_factor > 1)
    stop("scale_factor must be a single numeric value in (0, 1].", call. = FALSE)

  row        <- .get_sample_row(sample_sheet, sample_id)
  out_dir    <- file.path(row$output_dir, sample_id, "aligned")
  output_bam <- file.path(out_dir, paste0(sample_id, "_downsampled.bam"))

  if (!isTRUE(force) && file.exists(output_bam)) {
    message("Downsampled BAM already exists for: ", sample_id,
            " — skipping (use force=TRUE to redo)")
    return(invisible(output_bam))
  }

  input_bam <- .latest_bam(sample_sheet, sample_id)
  thr       <- as.integer(threads)

  if (scale_factor >= 1.0) {
    message("[", sample_id, "] scale_factor = 1 — no downsampling needed; ",
            "copying BAM for pipeline consistency")
    file.copy(input_bam, output_bam, overwrite = TRUE)
    bai_in <- paste0(input_bam, ".bai")
    if (file.exists(bai_in))
      file.copy(bai_in, paste0(output_bam, ".bai"), overwrite = TRUE)
  } else {
    # samtools -s flag: integer part = seed, decimal part = fraction
    # e.g. seed=42, fraction=0.75 → "-s 42.7500"
    frac_str <- sub("^0", "", formatC(scale_factor, digits = 4L, format = "f"))
    s_flag   <- paste0(as.integer(seed), frac_str)

    message("[", sample_id, "] Downsampling BAM (fraction = ",
            round(scale_factor, 4L), ", seed = ", seed, ")")
    exit <- .horizon_run_cli(
      "samtools",
      c("view", "-b", "--threads", thr, "-s", s_flag, "-o", output_bam, input_bam)
    )
    if (exit != 0)
      stop("samtools view (downsample) failed for: ", sample_id, call. = FALSE)
  }

  .horizon_run_cli("samtools", c("index", "--threads", thr, output_bam))

  if (isTRUE(remove_input)) {
    file.remove(input_bam)
    bai <- paste0(input_bam, ".bai")
    if (file.exists(bai)) file.remove(bai)
  }

  message("Downsampling complete for: ", sample_id,
          "\n  Output: ", output_bam)
  invisible(output_bam)
}


# ---------------------------------------------------------------------------
# .latest_bam — internal helper
#
# Returns the path to the most downstream processed BAM that exists for a
# given sample.  Priority (highest to lowest):
#   1. <sample_id>_downsampled.bam
#   2. <sample_id>_blacklist_filtered.bam
#   3. <sample_id>_processed.bam
#   4. <sample_id>_host.bam
#   5. <sample_id>_mapq_filtered.bam
# ---------------------------------------------------------------------------
.latest_bam <- function(sample_sheet, sample_id) {
  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")
  candidates <- c(
    file.path(out_dir, paste0(sample_id, "_downsampled.bam")),
    file.path(out_dir, paste0(sample_id, "_blacklist_filtered.bam")),
    file.path(out_dir, paste0(sample_id, "_processed.bam")),
    file.path(out_dir, paste0(sample_id, "_host.bam")),
    file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0)
    stop("No processed BAM found for sample: ", sample_id,
         "\n  Expected one of:\n  ",
         paste(candidates, collapse = "\n  "),
         call. = FALSE)
  found[1]
}
