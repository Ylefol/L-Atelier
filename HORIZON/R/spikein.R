# ==============================================================================
# HORIZON — Spike-in normalization helpers
#
# Functions for CUT&RUN / CUT&TAG spike-in calibration normalization.
# Workflow:
#   1. HORIZON_build_combined_index()   — concatenate host + spike-in FASTAs,
#                                         build Bowtie2 index
#   2. HORIZON_run_bowtie2()            — align to combined index (unchanged)
#   3. HORIZON_separate_spike_in()      — split BAM into host + spike-in BAMs
#   4. HORIZON_process_bam() onwards    — operate on host BAM only
#   5. HORIZON_compute_spike_in_factors() — count spike-in reads, compute
#                                           per-sample scale factors
#   6. HORIZON_bam_to_bigwig(..., scale_factor=) — spike-in normalised BigWig
# ==============================================================================


#' Build a Bowtie2 index from a combined host + spike-in genome
#'
#' Concatenates a host reference FASTA with a spike-in reference FASTA,
#' prefixing all spike-in chromosome names so that host and spike-in reads can
#' be cleanly separated after alignment.  The combined FASTA is written to disk
#' and reused on subsequent calls (set \code{overwrite = TRUE} to regenerate).
#' A Bowtie2 index is then built from the combined FASTA via
#' \code{\link{HORIZON_build_bowtie2_index}}.
#'
#' \strong{Why a combined index?}  Aligning to a single combined reference
#' forces reads to compete for the better mapping location.  Reads that
#' originate from the spike-in genome map uniquely to spike-in chromosomes,
#' while host reads align to host chromosomes.  This is more accurate than
#' aligning to each genome separately.
#'
#' @param host_fasta Character. Path to the host reference FASTA
#'   (e.g. \file{GRCh38.fa}).
#' @param spikein_fasta Character. Path to the spike-in reference FASTA
#'   (e.g. \file{sacCer3.fa} for yeast).
#' @param combined_dir Character. Directory where the combined FASTA is written.
#' @param index_dir Character. Directory for the Bowtie2 index files.
#'   Defaults to \code{combined_dir}.
#' @param index_name Character. Basename prefix for the Bowtie2 index and the
#'   combined FASTA file.  Default \code{"combined"}.
#' @param spikein_prefix Character. String prepended to every spike-in
#'   chromosome name (\code{>header} lines in the FASTA).  Must not be a
#'   prefix of any host chromosome name.  Default \code{"spikein_"}.
#'   Pass the same value to \code{\link{HORIZON_separate_spike_in}}.
#' @param threads Integer. Threads for \code{bowtie2-build}. Default 4.
#' @param memory Integer. Memory in MB for \code{bowtie2-build}. Default 8000.
#' @param overwrite Logical. Regenerate the combined FASTA even if it already
#'   exists.  Default \code{FALSE}.
#'
#' @return Character. Bowtie2 index basename (pass to
#'   \code{\link{HORIZON_run_bowtie2}}), invisibly.
#' @export
HORIZON_build_combined_index <- function(host_fasta,
                                          spikein_fasta,
                                          combined_dir,
                                          index_dir      = combined_dir,
                                          index_name     = "combined",
                                          spikein_prefix = "spikein_",
                                          threads        = 4L,
                                          memory         = 8000L,
                                          overwrite      = FALSE) {

  if (!file.exists(host_fasta))
    stop("Host FASTA not found: ", host_fasta, call. = FALSE)
  if (!file.exists(spikein_fasta))
    stop("Spike-in FASTA not found: ", spikein_fasta, call. = FALSE)
  if (!nzchar(spikein_prefix))
    stop("spikein_prefix must be a non-empty string.", call. = FALSE)

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(index_dir,    recursive = TRUE, showWarnings = FALSE)

  combined_fasta <- file.path(combined_dir, paste0(index_name, "_combined.fa"))

  if (file.exists(combined_fasta) && !isTRUE(overwrite)) {
    cat("[HORIZON] Combined FASTA already exists: ", combined_fasta,
            "\n  Set overwrite=TRUE to regenerate.")
  } else {
    cat("[HORIZON] Building combined FASTA...")
    cat("  Host:     ", host_fasta)
    cat("  Spike-in: ", spikein_fasta,
            "  [chromosome prefix = '", spikein_prefix, "']")

    # Copy host FASTA unchanged — avoids loading a large file into R memory
    file.copy(host_fasta, combined_fasta, overwrite = TRUE)

    # Stream spike-in FASTA, prefixing header lines, and append to combined
    con_in  <- file(spikein_fasta, "r")
    con_out <- file(combined_fasta, "a")
    on.exit({ try(close(con_in), silent = TRUE)
              try(close(con_out), silent = TRUE) }, add = TRUE)

    repeat {
      lines <- readLines(con_in, n = 10000L, warn = FALSE)
      if (length(lines) == 0L) break
      hdr <- startsWith(lines, ">")
      lines[hdr] <- sub("^>", paste0(">", spikein_prefix), lines[hdr])
      writeLines(lines, con_out)
    }

    close(con_in);  on.exit(NULL, add = FALSE)
    close(con_out); on.exit(NULL, add = FALSE)

    cat("[HORIZON] Combined FASTA written: ", combined_fasta)
  }

  # Build Bowtie2 index from the combined FASTA
  # Note: memory parameter is not passed — bowtie2-build has no memory control
  # (unlike Rsubread::buildindex). Memory usage during index build is determined
  # by bowtie2-build internally based on genome size.
  idx <- HORIZON_build_bowtie2_index(
    reference  = combined_fasta,
    index_dir  = index_dir,
    index_name = index_name,
    threads    = threads
  )

  invisible(idx)
}


#' Split a combined BAM into host and spike-in BAMs
#'
#' After aligning to a combined host + spike-in genome with
#' \code{\link{HORIZON_run_bowtie2}}, this function splits the
#' \code{_mapq_filtered.bam} into two BAM files by chromosome prefix:
#' \itemize{
#'   \item \code{<sample_id>_host.bam}    — host reads (feed to
#'         \code{\link{HORIZON_process_bam}} and downstream)
#'   \item \code{<sample_id>_spikein.bam} — spike-in reads (used only for
#'         counting in \code{\link{HORIZON_compute_spike_in_factors}})
#' }
#' The host BAM is automatically picked up by \code{.latest_bam()} and all
#' subsequent pipeline steps without any changes to those functions.
#'
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param spikein_prefix Character. Prefix that identifies spike-in chromosomes
#'   in the BAM header.  Must match the value used in
#'   \code{\link{HORIZON_build_combined_index}}.  Default \code{"spikein_"}.
#' @param threads Integer. Threads for samtools commands. Default 4.
#' @param remove_combined Logical. Remove the combined
#'   \code{_mapq_filtered.bam} after separation to save disk space.
#'   Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip separation when both
#'   \code{_host.bam} and \code{_spikein.bam} already exist.  Set \code{TRUE}
#'   to re-separate and overwrite.
#'
#' @return Named list with elements \code{host} and \code{spikein} giving the
#'   paths to the two output BAMs, invisibly.
#' @export
HORIZON_separate_spike_in <- function(sample_sheet,
                                       sample_id,
                                       spikein_prefix  = "spikein_",
                                       threads         = 4L,
                                       remove_combined = TRUE,
                                       force           = FALSE) {

  row     <- .get_sample_row(sample_sheet, sample_id)
  out_dir <- file.path(row$output_dir, sample_id, "aligned")

  combined_bam <- file.path(out_dir, paste0(sample_id, "_mapq_filtered.bam"))
  host_bam     <- file.path(out_dir, paste0(sample_id, "_host.bam"))
  spikein_bam  <- file.path(out_dir, paste0(sample_id, "_spikein.bam"))

  if (!isTRUE(force) && file.exists(host_bam) && file.exists(spikein_bam)) {
    cat("Separated BAMs already exist for: ", sample_id,
            " — skipping (use force=TRUE to re-separate)")
    return(invisible(list(host = host_bam, spikein = spikein_bam)))
  }

  if (!file.exists(combined_bam))
    stop("Combined BAM not found: ", combined_bam,
         ". Run HORIZON_run_bowtie2() first.", call. = FALSE)

  thr <- as.integer(threads)

  # Index combined BAM if needed
  if (!file.exists(paste0(combined_bam, ".bai"))) {
    cat("[", sample_id, "] Indexing combined BAM")
    .horizon_run_cli("samtools",
                      c("index", "--threads", thr, combined_bam))
  }

  # Read chromosome names from BAM header
  header_lines <- .horizon_run_cli(
    "samtools",
    c("view", "-H", combined_bam),
    stdout = TRUE
  )
  sq_lines  <- grep("^@SQ\t", header_lines, value = TRUE)
  chr_names <- sub(".*\tSN:([^\t]+).*", "\\1", sq_lines)

  host_chrs    <- chr_names[!startsWith(chr_names, spikein_prefix)]
  spikein_chrs <- chr_names[ startsWith(chr_names, spikein_prefix)]

  if (length(spikein_chrs) == 0L)
    stop("No chromosomes with prefix '", spikein_prefix,
         "' found in BAM header.\n",
         "  Found: ", paste(head(chr_names, 10L), collapse = ", "),
         if (length(chr_names) > 10L) "..." else "",
         "\n  Ensure HORIZON_build_combined_index() was called with the ",
         "same spikein_prefix.", call. = FALSE)

  cat("[", sample_id, "] Separating reads: ",
          length(host_chrs), " host chr, ",
          length(spikein_chrs), " spike-in chr")

  # Extract host reads
  exit <- .horizon_run_cli(
    "samtools",
    c("view", "-b", "--threads", thr, "-o", host_bam, combined_bam, host_chrs)
  )
  if (exit != 0) stop("samtools view (host) failed.", call. = FALSE)
  .horizon_run_cli("samtools",
                    c("index", "--threads", thr, host_bam))

  # Extract spike-in reads
  exit <- .horizon_run_cli(
    "samtools",
    c("view", "-b", "--threads", thr, "-o", spikein_bam, combined_bam, spikein_chrs)
  )
  if (exit != 0) stop("samtools view (spike-in) failed.", call. = FALSE)
  .horizon_run_cli("samtools",
                    c("index", "--threads", thr, spikein_bam))

  if (isTRUE(remove_combined)) {
    file.remove(combined_bam)
    bai <- paste0(combined_bam, ".bai")
    if (file.exists(bai)) file.remove(bai)
  }

  cat("Separation complete for: ", sample_id,
          "\n  Host BAM:     ", host_bam,
          "\n  Spike-in BAM: ", spikein_bam)
  invisible(list(host = host_bam, spikein = spikein_bam))
}


#' Compute per-sample scale factors from spike-in read counts
#'
#' Compute per-sample scale factors from spike-in read counts
#'
#' Counts aligned read pairs in each spike-in BAM (produced by
#' \code{\link{HORIZON_separate_spike_in}}) and returns a data frame containing
#' two normalization factors per sample derived from the same underlying spike-in
#' read counts, for use at different stages of the analysis:
#'
#' \describe{
#'   \item{\code{scale_factor_bw}}{For \code{\link{HORIZON_bam_to_bigwig}}.
#'     Formula: \eqn{\min(N) / N_i}.  Values \eqn{\leq 1}; the
#'     least-sequenced sample gets 1.0 and all others are scaled down.
#'     Passed to \code{bamCoverage --scaleFactor}.}
#'   \item{\code{size_factor_deseq2}}{For \code{DESeq2::sizeFactors()}.
#'     Formula: \eqn{N_i / \min(N)}.  Values \eqn{\geq 1}; DESeq2 divides
#'     raw counts by this value, so the least-sequenced sample (factor = 1.0)
#'     is unchanged and more deeply sequenced samples are scaled down.
#'     Assign via \code{sizeFactors(dds) <- factors$size_factor_deseq2} to
#'     bypass DESeq2's internal normalization and use spike-in calibration
#'     instead.  This is conceptually equivalent to RNA-seq spike-in
#'     normalization but computed from whole-genome read counts rather than
#'     spike-in gene rows in a count matrix.}
#' }
#'
#' @param spikein_bams Named character vector or list.  Names are sample IDs;
#'   values are paths to spike-in BAMs from
#'   \code{\link{HORIZON_separate_spike_in}}.
#' @param output_file Character or \code{NULL}.  Path to write a CSV of all
#'   factors (e.g. \file{output/spikein_factors.csv}).  The file can be read
#'   back later to supply factors to DESeq2 without re-running this function.
#'   Default \code{NULL} (no file written).
#' @param threads Integer. Threads for \code{samtools view}. Default 4.
#'
#' @return A data frame with columns \code{sample_id}, \code{spikein_reads},
#'   \code{scale_factor_bw}, \code{size_factor_deseq2}, invisibly.
#'   Samples with zero spike-in reads receive \code{NA} factors with a warning.
#' @export
HORIZON_compute_spike_in_factors <- function(spikein_bams,
                                              output_file = NULL,
                                              threads     = 4L) {

  if (is.null(names(spikein_bams)) || any(names(spikein_bams) == ""))
    stop("spikein_bams must be a fully named vector or list (names = sample IDs).",
         call. = FALSE)

  thr <- as.integer(threads)

  # Count aligned read pairs: -f 64 = read1, -F 4 = mapped → fragment count
  counts <- vapply(names(spikein_bams), function(sid) {
    bam <- spikein_bams[[sid]]
    if (!file.exists(bam))
      stop("Spike-in BAM not found for '", sid, "': ", bam, call. = FALSE)
    out <- .horizon_run_cli(
      "samtools",
      c("view", "-c", "-F", "4", "-f", "64", "--threads", thr, bam),
      stdout = TRUE
    )
    n <- suppressWarnings(as.integer(trimws(out[[1L]])))
    if (is.na(n))
      stop("Failed to count reads in spike-in BAM for: ", sid, call. = FALSE)
    n
  }, integer(1L))

  cat("Spike-in fragment counts:")
  for (sid in names(counts))
    cat("  ", sid, ": ", format(counts[[sid]], big.mark = ","))

  if (any(counts == 0L)) {
    zero_samps <- names(counts)[counts == 0L]
    warning("Zero spike-in reads in: ", paste(zero_samps, collapse = ", "),
            "\n  These samples will receive NA factors. ",
            "Check alignment to the combined genome.", call. = FALSE)
  }

  min_count <- min(counts[counts > 0L])

  # scale_factor_bw:    min / count  → multiply BigWig signal (≤ 1)
  # size_factor_deseq2: count / min  → divide DESeq2 counts   (≥ 1)
  factors <- data.frame(
    sample_id         = names(counts),
    spikein_reads     = as.integer(counts),
    scale_factor_bw   = ifelse(counts > 0L, min_count / counts, NA_real_),
    size_factor_deseq2 = ifelse(counts > 0L, counts / min_count, NA_real_),
    stringsAsFactors  = FALSE
  )

  cat("Normalization factors:")
  cat(sprintf("  %-20s  %12s  %14s  %18s",
                  "sample_id", "spikein_reads", "scale_factor_bw",
                  "size_factor_deseq2"))
  for (i in seq_len(nrow(factors))) {
    cat(sprintf("  %-20s  %12s  %14s  %18s",
                    factors$sample_id[i],
                    format(factors$spikein_reads[i], big.mark = ","),
                    if (is.na(factors$scale_factor_bw[i])) "NA"
                    else round(factors$scale_factor_bw[i], 5L),
                    if (is.na(factors$size_factor_deseq2[i])) "NA"
                    else round(factors$size_factor_deseq2[i], 5L)))
  }

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(factors, output_file, row.names = FALSE)
    cat("Spike-in factors written to: ", output_file)
  }

  invisible(factors)
}
