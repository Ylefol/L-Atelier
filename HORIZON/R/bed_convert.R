#' Convert a processed BAM to a fragment-level BED file
#'
#' Produces a 6-column fragment BED from a paired-end BAM.  Optionally applies
#' the standard Tn5 insertion-site correction (+4 bp on the + strand,
#' −5 bp on the − strand) to shift each fragment toward the actual cut site.
#'
#' The output BED is directly compatible with
#' \code{ELEUTHIA_load_bed()} and \code{ELEUTHIA_merge_fragments()} in GAIA.
#'
#' Internally:
#' \enumerate{
#'   \item Name-sorts the BAM (\code{samtools sort -n}).
#'   \item Converts to BEDPE (\code{bedtools bamtobed -bedpe}).
#'   \item Streams the BEDPE through \code{awk} to build fragment coordinates
#'         and optionally apply the Tn5 shift, then sorts with \code{sort}.
#' }
#'
#' \strong{When to use \code{shift_reads}:}
#' \itemize{
#'   \item ATAC-seq — \code{TRUE}: the Tn5 transposase cuts and ligates adapters,
#'         so the read start is offset from the actual insertion site by 4/5 bp.
#'   \item ChIP-seq / CUT&RUN / CUT&TAG — \code{FALSE}: no insertion-site
#'         correction is needed; fragment coordinates are used as-is.
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame.  Required when
#'   \code{bam_file} is \code{NULL}.
#' @param sample_id Character. Sample ID to process.
#' @param bam_file Character or \code{NULL}.  Path to a coordinate-sorted,
#'   indexed BAM file.  When supplied, \code{.latest_bam()} is bypassed and
#'   this BAM is used directly.  Output is written to \code{dirname(bam_file)}.
#'   Default \code{NULL}.
#' @param tag Character or \code{NULL}.  Optional suffix inserted before
#'   \code{_fragments.bed} in the output filename, producing
#'   \code{<sample_id>_<tag>_fragments.bed}.  Use this to avoid overwriting
#'   when calling the function twice for the same sample at different pipeline
#'   stages (e.g. \code{tag = "quantification"} for the pre-downsampled BED
#'   used by DESeq2, leaving the untagged file for peak calling).  Default
#'   \code{NULL} (no tag — output is \code{<sample_id>_fragments.bed}).
#' @param shift_reads Logical. Apply Tn5 insertion-site correction
#'   (+4/−5 bp).  Set \code{TRUE} for ATAC-seq, \code{FALSE} for all other
#'   assays.  Default \code{TRUE}.
#' @param threads Integer. Threads for \code{samtools sort}. Default 4.
#' @param remove_tmp Logical. Remove intermediate name-sorted BAM and BEDPE
#'   after conversion.  Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip conversion when the
#'   output BED already exists.  Set \code{TRUE} to regenerate.
#'
#' @return Character. Path to the fragment BED file, invisibly.
#' @export
HORIZON_bam_to_bed <- function(sample_sheet,
                                sample_id,
                                bam_file    = NULL,
                                tag         = NULL,
                                shift_reads = TRUE,
                                threads     = 4L,
                                remove_tmp  = TRUE,
                                force       = FALSE) {

  if (!is.null(bam_file)) {
    if (!file.exists(bam_file))
      stop("BAM file not found: ", bam_file, call. = FALSE)
    bam_path <- bam_file
    out_dir  <- dirname(bam_file)
  } else {
    bam_path <- .latest_bam(sample_sheet, sample_id)
    row      <- .get_sample_row(sample_sheet, sample_id)
    out_dir  <- file.path(row$output_dir, sample_id, "aligned")
  }

  bed_stem <- if (!is.null(tag)) paste0(sample_id, "_", tag) else sample_id
  ns_bam   <- file.path(out_dir, paste0(bed_stem, "_namesorted.bam"))
  bedpe_f  <- file.path(out_dir, paste0(bed_stem, "_raw.bedpe"))
  out_bed  <- file.path(out_dir, paste0(bed_stem, "_fragments.bed"))

  if (!isTRUE(force) && file.exists(out_bed)) {
    message("Fragment BED already exists for: ", sample_id,
            " — skipping (use force=TRUE to regenerate)")
    return(invisible(out_bed))
  }

  # Step 1: name-sort BAM -------------------------------------------------------
  message("[", sample_id, "] Name-sorting BAM for BEDPE conversion")
  exit <- .horizon_run_cli(
    "samtools",
    c("sort", "-n", "--threads", as.integer(threads),
      bam_path, "-o", ns_bam)
  )
  if (exit != 0) stop("samtools sort -n failed.", call. = FALSE)

  # Step 2: BAM → BEDPE --------------------------------------------------------
  message("[", sample_id, "] Converting BAM to BEDPE")
  exit <- .horizon_run_cli(
    "bedtools",
    c("bamtobed", "-bedpe", "-i", ns_bam),
    stdout = bedpe_f
  )
  if (exit != 0) stop("bedtools bamtobed failed.", call. = FALSE)

  if (isTRUE(remove_tmp)) file.remove(ns_bam)

  # Step 3: awk — BEDPE → fragment BED (streaming, no R memory) ---------------
  # Writes an awk script to a temp file to avoid shell-quoting issues, then
  # streams the BEDPE through awk into an unsorted temp BED.
  message("[", sample_id, "] Building fragment BED (streaming awk)")

  awk_script <- tempfile(fileext = ".awk")
  if (isTRUE(shift_reads)) {
    writeLines(c(
      'BEGIN { FS="\t"; OFS="\t" }',
      '$1 == $4 {',
      '  start = ($2 < $5) ? $2 : $5',
      '  end   = ($3 > $6) ? $3 : $6',
      '  if ($9 == "+") { start += 4; end += 4 }',
      '  else { start = (start >= 5) ? start - 5 : 0; end = (end >= 5) ? end - 5 : 0 }',
      '  if (end > start) print $1, start, end, NR, 0, $9',
      '}'
    ), awk_script)
  } else {
    writeLines(c(
      'BEGIN { FS="\t"; OFS="\t" }',
      '$1 == $4 {',
      '  start = ($2 < $5) ? $2 : $5',
      '  end   = ($3 > $6) ? $3 : $6',
      '  if (end > start) print $1, start, end, NR, 0, $9',
      '}'
    ), awk_script)
  }

  tmp_bed <- file.path(out_dir, paste0(bed_stem, "_unsorted.bed"))
  exit <- system2("awk", args = c("-f", awk_script, bedpe_f),
                  stdout = tmp_bed, stderr = "")
  file.remove(awk_script)
  if (exit != 0) stop("awk fragment conversion failed.", call. = FALSE)
  if (isTRUE(remove_tmp)) file.remove(bedpe_f)

  # Step 4: sort by chr, start -------------------------------------------------
  message("[", sample_id, "] Sorting fragment BED")
  exit <- system2("sort", args = c("-k1,1", "-k2,2n", tmp_bed),
                  stdout = out_bed, stderr = "")
  file.remove(tmp_bed)
  if (exit != 0) stop("sort failed for fragment BED.", call. = FALSE)

  message("Fragment BED complete for: ", sample_id,
          "\n  Output: ", out_bed)
  invisible(out_bed)
}

