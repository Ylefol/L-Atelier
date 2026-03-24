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
#'   \item Reads the BEDPE in R, builds fragment coordinates, optionally
#'         applies the Tn5 shift, then writes a sorted 6-column BED.
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
#' @param sample_sheet Validated sample sheet data.frame.
#' @param sample_id Character. Sample ID to process.
#' @param shift_reads Logical. Apply Tn5 insertion-site correction
#'   (+4/−5 bp).  Set \code{TRUE} for ATAC-seq, \code{FALSE} for all other
#'   assays.  Default \code{TRUE}.
#' @param threads Integer. Threads for \code{samtools sort}. Default 4.
#' @param remove_tmp Logical. Remove intermediate name-sorted BAM and BEDPE
#'   after conversion.  Default \code{TRUE}.
#' @param force Logical. If \code{FALSE} (default), skip conversion when
#'   \code{_fragments.bed} already exists.  Set \code{TRUE} to regenerate.
#'
#' @return Character. Path to the fragment BED file, invisibly.
#' @export
HORIZON_bam_to_bed <- function(sample_sheet,
                                sample_id,
                                shift_reads = TRUE,
                                threads     = 4L,
                                remove_tmp  = TRUE,
                                force       = FALSE) {

  bam_path <- .latest_bam(sample_sheet, sample_id)
  row      <- .get_sample_row(sample_sheet, sample_id)
  out_dir  <- file.path(row$output_dir, sample_id, "aligned")

  ns_bam   <- file.path(out_dir, paste0(sample_id, "_namesorted.bam"))
  bedpe_f  <- file.path(out_dir, paste0(sample_id, "_raw.bedpe"))
  out_bed  <- file.path(out_dir, paste0(sample_id, "_fragments.bed"))

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

  # Step 3: R: BEDPE → fragment BED -------------------------------------------
  message("[", sample_id, "] Building fragment BED")
  bedpe <- utils::read.table(bedpe_f, header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE,
                              col.names = c("chr1", "start1", "end1",
                                            "chr2", "start2", "end2",
                                            "name", "score",
                                            "strand1", "strand2"))
  if (isTRUE(remove_tmp)) file.remove(bedpe_f)

  # Keep only concordant pairs (same chromosome)
  bedpe <- bedpe[bedpe$chr1 == bedpe$chr2, ]

  # Fragment coordinates: span both mates
  frag <- data.frame(
    chr    = bedpe$chr1,
    start  = pmin(bedpe$start1, bedpe$start2),
    end    = pmax(bedpe$end1,   bedpe$end2),
    name   = seq_len(nrow(bedpe)),
    score  = 0L,
    strand = bedpe$strand1,
    stringsAsFactors = FALSE
  )

  # Step 4: Tn5 shift (optional) -----------------------------------------------
  if (isTRUE(shift_reads)) {
    plus_idx  <- frag$strand == "+"
    minus_idx <- frag$strand == "-"
    frag$start[plus_idx] <- frag$start[plus_idx] + 4L
    frag$end[plus_idx]   <- frag$end[plus_idx]   + 4L
    frag$start[minus_idx] <- pmax(0L, frag$start[minus_idx] - 5L)
    frag$end[minus_idx]   <- pmax(0L, frag$end[minus_idx]   - 5L)
  }

  # Remove any zero-width or negative-width fragments after shift
  frag <- frag[frag$end > frag$start, ]

  # Step 5: sort by chr, start -------------------------------------------------
  frag <- frag[order(frag$chr, frag$start), ]

  # Write 6-column BED ---------------------------------------------------------
  utils::write.table(frag, file = out_bed, quote = FALSE,
                     sep = "\t", row.names = FALSE, col.names = FALSE)

  message("Fragment BED complete for: ", sample_id,
          "\n  Output: ", out_bed,
          "\n  Fragments: ", nrow(frag))
  invisible(out_bed)
}

