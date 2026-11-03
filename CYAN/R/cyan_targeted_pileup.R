#' Targeted per-base pileup at a fixed list of genomic positions
#'
#' Computes read coverage and per-base (A/C/G/T) read counts at a short,
#' sparse list of already-known positions (e.g. RNA-editing candidate sites,
#' when confirming germline status against a matched WGS BAM). Unlike
#' \code{\link{CYAN_run_reditools}}, which scans every read across the whole
#' BAM one base at a time regardless of \code{bed_file} (that option only
#' filters what gets WRITTEN, not what gets SCANNED), this seeks directly
#' to each BED site via pysam's indexed pileup -- turning a whole-genome
#' scan into \code{O(sites)} indexed lookups.
#'
#' @details
#' This is purpose-built for the germline-confirmation use case
#' (\code{\link{CYAN_flag_germline_variants}}): it does not detect
#' substitutions, infer strand, or apply any of REDItools2's editing-specific
#' logic (omopolymeric/splice-site exclusion, frequency calculation). It only
#' answers "how many reads support each base at this exact position" --
#' every row of \code{bed_file} is returned exactly once, including
#' zero-coverage sites, which is the information germline-status calls
#' need (REDItools2 would simply omit positions it finds no coverage at,
#' and it would discard any analogue of a "no substitution at the DNA
#' position" row when restricted to the editing-sites output).
#'
#' Read/base-level filtering mirrors REDItools2's own
#' \code{filter_read()}/\code{filter_base()} logic (unmapped, QC-fail,
#' mapping quality, read length, secondary/supplementary, PCR duplicate,
#' non-proper-pair, chimeric/\code{SA}-tagged reads, base quality) so
#' coverage here stays comparable to the RNA-side REDItools2 numbers.
#'
#' As with \code{\link{CYAN_run_reditools}}, an existing raw output
#' checkpoint in \code{output_dir} for the same sample is reused rather than
#' re-run (same caveat: this does not check whether it matches the
#' parameters you are calling with now).
#'
#' @param bam             Character. Path to a sorted, indexed BAM file.
#' @param reference       Character or \code{NULL}. Path to the reference
#'   FASTA (must be indexed with \code{samtools faidx}). Only used to report
#'   the reference base per position; coverage/counts do not depend on it.
#' @param bed_file        Character. Path to a BED file of target positions
#'   (0-based, half-open -- as written by \code{\link{CYAN_check_bam}}-style
#'   tooling or \code{.cyan_write_bed}).
#' @param output_dir      Character. Directory for final result files
#'   (created if absent).
#' @param work_dir        Character. Directory for intermediate files
#'   (created if absent). Defaults to \code{tempdir()}.
#' @param sample_name     Character or \code{NULL}. Label used in output
#'   filenames. Derived from the BAM filename if \code{NULL}.
#' @param min_base_quality Integer. Minimum per-base quality score (default 30,
#'   matching \code{\link{CYAN_run_reditools}}'s native default).
#' @param min_read_quality Integer. Minimum read mapping quality (default 20).
#' @param min_read_length  Integer. Minimum read length; shorter reads are
#'   discarded (default 30).
#' @param require_dup_marked Logical. If \code{TRUE} (default), hard-stops
#'   before running unless the BAM has at least one read flagged as a
#'   duplicate (SAM FLAG 0x400) -- see
#'   \code{\link{CYAN_run_reditools}}'s parameter of the same name for the
#'   full rationale; the duplicate filter here relies on the same upstream
#'   marking.
#' @param keep_tmp        Logical. If \code{FALSE} (default), intermediate
#'   files in \code{work_dir} are removed after completion.
#' @param save_rds        Logical. If \code{TRUE} (default), the parsed
#'   result is saved as an RDS file in \code{output_dir}.
#' @param verbose         Logical. Print \code{[CYAN]} progress messages
#'   (default \code{TRUE}).
#' @param debug           Logical. If \code{TRUE}, passes \code{-V} and
#'   routes the Python script's stderr to the console (default \code{FALSE}).
#'
#' @return A \code{data.frame} with one row per \code{bed_file} position and
#'   columns \code{region}, \code{position} (1-based), \code{reference},
#'   \code{coverage}, \code{count_A}, \code{count_C}, \code{count_G},
#'   \code{count_T}.
#'
#' @seealso \code{\link{CYAN_run_reditools}}, \code{\link{CYAN_flag_germline_variants}}
#'
#' @export
CYAN_run_targeted_pileup <- function(
  bam,
  reference          = NULL,
  bed_file,
  output_dir,
  work_dir           = tempdir(),
  sample_name        = NULL,
  min_base_quality   = 30L,
  min_read_quality   = 20L,
  min_read_length    = 30L,
  require_dup_marked = TRUE,
  keep_tmp           = FALSE,
  save_rds           = TRUE,
  verbose            = TRUE,
  debug              = FALSE
) {
  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  bam <- CYAN_check_bam(bam)
  if (!is.null(reference)) reference <- CYAN_check_references(reference, bed_file)
  if (!file.exists(bed_file))
    stop("[CYAN] `bed_file` not found: ", bed_file, call. = FALSE)

  sample <- .cyan_sample_name(bam, sample_name)

  output_dir <- .cyan_ensure_dir(output_dir, "output directory")
  work_dir   <- .cyan_ensure_dir(work_dir,   "work directory")

  # ------------------------------------------------------------------
  # 2. Check for an existing raw checkpoint
  # ------------------------------------------------------------------
  raw_out <- file.path(output_dir, paste0(sample, "_targeted_pileup.txt"))
  rds_out <- file.path(output_dir, paste0(sample, "_targeted_pileup.rds"))

  if (file.exists(raw_out)) {
    finfo <- file.info(raw_out)
    if (verbose)
      message("[CYAN] Existing targeted-pileup output found (made ",
              format(finfo$mtime, "%Y-%m-%d %H:%M"), "), skipping re-run: ",
              raw_out, "\n  Delete the file to force a re-run.")
    result <- .cyan_parse_targeted_pileup_output(raw_out)
    if (save_rds) {
      saveRDS(result, rds_out)
      if (verbose) message("[CYAN] RDS saved: ", rds_out)
    }
    return(invisible(result))
  }

  # ------------------------------------------------------------------
  # 2b. Duplicate-marking check (same rationale as CYAN_run_reditools)
  # ------------------------------------------------------------------
  if (require_dup_marked && !.cyan_check_duplicates_marked(bam, verbose = verbose)) {
    stop("[CYAN] No duplicate-flagged reads found in this BAM (SAM FLAG 0x400, ",
         "scanned up to 5,000,000 reads).\n",
         "  This function's duplicate filter relies entirely on this flag already ",
         "being set upstream -- on a BAM that was never duplicate-marked, that filter ",
         "is a silent no-op.\n",
         "  Run 'samtools markdup' or Picard MarkDuplicates on this BAM first.\n",
         "  If you are confident this BAM genuinely has no duplicates to flag, pass ",
         "require_dup_marked = FALSE to skip this check.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 3. Build CLI arguments and run
  # ------------------------------------------------------------------
  script_path <- system.file("python", "targeted_pileup.py", package = "CYAN")
  if (!nzchar(script_path))
    stop("[CYAN] Cannot locate inst/python/targeted_pileup.py inside the CYAN package.",
         call. = FALSE)

  cli_args <- c(
    "-f",   bam,
    "-B",   bed_file,
    "-o",   raw_out,
    "-q",   as.integer(min_read_quality),
    "-bq",  as.integer(min_base_quality),
    "-mrl", as.integer(min_read_length)
  )
  if (!is.null(reference)) cli_args <- c(cli_args, "-r", reference)
  if (isTRUE(debug))       cli_args <- c(cli_args, "-V")

  if (verbose) message("[CYAN] Starting targeted pileup for sample: ", sample)

  .cyan_reditools_run(script_path, cli_args, verbose, debug)

  if (!file.exists(raw_out))
    stop("[CYAN] Targeted pileup completed but output file was not created: ",
         raw_out, call. = FALSE)

  if (verbose) message("[CYAN] Targeted pileup finished. Parsing output...")

  # ------------------------------------------------------------------
  # 4. Parse and save
  # ------------------------------------------------------------------
  result <- .cyan_parse_targeted_pileup_output(raw_out)

  if (save_rds) {
    saveRDS(result, rds_out)
    if (verbose) message("[CYAN] RDS saved: ", rds_out)
  }

  if (!keep_tmp) {
    tmp_files <- list.files(work_dir, full.names = TRUE)
    if (length(tmp_files) > 0) {
      file.remove(tmp_files)
      if (verbose) message("[CYAN] Temporary files removed from: ", work_dir)
    }
  }

  if (verbose)
    message("[CYAN] Done. ", nrow(result), " positions parsed into memory. Raw output: ", raw_out)

  invisible(result)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Parse targeted_pileup.py tab-delimited output into a data.frame
#' @keywords internal
.cyan_parse_targeted_pileup_output <- function(path) {
  raw <- data.table::fread(path, sep = "\t", header = TRUE, colClasses = "character")

  col_map <- c(
    "Region"              = "region",
    "Position"            = "position",
    "Reference"           = "reference",
    "BaseCount[A,C,G,T]"  = "base_count"
  )
  for (old in names(col_map)) {
    idx <- match(old, colnames(raw))
    if (!is.na(idx)) colnames(raw)[idx] <- col_map[[old]]
  }
  cov_idx <- grep("^Coverage-q", colnames(raw))
  if (length(cov_idx) == 1) colnames(raw)[cov_idx] <- "coverage"

  bc    <- gsub("[][ ]", "", raw$base_count)
  parts <- data.table::tstrsplit(bc, ",", fixed = TRUE)
  raw[, count_A := as.integer(parts[[1]])]
  raw[, count_C := as.integer(parts[[2]])]
  raw[, count_G := as.integer(parts[[3]])]
  raw[, count_T := as.integer(parts[[4]])]
  raw[, base_count := NULL]

  raw[, position := as.integer(position)]
  raw[, coverage := as.integer(coverage)]

  data.table::setcolorder(raw, c("region", "position", "reference", "coverage",
                                  "count_A", "count_C", "count_G", "count_T"))
  data.table::setDF(raw)
  raw
}
