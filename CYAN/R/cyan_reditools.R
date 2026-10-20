#' Run REDItools2 RNA editing detection
#'
#' Detects RNA editing events in a BAM file by running REDItools2, a
#' pysam-based pileup engine that identifies positions with nucleotide
#' substitutions relative to the reference genome. The canonical use case
#' is A-to-I (A>G) editing detection in RNA-seq data.
#'
#' @details
#' REDItools2 is run as a subprocess inside a managed basilisk Python
#' environment. The raw output (tab-delimited text) is written to
#' \code{output_dir} as a checkpoint before parsing, so the expensive
#' computation is never lost. The parsed \code{data.frame} is saved
#' alongside as an RDS file (opt-out with \code{save_rds = FALSE}).
#'
#' Before running, this checks \code{output_dir} for an existing raw
#' checkpoint for the same sample and reuses it if found, re-parsing rather
#' than re-running REDItools2 (reported with the file's modification time
#' and size). This does not check whether that file was generated with the
#' same parameters you are calling with now -- delete it first to force a
#' re-run with different settings; keep your own record of which parameters
#' produced it (e.g. via logs) if that matters.
#'
#' \strong{File outputs} (in \code{output_dir}):
#' \itemize{
#'   \item \code{{sample}_reditools.txt}      — full raw REDItools2 table
#'     (all covered positions)
#'   \item \code{{sample}_reditools_subs.txt} — substitution sites only
#'     (rows where \code{all_subs != "-"})
#'   \item \code{{sample}_reditools.rds}      — parsed \code{cyan_editing_result}
#'     (if \code{save_rds = TRUE})
#' }
#'
#' @section Parameter notes:
#' \describe{
#'   \item{\code{strict}}{When \code{TRUE} (default), only positions with at
#'     least one detected substitution are emitted. Set to \code{FALSE} to
#'     have REDItools2 report all covered positions in the raw \code{.txt}
#'     checkpoint (large output) -- note that the parsed/returned/RDS'd
#'     result is always substitution-rows-only regardless of this setting
#'     (see \code{\link{.cyan_parse_reditools_output}}); the non-strict rows
#'     are only ever available by reading the raw checkpoint directly.}
#'   \item{\code{min_coverage}}{Minimum number of reads covering a position
#'     (after base-quality filtering) for it to be reported.}
#'   \item{\code{min_edits}}{Minimum number of edited reads at a position.}
#'   \item{\code{strand_mode}}{0 = unstranded; 1 = second-strand (dUTP/fr-firststrand);
#'     2 = first-strand (fr-secondstrand). Match to your library protocol.}
#'   \item{\code{omopolymeric_file}}{Path to a file of homopolymeric regions to
#'     exclude — strongly recommended to reduce false positives. Can be generated
#'     by REDItools2 with the \code{-c} flag (not yet exposed here).}
#'   \item{\code{splicing_file}}{Path to a file of splice-site positions to
#'     exclude — recommended to avoid mismatch artefacts near splice junctions.}
#' }
#'
#' @param bam             Character. Path to a sorted, indexed BAM file.
#' @param reference       Character. Path to the reference FASTA (must be
#'   indexed with \code{samtools faidx}).
#' @param output_dir      Character. Directory for final result files
#'   (created if absent).
#' @param work_dir        Character. Directory for intermediate/temporary files
#'   (created if absent). Defaults to \code{tempdir()}. Can point to an
#'   external drive for large datasets.
#' @param sample_name     Character or \code{NULL}. Label used in output
#'   filenames. Derived from the BAM filename if \code{NULL}.
#' @param strand_mode     Integer. Library strandedness: 0 = unstranded (default),
#'   1 = second-strand, 2 = first-strand.
#' @param min_coverage    Integer. Minimum read coverage at a position (default 1,
#'   matching REDItools2 native default). Increase (e.g. 10) to suppress
#'   low-coverage noise post-hoc or here to reduce output size.
#' @param min_base_quality Integer. Minimum per-base quality score (default 30,
#'   matching REDItools2 native default).
#' @param min_read_quality Integer. Minimum read mapping quality (default 20,
#'   matching REDItools2 native default).
#' @param min_read_length  Integer. Minimum read length; shorter reads are
#'   discarded (default 30, matching REDItools2 native default).
#' @param min_edits       Integer. Minimum number of edited reads at a position
#'   (default 0, matching REDItools2 native default). Filtering on this value
#'   post-hoc is straightforward: \code{result[result$frequency * result$coverage >= 3, ]}.
#' @param strict          Logical. If \code{TRUE}, only positions with detected
#'   substitutions are written to the output. Default \code{FALSE} (matching
#'   REDItools2 native default) — all covered positions are reported. Set to
#'   \code{TRUE} to substantially reduce output size when only editing sites
#'   are of interest.
#' @param region          Character or \code{NULL}. Restrict analysis to a
#'   genomic region in samtools format (e.g. \code{"chr21:1-10000000"}).
#' @param omopolymeric_file Character or \code{NULL}. Path to homopolymeric
#'   regions file. Recommended to reduce false positives.
#' @param splicing_file   Character or \code{NULL}. Path to splice-site positions
#'   file. Recommended to reduce mismatch artefacts near junctions.
#' @param bed_file        Character or \code{NULL}. Path to a BED file of target
#'   regions. Analysis is restricted to these regions when supplied.
#' @param require_dup_marked Logical. If \code{TRUE} (default), hard-stops
#'   before running REDItools2 unless the BAM has at least one read flagged
#'   as a duplicate (SAM FLAG 0x400). REDItools2's own duplicate filter
#'   relies entirely on this flag already being set upstream (e.g. via
#'   \code{samtools markdup} or Picard \code{MarkDuplicates}) -- on a BAM
#'   that was never duplicate-marked, that filter is a silent no-op. Set to
#'   \code{FALSE} only if you are confident this BAM genuinely has no
#'   duplicates to flag.
#' @param keep_tmp        Logical. If \code{FALSE} (default), intermediate files
#'   in \code{work_dir} are removed after completion.
#' @param save_rds        Logical. If \code{TRUE} (default), the parsed result
#'   is saved as an RDS file in \code{output_dir}.
#' @param verbose         Logical. Print \code{[CYAN]} R-level progress messages
#'   (default \code{TRUE}). Does not affect REDItools2's own output.
#' @param debug           Logical. If \code{TRUE}, passes \code{-V} to REDItools2
#'   and routes its stdout/stderr to the console (default \code{FALSE}). Use this
#'   when a run fails and you need to see the raw Python output to diagnose the
#'   problem.
#'
#' @return An object of class \code{cyan_editing_result}: a \code{data.frame}
#'   with one row per \strong{substitution-carrying} position -- even when
#'   \code{strict = FALSE}, only rows with \code{all_subs != "-"} are parsed
#'   into this object (see \code{\link{.cyan_parse_reditools_output}}) --
#'   and columns:
#'   \describe{
#'     \item{region}{Chromosome / contig name}
#'     \item{position}{1-based genomic position}
#'     \item{reference}{Reference base}
#'     \item{strand}{Strand (0/1/2)}
#'     \item{coverage}{Read coverage after quality filtering}
#'     \item{mean_quality}{Mean base quality at this position}
#'     \item{count_A, count_C, count_G, count_T}{Per-base read counts}
#'     \item{all_subs}{Space-separated substitution codes (e.g. \code{"AG"},
#'       \code{"AG TC"}); \code{"-"} if no substitution detected}
#'     \item{frequency}{Fraction of reads supporting the most frequent
#'       substitution (0–1)}
#'   }
#'   The object additionally carries a \code{params} attribute listing the
#'   call parameters for reproducibility.
#'
#' @seealso \code{\link{CYAN_check_bam}}, \code{\link{CYAN_check_references}}
#'
#' @export
CYAN_run_reditools <- function(
  bam,
  reference,
  output_dir,
  work_dir          = tempdir(),
  sample_name       = NULL,
  strand_mode       = 0L,
  min_coverage      = 1L,
  min_base_quality  = 30L,
  min_read_quality  = 20L,
  min_read_length   = 30L,
  min_edits         = 0L,
  strict            = FALSE,
  region            = NULL,
  omopolymeric_file = NULL,
  splicing_file     = NULL,
  bed_file          = NULL,
  require_dup_marked = TRUE,
  keep_tmp          = FALSE,
  save_rds          = TRUE,
  verbose           = TRUE,
  debug             = FALSE
) {
  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  bam       <- CYAN_check_bam(bam)
  reference <- CYAN_check_references(reference, omopolymeric_file,
                                     splicing_file, bed_file)

  sample <- .cyan_sample_name(bam, sample_name)

  output_dir <- .cyan_ensure_dir(output_dir, "output directory")
  work_dir   <- .cyan_ensure_dir(work_dir,   "work directory")

  # ------------------------------------------------------------------
  # 2. Define output file paths and check for an existing raw checkpoint
  # ------------------------------------------------------------------
  raw_out  <- file.path(output_dir, paste0(sample, "_reditools.txt"))
  subs_out <- file.path(output_dir, paste0(sample, "_reditools_subs.txt"))
  rds_out  <- file.path(output_dir, paste0(sample, "_reditools.rds"))

  if (file.exists(raw_out)) {
    finfo   <- file.info(raw_out)
    size_gb <- round(finfo$size / 1024^3, 2)
    message("[CYAN] Existing raw output found (made ",
            format(finfo$mtime, "%Y-%m-%d %H:%M"), ", ", size_gb,
            " GB), skipping REDItools2 run: ", raw_out,
            "\n  Delete the file to force a re-run.")
    result <- .cyan_parse_reditools_output(raw_out)
    result <- .cyan_attach_params(result, as.list(match.call()[-1]))
    .cyan_save_subs(result, subs_out, verbose)
    if (save_rds) {
      saveRDS(result, rds_out)
      if (verbose) message("[CYAN] RDS saved: ", rds_out)
    }
    return(invisible(result))
  }

  # ------------------------------------------------------------------
  # 2b. Duplicate-marking check (only when actually about to run --
  #     skipped above on a cache hit, since no filtering happens then)
  # ------------------------------------------------------------------
  if (require_dup_marked && !.cyan_check_duplicates_marked(bam, verbose = verbose)) {
    stop("[CYAN] No duplicate-flagged reads found in this BAM (SAM FLAG 0x400, ",
         "scanned up to 5,000,000 reads).\n",
         "  REDItools2's own duplicate filter relies entirely on this flag already ",
         "being set upstream -- on a BAM that was never duplicate-marked, that filter ",
         "is a silent no-op and PCR duplicates will inflate editing-frequency estimates.\n",
         "  Run 'samtools markdup' or Picard MarkDuplicates on this BAM first (flagging ",
         "is enough -- REDItools2 excludes flagged reads itself, no need to remove them).\n",
         "  If you are confident this BAM genuinely has no duplicates to flag, pass ",
         "require_dup_marked = FALSE to skip this check.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 3. Build CLI arguments
  # ------------------------------------------------------------------
  script_path <- system.file("python", "reditools.py", package = "CYAN")
  if (!nzchar(script_path))
    stop("[CYAN] Cannot locate inst/python/reditools.py inside the CYAN package.",
         call. = FALSE)

  cli_args <- .cyan_reditools_build_args(
    bam               = bam,
    reference         = reference,
    output_file       = raw_out,
    strand_mode       = strand_mode,
    min_coverage      = min_coverage,
    min_base_quality  = min_base_quality,
    min_read_quality  = min_read_quality,
    min_read_length   = min_read_length,
    min_edits         = min_edits,
    strict            = strict,
    region            = region,
    omopolymeric_file = omopolymeric_file,
    splicing_file     = splicing_file,
    bed_file          = bed_file,
    debug             = debug
  )

  # ------------------------------------------------------------------
  # 4. Run REDItools2 via basilisk + system2
  # ------------------------------------------------------------------
  if (verbose) message("[CYAN] Starting REDItools2 for sample: ", sample)

  .cyan_reditools_run(script_path, cli_args, verbose, debug)

  if (!file.exists(raw_out))
    stop("[CYAN] REDItools2 completed but output file was not created: ", raw_out,
         call. = FALSE)

  if (verbose) message("[CYAN] REDItools2 finished. Parsing output...")

  # ------------------------------------------------------------------
  # 5. Parse output and save
  # ------------------------------------------------------------------
  result <- .cyan_parse_reditools_output(raw_out)
  result <- .cyan_attach_params(result, as.list(match.call()[-1]))
  .cyan_save_subs(result, subs_out, verbose)

  if (save_rds) {
    saveRDS(result, rds_out)
    if (verbose) message("[CYAN] RDS saved: ", rds_out)
  }

  # ------------------------------------------------------------------
  # 6. Clean up work_dir intermediates
  # ------------------------------------------------------------------
  if (!keep_tmp) {
    tmp_files <- list.files(work_dir, full.names = TRUE)
    if (length(tmp_files) > 0) {
      file.remove(tmp_files)
      if (verbose) message("[CYAN] Temporary files removed from: ", work_dir)
    }
  }

  if (verbose) {
    message("[CYAN] Done. ", nrow(result), " substitution sites parsed into memory.")
    message("[CYAN] Full checkpoint (all covered positions, on disk only): ", raw_out)
    message("[CYAN] Subs output : ", subs_out)
    if (save_rds) message("[CYAN] RDS output  : ", rds_out)
  }

  invisible(result)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build the REDItools2 CLI argument vector
#' @keywords internal
.cyan_reditools_build_args <- function(
  bam, reference, output_file,
  strand_mode, min_coverage, min_base_quality, min_read_quality,
  min_read_length, min_edits, strict, region,
  omopolymeric_file, splicing_file, bed_file, debug
) {
  args <- c(
    "-f",   bam,
    "-r",   reference,
    "-o",   output_file,
    "-s",   as.integer(strand_mode),
    "-l",   as.integer(min_coverage),     # min column length = coverage
    "-bq",  as.integer(min_base_quality),
    "-q",   as.integer(min_read_quality),
    "-mrl", as.integer(min_read_length),
    "-me",  as.integer(min_edits)
  )

  if (isTRUE(strict))              args <- c(args, "-S")
  if (!is.null(region))            args <- c(args, "-g", region)
  if (!is.null(omopolymeric_file)) args <- c(args, "-m", omopolymeric_file)
  if (!is.null(splicing_file))     args <- c(args, "-sf", splicing_file)
  if (!is.null(bed_file))          args <- c(args, "-B", bed_file)
  if (isTRUE(debug))               args <- c(args, "-V")  # REDItools2 verbose → debug only

  as.character(args)
}

#' Check a BAM for duplicate-flagged reads (SAM FLAG 0x400)
#'
#' REDItools2's own duplicate filter (\code{read.is_duplicate}, hardcoded in
#' \code{reditools.py}) only excludes reads that already carry this flag --
#' it does not mark duplicates itself. On a BAM that was never run through a
#' duplicate-marking tool (e.g. \code{samtools markdup}, Picard
#' \code{MarkDuplicates}), the flag is never set and that filter is a silent
#' no-op, letting PCR-duplicate reads inflate editing-frequency estimates
#' uncaught. This scans for at least one flagged read as a cheap upstream
#' sanity check, exiting immediately once one is found.
#'
#' @param bam Character. Path to a sorted, indexed BAM file.
#' @param max_reads Integer. Maximum number of reads to scan before
#'   concluding none are flagged (default 5,000,000). Only reached for a
#'   genuinely duplicate-free BAM -- any real duplicate-marked BAM exits on
#'   the first flagged read.
#' @param verbose Logical. Print \code{[CYAN]} progress messages.
#'
#' @return Logical. TRUE if at least one duplicate-flagged read was found.
#' @keywords internal
.cyan_check_duplicates_marked <- function(bam, max_reads = 5e6, verbose = TRUE) {
  script_path <- system.file("python", "check_duplicates.py", package = "CYAN")
  if (!nzchar(script_path))
    stop("[CYAN] Cannot locate inst/python/check_duplicates.py inside the CYAN package.",
         call. = FALSE)

  if (verbose) message("[CYAN] Checking BAM for duplicate-marking (SAM FLAG 0x400)...")

  out <- basilisk::basiliskRun(
    env = .cyan_reditools_env,
    fun = function(script, bam_path, cap) {
      python_bin <- reticulate::py_exe()
      system2(
        command = python_bin,
        args    = c(shQuote(script), shQuote(bam_path), as.character(cap)),
        stdout  = TRUE,
        stderr  = TRUE
      )
    },
    script   = script_path,
    bam_path = bam,
    cap      = as.integer(max_reads)
  )

  identical(trimws(out[length(out)]), "TRUE")
}

#' Execute REDItools2 via basilisk + system2
#'
#' Activates the managed Python environment, retrieves the Python binary,
#' and calls the REDItools2 script as a subprocess.
#' Python stdout/stderr are suppressed unless \code{debug = TRUE}.
#' @keywords internal
.cyan_reditools_run <- function(script_path, cli_args, verbose, debug) {
  basilisk::basiliskRun(
    env = .cyan_reditools_env,
    fun = function(script, args, verbose, debug) {
      python_bin <- reticulate::py_exe()
      if (verbose)
        message("[CYAN] Python: ", python_bin)
      exit_code <- system2(
        command = python_bin,
        args    = c(shQuote(script), args),
        stdout  = if (debug) "" else FALSE,
        stderr  = if (debug) "" else FALSE
      )
      if (exit_code != 0)
        stop("[CYAN] REDItools2 exited with code ", exit_code,
             ". Re-run with debug = TRUE to see Python output.", call. = FALSE)
      invisible(exit_code)
    },
    script  = script_path,
    args    = cli_args,
    verbose = verbose,
    debug   = debug
  )
}

#' Parse REDItools2 tab-delimited output into a data.frame
#'
#' REDItools2 emits a header line beginning with "Region" and then one row
#' per position. The BaseCount[ACGT] column (e.g. \code{"[0, 5, 100, 3]"}) is
#' split into four integer columns.
#'
#' Reads the file in fixed-size line chunks rather than loading it whole --
#' a \code{strict=FALSE} genome-wide run writes one row per COVERED position
#' (not just substitution-carrying ones), which can reach hundreds of
#' millions of rows and tens of GB of raw text. The previous
#' \code{read.table()} + \code{do.call(rbind, strsplit(...))} approach does
#' not scale to that -- confirmed: crashed (OOM) on a 37GB file on a machine
#' with 58GB available RAM, 2026-10-14. Each chunk is parsed with
#' \code{data.table::fread()} and its \code{BaseCount[ACGT]} column is split
#' with \code{data.table::tstrsplit()} (vectorised -- avoids the
#' intermediate list-of-vectors that \code{strsplit()} + \code{do.call(rbind, .)}
#' builds, which is what actually drove the memory blowup, not \code{fread()}
#' vs \code{read.table()} alone).
#'
#' Each chunk is also filtered down to substitution rows
#' (\code{all_subs != "-"}) \strong{before} being added to the accumulator --
#' the object this function returns (and that \code{\link{CYAN_run_reditools}}
#' retains/saves as RDS) is always at the scale a \code{strict=TRUE} run would
#' have produced, regardless of the \code{strict} value the raw file was
#' actually generated with (a no-op filter on an already-strict file, since
#' every row already qualifies). This is what keeps per-sample memory bounded
#' by the number of editing sites rather than the number of covered positions,
#' and is why \code{\link{CYAN_run_reditools_batch}} can accumulate many
#' samples' results without a separate rewrite.
#'
#' Non-substitution rows are not retained or written anywhere by this
#' function. The full covered-position table still exists on disk as the raw
#' REDItools2 \code{.txt} checkpoint written by \code{CYAN_run_reditools}
#' before parsing ever runs -- recovering the non-strict background (e.g. for
#' APOBEC signature work) means a separate, targeted read of that raw file,
#' not a slice of an in-memory object from this function.
#'
#' Chunks accumulate as a list and are combined once with
#' \code{data.table::rbindlist()} at the end.
#'
#' @param path Character. Path to the raw REDItools2 output file.
#' @param chunk_size Integer. Number of lines to read and parse per chunk.
#'   Default 10 million -- large enough to keep the number of chunks (and
#'   their per-chunk \code{fread()} overhead) small, small enough to keep
#'   any one chunk's memory footprint modest. Lower this if parsing still
#'   runs out of memory on a very large file or a more memory-constrained
#'   machine.
#' @keywords internal
.cyan_parse_reditools_output <- function(path, chunk_size = 10e6) {
  # REDItools2 column names vary slightly across versions; normalise them
  col_map <- c(
    "Region"            = "region",
    "Position"          = "position",
    "Reference"         = "reference",
    "Strand"            = "strand",
    "Coverage-q30"      = "coverage",
    "Coverage-q20"      = "coverage",
    "Coverage"          = "coverage",
    "MeanQ"             = "mean_quality",
    "BaseCount[ACGT]"   = "base_count",
    "AllSubs"           = "all_subs",
    "Frequency"         = "frequency"
  )
  core_cols <- c("region", "position", "reference", "strand",
                 "coverage", "mean_quality", "base_count", "all_subs", "frequency")

  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  header_raw <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]

  chunks   <- list()
  n_chunks <- 0L
  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break
    n_chunks <- n_chunks + 1L

    raw <- data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t",
                             header = FALSE, col.names = header_raw,
                             colClasses = "character")

    for (old in names(col_map)) {
      idx <- match(old, colnames(raw))
      if (!is.na(idx)) colnames(raw)[idx] <- col_map[[old]]
    }
    present <- intersect(core_cols, colnames(raw))
    dt <- raw[, ..present]

    # Keep only substitution rows -- see function doc for rationale. Filtered
    # before the BaseCount parse below so that parse never runs on rows we
    # are about to discard.
    if ("all_subs" %in% colnames(dt)) {
      dt <- dt[all_subs != "-" & !is.na(all_subs)]
    }

    # Parse BaseCount[ACGT] string "[A, C, G, T]" -> four integer columns
    # (tstrsplit is vectorised -- this is the fix for the memory blowup, see
    # function doc above)
    if ("base_count" %in% colnames(dt)) {
      bc    <- gsub("[\\[\\] ]", "", dt$base_count)
      parts <- data.table::tstrsplit(bc, ",", fixed = TRUE)
      dt[, count_A := as.integer(parts[[1]])]
      dt[, count_C := as.integer(parts[[2]])]
      dt[, count_G := as.integer(parts[[3]])]
      dt[, count_T := as.integer(parts[[4]])]
      dt[, base_count := NULL]
    }

    if ("position" %in% colnames(dt))  dt[, position  := as.integer(position)]
    if ("coverage" %in% colnames(dt))  dt[, coverage  := as.integer(coverage)]
    if ("frequency" %in% colnames(dt)) dt[, frequency := as.numeric(frequency)]

    chunks[[n_chunks]] <- dt
    cat("    Chunk", n_chunks, "parsed (", length(lines), "lines ->",
        nrow(dt), "substitution rows)\n")
  }

  df <- data.table::rbindlist(chunks, use.names = TRUE)

  # Reorder columns for clarity
  col_order <- c("region", "position", "reference", "strand",
                 "coverage", "mean_quality",
                 "count_A", "count_C", "count_G", "count_T",
                 "all_subs", "frequency")
  data.table::setcolorder(df, intersect(col_order, colnames(df)))
  data.table::setDF(df)  # match the original plain-data.frame return type

  class(df) <- c("cyan_editing_result", "data.frame")
  df
}

#' Write the substitution-rows checkpoint file
#'
#' \code{result} is already substitution-rows-only by the time this is
#' called -- \code{\link{.cyan_parse_reditools_output}} filters to them
#' during parsing -- so this just persists that table to \code{path}.
#' @keywords internal
.cyan_save_subs <- function(result, path, verbose) {
  utils::write.table(result, path, sep = "\t", row.names = FALSE, quote = FALSE)
  if (verbose)
    message("[CYAN] Substitution sites (", nrow(result), " rows) saved: ", path)
  invisible(result)
}

#' Attach call parameters as an attribute for reproducibility
#' @keywords internal
.cyan_attach_params <- function(result, params) {
  attr(result, "params") <- params
  result
}

#' Print method for cyan_editing_result
#'
#' \code{x} is always substitution-rows-only (see
#' \code{\link{.cyan_parse_reditools_output}}), regardless of the
#' \code{strict} value the underlying raw file was generated with.
#' @export
#' @method print cyan_editing_result
print.cyan_editing_result <- function(x, ...) {
  cat("cyan_editing_result\n")
  cat("  Substitution sites:", nrow(x), "\n")
  if ("frequency" %in% colnames(x))
    cat("  Freq range:", round(min(x$frequency, na.rm = TRUE), 3), "–",
        round(max(x$frequency, na.rm = TRUE), 3), "\n")
  cat("\n")
  print(utils::head(as.data.frame(x), 6), ...)
  if (nrow(x) > 6) cat("  ... [", nrow(x) - 6, " more rows]\n", sep = "")
  invisible(x)
}

# ===========================================================================
# Batch wrapper
# ===========================================================================

#' Run REDItools2 on multiple samples
#'
#' Convenience wrapper around \code{\link{CYAN_run_reditools}} that iterates
#' over a vector of BAM files. If the vector is named, each name is forwarded
#' as the \code{sample_name} argument, controlling output file prefixes. If
#' unnamed, sample names are derived from the BAM filenames as usual.
#'
#' Samples are processed sequentially. A failed sample issues a warning and
#' is recorded as \code{NULL} in the returned list; remaining samples continue.
#' Because \code{CYAN_run_reditools} checkpoints its raw output, re-running
#' after a partial failure skips already-completed samples automatically.
#'
#' Each retained \code{cyan_editing_result} is substitution-rows-only (see
#' \code{\link{.cyan_parse_reditools_output}}), so holding every sample's
#' result in this function's return list stays bounded by the number of
#' editing sites per sample, not the number of covered positions -- even for
#' a genome-wide \code{strict=FALSE} run across many samples.
#'
#' @param bam        Named or unnamed character vector of BAM file paths.
#'   Names, if present, are used as \code{sample_name} for each sample.
#' @param reference  Character. Path to the reference FASTA (shared across all
#'   samples).
#' @param output_dir Character. Directory for all result files (created if
#'   absent). Each sample writes its own prefixed files here.
#' @param ...        Additional arguments passed to \code{CYAN_run_reditools}
#'   (e.g. \code{strand_mode}, \code{min_coverage}, \code{keep_tmp}).
#' @param verbose    Logical. Print \code{[CYAN]} progress messages including
#'   a per-sample summary on completion (default \code{TRUE}).
#'
#' @return An object of class \code{cyan_batch_result}: a named list with one
#'   element per BAM file. Each element is either a \code{cyan_editing_result}
#'   (success) or \code{NULL} (failure). Failed samples are identified in the
#'   printed summary.
#'
#' @seealso \code{\link{CYAN_run_reditools}}
#'
#' @export
CYAN_run_reditools_batch <- function(bam, reference, output_dir, ...,
                                     verbose = TRUE) {
  sample_names <- names(bam)

  # Result list named by sample label (for easy downstream access)
  labels <- if (!is.null(sample_names)) sample_names else
    tools::file_path_sans_ext(basename(bam))
  results <- vector("list", length(bam))
  names(results) <- labels

  for (i in seq_along(bam)) {
    label <- labels[i]
    sname <- if (!is.null(sample_names)) sample_names[i] else NULL

    if (verbose)
      message("[CYAN] Batch [", i, "/", length(bam), "] Sample: ", label)

    tryCatch(
      results[[i]] <- CYAN_run_reditools(
        bam         = bam[i],
        reference   = reference,
        output_dir  = output_dir,
        sample_name = sname,
        verbose     = verbose,
        ...
      ),
      error = function(e) {
        warning("[CYAN] Sample '", label, "' failed: ", conditionMessage(e),
                call. = FALSE)
      }
    )
  }

  succeeded <- !vapply(results, is.null, logical(1))

  if (verbose) {
    message("\n[CYAN] Batch complete: ", sum(succeeded), "/", length(bam),
            " samples succeeded.")
    if (any(succeeded)) {
      n_sites <- vapply(results[succeeded], nrow, integer(1))
      for (nm in names(n_sites))
        message("[CYAN]   ", nm, ": ", n_sites[nm], " substitution sites")
    }
    if (any(!succeeded))
      message("[CYAN] Failed samples: ",
              paste(labels[!succeeded], collapse = ", "))
  }

  class(results) <- c("cyan_batch_result", "list")
  invisible(results)
}

#' Print method for cyan_batch_result
#' @export
#' @method print cyan_batch_result
print.cyan_batch_result <- function(x, ...) {
  succeeded <- !vapply(x, is.null, logical(1))
  cat("cyan_batch_result\n")
  cat("  Samples   :", length(x), "\n")
  cat("  Succeeded :", sum(succeeded), "\n")
  if (any(!succeeded))
    cat("  Failed    :", paste(names(x)[!succeeded], collapse = ", "), "\n")
  cat("\n")
  for (nm in names(x)[succeeded]) {
    r <- x[[nm]]
    cat("  ", nm, ": ", nrow(r), " substitution sites\n", sep = "")
  }
  invisible(x)
}
