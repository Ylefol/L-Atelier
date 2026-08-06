# ==============================================================================
# HORIZON - TEcount Transposable Element Quantification
# ==============================================================================
# TEcount (bundled with the bioconda `tetranscripts` package, alongside the
# 2-group `TEtranscripts` DE binary which is NOT used here) quantifies genes
# and TE subfamilies for a single BAM at a time -- architecturally like
# `salmon quant`, one sample in, one count table out, with no 2-group
# treatment/control framing required. This is what makes it usable for
# designs with more than two conditions: run it once per sample regardless of
# condition, aggregate into a matrix, and hand the TE half to GAIA's own
# ARTEMIS_normalize_counts()/ARTEMIS_differential_counts() for statistics.
#
# TEcount requires a BAM with multi-mapping alignments explicitly reported
# (see HORIZON_run_star_align() in star_align.R) -- it redistributes those
# reads across candidate TE loci via an EM algorithm rather than discarding
# them, which is the entire point of using it instead of featureCounts on a
# RepeatMasker GTF directly.
#
# Flag names below (-b, --GTF, --TE, --mode, --stranded, --sortByPos,
# --project, --outdir) confirmed against a real `TEcount --help` (2026-07-28).
# One correction from the original published-docs-only draft: --stranded
# takes "no"/"forward"/"reverse", not "no"/"yes"/"reverse" -- fixed in
# .strand_to_tecount() (internals.R). --mode defaults to "multi" in TEcount
# itself, matching this wrapper's own default.
#
# STAR and TEcount share the standard horizon_cli conda environment (see
# HORIZON_set_conda_env() / basilisk.R) -- both are called via
# .horizon_run_cli(), the same mechanism used for samtools/bedtools/macs3/
# deeptools/salmon. There is no dependency conflict requiring a dedicated
# environment (unlike split-pipe, which does need one -- see parse.R).
# ==============================================================================


#' Quantify genes + TE subfamilies for a single BAM with TEcount
#'
#' Wraps the \code{TEcount} CLI (part of the bioconda \code{tetranscripts}
#' package) to quantify both genes and transposable element (TE) subfamilies
#' from a single BAM, using an EM algorithm to redistribute multi-mapping
#' reads across candidate TE loci rather than discarding them.
#'
#' \code{TEcount} must be available in the conda environment registered with
#' \code{\link{HORIZON_set_conda_env}}. Install it alongside \code{STAR} (see
#' \code{\link{HORIZON_build_star_index}}) and the existing \code{horizon_cli}
#' tools:
#' \preformatted{
#'   conda install -n horizon_cli -c bioconda -c conda-forge star tetranscripts
#' }
#'
#' \strong{Only the TE rows of \code{TEcount}'s output are meant to be used
#' downstream} (see \code{\link{HORIZON_aggregate_tecount}}) -- its gene-count
#' half comes from a different aligner/counting method than HORIZON's
#' canonical gene-level pipeline (\code{\link{HORIZON_run_align}} +
#' \code{\link{HORIZON_run_count}}) and would be inconsistent with
#' already-computed differential expression results if mixed in. Canonical
#' gene counts should keep coming from the existing Rsubread/featureCounts
#' path; merge the two feature sets with
#' \code{\link{HORIZON_merge_gene_te_counts}}.
#'
#' \strong{Two usage modes:}
#' \enumerate{
#'   \item \strong{Sample-sheet mode} (default): provide \code{sample_sheet}
#'     and \code{sample_id}. \code{strandedness} is read from the sample
#'     sheet row, and the BAM is located at the standard HORIZON path
#'     \code{<output_dir>/<sample_id>/aligned_star/<sample_id>_sorted.bam}
#'     (i.e. \code{\link{HORIZON_run_star_align}}'s output -- \strong{not}
#'     the Rsubread BAM). Output is written to
#'     \code{<output_dir>/<sample_id>/tecount/}.
#'   \item \strong{Direct BAM mode}: provide \code{bam_file}, plus
#'     \code{strandedness} explicitly, since it cannot be reliably inferred
#'     without a sample sheet row. \code{sample_id} is inferred from the
#'     filename if omitted. Output goes to \code{output_dir} (defaults to the
#'     directory containing the BAM).
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}. Required in sample-sheet
#'   mode; must be \code{NULL} when \code{bam_file} is supplied.
#' @param sample_id Character. Sample ID to process (sample-sheet mode).
#'   Optional in direct BAM mode -- inferred from the filename if omitted.
#' @param bam_file Character or \code{NULL}. Path to a STAR-aligned,
#'   coordinate-sorted BAM with multi-mapping alignments reported (see
#'   \code{\link{HORIZON_run_star_align}}). When supplied, \code{sample_sheet}
#'   must be \code{NULL}. Default \code{NULL}.
#' @param output_dir Character or \code{NULL}. Output directory for direct BAM
#'   mode. Ignored in sample-sheet mode. Defaults to the directory containing
#'   \code{bam_file}.
#' @param strandedness Character or \code{NULL}. One of \code{"unstranded"},
#'   \code{"forward"}, \code{"reverse"}. Required in direct BAM mode; read
#'   from the sample sheet in sample-sheet mode. Mapped internally to
#'   \code{TEcount}'s own \code{--stranded} vocabulary
#'   (\code{"no"}/\code{"forward"}/\code{"reverse"}).
#' @param gene_gtf Character. Path to the gene GTF annotation -- should be the
#'   same GTF already used for the project's canonical gene-level counting,
#'   for consistency (even though its gene-count output is discarded here).
#' @param te_gtf Character. Path to the RepeatMasker-derived TE GTF
#'   (\code{gene_id}/\code{family_id}/\code{class_id} composite format).
#' @param mode Character. \code{TEcount}'s \code{--mode}: \code{"multi"}
#'   (default) applies EM-based redistribution of multi-mapping reads across
#'   TE loci -- the actual point of using TEcount over featureCounts.
#'   \code{"uniq"} counts unique-mapping reads only, useful as a robustness
#'   comparison against the \code{"multi"} result but not the primary mode.
#' @param force Logical. If \code{FALSE} (default), skip quantification when
#'   the output \code{.cntTable} already exists. Set \code{TRUE} to rerun.
#' @param verbose Logical. If \code{TRUE} (default), prints the full
#'   \code{TEcount} command before executing it.
#' @param ... Additional \code{TEcount} flags passed verbatim.
#'
#' @return Character. Path to the \code{.cntTable} file, invisibly.
#' @export
HORIZON_run_tecount <- function(sample_sheet  = NULL,
                                 sample_id     = NULL,
                                 bam_file      = NULL,
                                 output_dir    = NULL,
                                 strandedness  = NULL,
                                 gene_gtf,
                                 te_gtf,
                                 mode          = "multi",
                                 force         = FALSE,
                                 verbose       = TRUE,
                                 ...) {

  if (!file.exists(gene_gtf)) stop("Gene GTF not found: ", gene_gtf, call. = FALSE)
  if (!file.exists(te_gtf))   stop("TE GTF not found: ", te_gtf, call. = FALSE)

  valid_modes <- c("uniq", "multi")
  if (!mode %in% valid_modes)
    stop("mode must be one of: ", paste(valid_modes, collapse = ", "), call. = FALSE)

  if (!is.null(bam_file) && !is.null(sample_sheet))
    stop("Provide either 'bam_file' or 'sample_sheet', not both.", call. = FALSE)
  if (is.null(bam_file) && is.null(sample_sheet))
    stop("One of 'bam_file' or 'sample_sheet' must be provided.", call. = FALSE)

  if (!is.null(bam_file)) {
    # Direct BAM mode
    if (!file.exists(bam_file))
      stop("BAM file not found: ", bam_file, call. = FALSE)
    if (is.null(sample_id))
      sample_id <- tools::file_path_sans_ext(
                     tools::file_path_sans_ext(basename(bam_file)))
    if (is.null(strandedness))
      stop("'strandedness' is required in direct BAM mode.", call. = FALSE)
    bam     <- bam_file
    out_dir <- if (!is.null(output_dir)) output_dir else dirname(bam_file)
  } else {
    # Sample-sheet mode
    if (is.null(sample_id))
      stop("'sample_id' is required when using sample-sheet mode.", call. = FALSE)
    row          <- .get_sample_row(sample_sheet, sample_id)
    strandedness <- row$strandedness
    out_dir      <- file.path(row$output_dir, sample_id, "tecount")
    bam          <- file.path(row$output_dir, sample_id, "aligned_star",
                              paste0(sample_id, "_sorted.bam"))
    if (!file.exists(bam)) {
      stop("STAR-aligned BAM not found: ", bam,
           "\nRun HORIZON_run_star_align() first.", call. = FALSE)
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  strand_str <- .strand_to_tecount(strandedness)
  cnt_file   <- file.path(out_dir, paste0(sample_id, ".cntTable"))

  if (!isTRUE(force) && file.exists(cnt_file)) {
    cat("[TECOUNT] Output already exists for:", sample_id,
        "— skipping (use force=TRUE to rerun)\n")
    return(invisible(cnt_file))
  }

  args <- c(
    "-b",         bam,
    "--GTF",      gene_gtf,
    "--TE",       te_gtf,
    "--mode",     mode,
    "--stranded", strand_str,
    "--sortByPos",
    "--project",  sample_id,
    "--outdir",   out_dir,
    c(...)
  )

  cat("[TECOUNT] Quantifying genes + TEs for:", sample_id, "\n")
  cat("    BAM:          ", bam, "\n")
  cat("    Mode:         ", mode, "\n")
  cat("    Strandedness: ", strandedness, " (", strand_str, ")\n", sep = "")
  cat("    Output:       ", out_dir, "\n\n")

  if (isTRUE(verbose))
    cat("[TECOUNT] Command: TEcount", paste(args, collapse = " "), "\n\n")

  ret <- .horizon_run_cli("TEcount", args)

  if (ret != 0L)
    stop("TEcount failed for: ", sample_id, " (exit code ", ret, ")", call. = FALSE)

  if (!file.exists(cnt_file))
    stop("TEcount ran but expected output not found: ", cnt_file, call. = FALSE)

  cat("[TECOUNT] Quantification complete for:", sample_id, "\n")
  cat("    Output:", cnt_file, "\n")
  invisible(cnt_file)
}


# Internal: parse a RepeatMasker-derived TE GTF and return the set of unique
# composite TE IDs ("gene_id:family_id:class_id"), used to distinguish TE
# rows from gene rows in TEcount's combined cntTable output.
.tecount_te_ids <- function(te_gtf) {
  for (pkg in c("rtracklayer", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required to parse the TE GTF. ",
           "Install with: BiocManager::install('", pkg, "')", call. = FALSE)
  }

  gr   <- rtracklayer::import(te_gtf, format = "gtf")
  meta <- S4Vectors::mcols(gr)

  missing_cols <- setdiff(c("gene_id", "family_id", "class_id"), colnames(meta))
  if (length(missing_cols) > 0)
    stop("TE GTF is missing expected attribute(s): ",
         paste(missing_cols, collapse = ", "),
         ". Expected a RepeatMasker-derived GTF with gene_id/family_id/class_id.",
         call. = FALSE)

  unique(paste(meta$gene_id, meta$family_id, meta$class_id, sep = ":"))
}


#' Aggregate per-sample TEcount output into a TE count matrix
#'
#' Reads the per-sample \code{.cntTable} files produced by
#' \code{\link{HORIZON_run_tecount}} and combines the \strong{TE rows only}
#' (gene rows are dropped -- see \code{\link{HORIZON_run_tecount}} for why)
#' into a single TE-subfamily x sample count matrix.
#'
#' TE rows are identified by cross-referencing each \code{.cntTable}'s row
#' IDs against the ground-truth set of composite TE IDs parsed directly from
#' \code{te_gtf}, rather than assuming a fixed ID-format/delimiter convention.
#'
#' Output is written to \code{<output_root>/aggregated/tecount/}:
#' \itemize{
#'   \item \code{te_count_matrix.csv} / \code{.rds} -- TE subfamily x sample
#'     count matrix
#'   \item \code{sample_metadata.csv} -- sample sheet columns carried through
#'     (if \code{save_metadata = TRUE})
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param te_gtf Character. Path to the same TE GTF passed to
#'   \code{\link{HORIZON_run_tecount}} for every sample being aggregated.
#' @param output_root Character or NULL. Root output directory. If NULL,
#'   inferred from the first row's \code{output_dir}.
#' @param save_rds Logical. Also save the matrix as an RDS file. Default TRUE.
#' @param save_metadata Logical. Save the sample sheet as
#'   \code{sample_metadata.csv} alongside the matrix. Default TRUE.
#'
#' @return The aggregated TE count matrix (TE subfamily x sample), invisibly.
#' @export
HORIZON_aggregate_tecount <- function(sample_sheet,
                                       te_gtf,
                                       output_root   = NULL,
                                       save_rds      = TRUE,
                                       save_metadata = TRUE) {
  if (is.null(output_root)) {
    output_root <- sample_sheet$output_dir[1]
    cat("output_root inferred from first sample: ", output_root, "\n")
  }

  cat("[TECOUNT] Parsing TE GTF for ground-truth TE ID set...\n")
  te_ids <- .tecount_te_ids(te_gtf)
  cat("    TE subfamilies:", length(te_ids), "\n")

  agg_dir <- file.path(output_root, "aggregated", "tecount")
  dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)

  count_list <- list()

  for (sid in sample_sheet$sample_id) {
    sample_out_dir <- sample_sheet$output_dir[sample_sheet$sample_id == sid]
    cnt_file <- file.path(sample_out_dir, sid, "tecount", paste0(sid, ".cntTable"))

    if (!file.exists(cnt_file)) {
      warning("cntTable not found for sample '", sid, "' — skipped.")
      next
    }

    df <- read.table(cnt_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    ids    <- df[[1]]
    counts <- df[[2]]

    is_te <- ids %in% te_ids
    count_list[[sid]] <- setNames(counts[is_te], ids[is_te])
  }

  if (length(count_list) == 0) {
    stop("No cntTable files found. Run HORIZON_run_tecount() for at least one sample first.")
  }

  all_te       <- Reduce(union, lapply(count_list, names))
  count_matrix <- do.call(cbind, lapply(count_list, function(x) x[all_te]))
  rownames(count_matrix) <- all_te

  n_te      <- nrow(count_matrix)
  n_samples <- ncol(count_matrix)
  cat("[TECOUNT] Aggregated TE count matrix: ", n_te, " TE subfamilies x ",
      n_samples, " samples\n", sep = "")

  out_csv <- file.path(agg_dir, "te_count_matrix.csv")
  write.csv(count_matrix, out_csv)
  cat("    TE count matrix written to:", out_csv, "\n")

  if (save_rds) {
    saveRDS(count_matrix, file.path(agg_dir, "te_count_matrix.rds"))
  }

  if (save_metadata) {
    meta <- sample_sheet[sample_sheet$sample_id %in% colnames(count_matrix), ]
    write.csv(meta, file.path(agg_dir, "sample_metadata.csv"), row.names = FALSE)
    cat("    Sample metadata written to:", file.path(agg_dir, "sample_metadata.csv"), "\n")
  }

  invisible(count_matrix)
}


#' Merge gene-level and TE-level count matrices
#'
#' Combines a canonical gene-level count matrix (e.g. from
#' \code{\link{HORIZON_aggregate_counts}}) with a TE-subfamily count matrix
#' (from \code{\link{HORIZON_aggregate_tecount}}) into a single feature x
#' sample matrix suitable for \code{ARTEMIS_normalize_counts()}. There is no
#' prior HORIZON/GAIA convention for combining two distinct feature sets
#' (existing aggregators only ever combine samples for a single feature set),
#' so this fills that gap.
#'
#' @param gene_counts Matrix. Gene x sample count matrix, with gene IDs as
#'   rownames and sample IDs as colnames.
#' @param te_counts Matrix. TE subfamily x sample count matrix, with composite
#'   TE IDs as rownames and sample IDs as colnames.
#' @param verbose Logical. Print merge diagnostics (dimensions, any dropped
#'   samples). Default TRUE.
#'
#' @return The combined matrix (genes + TEs) x sample, restricted to samples
#'   present in both inputs.
#' @export
HORIZON_merge_gene_te_counts <- function(gene_counts, te_counts, verbose = TRUE) {
  if (is.null(rownames(gene_counts)) || is.null(rownames(te_counts)))
    stop("Both gene_counts and te_counts must have rownames (feature IDs).", call. = FALSE)
  if (is.null(colnames(gene_counts)) || is.null(colnames(te_counts)))
    stop("Both gene_counts and te_counts must have colnames (sample IDs).", call. = FALSE)

  collisions <- intersect(rownames(gene_counts), rownames(te_counts))
  if (length(collisions) > 0)
    stop("Rowname collision(s) between gene_counts and te_counts: ",
         paste(utils::head(collisions, 10), collapse = ", "),
         if (length(collisions) > 10) ", ..." else "",
         call. = FALSE)

  common_samples <- intersect(colnames(gene_counts), colnames(te_counts))
  if (length(common_samples) == 0)
    stop("No sample IDs in common between gene_counts and te_counts.", call. = FALSE)

  dropped <- union(setdiff(colnames(gene_counts), common_samples),
                    setdiff(colnames(te_counts), common_samples))
  if (length(dropped) > 0)
    warning("Sample(s) present in only one of gene_counts/te_counts, dropped: ",
            paste(dropped, collapse = ", "))

  combined <- rbind(gene_counts[, common_samples, drop = FALSE],
                     te_counts[, common_samples, drop = FALSE])

  if (isTRUE(verbose)) {
    cat("[TECOUNT] Merged gene + TE counts: ", nrow(combined), " features x ",
        length(common_samples), " samples\n", sep = "")
    cat("    Genes: ", nrow(gene_counts), "  TEs: ", nrow(te_counts), "\n", sep = "")
    if (length(dropped) > 0)
      cat("    Dropped sample(s): ", paste(dropped, collapse = ", "), "\n", sep = "")
  }

  combined
}
