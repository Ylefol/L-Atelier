# ==============================================================================
# HORIZON - Salmon Transcript Quantification
# ==============================================================================
# Decoy-aware transcriptome index build + per-sample selective-alignment
# quantification, for isoform-level (transcript) abundance estimation.
# Unlike HORIZON_run_count() (gene/exon-union counting via featureCounts on a
# genome-aligned BAM), Salmon maps directly against a transcriptome index and
# resolves reads shared between overlapping isoforms via its EM algorithm --
# the standard approach for transcript-level quantification.
# ==============================================================================


#' Build a decoy-aware Salmon transcriptome index
#'
#' @description Extracts transcript sequences from a genome FASTA + GTF
#' annotation (via \code{\link[GenomicFeatures]{extractTranscriptSeqs}}),
#' builds a decoy-aware Salmon index following Salmon's documented best
#' practice (the whole genome is included as "decoy" sequence, which prevents
#' reads originating from unspliced/intergenic genomic sequence from being
#' spuriously assigned to a similar transcript), and writes the index to
#' \code{index_dir}.
#'
#' This is a one-time operation per reference genome/annotation pair. The
#' \code{salmon} binary must be available in the conda environment registered
#' with \code{\link{HORIZON_set_conda_env}}.
#'
#' @param gtf Character. Path to the GTF annotation used to define transcripts.
#' @param genome_fasta Character. Path to the reference genome FASTA. Indexed
#'   in place (\code{.fai}) if not already indexed.
#' @param index_dir Character. Output directory for the Salmon index.
#' @param kmer_length Integer. Salmon \code{-k} (minimum acceptable match
#'   length for quasi-mapping). Default 31, appropriate for reads >= 75bp;
#'   reduce (e.g. to 23) for shorter reads. See the Salmon documentation for
#'   guidance.
#' @param keep_duplicates Logical. Passes \code{--keepDuplicates} to
#'   \code{salmon index} when \code{TRUE}. Without it, Salmon
#'   collapses transcripts with identical sequence into a single representative
#'   ID at index time, silently dropping the others from quantification --
#'   this matters for isoform-level analyses (e.g. downstream isoform-switch
#'   testing) where every transcript ID needs its own entry, and is generally
#'   recommended by Salmon/tximport's own tutorials regardless of downstream
#'   use. Set \code{FALSE} (default) to match Salmon's own out-of-the-box default.
#' @param threads Integer. Threads for \code{salmon index}. Default 4.
#' @param force Logical. Rebuild even if an index already exists at
#'   \code{index_dir}. Default \code{FALSE}.
#'
#' @return Character. \code{index_dir}, invisibly.
#'
#' @details
#' Intermediate files (the extracted transcript FASTA, the genome+transcript
#' "gentrome" FASTA used as Salmon's direct input, and the decoy sequence name
#' list) are written to \code{<index_dir>_prep/} alongside the index, for
#' provenance/debugging. They are not required after the index is built and
#' may be deleted to reclaim disk space.
#'
#' Requires: \code{GenomicFeatures}, \code{Biostrings}, \code{Rsamtools}.
#'
#' @export
HORIZON_build_salmon_index <- function(gtf,
                                        genome_fasta,
                                        index_dir,
                                        kmer_length     = 31L,
                                        keep_duplicates = FALSE,
                                        threads         = 4L,
                                        force           = FALSE) {

  if (!file.exists(gtf)) stop("GTF file not found: ", gtf, call. = FALSE)
  if (!file.exists(genome_fasta))
    stop("Genome FASTA not found: ", genome_fasta, call. = FALSE)

  for (pkg in c("GenomicFeatures", "Biostrings", "Rsamtools")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. ",
           "Install with: BiocManager::install('", pkg, "')", call. = FALSE)
  }

  info_file <- file.path(index_dir, "info.json")
  if (!isTRUE(force) && file.exists(info_file)) {
    cat("Salmon index already exists at: ", index_dir,
            " — skipping (use force=TRUE to rebuild)")
    return(invisible(index_dir))
  }

  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  prep_dir <- paste0(index_dir, "_prep")
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

  cat("[SALMON] Building TxDb from GTF (this may take a while)...")
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf, format = "gtf",
                                            dataSource = basename(gtf))

  fa_index <- paste0(genome_fasta, ".fai")
  if (!file.exists(fa_index)) {
    cat("[SALMON] Indexing genome FASTA...")
    Rsamtools::indexFa(genome_fasta)
  }
  genome <- Rsamtools::FaFile(genome_fasta)

  cat("[SALMON] Extracting transcript sequences...")
  exons_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  tx_seqs     <- GenomicFeatures::extractTranscriptSeqs(genome, exons_by_tx)
  cat("    Transcripts: ", length(tx_seqs))

  transcript_fasta <- file.path(prep_dir, "transcripts.fa")
  Biostrings::writeXStringSet(tx_seqs, filepath = transcript_fasta)

  cat("[SALMON] Building decoy-aware gentrome (transcripts + genome)...")
  decoys_path  <- file.path(prep_dir, "decoys.txt")
  genome_names <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(genome_fasta)))
  writeLines(genome_names, decoys_path)

  gentrome_fasta <- file.path(prep_dir, "gentrome.fa")
  file.copy(transcript_fasta, gentrome_fasta, overwrite = TRUE)
  file.append(gentrome_fasta, genome_fasta)

  cat("[SALMON] Building index (this may take a while for large genomes)...")
  args <- c(
    "index",
    "-t", gentrome_fasta,
    "-d", decoys_path,
    "-i", index_dir,
    "-k", as.integer(kmer_length),
    "-p", as.integer(threads)
  )
  if (isTRUE(keep_duplicates)) args <- c(args, "--keepDuplicates")
  exit_code <- .horizon_run_cli("salmon", args)
  if (exit_code != 0)
    stop("salmon index failed (exit code ", exit_code, ").", call. = FALSE)

  cat("Salmon index complete: ", index_dir)
  invisible(index_dir)
}


#' Quantify transcript abundance with Salmon
#'
#' @description Runs \code{salmon quant} in selective-alignment (mapping-based)
#' mode directly against FASTQ reads and a pre-built decoy-aware index (see
#' \code{\link{HORIZON_build_salmon_index}}). Unlike \code{\link{HORIZON_run_count}},
#' this does not use the genome-aligned BAM produced by
#' \code{\link{HORIZON_run_align}} at all -- Salmon performs its own
#' lightweight mapping against the transcriptome, which is what allows it to
#' resolve reads shared between overlapping isoforms via its EM algorithm.
#'
#' \strong{Two usage modes:}
#' \enumerate{
#'   \item \strong{Sample-sheet mode} (default): provide \code{sample_sheet}
#'     and \code{sample_id}. Reads are taken from the trimmed FASTQs produced
#'     by \code{\link{HORIZON_run_qc_trim}}
#'     (\code{<output_dir>/<sample_id>/qc/<sample_id>_trimmed_R1/R2.fastq.gz}).
#'     Output is written to \code{<output_dir>/<sample_id>/salmon/}.
#'   \item \strong{Direct FASTQ mode}: provide \code{fastq_r1} (and
#'     \code{fastq_r2} for paired-end data) directly. \code{sample_id} is
#'     inferred from the R1 filename if omitted. Output goes to
#'     \code{output_dir} (defaults to the directory containing \code{fastq_r1}).
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame. Required in
#'   sample-sheet mode; must be \code{NULL} when \code{fastq_r1} is supplied.
#' @param sample_id Character. Sample ID (sample-sheet mode). Optional in
#'   direct FASTQ mode -- inferred from the first element of \code{fastq_r1}
#'   if omitted.
#' @param fastq_r1 Character vector or \code{NULL}. Path(s) to (trimmed) R1
#'   FASTQ. Multiple files (e.g. multi-lane runs) are passed straight through
#'   to \code{salmon quant -1}, which treats them as one combined library --
#'   no pre-concatenation needed. When supplied, \code{sample_sheet} must be
#'   \code{NULL}. Default \code{NULL}.
#' @param fastq_r2 Character vector or \code{NULL}. Path(s) to R2 FASTQ for
#'   paired-end data (direct FASTQ mode only), positionally matched to
#'   \code{fastq_r1} (same length, same lane order). \code{NULL} for
#'   single-end.
#' @param output_dir Character or \code{NULL}. Output directory for direct
#'   FASTQ mode. Ignored in sample-sheet mode. Defaults to the directory
#'   containing \code{fastq_r1}.
#' @param index Character. Path to a Salmon index directory (from
#'   \code{\link{HORIZON_build_salmon_index}}).
#' @param lib_type Character. Salmon library type (\code{-l}). Default
#'   \code{"A"} (automatic detection -- recommended unless you have a specific
#'   reason to pin it; see the Salmon documentation for the ISR/ISF/IU codes
#'   corresponding to the sample sheet's \code{strandedness} values).
#' @param threads Integer. Threads for \code{salmon quant}. Default 4.
#' @param extra_flags Character vector of additional \code{salmon quant} flags
#'   appended verbatim (e.g. \code{c("--seqBias", "--gcBias")}). Default
#'   \code{character(0)} -- \code{--validateMappings} used to be the default
#'   here, but current Salmon versions warn that it has no effect (selective
#'   alignment, what \code{--validateMappings} used to enable, is now always
#'   the default mapping mode; \code{--sketch} is the flag to opt back into
#'   pseudoalignment).
#' @param force Logical. If \code{FALSE} (default), skip quantification when
#'   \code{quant.sf} already exists. Set \code{TRUE} to rerun.
#'
#' @return Character. Path to the \code{quant.sf} file, invisibly.
#' @export
HORIZON_run_salmon <- function(sample_sheet = NULL,
                                sample_id    = NULL,
                                fastq_r1     = NULL,
                                fastq_r2     = NULL,
                                output_dir   = NULL,
                                index,
                                lib_type     = "A",
                                threads      = 4L,
                                extra_flags  = character(0),
                                force        = FALSE) {

  if (!dir.exists(index))
    stop("Salmon index directory not found: ", index,
         "\nRun HORIZON_build_salmon_index() first.", call. = FALSE)

  if (!is.null(fastq_r1) && !is.null(sample_sheet))
    stop("Provide either 'fastq_r1' or 'sample_sheet', not both.", call. = FALSE)
  if (is.null(fastq_r1) && is.null(sample_sheet))
    stop("One of 'fastq_r1' or 'sample_sheet' must be provided.", call. = FALSE)

  if (!is.null(fastq_r1)) {
    # Direct FASTQ mode. fastq_r1/fastq_r2 may each be a vector of multiple
    # files (e.g. multi-lane runs) -- salmon accepts "-1 f1 f2 ... -2 g1 g2 ..."
    # directly and treats them as one combined library, so lanes are passed
    # through as-is rather than pre-concatenated.
    if (!all(file.exists(fastq_r1)))
      stop("FASTQ R1 not found: ", paste(fastq_r1[!file.exists(fastq_r1)], collapse = ", "),
           call. = FALSE)
    if (!is.null(fastq_r2) && !all(file.exists(fastq_r2)))
      stop("FASTQ R2 not found: ", paste(fastq_r2[!file.exists(fastq_r2)], collapse = ", "),
           call. = FALSE)
    if (!is.null(fastq_r2) && length(fastq_r1) != length(fastq_r2))
      stop("fastq_r1 and fastq_r2 must have the same length (", length(fastq_r1),
           " vs ", length(fastq_r2), ").", call. = FALSE)
    if (is.null(sample_id))
      sample_id <- sub("(_R1|_1)?\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_r1[1]))
    r1      <- fastq_r1
    r2      <- fastq_r2
    out_dir <- if (!is.null(output_dir)) output_dir else dirname(fastq_r1[1])
  } else {
    # Sample-sheet mode
    if (is.null(sample_id))
      stop("'sample_id' is required when using sample-sheet mode.", call. = FALSE)
    row    <- .get_sample_row(sample_sheet, sample_id)
    paired <- as.logical(row$paired_end)
    qc_dir <- file.path(row$output_dir, sample_id, "qc")
    r1 <- file.path(qc_dir, paste0(sample_id, "_trimmed_R1.fastq.gz"))
    r2 <- if (paired) file.path(qc_dir, paste0(sample_id, "_trimmed_R2.fastq.gz")) else NULL
    if (!file.exists(r1))
      stop("Trimmed FASTQ not found: ", r1,
           "\nRun HORIZON_run_qc_trim() first.", call. = FALSE)
    if (!is.null(r2) && !file.exists(r2))
      stop("Trimmed FASTQ not found: ", r2,
           "\nRun HORIZON_run_qc_trim() first.", call. = FALSE)
    out_dir <- file.path(row$output_dir, sample_id, "salmon")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  quant_file <- file.path(out_dir, "quant.sf")
  if (!isTRUE(force) && file.exists(quant_file)) {
    cat("Salmon quant already exists for: ", sample_id,
            " — skipping (use force=TRUE to rerun)")
    return(invisible(quant_file))
  }

  args <- c("quant", "-i", index, "-l", lib_type, "-p", as.integer(threads))
  if (!is.null(r2)) {
    args <- c(args, "-1", r1, "-2", r2)
  } else {
    args <- c(args, "-r", r1)
  }
  args <- c(args, "-o", out_dir, extra_flags)

  cat("[", sample_id, "] Running salmon quant")
  exit_code <- .horizon_run_cli("salmon", args)
  if (exit_code != 0)
    stop("salmon quant failed for '", sample_id, "' (exit code ", exit_code, ").",
         call. = FALSE)

  if (!file.exists(quant_file))
    stop("salmon quant ran but output file not found: ", quant_file, call. = FALSE)

  cat("Salmon quant complete for: ", sample_id, "\n  Output: ", quant_file)
  invisible(quant_file)
}


#' Aggregate per-sample Salmon quant.sf files into TPM / count matrices
#'
#' @description Reads the per-sample \code{quant.sf} files produced by
#' \code{\link{HORIZON_run_salmon}} and combines them into transcript x
#' sample matrices of TPM and estimated read counts (\code{NumReads}).
#' Samples missing a \code{quant.sf} file are skipped with a warning.
#'
#' Output is written to \code{<output_root>/aggregated/salmon/}:
#' \itemize{
#'   \item \code{tpm_matrix.csv} / \code{.rds} — transcript x sample TPM
#'   \item \code{counts_matrix.csv} / \code{.rds} — transcript x sample
#'     estimated counts (\code{NumReads})
#'   \item \code{sample_metadata.csv} — sample sheet columns carried through
#'     (if \code{save_metadata = TRUE})
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param output_root Character or NULL. Root output directory. If NULL,
#'   inferred from the first row's \code{output_dir}.
#' @param save_rds Logical. Also save matrices as RDS. Default TRUE.
#' @param save_metadata Logical. Save sample sheet as
#'   \code{sample_metadata.csv} alongside the matrices. Default TRUE.
#'
#' @return A list with elements \code{tpm} and \code{counts} (transcript x
#'   sample matrices), invisibly.
#' @export
HORIZON_aggregate_salmon <- function(sample_sheet,
                                      output_root   = NULL,
                                      save_rds      = TRUE,
                                      save_metadata = TRUE) {
  if (is.null(output_root)) {
    output_root <- sample_sheet$output_dir[1]
    cat("output_root inferred from first sample: ", output_root)
  }

  agg_dir <- file.path(output_root, "aggregated", "salmon")
  dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)

  tpm_list    <- list()
  counts_list <- list()

  for (sid in sample_sheet$sample_id) {
    sample_out_dir <- sample_sheet$output_dir[sample_sheet$sample_id == sid]
    quant_file <- file.path(sample_out_dir, sid, "salmon", "quant.sf")

    if (!file.exists(quant_file)) {
      warning("quant.sf not found for sample '", sid, "' — skipped.")
      next
    }

    df <- read.table(quant_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tpm_list[[sid]]    <- setNames(df$TPM,      df$Name)
    counts_list[[sid]] <- setNames(df$NumReads, df$Name)
  }

  if (length(tpm_list) == 0) {
    stop("No quant.sf files found. Run HORIZON_run_salmon() for at least one sample first.")
  }

  # Align all samples on the same transcript set (should be identical since
  # all samples share the same index, but union handles edge cases gracefully)
  all_tx        <- Reduce(union, lapply(tpm_list, names))
  tpm_matrix    <- do.call(cbind, lapply(tpm_list,    function(x) x[all_tx]))
  counts_matrix <- do.call(cbind, lapply(counts_list, function(x) x[all_tx]))
  rownames(tpm_matrix)    <- all_tx
  rownames(counts_matrix) <- all_tx

  n_tx      <- nrow(tpm_matrix)
  n_samples <- ncol(tpm_matrix)
  cat("Aggregated Salmon matrices: ", n_tx, " transcripts x ", n_samples, " samples")

  write.csv(tpm_matrix,    file.path(agg_dir, "tpm_matrix.csv"))
  write.csv(counts_matrix, file.path(agg_dir, "counts_matrix.csv"))
  cat("TPM matrix written to: ",    file.path(agg_dir, "tpm_matrix.csv"))
  cat("Counts matrix written to: ", file.path(agg_dir, "counts_matrix.csv"))

  if (save_rds) {
    saveRDS(tpm_matrix,    file.path(agg_dir, "tpm_matrix.rds"))
    saveRDS(counts_matrix, file.path(agg_dir, "counts_matrix.rds"))
  }

  if (save_metadata) {
    meta <- sample_sheet[sample_sheet$sample_id %in% colnames(tpm_matrix), ]
    write.csv(meta, file.path(agg_dir, "sample_metadata.csv"), row.names = FALSE)
    cat("Sample metadata written to: ", file.path(agg_dir, "sample_metadata.csv"))
  }

  invisible(list(tpm = tpm_matrix, counts = counts_matrix))
}
