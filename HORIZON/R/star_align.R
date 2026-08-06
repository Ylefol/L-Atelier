# ==============================================================================
# HORIZON - STAR Alignment (multi-mapping-aware, for TE quantification)
# ==============================================================================
# Rsubread::align (see align.R / HORIZON_run_align) is tuned for gene-level
# counting: it reports essentially one best location per read, which is the
# correct behavior for standard gene expression but discards exactly the
# signal transposable-element (TE) quantification depends on -- most
# individual TE copies are not unique sequence, so a large fraction of
# TE-derived reads genuinely multi-map across many near-identical loci.
# STAR, run with a relaxed multi-mapping filter, reports those alignments
# instead of dropping them, which is what HORIZON_run_tecount() (see
# tetranscripts.R) needs to redistribute multi-mapping reads across TE loci.
#
# This is a separate alignment path from HORIZON_run_align() / Rsubread --
# both can coexist per sample (see aligned/ vs aligned_star/ subfolders).
#
# STAR and TEcount share the standard horizon_cli conda environment (see
# HORIZON_set_conda_env() / basilisk.R) -- both are called via
# .horizon_run_cli(), the same mechanism used for samtools/bedtools/macs3/
# deeptools/salmon. There is no dependency conflict requiring a dedicated
# environment (unlike split-pipe, which does need one -- see parse.R).
# ==============================================================================


#' Build a STAR genome index
#'
#' Wraps \code{STAR --runMode genomeGenerate} to build a splice-aware genome
#' index for use with \code{\link{HORIZON_run_star_align}}. This is a one-time
#' operation per reference genome/annotation/read-length combination.
#'
#' \code{STAR} must be available in the conda environment registered with
#' \code{\link{HORIZON_set_conda_env}}. Install it alongside the existing
#' \code{horizon_cli} tools:
#' \preformatted{
#'   conda install -n horizon_cli -c bioconda -c conda-forge star tetranscripts
#' }
#'
#' @param genome_fasta Character. Path to the reference genome FASTA.
#' @param gtf Character. Path to the gene GTF annotation (used to build the
#'   splice junction database via \code{--sjdbGTFfile}).
#' @param index_dir Character. Output directory for the STAR index (STAR
#'   requires a real directory, unlike Rsubread/Salmon's basename-prefix
#'   index convention).
#' @param sjdb_overhang Integer. STAR's \code{--sjdbOverhang}, ideally
#'   \code{(read_length - 1)} for the actual FASTQ read length of the run this
#'   index will be used for -- \strong{not} a value to leave at the default
#'   without checking. Default \code{100L} is STAR's own documented
#'   general-purpose value (works reasonably for most Illumina read lengths
#'   via a few bp of junction-overhang slack) but should be verified against
#'   the real read length, e.g. via
#'   \code{zcat <fastq> | head -400 | awk 'NR\%4==2\{print length($0)\}' | sort -u}.
#' @param threads Integer. Threads for \code{--runThreadN}. Default 8.
#' @param extra_flags Character vector of additional STAR \code{genomeGenerate}
#'   flags appended verbatim (e.g. \code{c("--genomeSAindexNbases", "13")} for
#'   unusually small genomes). Default \code{character(0)} -- STAR's own
#'   default of 14 is already appropriate for a standard human/mouse-sized
#'   genome and does not need overriding here.
#' @param force Logical. Rebuild even if an index already exists at
#'   \code{index_dir}. Default \code{FALSE}.
#'
#' @return Character. \code{index_dir}, invisibly.
#' @export
HORIZON_build_star_index <- function(genome_fasta,
                                      gtf,
                                      index_dir,
                                      sjdb_overhang = 100L,
                                      threads       = 8L,
                                      extra_flags   = character(0),
                                      force         = FALSE) {

  if (!file.exists(genome_fasta))
    stop("Genome FASTA not found: ", genome_fasta, call. = FALSE)
  if (!file.exists(gtf)) stop("GTF file not found: ", gtf, call. = FALSE)

  index_marker <- file.path(index_dir, "SAindex")
  if (!isTRUE(force) && file.exists(index_marker)) {
    cat("[STAR] Index already exists at:", index_dir,
        "— skipping (use force=TRUE to rebuild)\n")
    return(invisible(index_dir))
  }

  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "--runMode",         "genomeGenerate",
    "--genomeDir",       index_dir,
    "--genomeFastaFiles", genome_fasta,
    "--sjdbGTFfile",     gtf,
    "--sjdbOverhang",    as.integer(sjdb_overhang),
    "--runThreadN",      as.integer(threads),
    extra_flags
  )

  cat("[STAR] Building genome index\n")
  cat("    FASTA:         ", genome_fasta, "\n")
  cat("    GTF:           ", gtf, "\n")
  cat("    sjdbOverhang:  ", sjdb_overhang, "\n")
  cat("    Output:        ", index_dir, "\n\n")

  ret <- .horizon_run_cli("STAR", args)

  if (ret != 0L)
    stop("STAR genomeGenerate failed (exit code ", ret, ")", call. = FALSE)

  cat("[STAR] Index build complete\n")
  cat("    Output:", index_dir, "\n")
  invisible(index_dir)
}


#' Align reads with STAR, reporting multi-mapping alignments
#'
#' Wraps \code{STAR --runMode alignReads} with a relaxed multi-mapping filter
#' suitable for downstream transposable-element (TE) quantification via
#' \code{\link{HORIZON_run_tecount}}. This is a separate alignment path from
#' \code{\link{HORIZON_run_align}} (Rsubread) -- output is written to a
#' distinct \code{aligned_star/} subfolder so both alignment methods can
#' coexist per sample.
#'
#' By default, uses the trimmed FASTQ files produced by
#' \code{\link{HORIZON_run_qc_trim}}. Raw FASTQs from the sample sheet can be
#' used instead by setting \code{use_trimmed = FALSE}.
#'
#' \strong{Why the multi-mapping filter matters here:} STAR's own default
#' (\code{--outFilterMultimapNmax 10}) silently discards any read mapping to
#' more than 10 genomic loci. For gene-level RNA-seq this is a reasonable
#' behavior, but it systematically discards reads from exactly the
#' high-copy-number TE families this workflow is meant to quantify. Relaxing
#' this filter (default here: 100) keeps those reads in the alignment so
#' \code{\link{HORIZON_run_tecount}} can redistribute them across TE loci via
#' its EM algorithm rather than losing them entirely.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param index Character. Path to a STAR genome index directory (from
#'   \code{\link{HORIZON_build_star_index}}).
#' @param use_trimmed Logical. Use trimmed FASTQs from the \code{qc/}
#'   subfolder (output of \code{\link{HORIZON_run_qc_trim}}). Default TRUE.
#'   Set FALSE to align the raw FASTQs listed in the sample sheet.
#' @param threads Integer. Number of alignment threads (\code{--runThreadN}).
#'   Default 8.
#' @param outfilter_multimap_nmax Integer. STAR's
#'   \code{--outFilterMultimapNmax} -- maximum number of loci a read is
#'   allowed to map to before being discarded entirely. Default 100, the value
#'   cited in the TEtranscripts tutorial/paper as standard practice for TE
#'   quantification; worth a second look if a TE family in this genome build
#'   has a pathologically higher copy number than that.
#' @param win_anchor_multimap_nmax Integer. STAR's
#'   \code{--winAnchorMultimapNmax} -- maximum number of loci an "anchor" seed
#'   is allowed to map to during seed search. Default 100, paired with
#'   \code{outfilter_multimap_nmax} per the same TEtranscripts-tutorial
#'   recipe.
#' @param tmp_dir Character or \code{NULL}. Parent directory for STAR's own
#'   working temp directory (\code{--outTmpDir}), where it creates a FIFO to
#'   feed \code{--readFilesCommand zcat} output through. STAR's own default
#'   (\code{<outFileNamePrefix>_STARtmp}, alongside the output BAM) breaks
#'   when the sample sheet's \code{output_dir} points at an exFAT/NTFS-mounted
#'   drive, since those filesystems can't hold POSIX FIFOs -- the alignment
#'   output BAM itself has no such restriction. Default \code{NULL} uses
#'   \code{tempdir()} (the R session's own temp directory), which is always on
#'   the local filesystem regardless of where \code{output_dir} points. Pass
#'   an explicit path if \code{tempdir()} doesn't have enough free space for
#'   STAR's temp files.
#' @param force Logical. If \code{FALSE} (default), skip alignment when the
#'   sorted output BAM already exists. Set \code{TRUE} to realign.
#' @param verbose Logical. If \code{TRUE} (default), prints the full STAR
#'   command before executing it.
#' @param ... Additional STAR flags passed verbatim (e.g.
#'   \code{c("--outSAMattributes", "NH", "HI", "AS", "nM")} to add fields
#'   beyond STAR's own \code{Standard} set). Note that \code{Standard} already
#'   includes the \code{NH} (number-of-hits) tag \code{HORIZON_run_tecount()}
#'   needs for multi-mapper redistribution -- do not strip it via a custom
#'   \code{--outSAMattributes} list here.
#'
#' @return Character. Path to the sorted, indexed BAM file, invisibly.
#' @export
HORIZON_run_star_align <- function(sample_sheet,
                                    sample_id,
                                    index,
                                    use_trimmed              = TRUE,
                                    threads                  = 8L,
                                    outfilter_multimap_nmax  = 100L,
                                    win_anchor_multimap_nmax = 100L,
                                    tmp_dir                  = NULL,
                                    force                    = FALSE,
                                    verbose                  = TRUE,
                                    ...) {

  if (!dir.exists(index))
    stop("STAR index directory not found: ", index,
         "\nRun HORIZON_build_star_index() first.", call. = FALSE)

  row    <- .get_sample_row(sample_sheet, sample_id)
  paired <- as.logical(row$paired_end)

  out_dir <- file.path(row$output_dir, sample_id, "aligned_star")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (use_trimmed) {
    qc_dir <- file.path(row$output_dir, sample_id, "qc")
    r1     <- file.path(qc_dir, paste0(sample_id, "_trimmed_R1.fastq.gz"))
    r2     <- file.path(qc_dir, paste0(sample_id, "_trimmed_R2.fastq.gz"))
  } else {
    r1 <- row$fastq_r1
    r2 <- row$fastq_r2
  }

  if (!file.exists(r1)) stop("R1 FASTQ not found: ", r1)
  if (paired && (is.na(r2) || !file.exists(r2))) stop("R2 FASTQ not found: ", r2)

  prefix     <- file.path(out_dir, paste0(sample_id, "_"))
  star_bam   <- paste0(prefix, "Aligned.sortedByCoord.out.bam")
  sorted_bam <- file.path(out_dir, paste0(sample_id, "_sorted.bam"))

  if (!isTRUE(force) && file.exists(sorted_bam)) {
    cat("[STAR] Sorted BAM already exists for:", sample_id,
        "— skipping (use force=TRUE to realign)\n")
    return(invisible(sorted_bam))
  }

  # STAR's own default tmp dir (<outFileNamePrefix>_STARtmp) sits next to the
  # output BAM -- fine on a local filesystem, but fails with "could not
  # create FIFO file" when output_dir is on an exFAT/NTFS-mounted drive
  # (confirmed via `STAR --help`: STAR fully clears/manages this directory
  # itself, so a fresh subdirectory under tempdir() is safe to hand it).
  star_tmp <- file.path(if (is.null(tmp_dir)) tempdir() else tmp_dir,
                         paste0("STARtmp_", sample_id))
  if (dir.exists(star_tmp)) unlink(star_tmp, recursive = TRUE)

  args <- c(
    "--runMode",               "alignReads",
    "--genomeDir",             index,
    "--readFilesIn",           r1, if (paired) r2,
    "--readFilesCommand",      "zcat",
    "--outFileNamePrefix",     prefix,
    "--outTmpDir",             star_tmp,
    "--outSAMtype",            "BAM", "SortedByCoordinate",
    "--outFilterMultimapNmax", as.integer(outfilter_multimap_nmax),
    "--winAnchorMultimapNmax", as.integer(win_anchor_multimap_nmax),
    "--runThreadN",            as.integer(threads),
    c(...)
  )

  cat("[STAR] Aligning reads for:", sample_id, "\n")
  cat("    Index   : ", index, "\n")
  cat("    R1      : ", r1, "\n")
  if (paired) cat("    R2      : ", r2, "\n")
  cat("    Output  : ", sorted_bam, "\n\n")

  if (isTRUE(verbose))
    cat("[STAR] Command: STAR", paste(args, collapse = " "), "\n\n")

  ret <- .horizon_run_cli("STAR", args)

  if (ret != 0L)
    stop("STAR alignment failed for: ", sample_id,
         " (exit code ", ret, ")", call. = FALSE)

  if (!file.exists(star_bam))
    stop("STAR ran but expected output BAM not found: ", star_bam, call. = FALSE)

  file.rename(star_bam, sorted_bam)

  cat("[STAR] Indexing BAM for:", sample_id, "\n")
  Rsamtools::indexBam(sorted_bam)

  cat("[STAR] Alignment complete for:", sample_id, "\n")
  cat("    Output:", sorted_bam, "\n")
  invisible(sorted_bam)
}
