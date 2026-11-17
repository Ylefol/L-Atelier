# ==============================================================================
# Motif-conditioned RNA-editing background (e.g. APOBEC TCW/UCC enrichment)
#
# Answers: "of every genomic occurrence of a given motif (read in TRANSCRIPT
# orientation), how many were actually covered by RNA-seq reads, and how many
# of those show the target edit?" -- compared against the same question for a
# control context (by default: the same reference base in any other flanking
# context).
#
# This REQUIRES the raw REDItools2 checkpoint file ({sample}_reditools.txt),
# not the parsed cyan_editing_result/RDS -- CYAN_run_reditools()'s own parser
# always drops non-substitution rows (see cyan_reditools.R), so the per-base
# "was this even covered" denominator only exists in that raw file, which is
# only genome-wide/complete when the sample was run with strict = FALSE.
#
# Strand handling (two independent, easily-conflated strand concepts):
#   1. REDItools2's own per-ROW strand call (0/1/2, the file's "Strand"
#      column). Empirically confirmed (2026-11-05, against the real GRCm39
#      FASTA): when strand == 1, the row's Reference/AllSubs are reported
#      as-is on the genomic (+) strand; when strand == 0, they are reported
#      complemented (i.e. relative to the genomic (-) strand); strand == 2
#      means REDItools couldn't determine this and the reported bases cannot
#      be trusted to resolve to a genomic strand at all.
#   2. The gene/transcript strand at that position, from the GTF -- this is
#      the orientation the motif itself is defined in (APOBEC's preference is
#      stated 5'->3' on the mRNA). Deliberately NOT inferred from #1: this
#      project already shipped one silent strand-direction bug
#      (strand_mode = 1 vs 2, see rna_editing_reditools.R history) caught only
#      indirectly via a REDIportal match-rate check, so the motif catalog's
#      orientation is pinned to an independent source (the GTF) rather than
#      trusting REDItools' own per-position call a second time.
# A row is therefore corrected in two steps: its own strand (0/1/2) first
# un-complements it back to genomic (+), then the GTF strand at that position
# re-complements it (if "-") into transcript orientation. strand == 2 rows
# keep their coverage (real read depth, strand-independent) but are excluded
# from the edited/not-edited call, since which base is genomically "alt" is
# exactly what strand == 2 means REDItools couldn't resolve.
#
# ADDED 2026-11-13: per-position alt_count/frequency (from REDItools2's own
# BaseCount[A,C,G,T] field, same two-step strand correction as above applied
# to all four counts together -- a strand flip relabels a count's base
# (A<->T, C<->G), it doesn't change its value). `edited` is unchanged
# (still the project's one boolean presence/absence definition); these are
# an additional, depth-normalized alternative for callers who want to
# threshold on read support/frequency rather than mere presence -- see
# CYAN_STATUS.txt (Mouse_A3B's `min_edits=0` depth-confound finding).
# ==============================================================================

#' Scan a motif's genomic occurrences for RNA-editing enrichment
#'
#' For a user-specified sequence motif (e.g. \code{"TCC"} for APOBEC's
#' preferred UCC context, tested at the middle base), finds every occurrence
#' of that motif in the exonic genome -- read in transcript (gene-strand)
#' orientation, not plain genomic orientation -- then cross-references the
#' raw REDItools2 checkpoint for one sample to determine, for each
#' occurrence: whether it was covered by RNA-seq reads at all, and whether it
#' shows the specified edit. The same is computed for a control context (by
#' default, the same reference base in any other flanking context) so the two
#' rates can be compared as an enrichment ratio.
#'
#' @details
#' This is a background/rate analysis layered on top of \emph{already-made}
#' REDItools2 substitution calls -- "edited" here means the row's own
#' \code{AllSubs} already contains the target substitution code (after strand
#' correction), not a newly-applied coverage/frequency threshold. This keeps
#' one definition of "edited" across the whole project rather than a second,
#' competing one.
#'
#' \code{edited} is a presence/absence call (REDItools2's own \code{AllSubs},
#' gated only by \code{CYAN_run_reditools()}'s \code{min_edits} -- a single
#' supporting read is enough under that function's default). This is
#' mechanically biased by sequencing depth: a higher-coverage position has
#' more chances for one erroneous/misaligned read to trip the same boolean,
#' independent of any real editing. \code{alt_count}/\code{frequency} (added
#' alongside \code{edited}, not replacing it) give the per-position
#' alt-base read count and its depth-normalized proportion
#' (\code{alt_count / coverage}), from REDItools2's own \code{BaseCount}
#' field -- letting a caller define their own, depth-robust edited/not
#' threshold (e.g. a minimum \code{alt_count} and/or \code{frequency})
#' downstream, without re-streaming the raw checkpoint a second time.
#'
#' Requires the raw REDItools2 \code{.txt} checkpoint
#' (\code{\link{CYAN_run_reditools}}'s \code{strict = FALSE} output), not the
#' parsed/returned \code{cyan_editing_result} -- see the file-level note.
#'
#' Motif matching only supports exact A/C/G/T strings (no IUPAC ambiguity
#' codes). \code{motif} must be given in DNA alphabet (\code{"TCC"}, not
#' \code{"UCC"}) since the reference FASTA and REDItools output never contain
#' \code{U}.
#'
#' Positions where an exon is annotated on both strands (overlapping
#' antisense genes) are excluded from the catalog entirely -- the motif's
#' reading direction can't be safely resolved there.
#'
#' @param raw_reditools_path Character. Path to a sample's raw
#'   \code{{sample}_reditools.txt} checkpoint file (one row per covered
#'   position; requires the sample to have been run with \code{strict =
#'   FALSE} in \code{\link{CYAN_run_reditools}}).
#' @param reference      Character. Reference FASTA path (must be
#'   \code{samtools faidx}-indexed), matching the genome build the sample was
#'   aligned/called against.
#' @param gtf             Character. Path to a GTF gene annotation matching
#'   \code{reference}, used only for its \code{exon} features' strand.
#' @param motif           Character. The motif to query, DNA alphabet, e.g.
#'   \code{"TCC"}. Required -- deliberately has no default, since no single
#'   motif is a safe assumption across use cases.
#' @param motif_position  Integer. 1-based index into \code{motif} of the
#'   base being tested for an edit (e.g. \code{2} for the middle base of
#'   \code{"TCC"}). Required, no default.
#' @param alt_base        Character. Single base (A/C/G/T): the edited
#'   product expected at \code{motif_position} (e.g. \code{"T"} for a C>T
#'   APOBEC edit). Required, no default.
#' @param output_dir      Character. Directory for result files (created if
#'   absent).
#' @param work_dir        Character. Directory for intermediate files
#'   (created if absent). Defaults to \code{tempdir()}.
#' @param sample_name     Character or \code{NULL}. Label for output
#'   filenames; derived from \code{raw_reditools_path} if \code{NULL}.
#' @param control_motif   Character or \code{NULL}. An alternate motif,
#'   same length as \code{motif} and sharing its base at \code{motif_position}
#'   (only the flanking context should differ), used as the control context
#'   instead of the default. Default \code{NULL}: control is every other
#'   occurrence of the same \code{motif_position} base, any flanking context.
#' @param min_coverage    Integer. Minimum read coverage for a catalog
#'   position to count as "surveilled" (toward the rate denominator). Default
#'   \code{10L} -- an unverified convenience default (matches
#'   \code{CYAN_flag_germline_variants()}'s \code{dna_min_coverage}), not an
#'   established standard; confirm before relying on it.
#' @param region          Character or \code{NULL}. Restrict the motif
#'   catalog to a genomic region in samtools format (e.g.
#'   \code{"chr21:1-10000000"}), or a bare chromosome name. \strong{Only
#'   restricts which catalog positions are searched for} -- the raw
#'   checkpoint is a flat text file with no per-chromosome index, so every
#'   line is still streamed regardless of \code{region}/\code{bed_file};
#'   this does not shrink the dominant (I/O) runtime cost. For a genuinely
#'   fast timing pilot, truncate the input file itself (e.g. \code{head -n
#'   <N>}) rather than relying on \code{region}.
#' @param bed_file        Character or \code{NULL}. Path to a BED file
#'   further restricting the catalog (same caveat as \code{region} above --
#'   does not reduce how much of the raw file is read). If both
#'   \code{region} and \code{bed_file} are given, their intersection is used.
#' @param chunk_size      Integer. Lines read and parsed per chunk while
#'   streaming the (potentially tens-of-GB) raw checkpoint file. Default
#'   \code{10e6}, matching \code{\link{CYAN_run_reditools}}'s own parser.
#' @param save_rds        Logical. If \code{TRUE} (default), save the full
#'   result as an RDS alongside the CSV outputs.
#' @param verbose         Logical. Print \code{[CYAN]} progress (default
#'   \code{TRUE}).
#'
#' @return An object of class \code{cyan_motif_result}: a list with
#'   \describe{
#'     \item{sites}{One row per catalog position (target + control): region,
#'       position, gtf_strand, group ("target"/"control"), coverage,
#'       surveilled, edited (\code{TRUE}/\code{FALSE}/\code{NA} -- \code{NA}
#'       when the row's own strand call was ambiguous), alt_count (reads
#'       supporting \code{alt_base} at \code{motif_position}, transcript-
#'       strand-corrected; \code{NA} under the same conditions as
#'       \code{edited}, or when the position had zero coverage), frequency
#'       (\code{alt_count / coverage}, \code{NA} where \code{alt_count} is
#'       \code{NA} or coverage is 0), transcript_subs, sample.}
#'     \item{summary}{One row per group: n_catalog, n_surveilled, n_edited,
#'       n_ambiguous_strand, edit_rate, sample.}
#'   }
#'   with the target/control enrichment ratio (\code{edit_rate} ratio)
#'   attached as the \code{"enrichment_ratio"} attribute.
#'
#' @seealso \code{\link{CYAN_scan_motif_editing_batch}},
#'   \code{\link{CYAN_run_reditools}}
#'
#' @export
CYAN_scan_motif_editing <- function(
  raw_reditools_path,
  reference,
  gtf,
  motif          = NULL,
  motif_position = NULL,
  alt_base       = NULL,
  output_dir,
  work_dir       = tempdir(),
  sample_name    = NULL,
  control_motif  = NULL,
  min_coverage   = 10L,
  region         = NULL,
  bed_file       = NULL,
  chunk_size     = 10e6,
  save_rds       = TRUE,
  verbose        = TRUE
) {
  .cyan_validate_motif_editing_args(motif, motif_position, alt_base)

  catalog_obj <- .cyan_build_motif_catalog(
    reference = reference, gtf = gtf, motif = motif,
    motif_position = motif_position, control_motif = control_motif,
    region = region, bed_file = bed_file, verbose = verbose
  )

  .cyan_scan_motif_editing_core(
    raw_reditools_path = raw_reditools_path, catalog_obj = catalog_obj,
    alt_base = toupper(alt_base), output_dir = output_dir, work_dir = work_dir,
    sample_name = sample_name, min_coverage = min_coverage,
    chunk_size = chunk_size, save_rds = save_rds, verbose = verbose
  )
}

#' Run \code{CYAN_scan_motif_editing()} across multiple samples
#'
#' Builds the motif catalog once (the expensive genome scan does not depend
#' on the sample) and reuses it across every raw REDItools2 checkpoint
#' supplied, rather than rebuilding it per sample.
#'
#' Samples are processed sequentially. A failed sample issues a warning and
#' is recorded as \code{NULL} in the returned list; remaining samples
#' continue.
#'
#' Each sample's full result (including its \code{sites} table, which can be
#' genome-catalog-sized -- tens of millions of rows) is always written to
#' disk per sample inside \code{\link{CYAN_scan_motif_editing}} (the
#' \code{_motif_sites.csv}/\code{_motif_summary.csv} outputs, plus
#' \code{_motif_result.rds} when \code{save_rds = TRUE}, the default). By
#' default this function does NOT also keep every sample's \code{sites}
#' table in the list it returns -- holding all of them simultaneously across
#' a real batch is what OOM-killed a 17-sample production run on its last
#' sample (confirmed 2026-11-09: a long-lived process accumulating one
#' catalog-sized data.frame per completed sample). Set \code{keep_sites =
#' TRUE} only for small batches (a handful of samples) where having every
#' \code{sites} table in-memory at once is actually wanted; for a real
#' cohort, read individual samples' \code{sites} back from disk as needed
#' instead.
#'
#' @param raw_reditools_paths Named or unnamed character vector of raw
#'   \code{{sample}_reditools.txt} checkpoint paths. Names (if present) are
#'   used as each sample's label; otherwise derived from the file names.
#' @param keep_sites   Logical. If \code{FALSE} (default), each sample's
#'   \code{sites} table is dropped (set to \code{NULL}) from the list this
#'   function returns, right after that sample's own on-disk outputs are
#'   written -- bounding this function's memory use to roughly one sample's
#'   \code{sites} table at a time instead of growing with every completed
#'   sample. If \code{TRUE}, the full per-sample result (including
#'   \code{sites}) is kept in the returned list, same as calling
#'   \code{\link{CYAN_scan_motif_editing}} directly in a loop and collecting
#'   the results yourself.
#' @param ... Passed to \code{CYAN_scan_motif_editing()} for each sample
#'   (e.g. \code{min_coverage}, \code{chunk_size}).
#' @inheritParams CYAN_scan_motif_editing
#'
#' @return A named list, one element per input path: each a
#'   \code{cyan_motif_result} (success) or \code{NULL} (failure). When
#'   \code{keep_sites = FALSE} (default), each non-NULL element's
#'   \code{sites} component is \code{NULL} -- only \code{summary} and the
#'   \code{"enrichment_ratio"}/\code{"sample_name"} attributes are retained;
#'   the full per-sample result is still on disk (see Details).
#'
#' @seealso \code{\link{CYAN_scan_motif_editing}}
#'
#' @export
CYAN_scan_motif_editing_batch <- function(
  raw_reditools_paths,
  reference,
  gtf,
  motif          = NULL,
  motif_position = NULL,
  alt_base       = NULL,
  output_dir,
  control_motif  = NULL,
  region         = NULL,
  bed_file       = NULL,
  keep_sites     = FALSE,
  verbose        = TRUE,
  ...
) {
  .cyan_validate_motif_editing_args(motif, motif_position, alt_base)

  catalog_obj <- .cyan_build_motif_catalog(
    reference = reference, gtf = gtf, motif = motif,
    motif_position = motif_position, control_motif = control_motif,
    region = region, bed_file = bed_file, verbose = verbose
  )

  labels <- names(raw_reditools_paths)
  if (is.null(labels))
    labels <- tools::file_path_sans_ext(sub("_reditools$", "",
                tools::file_path_sans_ext(basename(raw_reditools_paths))))

  results <- vector("list", length(raw_reditools_paths))
  names(results) <- labels

  for (i in seq_along(raw_reditools_paths)) {
    label <- labels[i]
    if (verbose) cat("\n[CYAN] ===== ", label, " =====\n", sep = "")
    res <- tryCatch(
      .cyan_scan_motif_editing_core(
        raw_reditools_path = raw_reditools_paths[[i]], catalog_obj = catalog_obj,
        alt_base = toupper(alt_base), output_dir = output_dir,
        work_dir = tempdir(), sample_name = label, verbose = verbose, ...
      ),
      error = function(e) {
        warning("[CYAN] Sample '", label, "' failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    # Drop the (potentially catalog-sized) sites table before accumulating
    # into `results` -- already safely on disk (see Details above); keeping
    # it here too is what grows this function's memory with every sample.
    if (!is.null(res) && !keep_sites) res$sites <- NULL
    results[[i]] <- res
  }

  n_ok <- sum(!vapply(results, is.null, logical(1)))
  if (verbose)
    cat("\n[CYAN] Motif-editing scan batch complete: ", n_ok, "/",
        length(results), " samples succeeded\n", sep = "")

  results
}

# ------------------------------------------------------------------------
# Internal: validation
# ------------------------------------------------------------------------

#' @keywords internal
.cyan_validate_motif_editing_args <- function(motif, motif_position, alt_base) {
  if (is.null(motif))
    stop("[CYAN] `motif` is required (no default) -- e.g. \"TCC\" for the ",
         "APOBEC-preferred motif. DNA alphabet only (\"T\", not \"U\").",
         call. = FALSE)
  if (is.null(motif_position))
    stop("[CYAN] `motif_position` is required (no default) -- the 1-based ",
         "index into `motif` of the base being tested for an edit.",
         call. = FALSE)
  if (is.null(alt_base))
    stop("[CYAN] `alt_base` is required (no default) -- the edited product ",
         "base expected at `motif_position` (e.g. \"T\" for a C>T edit).",
         call. = FALSE)
  if (!is.numeric(motif_position) || length(motif_position) != 1 ||
      motif_position < 1 || motif_position > nchar(motif))
    stop("[CYAN] `motif_position` must be a single integer between 1 and ",
         "nchar(motif).", call. = FALSE)
  alt_base <- toupper(alt_base)
  if (!grepl("^[ACGT]$", alt_base))
    stop("[CYAN] `alt_base` must be a single A/C/G/T base, got: ", alt_base,
         call. = FALSE)
  center_letter <- substr(toupper(motif), motif_position, motif_position)
  if (alt_base == center_letter)
    stop("[CYAN] `alt_base` ('", alt_base, "') is the same as the motif's ",
         "own base at `motif_position` ('", center_letter, "') -- that's not ",
         "an edit.", call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
.cyan_validate_motif <- function(motif, label = "motif") {
  motif <- toupper(motif)
  if (!grepl("^[ACGT]+$", motif))
    stop("[CYAN] `", label, "` must contain only A/C/G/T (DNA alphabet) -- got: ",
         motif, call. = FALSE)
  motif
}

# ------------------------------------------------------------------------
# Internal: base/substitution-code complementing
# ------------------------------------------------------------------------

#' @keywords internal
.cyan_complement_base <- function(x) chartr("ACGT", "TGCA", x)

#' Complement each code in a (possibly space-separated, possibly "-"/NA)
#' REDItools2 AllSubs string, preserving ref/alt pairing within each code.
#' @keywords internal
.cyan_complement_subs <- function(all_subs) {
  out <- all_subs
  not_na <- !is.na(all_subs) & all_subs != "-"
  if (any(not_na)) {
    codes_split <- strsplit(all_subs[not_na], " ", fixed = TRUE)
    out[not_na] <- vapply(
      codes_split,
      function(cs) paste(chartr("ACGT", "TGCA", cs), collapse = " "),
      character(1)
    )
  }
  out
}

# ------------------------------------------------------------------------
# Internal: region/BED restriction
# ------------------------------------------------------------------------

#' @keywords internal
.cyan_region_to_granges <- function(region = NULL, bed_file = NULL) {
  grs <- list()
  if (!is.null(region)) {
    if (grepl(":", region, fixed = TRUE)) {
      parts <- strsplit(region, "[:-]")[[1]]
      if (length(parts) != 3)
        stop("[CYAN] `region` must be \"chr\" or \"chr:start-end\", got: ",
             region, call. = FALSE)
      grs$region <- GenomicRanges::GRanges(
        parts[1], IRanges::IRanges(as.integer(parts[2]), as.integer(parts[3]))
      )
    } else {
      grs$region <- GenomicRanges::GRanges(region, IRanges::IRanges(1, .Machine$integer.max))
    }
  }
  if (!is.null(bed_file)) {
    if (!file.exists(bed_file))
      stop("[CYAN] `bed_file` not found: ", bed_file, call. = FALSE)
    grs$bed <- rtracklayer::import(bed_file)
  }
  if (length(grs) == 0) return(NULL)
  if (length(grs) == 1) return(grs[[1]])
  GenomicRanges::intersect(grs$region, grs$bed, ignore.strand = TRUE)
}

# ------------------------------------------------------------------------
# Internal: genome-wide, strand-aware motif scan
# ------------------------------------------------------------------------

#' Find every exonic, transcript-oriented occurrence of `motif`'s center
#' base, restricted to the (already ambiguous-excluded) per-strand exon
#' ranges. One pass per chromosome per strand over the full chromosome
#' sequence (Biostrings::matchPattern on ~100-250Mb DNAStrings is fast --
#' seconds, not minutes -- so this does not need chunking the way the raw
#' REDItools file does).
#' @keywords internal
.cyan_scan_motif_positions <- function(fa, exons_plus, exons_minus, motif,
                                        motif_position, seqlen, verbose = TRUE) {
  hits_list <- list()

  for (strnd in c("+", "-")) {
    exon_gr <- if (strnd == "+") exons_plus else exons_minus
    if (length(exon_gr) == 0) next
    chroms <- unique(as.character(GenomicRanges::seqnames(exon_gr)))

    for (chr in chroms) {
      len <- seqlen[[chr]]
      if (is.null(len)) next
      chr_gr  <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, len), strand = strnd)
      chr_seq <- Rsamtools::getSeq(fa, chr_gr)[[1]]

      m <- Biostrings::matchPattern(motif, chr_seq, fixed = TRUE)
      s <- BiocGenerics::start(m)
      if (length(s) == 0) next

      genomic_pos <- if (strnd == "+") {
        s + motif_position - 1L
      } else {
        len - (s + motif_position - 1L) + 1L
      }

      hit_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(genomic_pos, width = 1L), strand = strnd)
      hit_gr <- IRanges::subsetByOverlaps(hit_gr, exon_gr, ignore.strand = TRUE)

      if (length(hit_gr) > 0) hits_list[[length(hits_list) + 1L]] <- hit_gr
      if (verbose && length(hit_gr) > 0)
        cat("        ", chr, " (", strnd, "): ", length(hit_gr), " occurrences\n", sep = "")
    }
  }

  if (length(hits_list) == 0) return(GenomicRanges::GRanges())
  # Each element carries its own single-chromosome seqinfo (built from one
  # chromosome's length only) -- combining ~40+ of these raises a harmless
  # "no sequence levels in common" warning per pair; seqinfo isn't used
  # downstream (only seqnames/ranges/strand are), so this is cosmetic.
  suppressWarnings(do.call(c, hits_list))
}

# ------------------------------------------------------------------------
# Internal: catalog builder (target + control positions, strand-resolved)
# ------------------------------------------------------------------------

#' @keywords internal
.cyan_build_motif_catalog <- function(reference, gtf, motif, motif_position,
                                       control_motif = NULL, region = NULL,
                                       bed_file = NULL, verbose = TRUE) {
  reference <- CYAN_check_references(reference)
  if (!file.exists(gtf))
    stop("[CYAN] GTF file not found: ", gtf, call. = FALSE)

  motif <- .cyan_validate_motif(motif, "motif")
  motif_position <- as.integer(motif_position)
  center_letter <- substr(motif, motif_position, motif_position)

  if (!is.null(control_motif)) {
    control_motif <- .cyan_validate_motif(control_motif, "control_motif")
    if (nchar(control_motif) != nchar(motif))
      stop("[CYAN] `control_motif` must be the same length as `motif`.", call. = FALSE)
    if (substr(control_motif, motif_position, motif_position) != center_letter)
      stop("[CYAN] `control_motif` must share `motif`'s own base at ",
           "`motif_position` ('", center_letter, "') -- it should test the ",
           "same edit type in a different flanking context.", call. = FALSE)
  }

  if (verbose)
    cat("[CYAN] Building motif catalog: motif=", motif, " (position ",
        motif_position, " = '", center_letter, "'), gtf=", basename(gtf), "\n", sep = "")

  fa       <- Rsamtools::FaFile(reference)
  fa_index <- Rsamtools::scanFaIndex(fa)
  seqlen   <- setNames(BiocGenerics::width(fa_index), as.character(GenomicRanges::seqnames(fa_index)))

  restrict_gr <- .cyan_region_to_granges(region, bed_file)

  exons <- rtracklayer::import(gtf)
  exons <- exons[exons$type == "exon" & as.character(GenomicRanges::strand(exons)) %in% c("+", "-")]

  exons_plus  <- GenomicRanges::reduce(exons[GenomicRanges::strand(exons) == "+"], ignore.strand = TRUE)
  exons_minus <- GenomicRanges::reduce(exons[GenomicRanges::strand(exons) == "-"], ignore.strand = TRUE)
  GenomicRanges::strand(exons_plus)  <- "+"
  GenomicRanges::strand(exons_minus) <- "-"

  # Exon sequence annotated on both strands (overlapping antisense genes) --
  # the motif's reading direction can't be safely resolved there, so this
  # sequence is excluded from BOTH strand sets entirely rather than guessed.
  overlap_both   <- suppressWarnings(GenomicRanges::intersect(exons_plus, exons_minus, ignore.strand = TRUE))
  n_ambiguous_bp <- sum(IRanges::width(overlap_both))
  exons_plus  <- suppressWarnings(GenomicRanges::setdiff(exons_plus,  overlap_both, ignore.strand = TRUE))
  exons_minus <- suppressWarnings(GenomicRanges::setdiff(exons_minus, overlap_both, ignore.strand = TRUE))
  GenomicRanges::strand(exons_plus)  <- "+"
  GenomicRanges::strand(exons_minus) <- "-"

  if (n_ambiguous_bp > 0 && verbose)
    cat("    Excluded ", n_ambiguous_bp, " bp of dual-strand-annotated exon ",
        "sequence (ambiguous transcript orientation)\n", sep = "")

  if (!is.null(restrict_gr)) {
    exons_plus  <- suppressWarnings(GenomicRanges::intersect(exons_plus,  restrict_gr, ignore.strand = TRUE))
    exons_minus <- suppressWarnings(GenomicRanges::intersect(exons_minus, restrict_gr, ignore.strand = TRUE))
    GenomicRanges::strand(exons_plus)  <- "+"
    GenomicRanges::strand(exons_minus) <- "-"
  }

  if (verbose) cat("    Scanning for target motif '", motif, "'...\n", sep = "")
  target_gr <- .cyan_scan_motif_positions(fa, exons_plus, exons_minus, motif,
                                           motif_position, seqlen, verbose)

  if (!is.null(control_motif)) {
    if (verbose) cat("    Scanning for control motif '", control_motif, "'...\n", sep = "")
    control_gr <- .cyan_scan_motif_positions(fa, exons_plus, exons_minus, control_motif,
                                              motif_position, seqlen, verbose)
  } else {
    if (verbose) cat("    Scanning for all '", center_letter, "' occurrences ",
                      "(generic control)...\n", sep = "")
    all_center_gr <- .cyan_scan_motif_positions(fa, exons_plus, exons_minus, center_letter,
                                                 1L, seqlen, verbose)
    control_gr <- suppressWarnings(GenomicRanges::setdiff(all_center_gr, target_gr, ignore.strand = FALSE))
  }
  # Defensive: guarantees disjoint target/control sets even for a
  # user-supplied control_motif (should already be disjoint by construction,
  # since a genomic position's flanking context is fixed).
  control_gr <- suppressWarnings(GenomicRanges::setdiff(control_gr, target_gr, ignore.strand = FALSE))

  catalog <- data.table::data.table(
    region     = c(as.character(GenomicRanges::seqnames(target_gr)),
                   as.character(GenomicRanges::seqnames(control_gr))),
    position   = c(GenomicRanges::start(target_gr), GenomicRanges::start(control_gr)),
    gtf_strand = c(as.character(GenomicRanges::strand(target_gr)),
                   as.character(GenomicRanges::strand(control_gr))),
    group      = c(rep("target", length(target_gr)), rep("control", length(control_gr)))
  )
  data.table::setkey(catalog, region, position)

  if (verbose)
    cat("[CYAN] Catalog built: ", length(target_gr), " target / ",
        length(control_gr), " control positions\n", sep = "")

  list(
    catalog = catalog, center_letter = center_letter, motif = motif,
    motif_position = motif_position, control_motif = control_motif,
    n_ambiguous_bp = n_ambiguous_bp
  )
}

# ------------------------------------------------------------------------
# Internal: per-sample streaming scan against a prebuilt catalog
# ------------------------------------------------------------------------

#' @keywords internal
.cyan_scan_motif_editing_core <- function(raw_reditools_path, catalog_obj, alt_base,
                                           output_dir, work_dir = tempdir(),
                                           sample_name = NULL, min_coverage = 10L,
                                           chunk_size = 10e6, save_rds = TRUE,
                                           verbose = TRUE) {
  if (!file.exists(raw_reditools_path))
    stop("[CYAN] Raw REDItools2 checkpoint not found: ", raw_reditools_path,
         call. = FALSE)

  if (is.null(sample_name))
    sample_name <- tools::file_path_sans_ext(
      sub("_reditools$", "", tools::file_path_sans_ext(basename(raw_reditools_path)))
    )
  output_dir <- .cyan_ensure_dir(output_dir, "output directory")
  work_dir   <- .cyan_ensure_dir(work_dir, "work directory")

  catalog       <- catalog_obj$catalog
  center_letter <- catalog_obj$center_letter
  target_code   <- paste0(center_letter, alt_base)

  if (verbose)
    cat("[CYAN] Motif-editing scan for sample: ", sample_name, "\n",
        "    Catalog: ", sum(catalog$group == "target"), " target / ",
        sum(catalog$group == "control"), " control positions\n", sep = "")

  col_map <- c(
    "Region" = "region", "Position" = "position", "Reference" = "reference",
    "Strand" = "redi_strand",
    "Coverage-q30" = "coverage", "Coverage-q20" = "coverage", "Coverage" = "coverage",
    "AllSubs" = "all_subs",
    "BaseCount[ACGT]" = "base_count", "BaseCount[A,C,G,T]" = "base_count"
  )
  need_cols <- c("region", "position", "reference", "redi_strand", "coverage",
                 "all_subs", "base_count")

  con <- file(raw_reditools_path, open = "r")
  on.exit(close(con), add = TRUE)
  header_raw <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]

  hit_chunks  <- list()
  n_chunk     <- 0L
  n_rows_seen <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break
    n_chunk     <- n_chunk + 1L
    n_rows_seen <- n_rows_seen + length(lines)

    # fread's `text=` accepts the character vector directly (one line per
    # element) -- passing it straight through avoids an expensive
    # paste(collapse="\n") + re-tokenize round-trip. Benchmarked 2026-11-05
    # on a real 10M-line chunk of this project's data: ~33s (paste + fread
    # on the collapsed string) vs. ~8s (fread on the vector directly).
    raw <- data.table::fread(text = lines, sep = "\t",
                              header = FALSE, col.names = header_raw, colClasses = "character")
    for (old in names(col_map)) {
      idx <- match(old, colnames(raw))
      if (!is.na(idx)) colnames(raw)[idx] <- col_map[[old]]
    }
    present <- intersect(need_cols, colnames(raw))
    dt <- raw[, ..present]
    dt[, position := as.integer(position)]
    if ("coverage" %in% colnames(dt))    dt[, coverage := as.integer(coverage)]
    if ("redi_strand" %in% colnames(dt)) dt[, redi_strand := as.integer(redi_strand)]

    # Parse "BaseCount[A,C,G,T]" ("[n, n, n, n]") into four integer columns --
    # same vectorised tstrsplit pattern as .cyan_parse_reditools_output(), so
    # per-base read counts are available for a frequency (alt_count/coverage),
    # not just REDItools2's own presence/absence AllSubs call.
    if ("base_count" %in% colnames(dt)) {
      bc    <- gsub("[][ ]", "", dt$base_count)
      parts <- data.table::tstrsplit(bc, ",", fixed = TRUE)
      dt[, count_A := as.integer(parts[[1]])]
      dt[, count_C := as.integer(parts[[2]])]
      dt[, count_G := as.integer(parts[[3]])]
      dt[, count_T := as.integer(parts[[4]])]
      dt[, base_count := NULL]
    }

    hit <- merge(dt, catalog, by = c("region", "position"))
    if (nrow(hit) > 0) hit_chunks[[length(hit_chunks) + 1L]] <- hit

    if (verbose)
      cat("    Chunk", n_chunk, "(", length(lines), "lines) --", nrow(hit), "catalog matches\n")
  }

  if (verbose)
    cat("[CYAN] Finished streaming: ", n_rows_seen, " covered positions read, ",
        sum(vapply(hit_chunks, nrow, integer(1))), " catalog matches\n", sep = "")

  hits <- if (length(hit_chunks) > 0) {
    data.table::rbindlist(hit_chunks)
  } else {
    catalog[0][, `:=`(reference = character(0), redi_strand = integer(0),
                       coverage = integer(0), all_subs = character(0),
                       count_A = integer(0), count_C = integer(0),
                       count_G = integer(0), count_T = integer(0))]
  }

  # ---- Step 1: undo REDItools' own per-row strand complementing (0/1/2) to
  # get the genomic (+) orientation; strand == 2 -> unresolvable, kept as NA.
  hits[, genomic_ref := ifelse(
    redi_strand == 1L, reference,
    ifelse(redi_strand == 0L, .cyan_complement_base(reference), NA_character_)
  )]
  hits[, genomic_subs := ifelse(
    redi_strand == 1L, all_subs,
    ifelse(redi_strand == 0L, .cyan_complement_subs(all_subs), NA_character_)
  )]

  # ---- Step 2: apply the GTF (transcript) strand -- independent source from
  # step 1 -- to get the orientation the motif was defined in.
  hits[, transcript_ref := ifelse(gtf_strand == "+", genomic_ref, .cyan_complement_base(genomic_ref))]
  hits[, transcript_subs := ifelse(gtf_strand == "+", genomic_subs, .cyan_complement_subs(genomic_subs))]

  mismatch <- !is.na(hits$transcript_ref) & hits$transcript_ref != center_letter
  n_mismatch <- sum(mismatch)
  if (n_mismatch > 0 && verbose)
    cat("    NOTE: ", n_mismatch, " catalog matches have a reference-base ",
        "mismatch vs. the expected motif letter after strand correction -- ",
        "excluded from edited-calling (see function docs)\n", sep = "")

  hits[, edited := ifelse(
    redi_strand == 2L | mismatch, NA,
    vapply(strsplit(transcript_subs, " ", fixed = TRUE), function(codes) {
      !is.na(codes[1]) && target_code %in% codes
    }, logical(1))
  )]
  hits[, surveilled := coverage >= min_coverage]

  # ---- Per-base read counts, same two-step strand correction as
  # genomic_ref/transcript_ref above, applied to all four counts together --
  # relabeling a count under a strand flip means swapping which base it
  # belongs to (A<->T, C<->G), not changing its value.
  hits[, genomic_count_A := ifelse(redi_strand == 1L, count_A, ifelse(redi_strand == 0L, count_T, NA_integer_))]
  hits[, genomic_count_C := ifelse(redi_strand == 1L, count_C, ifelse(redi_strand == 0L, count_G, NA_integer_))]
  hits[, genomic_count_G := ifelse(redi_strand == 1L, count_G, ifelse(redi_strand == 0L, count_C, NA_integer_))]
  hits[, genomic_count_T := ifelse(redi_strand == 1L, count_T, ifelse(redi_strand == 0L, count_A, NA_integer_))]

  hits[, transcript_count_A := ifelse(gtf_strand == "+", genomic_count_A, genomic_count_T)]
  hits[, transcript_count_C := ifelse(gtf_strand == "+", genomic_count_C, genomic_count_G)]
  hits[, transcript_count_G := ifelse(gtf_strand == "+", genomic_count_G, genomic_count_C)]
  hits[, transcript_count_T := ifelse(gtf_strand == "+", genomic_count_T, genomic_count_A)]

  # alt_count/frequency: a depth-normalized alternative to the boolean
  # `edited` above (edited = "at least one supporting read", mechanically
  # biased toward higher-coverage positions). NA under the same conditions
  # as `edited` (unresolvable strand, or reference mismatch after
  # correction) -- a 0-count at those positions would be a fabricated value,
  # not a real "no reads support the edit" observation.
  alt_col <- paste0("transcript_count_", alt_base)
  hits[, alt_count := hits[[alt_col]]]
  hits[redi_strand == 2L | mismatch, alt_count := NA_integer_]
  hits[, frequency := ifelse(coverage > 0, alt_count / coverage, NA_real_)]

  # ---- Catalog positions never seen in the raw file at all: zero coverage
  # in this sample (not ambiguous -- genuinely not surveilled).
  if (nrow(hits) > 0) {
    matched_keys <- unique(hits[, .(region, position)])
    never_seen   <- catalog[!matched_keys, on = c("region", "position")]
  } else {
    never_seen <- data.table::copy(catalog)
  }
  if (nrow(never_seen) > 0)
    never_seen[, `:=`(
      reference = NA_character_, redi_strand = NA_integer_, coverage = 0L,
      all_subs = "-", genomic_ref = NA_character_, genomic_subs = NA_character_,
      transcript_ref = NA_character_, transcript_subs = NA_character_,
      edited = FALSE, surveilled = FALSE,
      alt_count = NA_integer_, frequency = NA_real_
    )]

  sites <- data.table::rbindlist(list(hits, never_seen), use.names = TRUE, fill = TRUE)
  sites[, sample := sample_name]

  summary_dt <- sites[, .(
    n_catalog          = .N,
    n_surveilled        = sum(surveilled, na.rm = TRUE),
    n_edited            = sum(edited & surveilled, na.rm = TRUE),
    n_ambiguous_strand  = sum(surveilled & is.na(edited), na.rm = TRUE)
  ), by = group]
  summary_dt[, edit_rate := n_edited / n_surveilled]
  summary_dt[, sample := sample_name]

  tgt_rate <- summary_dt[group == "target", edit_rate]
  ctl_rate <- summary_dt[group == "control", edit_rate]
  enrichment <- if (length(tgt_rate) == 1 && length(ctl_rate) == 1 &&
                     !is.na(ctl_rate) && ctl_rate > 0) tgt_rate / ctl_rate else NA_real_

  # fwrite(), not write.csv() -- benchmarked 2026-11-05 on a real ~29M-row
  # catalog: write.csv() ~49s vs. fwrite() ~2.8s for the same table. `sites`
  # can be genome-catalog-sized (every control-context position, covered or
  # not), so this is not a cosmetic difference at real scale.
  data.table::fwrite(sites, file.path(output_dir, paste0(sample_name, "_motif_sites.csv")))
  data.table::fwrite(summary_dt, file.path(output_dir, paste0(sample_name, "_motif_summary.csv")))

  # setDF() converts in place (no second copy) -- `sites` can be
  # genome-catalog-sized (see fwrite() comment above), and the prior
  # as.data.frame(sites) copy was confirmed (2026-11-09, Mouse_A3B
  # 17-sample batch) to OOM-kill the process on its last/largest sample,
  # right between the fwrite() calls above and the saveRDS() below.
  data.table::setDF(sites)
  result <- list(sites = sites, summary = as.data.frame(summary_dt))
  attr(result, "enrichment_ratio") <- enrichment
  attr(result, "sample_name") <- sample_name
  class(result) <- "cyan_motif_result"

  if (save_rds)
    saveRDS(result, file.path(output_dir, paste0(sample_name, "_motif_result.rds")))

  if (verbose) {
    tgt <- summary_dt[group == "target", ]
    ctl <- summary_dt[group == "control", ]
    cat("[CYAN] ", sample_name,
        " -- target: ", signif(tgt$edit_rate, 4), " (", tgt$n_edited, "/", tgt$n_surveilled, ")",
        ", control: ", signif(ctl$edit_rate, 4), " (", ctl$n_edited, "/", ctl$n_surveilled, ")",
        ", enrichment: ", signif(enrichment, 4), "\n", sep = "")
  }

  result
}
