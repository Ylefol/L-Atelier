#' Flag RNA-editing candidates that are actually germline variants
#'
#' Cross-references confident RNA-editing candidates (substitution-carrying
#' positions in a \code{cyan_editing_result} from \code{\link{CYAN_run_reditools}})
#' against matched whole-genome sequencing (WGS) data from the same
#' individual, to distinguish true post-transcriptional editing from
#' positions where the "edit" is really just a germline variant relative to
#' the reference genome.
#'
#' @details
#' This does not re-implement RNA-editing detection for DNA-seq data from
#' scratch. The WGS side is looked up with
#' \code{\link{CYAN_run_targeted_pileup}}, which seeks directly to exactly
#' the genomic positions flagged as RNA-editing candidates (via
#' \code{bed_file}) instead of scanning the whole BAM -- every candidate
#' position is reported, including reference-only/zero-coverage ones, so
#' "no row" is never ambiguous between "not covered" and "no variant".
#' (An earlier version of this function used \code{\link{CYAN_run_reditools}}
#' itself as the WGS-side pileup engine; REDItools2's \code{bed_file} option
#' only filters what gets written, not what gets scanned, so that approach
#' scanned the entire BAM regardless -- confirmed to take hours even
#' restricted to ~150 positions. \code{\link{CYAN_run_targeted_pileup}}
#' answers the same "coverage and base counts at these exact positions"
#' question without the whole-genome scan.)
#' The RNA-side call is left completely untouched by this function; nothing
#' about \code{\link{CYAN_run_reditools}} or \code{\link{CYAN_run_reditools_batch}}
#' changes, and this step is entirely optional -- an RNA-only workflow with
#' no matched WGS data is unaffected.
#'
#' For each RNA-editing candidate, the primary alt base (the substitution
#' with the most supporting RNA reads at that position) is looked up in the
#' WGS pileup and classified as:
#' \itemize{
#'   \item \code{"germline"} -- WGS shows the same alt allele at or above
#'     \code{dna_min_alt_reads} reads and \code{dna_min_alt_freq} frequency
#'   \item \code{"confirmed_somatic"} -- WGS is covered at
#'     \code{dna_min_coverage} or better and shows reference only
#'   \item \code{"undetermined"} -- WGS coverage at that position is below
#'     \code{dna_min_coverage}, too thin to call either way
#' }
#' Positions are tagged, never silently dropped -- callers decide downstream
#' whether to exclude \code{"undetermined"} sites or keep them with a caveat.
#'
#' \code{dna_min_coverage}, \code{dna_min_alt_reads}, and \code{dna_min_alt_freq}
#' default to the values used by the original REDItools' \code{REDItoolDnaRna.py}
#' (\code{-c 10}, \code{-v 3}, \code{-N 0.1} respectively) for exactly this
#' DNA/RNA comparison -- carried over as a reasonable starting point, not
#' independently re-derived. Verify against your own data before treating
#' them as settled.
#'
#' @param rna_result   A \code{cyan_editing_result} from
#'   \code{\link{CYAN_run_reditools}} (or one element of a
#'   \code{\link{CYAN_run_reditools_batch}} result), for the RNA-seq sample
#'   from the SAME individual as \code{dna_bam}.
#' @param dna_bam      Character. Path to a sorted, indexed, duplicate-marked
#'   WGS BAM from the same individual as \code{rna_result}, aligned to the
#'   same reference build.
#' @param reference    Character. Path to the reference FASTA (must match
#'   the one used for \code{rna_result}).
#' @param output_dir   Character. Directory for output files (created if
#'   absent).
#' @param work_dir     Character. Directory for intermediate/temporary files
#'   (created if absent). Defaults to \code{tempdir()}.
#' @param sample_name  Character or \code{NULL}. Label used in output
#'   filenames and as the DNA-side targeted-pileup sample prefix. Derived
#'   from \code{dna_bam}'s filename if \code{NULL}.
#' @param dna_min_coverage  Integer. Minimum WGS coverage at a position to
#'   call it confidently (either germline or confirmed-somatic) rather than
#'   \code{"undetermined"}. Default \code{10L}.
#' @param dna_min_alt_reads Integer. Minimum WGS reads supporting the alt
#'   allele to call a position germline. Default \code{3L}.
#' @param dna_min_alt_freq  Numeric. Minimum WGS alt-allele frequency to call
#'   a position germline. Default \code{0.1}.
#' @param require_dup_marked Logical. Passed through to the internal
#'   \code{\link{CYAN_run_targeted_pileup}} call on \code{dna_bam} -- see its
#'   own documentation. Default \code{TRUE}.
#' @param keep_tmp     Logical. If \code{FALSE} (default), the temporary
#'   candidate-site BED file (and the targeted-pileup call's own
#'   intermediates in \code{work_dir}) are removed after completion.
#' @param save_rds     Logical. If \code{TRUE} (default), the flagged result
#'   is saved as an RDS file in \code{output_dir}.
#' @param verbose      Logical. Print \code{[CYAN]} progress messages
#'   (default \code{TRUE}).
#' @param debug        Logical. Passed through to the internal
#'   \code{\link{CYAN_run_targeted_pileup}} call on \code{dna_bam} (default
#'   \code{FALSE}).
#' @param ...          Additional arguments forwarded to the internal
#'   \code{\link{CYAN_run_targeted_pileup}} call on \code{dna_bam} (e.g.
#'   \code{min_base_quality}, \code{min_read_quality}, \code{min_read_length}).
#'
#' @return \code{rna_result} with four columns appended:
#'   \code{dna_coverage}, \code{dna_alt_reads}, \code{dna_alt_freq},
#'   \code{germline_status}. Rows that were never RNA-editing candidates
#'   (\code{all_subs == "-"}) get \code{NA} in all four. The germline-flagging
#'   call parameters are attached as a \code{germline_params} attribute
#'   (the original \code{params} attribute from \code{rna_result} is left
#'   untouched).
#'
#' @seealso \code{\link{CYAN_run_reditools}}, \code{\link{CYAN_run_targeted_pileup}},
#'   \code{\link{CYAN_flag_germline_variants_batch}}
#'
#' @export
CYAN_flag_germline_variants <- function(
  rna_result,
  dna_bam,
  reference,
  output_dir,
  work_dir           = tempdir(),
  sample_name        = NULL,
  dna_min_coverage   = 10L,
  dna_min_alt_reads  = 3L,
  dna_min_alt_freq   = 0.1,
  require_dup_marked = TRUE,
  keep_tmp           = FALSE,
  save_rds           = TRUE,
  verbose            = TRUE,
  debug              = FALSE,
  ...
) {
  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  if (!inherits(rna_result, "cyan_editing_result"))
    stop("[CYAN] `rna_result` must be a cyan_editing_result object, ",
         "as returned by CYAN_run_reditools().", call. = FALSE)

  req_cols <- c("region", "position", "reference", "all_subs",
                "count_A", "count_C", "count_G", "count_T")
  missing_cols <- setdiff(req_cols, colnames(rna_result))
  if (length(missing_cols) > 0)
    stop("[CYAN] `rna_result` is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  dna_bam    <- CYAN_check_bam(dna_bam)
  reference  <- CYAN_check_references(reference)

  sample     <- .cyan_sample_name(dna_bam, sample_name)
  output_dir <- .cyan_ensure_dir(output_dir, "output directory")
  work_dir   <- .cyan_ensure_dir(work_dir,   "work directory")

  # ------------------------------------------------------------------
  # 2. Extract RNA-editing candidates and their primary alt base
  # ------------------------------------------------------------------
  is_candidate <- rna_result$all_subs != "-" & !is.na(rna_result$all_subs)
  if (!any(is_candidate))
    stop("[CYAN] `rna_result` contains no substitution-carrying positions ",
         "(all_subs == \"-\" for every row) -- nothing to compare against DNA.",
         call. = FALSE)

  cand <- rna_result[is_candidate, c("region", "position", "reference",
                                      "count_A", "count_C", "count_G", "count_T",
                                      "all_subs")]
  cand$primary_alt <- .cyan_primary_alt_base(
    cand$reference, cand$all_subs,
    cand$count_A, cand$count_C, cand$count_G, cand$count_T
  )
  n_unresolved <- sum(is.na(cand$primary_alt))
  cand <- cand[!is.na(cand$primary_alt), ]

  if (verbose) {
    cat("[CYAN] Germline-variant flagging for sample:", sample, "\n")
    cat("    RNA-editing candidates:", nrow(cand), "\n")
    if (n_unresolved > 0)
      cat("    Skipped (AllSubs didn't resolve to the reference base):",
          n_unresolved, "\n")
    cat("    DNA BAM:", dna_bam, "\n\n")
  }

  # ------------------------------------------------------------------
  # 3. Write candidate-site BED and run a targeted pileup on the WGS BAM
  # ------------------------------------------------------------------
  bed_path <- file.path(work_dir, paste0(sample, "_candidates.bed"))
  .cyan_write_bed(cand$region, cand$position, bed_path)

  dna_result <- CYAN_run_targeted_pileup(
    bam                 = dna_bam,
    reference           = reference,
    bed_file            = bed_path,
    output_dir          = output_dir,
    work_dir            = work_dir,
    sample_name         = paste0(sample, "_dna"),
    require_dup_marked  = require_dup_marked,
    keep_tmp            = keep_tmp,
    save_rds            = FALSE,
    verbose             = verbose,
    debug               = debug,
    ...
  )

  # ------------------------------------------------------------------
  # 4. Join and classify
  # ------------------------------------------------------------------
  dna_cols <- dna_result[, c("region", "position", "coverage",
                              "count_A", "count_C", "count_G", "count_T")]
  colnames(dna_cols) <- c("region", "position", "dna_coverage",
                           "dna_count_A", "dna_count_C", "dna_count_G", "dna_count_T")

  merged <- merge(cand, dna_cols, by = c("region", "position"), all.x = TRUE)

  alt_counts <- as.matrix(merged[, c("dna_count_A", "dna_count_C",
                                      "dna_count_G", "dna_count_T")])
  colnames(alt_counts) <- c("A", "C", "G", "T")
  merged$dna_alt_reads <- alt_counts[cbind(seq_len(nrow(merged)),
                                            match(merged$primary_alt, colnames(alt_counts)))]
  merged$dna_alt_freq <- ifelse(
    !is.na(merged$dna_coverage) & merged$dna_coverage > 0,
    merged$dna_alt_reads / merged$dna_coverage, NA_real_
  )

  merged$germline_status <- ifelse(
    is.na(merged$dna_coverage) | merged$dna_coverage < dna_min_coverage,
    "undetermined",
    ifelse(
      !is.na(merged$dna_alt_reads) &
        merged$dna_alt_reads >= dna_min_alt_reads &
        merged$dna_alt_freq  >= dna_min_alt_freq,
      "germline", "confirmed_somatic"
    )
  )

  # ------------------------------------------------------------------
  # 5. Merge classification back onto the full rna_result
  # ------------------------------------------------------------------
  out <- rna_result
  out$dna_coverage    <- NA_integer_
  out$dna_alt_reads    <- NA_integer_
  out$dna_alt_freq     <- NA_real_
  out$germline_status  <- NA_character_

  key_out <- paste(out$region, out$position)
  key_mrg <- paste(merged$region, merged$position)
  idx <- match(key_mrg, key_out)
  out$dna_coverage[idx]   <- merged$dna_coverage
  out$dna_alt_reads[idx]  <- merged$dna_alt_reads
  out$dna_alt_freq[idx]   <- merged$dna_alt_freq
  out$germline_status[idx] <- merged$germline_status

  class(out) <- c("cyan_editing_result", "data.frame")
  attr(out, "germline_params") <- as.list(match.call()[-1])

  # ------------------------------------------------------------------
  # 6. Save outputs
  # ------------------------------------------------------------------
  txt_out <- file.path(output_dir, paste0(sample, "_germline_flagged.txt"))
  utils::write.table(
    merged[, c("region", "position", "reference", "primary_alt",
               "dna_coverage", "dna_alt_reads", "dna_alt_freq", "germline_status")],
    txt_out, sep = "\t", row.names = FALSE, quote = FALSE
  )

  if (save_rds) {
    rds_out <- file.path(output_dir, paste0(sample, "_germline_flagged.rds"))
    saveRDS(out, rds_out)
  }

  if (!keep_tmp && file.exists(bed_path)) file.remove(bed_path)

  if (verbose) {
    tab <- table(merged$germline_status)
    cat("[CYAN] Classification complete for", sample, "\n")
    for (nm in names(tab))
      cat("    ", nm, ": ", tab[[nm]], "\n", sep = "")
    cat("    Candidates table:", txt_out, "\n")
    if (save_rds) cat("    RDS saved:", file.path(output_dir, paste0(sample, "_germline_flagged.rds")), "\n")
  }

  invisible(out)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Determine the primary (highest-supported) alt base from an AllSubs code
#'
#' \code{AllSubs} may list more than one substitution code at a position
#' (e.g. \code{"AG AC"}). Only codes whose first character matches the
#' reference base are meaningful; among those, the alt base with the most
#' supporting reads is taken as primary.
#' @keywords internal
.cyan_primary_alt_base <- function(reference, all_subs, count_A, count_C, count_G, count_T) {
  counts <- cbind(A = count_A, C = count_C, G = count_G, T = count_T)
  vapply(seq_along(all_subs), function(i) {
    codes <- strsplit(all_subs[i], "\\s+")[[1]]
    codes <- codes[nchar(codes) == 2 & substr(codes, 1, 1) == reference[i]]
    if (length(codes) == 0) return(NA_character_)
    alts <- substr(codes, 2, 2)
    alts[which.max(counts[i, alts])]
  }, character(1))
}

#' Write genomic positions to a BED file (0-based, half-open)
#' @keywords internal
.cyan_write_bed <- function(region, position, path) {
  bed <- unique(data.frame(region = region, start = position - 1L,
                            end = position, stringsAsFactors = FALSE))
  bed <- bed[order(bed$region, bed$start), ]
  utils::write.table(bed, path, sep = "\t", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
  invisible(path)
}


# ===========================================================================
# Batch wrapper
# ===========================================================================

#' Flag germline variants across multiple matched RNA/WGS sample pairs
#'
#' Convenience wrapper around \code{\link{CYAN_flag_germline_variants}} for
#' the (expected common) case of only some RNA-seq samples having matched
#' WGS data available. \code{dna_bam} only needs entries for animals that
#' actually have a WGS BAM -- samples with no match are simply absent from
#' \code{dna_bam} and are skipped, with no special-casing needed.
#'
#' @param rna_results Named list of \code{cyan_editing_result} objects (e.g.
#'   a \code{\link{CYAN_run_reditools_batch}} result), names = sample/animal
#'   IDs.
#' @param dna_bam     Named character vector of WGS BAM paths. Names must be
#'   a subset of \code{names(rna_results)} -- one entry per animal with
#'   matched WGS data.
#' @param reference   Character. Path to the reference FASTA (shared across
#'   all samples).
#' @param output_dir  Character. Directory for all output files (created if
#'   absent).
#' @param ...         Additional arguments passed to
#'   \code{\link{CYAN_flag_germline_variants}} (e.g. \code{dna_min_coverage},
#'   \code{work_dir}).
#' @param verbose     Logical. Print \code{[CYAN]} progress messages
#'   including a per-sample summary on completion (default \code{TRUE}).
#'
#' @return An object of class \code{cyan_batch_result}: a named list with one
#'   element per entry in \code{dna_bam}, each either a flagged
#'   \code{cyan_editing_result} (success) or \code{NULL} (failure).
#'
#' @seealso \code{\link{CYAN_flag_germline_variants}}
#'
#' @export
CYAN_flag_germline_variants_batch <- function(rna_results, dna_bam, reference,
                                              output_dir, ..., verbose = TRUE) {
  labels <- names(dna_bam)
  if (is.null(labels) || any(!nzchar(labels)))
    stop("[CYAN] `dna_bam` must be a fully named vector ",
         "(names = animal IDs matching `rna_results`).", call. = FALSE)

  missing <- setdiff(labels, names(rna_results))
  if (length(missing) > 0)
    stop("[CYAN] No matching rna_result for: ", paste(missing, collapse = ", "),
         call. = FALSE)

  results <- vector("list", length(labels))
  names(results) <- labels

  for (i in seq_along(labels)) {
    label <- labels[i]

    if (verbose)
      cat("[CYAN] Batch [", i, "/", length(labels), "] Sample: ", label, "\n", sep = "")

    tryCatch(
      results[[i]] <- CYAN_flag_germline_variants(
        rna_result  = rna_results[[label]],
        dna_bam     = dna_bam[[label]],
        reference   = reference,
        output_dir  = output_dir,
        sample_name = label,
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
    cat("\n[CYAN] Batch complete: ", sum(succeeded), "/", length(labels),
        " samples succeeded.\n", sep = "")
    if (any(succeeded)) {
      for (nm in names(results)[succeeded]) {
        tab <- table(results[[nm]]$germline_status)
        cat("    ", nm, ": ", paste(names(tab), tab, sep = "=", collapse = ", "),
            "\n", sep = "")
      }
    }
    if (any(!succeeded))
      cat("    Failed samples: ", paste(labels[!succeeded], collapse = ", "), "\n", sep = "")
  }

  class(results) <- c("cyan_batch_result", "list")
  invisible(results)
}
