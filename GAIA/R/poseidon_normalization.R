###############################################################################
########################### Normalization #####################################
###############################################################################
# Omics-specific normalization methods. Currently mass spec normalization;
# this is where CPM/TPM/RPKM/quantile/TMM normalization would also live.

#' Normalize a Log2 Mass Spectrometry Intensity Matrix
#'
#' @description Normalizes a log2-scale DIA mass spectrometry intensity matrix
#' using either Variance Stabilizing Normalization (VSN, default) or per-sample
#' median centering.
#'
#' @param wide Numeric matrix, proteins (rows) × samples (cols). Values must be
#'   log2-transformed (as produced by \code{ELEUTHIA_load_massspec()}). Both
#'   methods expect log2 input; VSN back-transforms to raw intensities internally
#'   before fitting.
#' @param method Character. Normalization method: \code{"vsn"} (default) or
#'   \code{"median_centering"}.
#' @param verbose Logical. Print normalization summary. Default: TRUE.
#'
#' @return A numeric matrix of the same dimensions, normalized. NA values are
#'   preserved in their original positions.
#'
#' @details
#' \strong{VSN} (Huber et al. 2002, Bioinformatics): fits a generalized
#' log (glog/arsinh) calibration per sample via robust iterative regression on
#' the raw intensity scale. Simultaneously normalizes cross-sample shifts AND
#' stabilizes variance across the full intensity range (heteroscedasticity
#' correction) — proteins at low and high abundance are treated equally. Output
#' is on a log2-like scale (the "h" scale). Because the pipeline stores log2
#' data in \code{ms$wide}, VSN mode back-transforms (\code{2^wide}) before
#' fitting and re-applies the original NA mask after, since \code{vsn::justvsn()}
#' fitting is driven by observed values only.
#' Requires the \code{vsn} Bioconductor package:
#' \code{BiocManager::install("vsn")}.
#'
#' \strong{Median centering}: subtracts the per-sample median (non-NA values)
#' from each column. Removes run-to-run location shifts on the log2 scale but
#' does not correct heteroscedasticity. Retained as a lightweight fallback or
#' for comparison.
#'
#' For the analysis pipeline:
#' \preformatted{
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide, method = "vsn")
#' }
POSEIDON_normalize_massspec <- function(wide,
                                         method  = c("vsn", "median_centering"),
                                         verbose = TRUE) {

  if (!is.matrix(wide))
    stop("wide must be a numeric matrix (proteins x samples).")

  method <- match.arg(method)

  if (method == "vsn") {

    if (!requireNamespace("vsn", quietly = TRUE))
      stop("Package 'vsn' is required for VSN normalization.\n",
           "Install via: BiocManager::install('vsn')", call. = FALSE)

    na_mask <- is.na(wide)
    raw     <- 2^wide

    vsn_fit  <- vsn::justvsn(raw)
    norm_mat <- as.matrix(vsn_fit)
    dimnames(norm_mat) <- dimnames(wide)
    norm_mat[na_mask]  <- NA

    if (verbose) {
      cat("[POSEIDON] VSN normalization (mass spec):\n")
      cat("    Samples normalized:", ncol(norm_mat), "\n")
      cat("    Proteins:          ", nrow(norm_mat), "\n")
      cat("    NA values preserved:", sum(na_mask), "\n")
    }

  } else {

    col_medians <- apply(wide, 2L, stats::median, na.rm = TRUE)

    if (any(is.na(col_medians)))
      warning("Some samples have all-NA values; their medians are NA and they will not be centered.")

    norm_mat <- sweep(wide, 2L, col_medians, FUN = "-")

    if (verbose) {
      cat("[POSEIDON] Median centering (mass spec):\n")
      cat("    Samples normalized:", sum(!is.na(col_medians)), "\n")
      cat("    Median shift range: [",
          round(min(col_medians, na.rm = TRUE), 2), ", ",
          round(max(col_medians, na.rm = TRUE), 2), "]\n", sep = "")
    }

  }

  return(norm_mat)
}


#' Aggregate Transcript-Level Counts to Gene-Level Counts
#'
#' @description Collapses a transcript x sample count matrix (e.g. from
#' Nanopore/long-read quantification) into a gene x sample matrix by summing
#' transcript counts per parent gene, using the transcript-to-gene mapping
#' from a GTF annotation.
#'
#' @param counts Numeric matrix, transcripts (rows) x samples (cols).
#'   Row names must match the \code{transcript_id}(.\code{transcript_version})
#'   values in \code{gtf_file}.
#' @param gtf_file Path to the GTF annotation used to quantify \code{counts}.
#' @param verbose Logical. Print a summary. Default: TRUE.
#'
#' @return A numeric matrix, genes (rows) x samples (cols), with row names
#'   equal to \code{gene_id}(.\code{gene_version}) as found in \code{gtf_file}.
#'
#' @details
#' Transcript and gene IDs are versioned (e.g. \code{"ENST00000511072.5"})
#' when the GTF provides \code{transcript_version}/\code{gene_version}
#' attributes, and left unversioned otherwise -- matching whichever
#' convention the counts were generated with. Transcripts in \code{counts}
#' with no matching entry in the GTF are dropped (reported when
#' \code{verbose = TRUE}); this is expected for e.g. spike-ins or
#' annotation-version mismatches and should be checked if the dropped
#' fraction is large.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gene_counts <- POSEIDON_aggregate_transcript_to_gene(
#'   counts   = rna_data_Nano$counts,
#'   gtf_file = "data/Nanopore/Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf"
#' )
#' }
POSEIDON_aggregate_transcript_to_gene <- function(counts,
                                                    gtf_file,
                                                    verbose = TRUE) {

  if (!is.matrix(counts))
    stop("counts must be a matrix (transcripts x samples).", call. = FALSE)

  if (!file.exists(gtf_file))
    stop("GTF file not found: ", gtf_file, call. = FALSE)

  gtf <- rtracklayer::readGFF(gtf_file)
  gtf <- gtf[gtf$type == "transcript", ]

  transcript_id <- if ("transcript_version" %in% colnames(gtf)) {
    paste(gtf$transcript_id, gtf$transcript_version, sep = ".")
  } else {
    as.character(gtf$transcript_id)
  }
  gene_id <- if ("gene_version" %in% colnames(gtf)) {
    paste(gtf$gene_id, gtf$gene_version, sep = ".")
  } else {
    as.character(gtf$gene_id)
  }

  tx2gene <- stats::setNames(gene_id, transcript_id)

  matched <- rownames(counts) %in% names(tx2gene)
  n_total <- nrow(counts)
  n_matched <- sum(matched)

  if (n_matched == 0)
    stop("None of the transcript IDs in 'counts' were found in 'gtf_file'. ",
         "Check that both use the same ID format (versioned vs unversioned) ",
         "and annotation release.", call. = FALSE)

  gene_counts <- rowsum(counts[matched, , drop = FALSE],
                         group = tx2gene[rownames(counts)[matched]])

  if (verbose) {
    cat("[POSEIDON] Aggregated transcript counts to gene level:\n")
    cat("    Transcripts:", n_total, "\n")
    cat("    Matched to GTF:", n_matched, "\n")
    if (n_matched < n_total)
      cat("    Dropped (no GTF match):", n_total - n_matched, "\n")
    cat("    Genes:", nrow(gene_counts), "\n")
  }

  gene_counts
}
