#' Artemis - Isoform Switch Analysis Functions
#'
#' @description Differential isoform usage analysis using IsoformSwitchAnalyzeR.
#' Identifies genes where the relative usage of transcript isoforms changes
#' between conditions, independently of overall gene expression changes.


# ==============================================================================
# ISOFORM SWITCH TESTING
# ==============================================================================

#' Differential Isoform Usage (Isoform Switch) Analysis
#'
#' @description Tests for differential isoform usage between two conditions
#' using \pkg{IsoformSwitchAnalyzeR}. Identifies genes where the relative
#' proportion of transcripts changes between conditions, independently of total
#' gene expression changes.
#'
#' Complements \code{ARTEMIS_differential_counts()}: DE analysis detects
#' changes in total gene output; isoform switch analysis detects changes in
#' which form of the gene is produced. A gene can show switching without
#' overall DE, and vice versa.
#'
#' @param counts Numeric matrix of transcript counts (transcripts x samples).
#'   Rownames must be transcript IDs matching those in \code{gtf_path}. Raw
#'   integer counts (e.g. from IsoQuant or FLAMES) and estimated counts (e.g.
#'   from Salmon/kallisto) are both accepted.
#' @param targets Data.frame with sample metadata. Rownames must match
#'   colnames of \code{counts}. Must include a group column.
#' @param reference Character. The reference/baseline group label (denominator).
#' @param experiment Character. The experimental/comparison group label
#'   (numerator). Positive dIF means higher isoform usage in experiment.
#' @param group_col Character. Column in targets containing group labels.
#'   Default: \code{"group"}.
#' @param gtf_path Character. Path to GTF annotation file. Transcript IDs in
#'   the GTF must match rownames of \code{counts}.
#' @param method Character. Statistical testing method: \code{"DEXSeq"}
#'   (default, negative binomial) or \code{"satuRn"} (quasi-binomial;
#'   replaces DRIMSeq in IsoformSwitchAnalyzeR >= 2.x).
#' @param alpha Numeric. FDR significance threshold. Default: \code{0.05}.
#' @param dIF_cutoff Numeric. Minimum absolute delta Isoform Fraction required.
#'   A switch must change isoform usage by at least this amount to be reported.
#'   Default: \code{0.1} (10\%). This is the IsoformSwitchAnalyzeR default and
#'   is not an established universal standard — verify for your context.
#' @param min_gene_expression Numeric. Minimum mean expression per gene.
#'   Default: \code{1}.
#' @param min_transcript_expression Numeric. Minimum mean expression per
#'   transcript. Default: \code{1}.
#' @param fix_stringtie_annotation Logical. Passed to \code{importRdata()}'s
#'   \code{fixStringTieAnnotationProblem}. This ISAS feature re-derives
#'   \code{gene_id} for de novo/StringTie-style assemblies (which assign
#'   synthetic gene IDs) by mapping transcripts back to a reference gene via
#'   genomic overlap — but it does so by substituting a gene \emph{name}, not
#'   a stable gene ID, and gene names are not guaranteed unique across the
#'   genome (paralogs, pseudogenes, readthrough/antisense genes can share a
#'   base symbol). Confirmed empirically: enabling this on a clean, directly
#'   annotated reference GTF replaced the majority of \code{gene_id} values
#'   with gene symbols and produced spurious "same gene_id on different
#'   chromosomes" collisions, causing genes to be dropped entirely. Default:
#'   \code{FALSE} — only set \code{TRUE} if \code{gtf_path} is itself a
#'   StringTie/discovery-mode annotation with synthetic gene IDs that need
#'   rescuing.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @details
#' One warning from \pkg{IsoformSwitchAnalyzeR} is expected and not a bug:
#' "Using row.names as isoform_id" simply confirms it used the counts matrix
#' rownames as transcript IDs (this function's design).
#'
#' If you see "gene_ids or isoform_ids were not unique... removed N gene_id"
#' with \code{fix_stringtie_annotation = FALSE}, that reflects genuine
#' gene_ids spanning more than one chromosome in \code{gtf_path} — inspect
#' your annotation. If it appears with \code{fix_stringtie_annotation = TRUE},
#' it is likely the gene-symbol substitution artifact described above rather
#' than a real annotation problem.
#'
#' When \code{method = "satuRn"}, this function attaches
#' \pkg{SummarizedExperiment} internally (no-op if already attached). This
#' works around a confirmed upstream bug: \code{isoformSwitchTestSatuRn()}
#' calls the unqualified \code{rowData()} without importing it in
#' IsoformSwitchAnalyzeR's NAMESPACE (present in 2.2.0, the current
#' Bioconductor release). Since this function calls IsoformSwitchAnalyzeR via
#' \code{::} rather than \code{library()}, its Depends chain is never
#' attached, so \code{rowData()} would otherwise error with "could not find
#' function".
#'
#' \code{method = "satuRn"} also disables \code{satuRn::testDTU()}'s built-in
#' p-value calibration diagplots (a Frequency histogram and a Density plot,
#' one pair per contrast). Those are satuRn's own diagnostics, drawn as a side
#' effect, not part of this package's plotting layer (\pkg{AETHER}).
#'
#' This function does \emph{not} support confounder/covariate adjustment
#' (e.g. batch, sex, age). The design matrix built internally only contains
#' \code{sampleID} and \code{condition} columns. IsoformSwitchAnalyzeR itself
#' supports additional cofactor columns on the design matrix (taken into
#' account automatically by both \code{isoformSwitchTestDEXSeq()} and
#' \code{isoformSwitchTestSatuRn()}; see its vignette section "How to handle
#' confounding effects (including batches)"), but that path is not currently
#' exposed by this wrapper. If your data has known confounders, correct for
#' them prior to calling this function rather than relying on it here.
#'
#' @return An S3 object of class \code{"artemis_isoform_switch"} containing:
#'   \describe{
#'     \item{switch_list}{Full \code{switchAnalyzeRlist} (native ISAS format;
#'       pass directly to IsoformSwitchAnalyzeR functions for advanced use)}
#'     \item{isoform_results}{Data.frame of per-isoform results: gene_id,
#'       isoform_id, condition_1, condition_2, dIF, isoform_switch_q_value,
#'       gene_switch_q_value, iso_significant. Sorted by q-value.}
#'     \item{summary}{Top-level switch count summary from
#'       \code{extractSwitchSummary()}}
#'     \item{n_switches}{Number of genes with at least one significant isoform
#'       switch}
#'     \item{params}{List of parameters used}
#'   }
#'
#' @export
ARTEMIS_isoform_switch <- function(counts,
                                    targets,
                                    reference,
                                    experiment,
                                    group_col                 = "group",
                                    gtf_path,
                                    method                    = c("DEXSeq", "satuRn"),
                                    alpha                     = 0.05,
                                    dIF_cutoff                = 0.1,
                                    min_gene_expression       = 1,
                                    min_transcript_expression = 1,
                                    fix_stringtie_annotation  = FALSE,
                                    verbose                   = TRUE) {

  method <- match.arg(method)

  # --- Validate inputs --------------------------------------------------------
  if (!is.matrix(counts)) counts <- as.matrix(counts)

  if (is.null(rownames(counts)))
    stop("counts matrix must have rownames (transcript IDs)")

  if (!group_col %in% colnames(targets))
    stop("group_col '", group_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))

  common_samples <- intersect(colnames(counts), rownames(targets))
  if (length(common_samples) == 0)
    stop("No matching samples between counts colnames and targets rownames.\n",
         "  Counts columns : ", paste(head(colnames(counts), 3), collapse = ", "), "...\n",
         "  Targets rownames: ", paste(head(rownames(targets), 3), collapse = ", "), "...")

  counts  <- counts[, common_samples, drop = FALSE]
  targets <- targets[common_samples, , drop = FALSE]

  groups        <- as.character(targets[[group_col]])
  unique_groups <- unique(groups)

  if (!reference %in% unique_groups)
    stop("reference '", reference, "' not found in targets. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))
  if (!experiment %in% unique_groups)
    stop("experiment '", experiment, "' not found in targets. ",
         "Available groups: ", paste(unique_groups, collapse = ", "))

  # Subset to the two comparison groups
  keep    <- groups %in% c(reference, experiment)
  counts  <- counts[, keep, drop = FALSE]
  targets <- targets[keep, , drop = FALSE]
  groups  <- as.character(targets[[group_col]])

  n_reference  <- sum(groups == reference)
  n_experiment <- sum(groups == experiment)

  if (n_reference < 2)
    stop("Reference group '", reference, "' has only ", n_reference,
         " sample(s). At least 2 replicates required.")
  if (n_experiment < 2)
    stop("Experiment group '", experiment, "' has only ", n_experiment,
         " sample(s). At least 2 replicates required.")

  if (!file.exists(gtf_path))
    stop("gtf_path not found: ", gtf_path)

  if (verbose) {
    cat("[ARTEMIS] Isoform Switch Analysis (IsoformSwitchAnalyzeR)\n")
    cat("    Reference   : ", reference, " (n = ", n_reference, ")\n", sep = "")
    cat("    Experiment  : ", experiment, " (n = ", n_experiment, ")\n", sep = "")
    cat("    Transcripts : ", nrow(counts), "\n", sep = "")
    cat("    Method      : ", method, "\n", sep = "")
    cat("    alpha       : ", alpha, " | dIF cutoff: ", dIF_cutoff, "\n\n", sep = "")
  }

  # --- Build ISAS design matrix -----------------------------------------------
  design_mat <- data.frame(
    sampleID  = rownames(targets),
    condition = groups,
    stringsAsFactors = FALSE
  )

  comparisons_df <- data.frame(
    condition_1      = reference,
    condition_2      = experiment,
    stringsAsFactors = FALSE
  )

  # --- Import into switchAnalyzeRlist -----------------------------------------
  if (verbose) cat("[ARTEMIS] Importing data into switchAnalyzeRlist...\n")

  sar <- IsoformSwitchAnalyzeR::importRdata(
    isoformCountMatrix  = counts,
    designMatrix        = design_mat,
    isoformExonAnnoation = gtf_path,
    comparisonsToMake   = comparisons_df,
    addAnnotatedORFs    = TRUE,
    fixStringTieAnnotationProblem = fix_stringtie_annotation,
    quiet               = !verbose
  )

  # --- Pre-filter low-expression features -------------------------------------
  if (verbose) cat("[ARTEMIS] Pre-filtering low-expression transcripts...\n")

  sar <- IsoformSwitchAnalyzeR::preFilter(
    switchAnalyzeRlist       = sar,
    geneExpressionCutoff     = min_gene_expression,
    isoformExpressionCutoff  = min_transcript_expression,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes   = FALSE,
    quiet                    = !verbose
  )

  # --- Statistical testing ----------------------------------------------------
  if (verbose) cat("[ARTEMIS] Running", method, "switch test...\n")

  if (method == "DEXSeq") {
    sar <- IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(
      switchAnalyzeRlist  = sar,
      alpha               = alpha,
      dIFcutoff           = dIF_cutoff,
      reduceToSwitchingGenes = FALSE,
      quiet               = !verbose
    )
  } else {
    # Workaround for a confirmed upstream bug: isoformSwitchTestSatuRn() calls
    # the unqualified rowData() without importing it in IsoformSwitchAnalyzeR's
    # NAMESPACE (still present in 2.2.0, the current Bioconductor release).
    # Because we call IsoformSwitchAnalyzeR via `::` rather than library(), its
    # Depends chain (which would normally attach SummarizedExperiment) never
    # gets attached, so rowData() is otherwise unresolvable at call time.
    # Attaching it here makes it reachable through R's namespace fallback
    # lookup, with no effect if it is already attached.
    if (!"package:SummarizedExperiment" %in% search()) {
      suppressPackageStartupMessages(attachNamespace("SummarizedExperiment"))
    }

    # satuRn::testDTU() (called internally) draws p-value calibration diagplots
    # (a Frequency histogram and a Density plot per contrast) as a side effect.
    # These are satuRn's own diagnostics, not part of this package's plotting
    # layer (AETHER) -- disabled here rather than left to print unexpectedly.
    sar <- IsoformSwitchAnalyzeR::isoformSwitchTestSatuRn(
      switchAnalyzeRlist  = sar,
      alpha               = alpha,
      dIFcutoff           = dIF_cutoff,
      reduceToSwitchingGenes = FALSE,
      diagplots           = FALSE,
      quiet               = !verbose
    )
  }

  # --- Extract per-isoform results table --------------------------------------
  iso_feat    <- sar$isoformFeatures
  result_cols <- intersect(
    c("gene_id", "isoform_id", "condition_1", "condition_2",
      "dIF", "isoform_switch_q_value", "gene_switch_q_value",
      "iso_significant", "gene_significant"),
    colnames(iso_feat)
  )
  isoform_results <- iso_feat[, result_cols, drop = FALSE]
  isoform_results <- isoform_results[
    order(isoform_results$isoform_switch_q_value, na.last = TRUE), ]
  rownames(isoform_results) <- NULL

  # --- Top-level summary ------------------------------------------------------
  switch_summary <- tryCatch(
    IsoformSwitchAnalyzeR::extractSwitchSummary(
      switchAnalyzeRlist = sar,
      alpha              = alpha,
      dIFcutoff          = dIF_cutoff,
      onlySigIsoforms    = FALSE
    ),
    error = function(e) NULL
  )

  # Count significant switching genes
  if (!is.null(switch_summary) && "Genes" %in% colnames(switch_summary)) {
    n_switches <- sum(switch_summary$Genes, na.rm = TRUE)
  } else if ("gene_significant" %in% colnames(isoform_results)) {
    n_switches <- length(unique(
      isoform_results$gene_id[isoform_results$gene_significant %in% TRUE]
    ))
  } else {
    n_switches <- sum(
      !is.na(isoform_results$isoform_switch_q_value) &
        isoform_results$isoform_switch_q_value < alpha, na.rm = TRUE
    )
  }

  if (verbose) {
    cat("\n[ARTEMIS] Results (FDR < ", alpha, ", |dIF| >= ", dIF_cutoff, "):\n", sep = "")
    cat("    Significant switching genes: ", n_switches, "\n", sep = "")
    if (n_switches > 0) {
      sig_rows <- isoform_results[
        !is.na(isoform_results$isoform_switch_q_value) &
          isoform_results$isoform_switch_q_value < alpha, ]
      top <- head(sig_rows, 5)
      cat("    Top switches:\n")
      for (i in seq_len(nrow(top))) {
        cat("        ", top$gene_id[i], " | ", top$isoform_id[i],
            " | dIF = ", round(top$dIF[i], 3),
            " | q = ", format.pval(top$isoform_switch_q_value[i], digits = 2),
            "\n", sep = "")
      }
    }
    cat("\n")
  }

  result <- list(
    switch_list     = sar,
    isoform_results = isoform_results,
    summary         = switch_summary,
    n_switches      = n_switches,
    params          = list(
      reference                 = reference,
      experiment                = experiment,
      group_col                 = group_col,
      method                    = method,
      alpha                     = alpha,
      dIF_cutoff                = dIF_cutoff,
      gtf_path                  = gtf_path,
      min_gene_expression       = min_gene_expression,
      min_transcript_expression = min_transcript_expression,
      fix_stringtie_annotation  = fix_stringtie_annotation
    )
  )
  class(result) <- c("artemis_isoform_switch", "list")
  return(result)
}


#' @method print artemis_isoform_switch
#' @export
print.artemis_isoform_switch <- function(x, ...) {
  cat("Isoform Switch Analysis\n")
  cat("------------------------------\n")
  cat("Comparison  :", x$params$experiment, "vs", x$params$reference, "\n")
  cat("Method      :", x$params$method, "\n")
  cat("alpha       :", x$params$alpha, " | dIF cutoff:", x$params$dIF_cutoff, "\n")
  cat("Sig. genes  :", x$n_switches, "\n")
  if (!is.null(x$consequence_summary))
    cat("Consequences: annotated\n")
  invisible(x)
}


# ==============================================================================
# CONSEQUENCE ANALYSIS
# ==============================================================================

#' Analyze Functional Consequences of Isoform Switches
#'
#' @description Adds biological consequence annotations to isoform switch
#' results. Integrates outputs from external tools (CPAT, SignalP, Pfam/HMMER)
#' and identifies functional differences between switching isoform pairs:
#' changes in coding potential, NMD sensitivity, protein domains, signal
#' peptides, and intron retention.
#'
#' @param switch_result An \code{artemis_isoform_switch} object from
#'   \code{ARTEMIS_isoform_switch()}.
#' @param cpat_file Character or NULL. Path to CPAT output file for coding
#'   potential annotation. Default: NULL (skipped).
#' @param signalp_file Character or NULL. Path to SignalP output file for
#'   signal peptide annotation. Default: NULL (skipped).
#' @param pfam_file Character or NULL. Path to Pfam/HMMER output file for
#'   protein domain annotation. Default: NULL (skipped).
#' @param cpat_cutoff Numeric. Coding probability cutoff for CPAT. Human
#'   default is 0.725, mouse is 0.44. No universal standard — verify for your
#'   organism. Only used if \code{cpat_file} is provided. Default: \code{0.725}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return The input \code{artemis_isoform_switch} object with two additional
#'   fields:
#'   \describe{
#'     \item{consequence_summary}{Data.frame from
#'       \code{sar$switchConsequence} — one row per isoform-switch pair with
#'       annotated consequence types}
#'     \item{switch_list}{Updated switchAnalyzeRlist with consequence
#'       annotations}
#'   }
#'
#' @details
#' At least one of \code{cpat_file}, \code{signalp_file}, or \code{pfam_file}
#' must be provided. Structural consequences (intron retention, ORF sequence
#' similarity, NMD sensitivity) are always included when ORF annotations are
#' present from the initial \code{ARTEMIS_isoform_switch()} call (they require
#' no external tools).
#'
#' External tools must be run independently before calling this function.
#' See the IsoformSwitchAnalyzeR vignette for expected file formats and
#' instructions for running each tool.
#'
#' @export
ARTEMIS_isoform_switch_consequences <- function(switch_result,
                                                  cpat_file    = NULL,
                                                  signalp_file = NULL,
                                                  pfam_file    = NULL,
                                                  cpat_cutoff  = 0.725,
                                                  verbose      = TRUE) {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()")

  if (is.null(cpat_file) && is.null(signalp_file) && is.null(pfam_file))
    stop("At least one of cpat_file, signalp_file, or pfam_file must be provided.")

  sar <- switch_result$switch_list

  if (verbose) cat("[ARTEMIS] Analyzing isoform switch consequences...\n\n")

  # --- Integrate external tool outputs ----------------------------------------
  if (!is.null(cpat_file)) {
    if (!file.exists(cpat_file))
      stop("cpat_file not found: ", cpat_file)
    if (verbose) cat("[ARTEMIS] Adding CPAT coding potential...\n")
    sar <- IsoformSwitchAnalyzeR::analyzeCPAT(
      switchAnalyzeRlist   = sar,
      pathToCPATresultFile = cpat_file,
      codingCutoff         = cpat_cutoff,
      removeNoncodinORFs   = TRUE,
      quiet                = !verbose
    )
  }

  if (!is.null(signalp_file)) {
    if (!file.exists(signalp_file))
      stop("signalp_file not found: ", signalp_file)
    if (verbose) cat("[ARTEMIS] Adding SignalP signal peptides...\n")
    sar <- IsoformSwitchAnalyzeR::analyzeSignalP(
      switchAnalyzeRlist      = sar,
      pathToSignalPresultFile = signalp_file,
      quiet                   = !verbose
    )
  }

  if (!is.null(pfam_file)) {
    if (!file.exists(pfam_file))
      stop("pfam_file not found: ", pfam_file)
    if (verbose) cat("[ARTEMIS] Adding Pfam protein domains...\n")
    sar <- IsoformSwitchAnalyzeR::analyzePFAM(
      switchAnalyzeRlist   = sar,
      pathToPFAMresultFile = pfam_file,
      quiet                = !verbose
    )
  }

  # --- Build consequence set dynamically --------------------------------------
  # Structural consequences are always included (require only ORF annotations
  # from addAnnotatedORFs=TRUE in importRdata, no external tool needed)
  consequences <- c("intron_retention", "ORF_seq_similarity", "NMD_status")
  if (!is.null(cpat_file))    consequences <- c(consequences, "coding_potential")
  if (!is.null(pfam_file))    consequences <- c(consequences, "domains_identified")
  if (!is.null(signalp_file)) consequences <- c(consequences, "signal_peptide_identified")

  if (verbose) {
    cat("[ARTEMIS] Consequences to analyze:\n")
    cat("    ", paste(consequences, collapse = ", "), "\n\n")
  }

  sar <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
    switchAnalyzeRlist    = sar,
    consequencesToAnalyze = consequences,
    quiet                 = !verbose
  )

  # Extract consequence table from the switchAnalyzeRlist slot
  consequence_summary <- sar$switchConsequence

  if (verbose) {
    n_genes <- if (!is.null(consequence_summary) && "gene_id" %in% colnames(consequence_summary))
      length(unique(consequence_summary$gene_id)) else 0L
    cat("[ARTEMIS] Consequences annotated for", n_genes, "genes.\n\n")
  }

  switch_result$switch_list         <- sar
  switch_result$consequence_summary <- consequence_summary
  return(switch_result)
}
