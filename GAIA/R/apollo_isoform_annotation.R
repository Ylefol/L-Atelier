#' Apollo - Isoform Switch Consequence Tool Wrappers
#'
#' @description Runs external protein/transcript annotation tools (CPAT,
#' SignalP, Pfam via \code{pfam_scan.pl}) against isoforms from an
#' \code{artemis_isoform_switch} object, producing result files consumable by
#' \code{ARTEMIS_isoform_switch_consequences()}. All three run inside
#' isolated \pkg{basilisk} environments (see \code{GAIA/R/basilisk.R}) --
#' CPAT and SignalP are Python tools; Pfam is packaged via bioconda's
#' \code{pfam_scan} recipe, which basilisk installs the same way regardless
#' of it being a Perl/native-binary tool rather than a Python one.


# ==============================================================================
# SEQUENCE EXTRACTION
# ==============================================================================

#' Extract Isoform Nucleotide/Amino Acid Sequences
#'
#' @description Extracts nucleotide and amino acid FASTA sequences for the
#' isoforms in an \code{artemis_isoform_switch} object, via
#' \code{IsoformSwitchAnalyzeR::extractSequence()}. This is the required
#' first step before running any of \code{APOLLO_run_cpat()},
#' \code{APOLLO_run_signalp()}, or \code{APOLLO_run_pfam()} -- the nucleotide
#' FASTA feeds CPAT, the amino acid FASTA feeds SignalP and Pfam.
#'
#' @param switch_result An \code{artemis_isoform_switch} object from
#'   \code{ARTEMIS_isoform_switch()}.
#' @param genome_object A \code{BSgenome} object matching the genome used to
#'   build \code{switch_result} (e.g. \code{BSgenome.Hsapiens.UCSC.hg38::Hsapiens}
#'   after loading that annotation package). This is a hard requirement of
#'   \code{IsoformSwitchAnalyzeR::extractSequence()} -- a plain genome FASTA
#'   path is not accepted.
#' @param output_dir Character. Directory the FASTA files are written to
#'   (created if absent).
#' @param only_switching_genes Logical. If \code{TRUE} (default), only extract
#'   sequences for genes already called significant in \code{switch_result}
#'   (using its own \code{alpha}/\code{dIF_cutoff}). Set \code{FALSE} to
#'   extract for every tested gene.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return The input \code{artemis_isoform_switch} object with two additional
#'   fields:
#'   \describe{
#'     \item{sequence_files}{List with \code{nt_fasta} and \code{aa_fasta}
#'       paths, used as defaults by \code{APOLLO_run_cpat()},
#'       \code{APOLLO_run_signalp()}, and \code{APOLLO_run_pfam()}.}
#'     \item{switch_list}{Updated switchAnalyzeRlist with the sequences also
#'       cached internally (\code{ntSequence}/\code{aaSequence}).}
#'   }
#'
#' @export
APOLLO_extract_isoform_sequences <- function(switch_result,
                                              genome_object,
                                              output_dir,
                                              only_switching_genes = TRUE,
                                              verbose = TRUE) {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (!methods::is(genome_object, "BSgenome"))
    stop("genome_object must be a BSgenome object (e.g. load the appropriate ",
         "BSgenome.* annotation package and pass its genome object, such as ",
         "BSgenome.Hsapiens.UCSC.hg38::Hsapiens). A plain genome FASTA path is ",
         "not accepted by IsoformSwitchAnalyzeR::extractSequence().", call. = FALSE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (verbose) {
    cat("[APOLLO] Extracting isoform sequences (IsoformSwitchAnalyzeR)\n")
    cat("    Output dir         : ", output_dir, "\n", sep = "")
    cat("    Only switching genes: ", only_switching_genes, "\n\n", sep = "")
  }

  sar <- IsoformSwitchAnalyzeR::extractSequence(
    switchAnalyzeRlist = switch_result$switch_list,
    genomeObject        = genome_object,
    onlySwitchingGenes  = only_switching_genes,
    alpha               = switch_result$params$alpha,
    dIFcutoff           = switch_result$params$dIF_cutoff,
    writeToFile         = TRUE,
    pathToOutput        = output_dir,
    outputPrefix        = "isoform_switch",
    quiet               = !verbose
  )

  nt_fasta <- file.path(output_dir, "isoform_switch_nt.fasta")
  aa_fasta <- file.path(output_dir, "isoform_switch_AA.fasta")

  if (!file.exists(nt_fasta))
    stop("extractSequence() completed but nucleotide FASTA was not created: ",
         nt_fasta, call. = FALSE)
  if (!file.exists(aa_fasta))
    stop("extractSequence() completed but amino acid FASTA was not created: ",
         aa_fasta, ". This usually means no ORFs are annotated -- check that ",
         "ARTEMIS_isoform_switch() was run with its default addAnnotatedORFs=TRUE.",
         call. = FALSE)

  if (verbose) {
    cat("[APOLLO] Sequences extracted.\n")
    cat("    Nucleotide FASTA: ", nt_fasta, "\n", sep = "")
    cat("    Amino acid FASTA: ", aa_fasta, "\n\n", sep = "")
  }

  switch_result$switch_list    <- sar
  switch_result$sequence_files <- list(nt_fasta = nt_fasta, aa_fasta = aa_fasta)
  switch_result
}


# ==============================================================================
# CPAT -- coding potential
# ==============================================================================

#' Run CPAT Coding Potential Prediction
#'
#' @description Runs CPAT (Coding-Potential Assessment Tool) against the
#' nucleotide FASTA from \code{APOLLO_extract_isoform_sequences()}, inside an
#' isolated \pkg{basilisk} Python environment (\code{.gaia_cpat_env}). No
#' manual Python/conda setup is required beyond having internet access on
#' first use (basilisk builds the environment then).
#'
#' @param switch_result An \code{artemis_isoform_switch} object, normally one
#'   that has already been through \code{APOLLO_extract_isoform_sequences()}.
#' @param hexamer_file Character. Path to a CPAT hexamer frequency table
#'   (species-specific; CPAT publishes prebuilt tables for human/mouse/fly/
#'   zebrafish, or generate one with CPAT's \code{make_hexamer_tab.py}).
#' @param logit_model_file Character. Path to the matching CPAT logit model
#'   file (\code{.RData}).
#' @param output_dir Character. Directory for CPAT's raw output and the
#'   reformatted result file (created if absent).
#' @param nt_fasta Character or \code{NULL}. Nucleotide FASTA path. Default
#'   \code{NULL} uses \code{switch_result$sequence_files$nt_fasta}.
#' @param verbose Logical. Print progress and CPAT's own console output.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly, the path to a reformatted result file compatible with
#'   \code{ARTEMIS_isoform_switch_consequences(cpat_file = ...)}.
#'
#' @details
#' CPAT's own output columns (\code{ID, mRNA, ORF_strand, ORF_frame,
#' ORF_start, ORF_end, ORF, Fickett, Hexamer, Coding_prob}) do not match what
#' \code{IsoformSwitchAnalyzeR::analyzeCPAT()} expects to parse (a 5- or
#' 8-column legacy format) -- confirmed by inspecting \code{analyzeCPAT()}'s
#' source directly, and a known point of confusion in the community (e.g.
#' \url{https://www.biostars.org/p/9531584/}). This function reformats CPAT's
#' raw output into the exact column set/order \code{analyzeCPAT()} requires
#' before returning the path, so no manual reformatting is needed downstream.
#'
#' @export
APOLLO_run_cpat <- function(switch_result,
                             hexamer_file,
                             logit_model_file,
                             output_dir,
                             nt_fasta = NULL,
                             verbose  = TRUE) {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (is.null(nt_fasta)) nt_fasta <- switch_result$sequence_files$nt_fasta
  if (is.null(nt_fasta))
    stop("No nucleotide FASTA available. Run APOLLO_extract_isoform_sequences() ",
         "first, or supply nt_fasta directly.", call. = FALSE)
  if (!file.exists(nt_fasta))
    stop("nt_fasta not found: ", nt_fasta, call. = FALSE)
  if (!file.exists(hexamer_file))
    stop("hexamer_file not found: ", hexamer_file, call. = FALSE)
  if (!file.exists(logit_model_file))
    stop("logit_model_file not found: ", logit_model_file, call. = FALSE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_prefix <- file.path(output_dir, "cpat_out")
  cli_args   <- c("-g", nt_fasta, "-x", hexamer_file, "-d", logit_model_file,
                   "-o", out_prefix)

  if (verbose) cat("[APOLLO] Running CPAT (basilisk Python environment)...\n")

  start_time <- proc.time()
  .apollo_cpat_run(cli_args, verbose)
  elapsed <- (proc.time() - start_time)["elapsed"]

  raw_file <- paste0(out_prefix, ".ORF_prob.best.tsv")
  if (!file.exists(raw_file)) raw_file <- paste0(out_prefix, ".ORF_prob.tsv")
  if (!file.exists(raw_file))
    stop("CPAT completed but no output file was found at prefix: ", out_prefix,
         call. = FALSE)

  reformatted <- file.path(output_dir, "cpat_results_isas.tsv")
  .apollo_reformat_cpat(raw_file, reformatted)

  if (verbose) {
    cat(" CPAT completed in", round(elapsed, 1), "seconds.\n")
    cat("[APOLLO] ISAS-compatible result: ", reformatted, "\n\n", sep = "")
  }

  invisible(reformatted)
}


#' Execute CPAT inside the basilisk environment
#' @keywords internal
.apollo_cpat_run <- function(cli_args, verbose) {
  basilisk::basiliskRun(
    env = .gaia_cpat_env,
    fun = function(args, verbose) {
      python_bin <- reticulate::py_exe()
      cpat_bin   <- file.path(dirname(python_bin), "cpat")
      if (!file.exists(cpat_bin))
        stop("cpat executable not found in basilisk environment bin directory: ",
             dirname(python_bin), call. = FALSE)

      exit_code <- system2(
        command = cpat_bin,
        args    = args,
        stdout  = if (verbose) "" else FALSE,
        stderr  = if (verbose) "" else FALSE
      )
      if (exit_code != 0)
        stop("CPAT exited with code ", exit_code, ". Re-run with verbose = TRUE ",
             "to see CPAT's own output.", call. = FALSE)
      invisible(exit_code)
    },
    args = cli_args, verbose = verbose
  )
}


#' Reformat CPAT output to IsoformSwitchAnalyzeR's expected column layout
#' @keywords internal
.apollo_reformat_cpat <- function(raw_file, out_file) {
  raw <- utils::read.table(raw_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE)

  needed <- c("seq_ID", "mRNA", "ORF", "Fickett", "Hexamer", "Coding_prob")
  missing_cols <- setdiff(needed, colnames(raw))
  if (length(missing_cols) > 0)
    stop("CPAT output is missing expected column(s): ",
         paste(missing_cols, collapse = ", "),
         ". CPAT's output format may differ from the version this wrapper was ",
         "built against (CPAT 3.0.5) -- check ", raw_file, " and adjust ",
         ".apollo_reformat_cpat() if needed.", call. = FALSE)

  # analyzeCPAT()'s 5-column branch requires exactly 5 data columns (checked
  # with `==`, so order-sensitive) and takes the transcript id from the row
  # names, not from a 6th "id" column -- a 6-column file falls through to
  # its generic "problem with the CPAT result file" error instead. Use
  # seq_ID (the original transcript id from the input FASTA header), not
  # ID (CPAT's own per-ORF id, suffixed "_ORF_1") -- ID never matches the
  # switchAnalyzeRlist's isoform_id, which triggers analyzeCPAT()'s
  # "transcript ids ... does not match" error.
  out <- data.frame(
    mRNA_size     = raw$mRNA,
    ORF_size      = raw$ORF,
    Fickett_score = raw$Fickett,
    Hexamer_score = raw$Hexamer,
    coding_prob   = raw$Coding_prob,
    stringsAsFactors = FALSE
  )
  rownames(out) <- raw$seq_ID

  utils::write.table(out, out_file, sep = "\t", row.names = TRUE, quote = FALSE)
  invisible(out_file)
}


# ==============================================================================
# SignalP 6 -- signal peptides
# ==============================================================================

#' Run SignalP 6 Signal Peptide Prediction
#'
#' @description Runs SignalP 6 against the amino acid FASTA from
#' \code{APOLLO_extract_isoform_sequences()}, inside an isolated
#' \pkg{basilisk} Python environment (\code{.gaia_signalp_env}).
#'
#' @param switch_result An \code{artemis_isoform_switch} object, normally one
#'   that has already been through \code{APOLLO_extract_isoform_sequences()}.
#' @param output_dir Character. Directory for SignalP's output (created if
#'   absent). SignalP writes \code{prediction_results.txt} here.
#' @param aa_fasta Character or \code{NULL}. Amino acid FASTA path. Default
#'   \code{NULL} uses \code{switch_result$sequence_files$aa_fasta}.
#' @param organism Character. \code{"euk"} (eukaryote, default) or
#'   \code{"other"}, passed to SignalP's \code{--organism}.
#' @param mode Character. SignalP 6 prediction mode: \code{"fast"} (default),
#'   \code{"slow"}, or \code{"slow-sequential"} (higher accuracy, much slower).
#' @param verbose Logical. Print progress and SignalP's own console output.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly, the path to \code{prediction_results.txt}, directly
#'   compatible with \code{ARTEMIS_isoform_switch_consequences(signalp_file = ...)}
#'   -- confirmed by inspecting \code{analyzeSignalP()}'s source, which
#'   natively parses SignalP 6's output format (detected via a
#'   \code{"SignalP-6"} header line), no reformatting needed.
#'
#' @details
#' SignalP 6 is license-gated (DTU Health Tech) and not on public PyPI, so
#' \code{.gaia_signalp_env} only provisions a base Python interpreter --
#' SignalP 6 itself must be installed once, manually, into that environment
#' (download the licensed wheel, then install with the environment's own
#' pip; see \code{basilisk::obtainEnvironmentPath(.gaia_signalp_env)} to
#' locate it). Follow DTU's installation instructions in full, including the
#' separate model-weights copy step. This function checks for the installed
#' \code{signalp6} executable and stops with setup guidance if it's missing,
#' rather than attempting to install it automatically.
#'
#' @export
APOLLO_run_signalp <- function(switch_result,
                                output_dir,
                                aa_fasta = NULL,
                                organism = c("euk", "other"),
                                mode     = c("fast", "slow", "slow-sequential"),
                                verbose  = TRUE) {

  organism <- match.arg(organism)
  mode     <- match.arg(mode)

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (is.null(aa_fasta)) aa_fasta <- switch_result$sequence_files$aa_fasta
  if (is.null(aa_fasta))
    stop("No amino acid FASTA available. Run APOLLO_extract_isoform_sequences() ",
         "first, or supply aa_fasta directly.", call. = FALSE)
  if (!file.exists(aa_fasta))
    stop("aa_fasta not found: ", aa_fasta, call. = FALSE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (verbose) cat("[APOLLO] Running SignalP 6 (basilisk Python environment)...\n")

  start_time <- proc.time()
  .apollo_signalp_run(aa_fasta, output_dir, organism, mode, verbose)
  elapsed <- (proc.time() - start_time)["elapsed"]

  result_file <- file.path(output_dir, "prediction_results.txt")
  if (!file.exists(result_file))
    stop("SignalP completed but prediction_results.txt was not found in: ",
         output_dir, call. = FALSE)

  if (verbose) {
    cat(" SignalP completed in", round(elapsed, 1), "seconds.\n")
    cat("[APOLLO] Result: ", result_file, "\n\n", sep = "")
  }

  invisible(result_file)
}


#' Execute SignalP 6 inside the basilisk environment
#' @keywords internal
.apollo_signalp_run <- function(fasta, output_dir, organism, mode, verbose) {
  basilisk::basiliskRun(
    env = .gaia_signalp_env,
    fun = function(fasta, output_dir, organism, mode, verbose) {
      python_bin  <- reticulate::py_exe()
      signalp_bin <- file.path(dirname(python_bin), "signalp6")
      if (!file.exists(signalp_bin))
        stop("signalp6 executable not found in the basilisk environment (",
             dirname(python_bin), "). SignalP 6 is license-gated and cannot be ",
             "auto-installed -- download it from DTU Health Tech ",
             "(https://services.healthtech.dtu.dk/) and pip-install the wheel ",
             "into this environment's own pip (see ",
             "basilisk::obtainEnvironmentPath(.gaia_signalp_env)) following ",
             "DTU's installation instructions in full, including the separate ",
             "model-weights step, before calling this function.", call. = FALSE)

      exit_code <- system2(
        command = signalp_bin,
        args    = c("--fastafile", fasta, "--output_dir", output_dir,
                     "--format", "txt", "--organism", organism, "--mode", mode),
        stdout  = if (verbose) "" else FALSE,
        stderr  = if (verbose) "" else FALSE
      )
      if (exit_code != 0)
        stop("SignalP 6 exited with code ", exit_code, ". Re-run with ",
             "verbose = TRUE to see SignalP's own output.", call. = FALSE)
      invisible(exit_code)
    },
    fasta = fasta, output_dir = output_dir, organism = organism, mode = mode,
    verbose = verbose
  )
}


# ==============================================================================
# Pfam -- protein domains (via pfam_scan.pl / HMMER)
# ==============================================================================

#' Run Pfam Protein Domain Annotation
#'
#' @description Runs \code{pfam_scan.pl} (a Perl wrapper around HMMER's
#' \code{hmmscan}) against the amino acid FASTA from
#' \code{APOLLO_extract_isoform_sequences()}, inside an isolated
#' \pkg{basilisk} environment (\code{.gaia_pfam_env}). \code{pfam_scan.pl} is
#' packaged on bioconda as \code{pfam_scan}, which pulls in HMMER and its own
#' Perl module dependencies as part of the conda recipe -- so no manual
#' PfamScan.tar.gz/CPAN setup is needed.
#'
#' @param switch_result An \code{artemis_isoform_switch} object, normally one
#'   that has already been through \code{APOLLO_extract_isoform_sequences()}.
#' @param pfam_dir Character. Path to a directory containing the Pfam-A HMM
#'   database, pre-pressed with \code{hmmpress} (i.e. containing
#'   \code{Pfam-A.hmm}, \code{Pfam-A.hmm.dat}, and the \code{hmmpress}-
#'   generated index files).
#' @param output_dir Character. Directory for Pfam's output (created if
#'   absent).
#' @param aa_fasta Character or \code{NULL}. Amino acid FASTA path. Default
#'   \code{NULL} uses \code{switch_result$sequence_files$aa_fasta}.
#' @param nproc Integer. Number of CPUs for \code{pfam_scan.pl}'s
#'   \code{-cpu} option. Default: \code{4}.
#' @param verbose Logical. Print progress and \code{pfam_scan.pl}'s own
#'   console output. Default: \code{TRUE}.
#'
#' @return Invisibly, the path to the Pfam result file, directly compatible
#'   with \code{ARTEMIS_isoform_switch_consequences(pfam_file = ...)}.
#'
#' @details
#' \code{IsoformSwitchAnalyzeR::analyzePFAM()} expects \code{pfam_scan.pl}'s
#' specific output format (confirmed against the package's own bundled
#' example result file), not raw \code{hmmscan --domtblout} output -- this is
#' why the wrapper runs \code{pfam_scan.pl} rather than calling
#' \code{hmmscan} directly, even though the latter would be sufficient to
#' scan against Pfam-A on its own.
#'
#' @export
APOLLO_run_pfam <- function(switch_result,
                             pfam_dir,
                             output_dir,
                             aa_fasta = NULL,
                             nproc    = 4,
                             verbose  = TRUE) {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (is.null(aa_fasta)) aa_fasta <- switch_result$sequence_files$aa_fasta
  if (is.null(aa_fasta))
    stop("No amino acid FASTA available. Run APOLLO_extract_isoform_sequences() ",
         "first, or supply aa_fasta directly.", call. = FALSE)
  if (!file.exists(aa_fasta))
    stop("aa_fasta not found: ", aa_fasta, call. = FALSE)
  if (!dir.exists(pfam_dir))
    stop("pfam_dir not found: ", pfam_dir, call. = FALSE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_file <- file.path(output_dir, "pfam_results.txt")
  cli_args <- c("-fasta", aa_fasta, "-dir", pfam_dir, "-outfile", out_file,
                "-cpu", as.character(nproc))

  if (verbose) {
    cat("[APOLLO] Running pfam_scan.pl (basilisk environment)...\n")
    cat("    AA fasta: ", aa_fasta, "\n", sep = "")
    cat("    Pfam DB : ", pfam_dir, "\n", sep = "")
  }

  start_time <- proc.time()
  .apollo_pfam_run(cli_args, verbose)
  elapsed <- (proc.time() - start_time)["elapsed"]

  if (!file.exists(out_file))
    stop("pfam_scan.pl completed but output file was not created: ", out_file,
         call. = FALSE)

  if (verbose) {
    cat(" pfam_scan.pl completed in", round(elapsed, 1), "seconds.\n")
    cat("[APOLLO] Result: ", out_file, "\n\n", sep = "")
  }

  invisible(out_file)
}


#' Execute pfam_scan.pl inside the basilisk environment
#' @keywords internal
.apollo_pfam_run <- function(cli_args, verbose) {
  basilisk::basiliskRun(
    env = .gaia_pfam_env,
    fun = function(args, verbose) {
      python_bin    <- reticulate::py_exe()
      pfam_scan_bin <- file.path(dirname(python_bin), "pfam_scan.pl")
      if (!file.exists(pfam_scan_bin))
        stop("pfam_scan.pl executable not found in basilisk environment bin ",
             "directory: ", dirname(python_bin), call. = FALSE)

      exit_code <- system2(
        command = pfam_scan_bin,
        args    = args,
        stdout  = if (verbose) "" else FALSE,
        stderr  = if (verbose) "" else FALSE
      )
      if (exit_code != 0)
        stop("pfam_scan.pl exited with code ", exit_code, ". Re-run with ",
             "verbose = TRUE to see its own output.", call. = FALSE)
      invisible(exit_code)
    },
    args = cli_args, verbose = verbose
  )
}
