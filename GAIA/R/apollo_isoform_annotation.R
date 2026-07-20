#' Apollo - Isoform Switch Consequence Tool Wrappers
#'
#' @description Runs external protein/transcript annotation tools (CPAT,
#' SignalP, Pfam via \code{pfam_scan.pl}, DeepTMHMM) against isoforms from an
#' \code{artemis_isoform_switch} object, producing result files consumable by
#' \code{ARTEMIS_isoform_switch_consequences()}. All four run inside
#' isolated \pkg{basilisk} environments (see \code{GAIA/R/basilisk.R}) --
#' CPAT and SignalP are Python tools; Pfam is packaged via bioconda's
#' \code{pfam_scan} recipe, which basilisk installs the same way regardless
#' of it being a Perl/native-binary tool rather than a Python one; DeepTMHMM
#' is an academically licensed PyTorch model run from a user-supplied local
#' checkout rather than an installed executable.


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


# ==============================================================================
# DeepTMHMM -- topology prediction
# ==============================================================================

#' Run DeepTMHMM Topology Prediction
#'
#' @description Runs DeepTMHMM (transmembrane helix / signal peptide /
#' cell-membrane topology prediction) against the amino acid FASTA from
#' \code{APOLLO_extract_isoform_sequences()}, inside an isolated
#' \pkg{basilisk} Python environment (\code{.gaia_deeptmhmm_env}).
#'
#' @param switch_result An \code{artemis_isoform_switch} object, normally one
#'   that has already been through \code{APOLLO_extract_isoform_sequences()}.
#' @param deeptmhmm_dir Character. Path to the locally installed, academically
#'   licensed DeepTMHMM package directory -- a flat checkout containing
#'   \code{predict.py}, \code{utils.py}, the \code{experiments/} folder, and
#'   the five \code{deeptmhmm_cv_*.model} weight files.
#' @param output_dir Character. Directory \code{predict.py} writes its output
#'   to. Must \strong{not} already exist -- \code{predict.py} itself checks
#'   for this and exits with an error if \code{output_dir} is already
#'   present (unlike \code{APOLLO_run_cpat()}/\code{APOLLO_run_signalp()}/
#'   \code{APOLLO_run_pfam()}, this function does not pre-create it).
#' @param aa_fasta Character or \code{NULL}. Amino acid FASTA path. Default
#'   \code{NULL} uses \code{switch_result$sequence_files$aa_fasta}.
#' @param verbose Logical. Print progress and DeepTMHMM's own console output.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly, the path to the \code{TMRs.gff3} result file, directly
#'   compatible with
#'   \code{IsoformSwitchAnalyzeR::analyzeDeepTMHMM(pathToDeepTMHMMresultFile = ...)}
#'   / \code{ARTEMIS_isoform_switch_consequences(deeptmhmm_file = ...)}.
#'
#' @details
#' DeepTMHMM is academically licensed (DTU Health Tech / BioLib) and not on
#' public PyPI/bioconda channels, so \code{.gaia_deeptmhmm_env} only
#' provisions a base Python 3.8 interpreter -- \code{predict.py}'s own
#' dependencies (PyTorch, per its bundled \code{requirements.txt}) must be
#' installed once, manually, into that environment following the license
#' holder's own \code{README.txt} before calling this function. This
#' function checks for \code{predict.py} in \code{deeptmhmm_dir} and stops
#' with setup guidance if it's missing, rather than attempting to install it
#' automatically (same pattern as \code{APOLLO_run_signalp()}).
#'
#' Unlike CPAT/SignalP/Pfam, \code{predict.py} is invoked from within its own
#' package directory rather than as an installed executable on \code{PATH}:
#' it loads its five \code{deeptmhmm_cv_*.model} weight files and its own
#' \code{utils}/\code{experiments.tmhmm3.tm_util} modules via bare relative
#' paths/imports (confirmed by reading \code{predict.py} directly), so this
#' function temporarily changes the working directory to
#' \code{deeptmhmm_dir} for the duration of the call (restored via
#' \code{on.exit()}, including on error).
#'
#' \code{predict.py}'s own \code{requirements.txt} pins
#' \code{torch==1.5.0+cu92} (a CUDA 9.2-era wheel) -- confirm this actually
#' resolves on your hardware/driver stack before assuming it as fixed; a
#' newer CPU-only or different CUDA-version torch build may be needed
#' instead. GPU vs CPU execution is auto-detected by \code{predict.py} itself
#' (\code{torch.cuda.is_available()}) and is not controlled by this wrapper.
#'
#' @export
APOLLO_run_deeptmhmm <- function(switch_result,
                                  deeptmhmm_dir,
                                  output_dir,
                                  aa_fasta = NULL,
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

  if (!dir.exists(deeptmhmm_dir))
    stop("deeptmhmm_dir not found: ", deeptmhmm_dir, call. = FALSE)

  predict_script <- file.path(deeptmhmm_dir, "predict.py")
  if (!file.exists(predict_script))
    stop("predict.py not found in deeptmhmm_dir: ", deeptmhmm_dir, ". DeepTMHMM ",
         "is academically licensed and not auto-installable -- request a copy ",
         "from DTU Health Tech/BioLib (licensing@biolib.com, see ",
         "https://biolib.com/DTU/DeepTMHMM/), extract it, and install its ",
         "PyTorch dependencies into .gaia_deeptmhmm_env following its own ",
         "README.txt before calling this function.", call. = FALSE)

  if (dir.exists(output_dir))
    stop("output_dir already exists: ", output_dir, ". predict.py refuses to ",
         "run if its output directory is already present -- remove it or ",
         "choose a new path.", call. = FALSE)

  aa_fasta_abs <- normalizePath(aa_fasta, mustWork = TRUE)
  # predict.py creates output_dir itself -- it doesn't exist yet, so it can't
  # be normalizePath(..., mustWork = TRUE) directly. On a non-existent path,
  # normalizePath(mustWork = FALSE) is unreliable -- on Linux it just returns
  # relative paths unchanged instead of resolving them. Since .apollo_deeptmhmm_run()
  # setwd()s into deeptmhmm_dir before invoking predict.py, a still-relative
  # output_dir_abs would resolve inside deeptmhmm_dir instead of the caller's
  # intended location. Resolve via the parent directory instead, which does
  # exist (or can be created without touching output_dir itself).
  output_parent <- dirname(output_dir)
  if (!dir.exists(output_parent)) dir.create(output_parent, recursive = TRUE)
  output_dir_abs <- file.path(normalizePath(output_parent, mustWork = TRUE),
                               basename(output_dir))

  if (verbose) {
    cat("[APOLLO] Running DeepTMHMM (basilisk Python environment)...\n")
    cat("    AA fasta : ", aa_fasta_abs, "\n", sep = "")
    cat("    DeepTMHMM: ", deeptmhmm_dir, "\n", sep = "")
  }

  start_time <- proc.time()
  .apollo_deeptmhmm_run(deeptmhmm_dir, aa_fasta_abs, output_dir_abs, verbose)
  elapsed <- (proc.time() - start_time)["elapsed"]

  gff3_file <- file.path(output_dir, "TMRs.gff3")
  if (!file.exists(gff3_file))
    stop("DeepTMHMM completed but TMRs.gff3 was not found in: ", output_dir,
         call. = FALSE)

  if (verbose) {
    cat(" DeepTMHMM completed in", round(elapsed, 1), "seconds.\n")
    cat("[APOLLO] Result: ", gff3_file, "\n\n", sep = "")
  }

  invisible(gff3_file)
}


#' Execute DeepTMHMM's predict.py inside the basilisk environment
#' @keywords internal
.apollo_deeptmhmm_run <- function(deeptmhmm_dir, fasta, output_dir, verbose) {
  basilisk::basiliskRun(
    env = .gaia_deeptmhmm_env,
    fun = function(deeptmhmm_dir, fasta, output_dir, verbose) {
      python_bin <- reticulate::py_exe()

      # predict.py loads its five deeptmhmm_cv_*.model weight files and its
      # own utils/experiments.tmhmm3.tm_util modules via bare relative
      # paths/imports -- it must be run with this directory as cwd.
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(deeptmhmm_dir)

      exit_code <- system2(
        command = python_bin,
        args    = c("predict.py", "--fasta", fasta, "--output-dir", output_dir),
        stdout  = if (verbose) "" else FALSE,
        stderr  = if (verbose) "" else FALSE
      )
      if (exit_code != 0)
        stop("DeepTMHMM (predict.py) exited with code ", exit_code, ". Re-run ",
             "with verbose = TRUE to see its own output.", call. = FALSE)
      invisible(exit_code)
    },
    deeptmhmm_dir = deeptmhmm_dir, fasta = fasta, output_dir = output_dir,
    verbose = verbose
  )
}
