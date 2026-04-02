# GAIA/Apollo/motif_enrichment.R
# HOMER motif enrichment analysis for ATAC-seq peaks
#
# Wraps HOMER's findMotifsGenome.pl for known motif enrichment and optional
# de novo motif discovery. Supports single and batch analysis of BED files.
#
# Requires HOMER to be installed and on PATH.
# See: http://homer.ucsd.edu/homer/


# ==============================================================================
# Peak preprocessing
# ==============================================================================

#' Extract summit-centered BED files from narrowPeak files
#'
#' Reads MACS2/MACS3 narrowPeak files and writes fixed-width BED files centered
#' on each peak summit (column 10: 0-based offset from peak start). The output
#' is a named vector of BED paths that feeds directly into
#' \code{APOLLO_homer_motif_enrichment_batch()}.
#'
#' @param peak_files Named character vector of narrowPeak file paths. Names are
#'   used as set identifiers and as the base of output BED filenames.
#' @param output_dir Directory to save summit BED files. Created if it doesn't
#'   exist.
#' @param half_window Integer. Half-width (bp) around each summit. Total window
#'   = \code{2 * half_window}. Default: 100 (200 bp total).
#' @param chr_mapping Optional chromosome name conversion applied to the chr
#'   column before writing the BED file. Useful when peak files use UCSC-style
#'   names (chr1, chr2) but the HOMER genome was indexed with NCBI accessions
#'   (e.g., NC_060925.1 for T2T). Accepts:
#'   \itemize{
#'     \item A genome name string (e.g., \code{"T2T"}) — uses the built-in
#'       mapping from \code{APOLLO_get_chr_mapping()}, automatically reversed
#'       to UCSC → NCBI direction.
#'     \item A named character vector where names are the source chr names and
#'       values are the target chr names.
#'     \item \code{NULL} (default) — no conversion.
#'   }
#'   Peaks whose chromosome is not found in the mapping are dropped with a
#'   warning.
#' @param skip_missing Logical. If TRUE, missing peak files are skipped with a
#'   warning rather than causing an error. Default: TRUE.
#' @param verbose Logical. Print per-set progress. Default: TRUE.
#'
#' @return Named character vector of output BED file paths (one per input set).
#'   Silently drops entries for files that were skipped or contained no valid
#'   summits.
#'
#' @details
#' NarrowPeak column 10 is the 0-based distance from the peak start to the
#' summit. Rows where this value is \code{-1} (as in broadPeak-style output with
#' no called summit) are skipped. Coordinates that would extend below 0 are
#' clamped to 0; no upper-boundary clamping is applied (chromosome sizes are
#' not required).
#'
#' The output BED is 6-column (chr, start, end, name, score, strand). Name and
#' score are taken from the narrowPeak where available.
#'
#' @examples
#' \dontrun{
#' peak_files <- c(
#'   KO1a = "data/ATACseq/KO1a_sorted_peaks.narrowPeak",
#'   WT1a = "data/ATACseq/WT1a_sorted_peaks.narrowPeak"
#' )
#'
#' # Standard hg38
#' summit_beds <- APOLLO_extract_summits(peak_files, "results/summits/")
#'
#' # T2T: convert chr1/chr2/... to NC_060925.1/NC_060926.1/...
#' summit_beds <- APOLLO_extract_summits(peak_files, "results/summits/", chr_mapping = "T2T")
#'
#' # Feed directly into HOMER batch
#' batch <- APOLLO_homer_motif_enrichment_batch(summit_beds, "T2T", "results/motifs/")
#' }
#' @export
APOLLO_extract_summits <- function(peak_files,
                                    output_dir,
                                    half_window = 100L,
                                    chr_mapping = NULL,
                                    skip_missing = TRUE,
                                    verbose = TRUE) {

  if (is.null(names(peak_files)) || any(names(peak_files) == "")) {
    stop("peak_files must be a named character vector. Names are used as set identifiers.",
         call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("[APOLLO] Created directory:", output_dir, "\n")
  }

  half_window <- as.integer(half_window)

  # Resolve chr_mapping to a source -> target named vector
  chr_map_vec <- NULL
  if (!is.null(chr_mapping)) {
    if (is.character(chr_mapping) && length(chr_mapping) == 1 &&
        (is.null(names(chr_mapping)) || names(chr_mapping) == "")) {
      # Built-in genome name (e.g., "T2T"): APOLLO_get_chr_mapping returns
      # NCBI -> UCSC; reverse it to get UCSC -> NCBI for BED chr renaming
      ncbi_to_ucsc <- APOLLO_get_chr_mapping(chr_mapping)
      chr_map_vec  <- stats::setNames(names(ncbi_to_ucsc), ncbi_to_ucsc)
      if (verbose) cat("    Chromosome mapping: UCSC -> NCBI (", chr_mapping, ")\n", sep = "")
    } else if (is.character(chr_mapping) && !is.null(names(chr_mapping))) {
      chr_map_vec <- chr_mapping
    } else {
      stop("chr_mapping must be a genome name string (e.g., 'T2T') or a named character vector.",
           call. = FALSE)
    }
  }

  out_paths <- stats::setNames(rep(NA_character_, length(peak_files)), names(peak_files))

  for (nm in names(peak_files)) {
    peak_file <- peak_files[[nm]]

    # Check file exists
    if (!file.exists(peak_file)) {
      if (skip_missing) {
        if (verbose) cat("    Skipping (file not found):", peak_file, "\n")
        next
      } else {
        stop("Peak file not found: ", peak_file, call. = FALSE)
      }
    }

    # Read narrowPeak (tab-separated, no header)
    np <- tryCatch(
      utils::read.table(peak_file, header = FALSE, stringsAsFactors = FALSE,
                        fill = TRUE, comment.char = "#"),
      error = function(e) {
        stop("Failed to read narrowPeak '", peak_file, "': ", e$message, call. = FALSE)
      }
    )

    if (ncol(np) < 10) {
      stop("File does not appear to be narrowPeak format (expected >= 10 columns): ",
           peak_file, call. = FALSE)
    }

    # Drop rows with no summit (col10 = -1, used in broadPeak-style output)
    valid <- !is.na(np[[10]]) & np[[10]] >= 0L
    if (any(!valid)) {
      if (verbose) {
        cat("    ", nm, ": dropping", sum(!valid),
            "peak(s) with no summit (col10 = -1)\n", sep = "")
      }
      np <- np[valid, , drop = FALSE]
    }

    if (nrow(np) == 0) {
      warning("No peaks with a valid summit in: ", peak_file)
      next
    }

    # Summit position (0-based): peak start + col10 offset
    summit    <- np[[2]] + np[[10]]
    start_new <- pmax(0L, as.integer(summit) - half_window)
    end_new   <- as.integer(summit) + half_window

    # 6-column BED
    has_name  <- ncol(np) >= 4 && !all(np[[4]] == ".")
    has_score <- ncol(np) >= 5
    has_strand <- ncol(np) >= 6 && !all(np[[6]] == ".")

    bed <- data.frame(
      chr    = np[[1]],
      start  = start_new,
      end    = end_new,
      name   = if (has_name) np[[4]] else paste0(nm, "_peak_", seq_len(nrow(np))),
      score  = if (has_score) np[[5]] else 0L,
      strand = if (has_strand) np[[6]] else ".",
      stringsAsFactors = FALSE
    )

    # Apply chromosome name conversion (e.g., chr1 -> NC_060925.1 for T2T)
    if (!is.null(chr_map_vec)) {
      mapped_chr <- chr_map_vec[bed$chr]
      n_unmapped <- sum(is.na(mapped_chr))
      if (n_unmapped > 0) {
        warning(nm, ": ", n_unmapped, " peak(s) on chromosomes not in mapping — dropped.")
        bed <- bed[!is.na(mapped_chr), , drop = FALSE]
        mapped_chr <- mapped_chr[!is.na(mapped_chr)]
      }
      bed$chr <- mapped_chr
    }

    out_file <- file.path(output_dir, paste0(nm, "_summits.bed"))
    utils::write.table(bed, out_file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)

    out_paths[[nm]] <- out_file

    if (verbose) {
      cat("    ", nm, ": ", nrow(bed), " summit regions (", 2L * half_window,
          " bp) -> ", basename(out_file), "\n", sep = "")
    }
  }

  # Drop entries for sets that were skipped
  out_paths <- out_paths[!is.na(out_paths)]

  if (verbose) {
    cat("    Summit BED files written:", length(out_paths), "/",
        length(peak_files), "\n")
  }

  return(out_paths)
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Check HOMER installation
#' @return Path to findMotifsGenome.pl or stops with error
#' @keywords internal
.homer_check_installation <- function() {
  homer_path <- Sys.which("findMotifsGenome.pl")
  if (homer_path == "") {
    stop(
      "HOMER not found on PATH.\n",
      "  Install HOMER: http://homer.ucsd.edu/homer/introduction/install.html\n",
      "  Ensure findMotifsGenome.pl is on your PATH.",
      call. = FALSE
    )
  }
  return(homer_path)
}


#' Build HOMER findMotifsGenome.pl command arguments
#' @keywords internal
.homer_build_command <- function(bed_file, genome, output_dir,
                                 size, mask, denovo, denovo_n,
                                 bg, nproc, extra_args) {

  args <- c(bed_file, genome, output_dir)

  # Region size
  args <- c(args, "-size", as.character(size))

  # Mask repeats
  if (mask) args <- c(args, "-mask")

  # De novo control
  if (!denovo) {
    args <- c(args, "-nomotif")
  } else {
    args <- c(args, "-S", as.character(denovo_n))
  }

  # Background
  if (!is.null(bg)) {
    if (!file.exists(bg)) stop("Background BED file not found: ", bg, call. = FALSE)
    args <- c(args, "-bg", bg)
  }

  # Processors
  args <- c(args, "-p", as.character(nproc))

  # Extra arguments
  if (!is.null(extra_args)) {
    args <- c(args, extra_args)
  }

  return(args)
}


#' Parse HOMER knownResults.txt
#'
#' @param known_file Path to knownResults.txt
#' @return data.frame with cleaned columns
#' @keywords internal
.homer_parse_known <- function(known_file) {

  if (!file.exists(known_file)) return(NULL)

  # Read raw - HOMER uses tab-separated with a header line
  raw <- tryCatch(
    utils::read.delim(known_file, header = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = ""),
    error = function(e) {
      warning("Failed to parse knownResults.txt: ", e$message)
      return(NULL)
    }
  )

  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # HOMER column names (may vary slightly):
  # "Motif Name" "Consensus" "P-value" "Log P-value" "q-value (Benjamini)"
  # "# of Target Sequences with Motif(of X)" "% of Target Sequences with Motif"
  # "# of Background Sequences with Motif(of Y)" "% of Background Sequences with Motif"
  colnames(raw) <- c("motif_name", "consensus", "p_value", "log_p_value",
                      "q_value", "n_target_raw", "pct_target_raw",
                      "n_background_raw", "pct_background_raw")[seq_len(ncol(raw))]

  # Parse count columns: "123.0(of 456)" -> n_target=123, total_target=456
  .parse_count <- function(x) {
    # Handle formats: "123.0(of 456)" or just "123"
    n <- as.numeric(sub("\\(.*", "", x))
    total <- as.numeric(sub(".*of\\s+", "", sub("\\).*", "", x)))
    list(n = n, total = total)
  }

  target_parsed <- .parse_count(raw$n_target_raw)
  bg_parsed <- .parse_count(raw$n_background_raw)

  # Parse percentage columns: "45.67%" -> 45.67
  .parse_pct <- function(x) {
    as.numeric(sub("%", "", x))
  }

  result <- data.frame(
    motif_name   = raw$motif_name,
    consensus    = raw$consensus,
    p_value      = as.numeric(raw$p_value),
    log_p_value  = as.numeric(raw$log_p_value),
    q_value      = as.numeric(raw$q_value),
    n_target     = target_parsed$n,
    total_target = target_parsed$total,
    pct_target   = .parse_pct(raw$pct_target_raw),
    n_background = bg_parsed$n,
    total_background = bg_parsed$total,
    pct_background   = .parse_pct(raw$pct_background_raw),
    stringsAsFactors = FALSE
  )

  # Clean motif names: extract TF family from complex HOMER names
  # e.g., "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer" -> keep full name but add clean column
  result$motif_family <- sub("/.*", "", result$motif_name)

  return(result)
}


#' Parse HOMER de novo motif file
#'
#' @param motif_file Path to homerMotifs.all.motifs
#' @return data.frame with de novo motif summary
#' @keywords internal
.homer_parse_denovo <- function(motif_file) {

  if (!file.exists(motif_file)) return(NULL)

  lines <- readLines(motif_file, warn = FALSE)
  header_lines <- lines[grepl("^>", lines)]

  if (length(header_lines) == 0) return(NULL)

  # Parse each header line
  # Format: >consensus \t name,BestGuess:match \t log_odds \t log_p_value \t 0 \t T:count(pct%),B:count(pct%),P:pvalue
  results <- lapply(header_lines, function(h) {
    # Remove leading >
    h <- sub("^>", "", h)
    fields <- strsplit(h, "\t")[[1]]

    consensus <- fields[1]

    # Parse name/best guess from field 2
    info <- fields[2]
    best_match <- ""
    if (grepl("BestGuess:", info)) {
      best_match <- sub(".*BestGuess:", "", info)
      best_match <- sub("\\(.*$", "", best_match)  # Remove trailing score
    }

    log_odds <- suppressWarnings(as.numeric(fields[3]))
    log_p_value <- suppressWarnings(as.numeric(fields[4]))

    # Parse T:count(pct%),B:count(pct%),P:pvalue from last field
    stats_field <- fields[length(fields)]
    n_target <- NA_real_
    pct_target <- NA_real_
    n_background <- NA_real_
    pct_background <- NA_real_
    p_value <- NA_real_

    if (grepl("T:", stats_field)) {
      t_match <- regmatches(stats_field, regexpr("T:[0-9.]+\\([0-9.]+%\\)", stats_field))
      if (length(t_match) > 0) {
        n_target <- as.numeric(sub("T:([0-9.]+)\\(.*", "\\1", t_match))
        pct_target <- as.numeric(sub(".*\\(([0-9.]+)%\\)", "\\1", t_match))
      }
      b_match <- regmatches(stats_field, regexpr("B:[0-9.]+\\([0-9.]+%\\)", stats_field))
      if (length(b_match) > 0) {
        n_background <- as.numeric(sub("B:([0-9.]+)\\(.*", "\\1", b_match))
        pct_background <- as.numeric(sub(".*\\(([0-9.]+)%\\)", "\\1", b_match))
      }
      p_match <- regmatches(stats_field, regexpr("P:[0-9eE.+-]+", stats_field))
      if (length(p_match) > 0) {
        p_value <- as.numeric(sub("P:", "", p_match))
      }
    }

    data.frame(
      consensus      = consensus,
      best_match     = best_match,
      log_odds       = log_odds,
      log_p_value    = log_p_value,
      p_value        = p_value,
      n_target       = n_target,
      pct_target     = pct_target,
      n_background   = n_background,
      pct_background = pct_background,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


# ==============================================================================
# Public functions
# ==============================================================================

#' Run HOMER motif enrichment on a single BED file
#'
#' Wraps findMotifsGenome.pl to perform known motif enrichment and optional
#' de novo motif discovery on a set of genomic regions.
#'
#' @param bed_file Path to BED file with peak coordinates.
#' @param genome HOMER genome string (e.g., "hg38", "mm10", "hg19").
#' @param output_dir Directory for HOMER output. Created if it doesn't exist.
#' @param size Region size around peak center for motif scanning.
#'   Default 200. Use "given" to use actual peak widths from the BED file.
#' @param mask Logical. Mask repeat sequences. Default TRUE (recommended for ATAC-seq).
#' @param denovo Logical. Enable de novo motif discovery. Default FALSE (known only).
#'   When FALSE, uses -nomotif flag for faster execution.
#' @param denovo_n Integer. Number of de novo motifs to find. Default 10. Maps to -S.
#' @param bg Optional path to background BED file for custom background comparison.
#'   E.g., all peaks as background when testing enrichment in differential subset.
#' @param nproc Number of processors. Default 4. Maps to -p.
#' @param extra_args Character vector of additional HOMER arguments passed directly.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return S3 object of class "homer_motif" with components:
#'   \describe{
#'     \item{known}{Data.frame of known motif enrichment results}
#'     \item{denovo}{Data.frame of de novo motifs (NULL if denovo=FALSE)}
#'     \item{output_dir}{Path to HOMER output directory}
#'     \item{metadata}{List with genome, parameters, run info}
#'   }
#'
#' @examples
#' \dontrun{
#' # Known motifs only (fast)
#' result <- APOLLO_homer_motif_enrichment(
#'   "peaks.bed", "hg38", "results/motifs/"
#' )
#'
#' # With custom background and de novo discovery
#' result <- APOLLO_homer_motif_enrichment(
#'   "diff_peaks.bed", "hg38", "results/motifs/",
#'   bg = "all_peaks.bed", denovo = TRUE
#' )
#'
#' }
#' @export
APOLLO_homer_motif_enrichment <- function(bed_file,
                                           genome,
                                           output_dir,
                                           size = 200,
                                           mask = TRUE,
                                           denovo = FALSE,
                                           denovo_n = 10,
                                           bg = NULL,
                                           nproc = 4,
                                           extra_args = NULL,
                                           verbose = TRUE) {

  # --- Validate inputs ---
  homer_path <- .homer_check_installation()
  if (verbose) cat("[APOLLO] HOMER found:", homer_path, "\n")

  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file, call. = FALSE)
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Build and run command ---
  args <- .homer_build_command(
    bed_file  = bed_file,
    genome    = genome,
    output_dir = output_dir,
    size      = size,
    mask      = mask,
    denovo    = denovo,
    denovo_n  = denovo_n,
    bg        = bg,
    nproc     = nproc,
    extra_args = extra_args
  )

  if (verbose) {
    cat("[APOLLO] Running HOMER findMotifsGenome.pl...\n")
    cat("    BED file:", bed_file, "\n")
    cat("    Genome:", genome, "\n")
    cat("    Output:", output_dir, "\n")
    cat("    De novo:", if (denovo) paste("yes (n =", denovo_n, ")") else "no", "\n")
    if (!is.null(bg)) cat("    Background:", bg, "\n")
  }

  start_time <- proc.time()

  exit_code <- system2(
    command = "findMotifsGenome.pl",
    args    = args,
    stdout  = file.path(output_dir, "homer_stdout.log"),
    stderr  = file.path(output_dir, "homer_stderr.log")
  )

  elapsed <- (proc.time() - start_time)["elapsed"]

  if (exit_code != 0) {
    stderr_log <- file.path(output_dir, "homer_stderr.log")
    stderr_content <- ""
    if (file.exists(stderr_log)) {
      stderr_content <- paste(readLines(stderr_log, n = 20), collapse = "\n")
    }
    stop(
      "HOMER exited with code ", exit_code, ".\n",
      "Check log: ", stderr_log, "\n",
      if (nchar(stderr_content) > 0) paste0("Last output:\n", stderr_content),
      call. = FALSE
    )
  }

  if (verbose) cat(" HOMER completed in", round(elapsed, 1), "seconds.\n")

  # --- Parse results ---
  result <- APOLLO_load_homer_results(output_dir, verbose = verbose)

  # Enrich metadata with run parameters
  result$metadata$genome <- genome
  result$metadata$bed_file <- bed_file
  result$metadata$bg_file <- bg
  result$metadata$size <- size
  result$metadata$mask <- mask
  result$metadata$denovo <- denovo
  result$metadata$denovo_n <- denovo_n
  result$metadata$nproc <- nproc
  result$metadata$elapsed_seconds <- as.numeric(elapsed)

  return(result)
}


#' Run HOMER motif enrichment on multiple BED files
#'
#' Iterates over a named vector of BED file paths, running
#' APOLLO_homer_motif_enrichment() for each and combining results.
#'
#' @param bed_files Named character vector of BED file paths.
#'   Names are used as peak set identifiers.
#' @param genome HOMER genome string (e.g., "hg38").
#' @param output_dir Base output directory. Subdirectories are created per peak set.
#' @param ... Additional arguments passed to APOLLO_homer_motif_enrichment()
#'   (size, mask, denovo, denovo_n, bg, nproc, extra_args).
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return S3 object of class "homer_motif_batch" with components:
#'   \describe{
#'     \item{results}{Named list of homer_motif objects}
#'     \item{combined}{Data.frame of all known motif results with peak_set column}
#'     \item{summary}{Data.frame with per-set summary statistics}
#'     \item{metadata}{List with genome, parameters, bed_files}
#'   }
#'
#' @examples
#' \dontrun{
#' bed_files <- c(
#'   up   = "results/peaks_up.bed",
#'   down = "results/peaks_down.bed"
#' )
#' batch <- APOLLO_homer_motif_enrichment_batch(
#'   bed_files, "hg38", "results/motifs/"
#' )
#'
#' }
#' @export
APOLLO_homer_motif_enrichment_batch <- function(bed_files,
                                                 genome,
                                                 output_dir,
                                                 ...,
                                                 verbose = TRUE) {

  # Validate
  if (is.null(names(bed_files)) || any(names(bed_files) == "")) {
    stop("bed_files must be a named character vector. Names are used as set identifiers.",
         call. = FALSE)
  }

  set_names <- names(bed_files)
  n_sets <- length(bed_files)

  if (verbose) {
    cat("[APOLLO] Running HOMER batch motif enrichment\n")
    cat("     Peak sets:", n_sets, "\n")
    cat("    Sets:", paste(set_names, collapse = ", "), "\n")
  }

  # Run each set
  results <- list()
  for (i in seq_along(bed_files)) {
    set_name <- set_names[i]
    set_dir <- file.path(output_dir, set_name)

    if (verbose) cat("[APOLLO] --- [", i, "/", n_sets, "] ", set_name, " ---\n", sep = "")

    results[[set_name]] <- APOLLO_homer_motif_enrichment(
      bed_file   = bed_files[i],
      genome     = genome,
      output_dir = set_dir,
      ...,
      verbose    = verbose
    )
  }

  # Build batch object
  batch <- .homer_build_batch(results, genome, bed_files)

  if (verbose) {
    cat("    Batch complete. Summary:\n")
    print(batch$summary)
  }

  return(batch)
}


#' Load HOMER results from an existing output directory
#'
#' Parses HOMER findMotifsGenome.pl output without re-running. Useful for
#' loading previously computed results.
#'
#' @param homer_dir Path to HOMER output directory (containing knownResults.txt).
#' @param verbose Logical. Print loading messages. Default TRUE.
#'
#' @return S3 object of class "homer_motif" with components:
#'   \describe{
#'     \item{known}{Data.frame of known motif enrichment results}
#'     \item{denovo}{Data.frame of de novo motifs (NULL if not available)}
#'     \item{output_dir}{Path to HOMER output directory}
#'     \item{metadata}{List with n_target_seqs, n_background_seqs}
#'   }
#'
#' @export
APOLLO_load_homer_results <- function(homer_dir, verbose = TRUE) {

  if (!dir.exists(homer_dir)) {
    stop("HOMER output directory not found: ", homer_dir, call. = FALSE)
  }

  if (verbose) cat("[APOLLO] Loading HOMER results from:", homer_dir, "\n")

  # Parse known motif results
  known_file <- file.path(homer_dir, "knownResults.txt")
  known <- .homer_parse_known(known_file)

  if (is.null(known)) {
    warning("No known motif results found in: ", homer_dir)
  } else if (verbose) {
    cat("    Known motifs:", nrow(known), "\n")
    n_sig <- sum(known$q_value < 0.05, na.rm = TRUE)
    cat("    Significant (q < 0.05):", n_sig, "\n")
  }

  # Parse de novo motif results (if available)
  denovo_file <- file.path(homer_dir, "homerMotifs.all.motifs")
  denovo <- .homer_parse_denovo(denovo_file)

  if (!is.null(denovo) && verbose) {
    cat("    De novo motifs:", nrow(denovo), "\n")
  }

  # Extract target/background counts from known results
  n_target <- NA_integer_
  n_background <- NA_integer_
  if (!is.null(known) && nrow(known) > 0) {
    n_target <- known$total_target[1]
    n_background <- known$total_background[1]
  }

  result <- list(
    known      = known,
    denovo     = denovo,
    output_dir = normalizePath(homer_dir, mustWork = FALSE),
    metadata   = list(
      n_target_seqs     = n_target,
      n_background_seqs = n_background
    )
  )
  class(result) <- c("homer_motif", "list")

  return(result)
}


#' Load multiple HOMER result directories into a batch object
#'
#' Scans a base directory for HOMER output subdirectories and loads them
#' into a combined homer_motif_batch object.
#'
#' @param output_dir Base directory containing subdirectories per peak set.
#' @param set_names Optional character vector of subdirectory names to load.
#'   If NULL, auto-detects from subdirectories containing knownResults.txt.
#' @param verbose Logical. Print loading messages. Default TRUE.
#'
#' @return S3 object of class "homer_motif_batch" (same as batch enrichment).
#'
#' @export
APOLLO_load_homer_batch <- function(output_dir, set_names = NULL, verbose = TRUE) {

  if (!dir.exists(output_dir)) {
    stop("Output directory not found: ", output_dir, call. = FALSE)
  }

  # Auto-detect set directories
  if (is.null(set_names)) {
    all_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = FALSE)
    set_names <- all_dirs[file.exists(file.path(output_dir, all_dirs, "knownResults.txt"))]
    if (length(set_names) == 0) {
      stop("No HOMER result subdirectories found in: ", output_dir, call. = FALSE)
    }
  }

  if (verbose) {
    cat("[APOLLO] Loading HOMER batch results from:", output_dir, "\n")
    cat("    Sets found:", length(set_names), "\n")
  }

  # Load each set
  results <- list()
  for (set_name in set_names) {
    set_dir <- file.path(output_dir, set_name)
    if (verbose) cat("    Loading:", set_name, "\n")
    results[[set_name]] <- APOLLO_load_homer_results(set_dir, verbose = FALSE)
  }

  batch <- .homer_build_batch(results, genome = NA_character_, bed_files = NULL)

  if (verbose) {
    cat("    Batch loaded. Summary:\n")
    print(batch$summary)
  }

  return(batch)
}


#' Filter motif enrichment results
#'
#' Filter known motif results by significance thresholds and/or top N.
#' Works with both single (homer_motif) and batch (homer_motif_batch) objects.
#'
#' @param homer_result A homer_motif or homer_motif_batch object.
#' @param q_thresh Maximum q-value (Benjamini). Default 0.05.
#' @param p_thresh Maximum raw p-value. If set, overrides q_thresh.
#' @param min_target_pct Minimum percentage of target sequences with motif.
#' @param top_n Keep top N motifs by q-value (per peak set for batch).
#' @param verbose Logical. Print filter statistics. Default TRUE.
#'
#' @return Same class as input, with filtered results.
#'
#' @export
APOLLO_filter_motifs <- function(homer_result,
                                  q_thresh = 0.05,
                                  p_thresh = NULL,
                                  min_target_pct = NULL,
                                  top_n = NULL,
                                  verbose = TRUE) {

  # --- Single result ---
  if (inherits(homer_result, "homer_motif")) {
    if (is.null(homer_result$known) || nrow(homer_result$known) == 0) {
      if (verbose) cat("[APOLLO] No known motifs to filter.\n")
      return(homer_result)
    }

    n_before <- nrow(homer_result$known)
    df <- homer_result$known

    # Apply filters
    if (!is.null(p_thresh)) {
      df <- df[!is.na(df$p_value) & df$p_value <= p_thresh, ]
    } else {
      df <- df[!is.na(df$q_value) & df$q_value <= q_thresh, ]
    }

    if (!is.null(min_target_pct)) {
      df <- df[!is.na(df$pct_target) & df$pct_target >= min_target_pct, ]
    }

    if (!is.null(top_n) && nrow(df) > top_n) {
      df <- df[order(df$q_value), ]
      df <- df[seq_len(top_n), ]
    }

    homer_result$known <- df

    if (verbose) cat("[APOLLO] Filtered:", n_before, "->", nrow(df), "known motifs\n")
    return(homer_result)
  }

  # --- Batch result ---
  if (inherits(homer_result, "homer_motif_batch")) {
    # Filter individual results
    for (nm in names(homer_result$results)) {
      homer_result$results[[nm]] <- APOLLO_filter_motifs(
        homer_result$results[[nm]],
        q_thresh = q_thresh, p_thresh = p_thresh,
        min_target_pct = min_target_pct, top_n = top_n,
        verbose = FALSE
      )
    }

    # Rebuild combined
    n_before <- nrow(homer_result$combined)
    homer_result <- .homer_rebuild_batch(homer_result)

    if (verbose) {
      cat("[APOLLO] Filtered combined:", n_before, "->", nrow(homer_result$combined), "rows\n")
    }
    return(homer_result)
  }

  stop("Input must be a homer_motif or homer_motif_batch object.", call. = FALSE)
}


# ==============================================================================
# Batch helpers
# ==============================================================================

#' Build homer_motif_batch from a list of homer_motif results
#' @keywords internal
.homer_build_batch <- function(results, genome, bed_files) {

  set_names <- names(results)

  # Combine known motif results
  combined_list <- lapply(set_names, function(nm) {
    df <- results[[nm]]$known
    if (!is.null(df) && nrow(df) > 0) {
      df$peak_set <- nm
      df
    } else {
      NULL
    }
  })
  combined <- do.call(rbind, combined_list)

  # Build summary
  summary_df <- data.frame(
    peak_set = set_names,
    n_known_motifs = sapply(results, function(r) {
      if (!is.null(r$known)) nrow(r$known) else 0L
    }),
    n_significant = sapply(results, function(r) {
      if (!is.null(r$known)) sum(r$known$q_value < 0.05, na.rm = TRUE) else 0L
    }),
    top_motif = sapply(results, function(r) {
      if (!is.null(r$known) && nrow(r$known) > 0) {
        r$known$motif_family[which.min(r$known$q_value)]
      } else {
        NA_character_
      }
    }),
    n_target_seqs = sapply(results, function(r) {
      r$metadata$n_target_seqs
    }),
    has_denovo = sapply(results, function(r) !is.null(r$denovo)),
    stringsAsFactors = FALSE
  )

  batch <- list(
    results  = results,
    combined = combined,
    summary  = summary_df,
    metadata = list(
      genome    = genome,
      bed_files = bed_files,
      n_sets    = length(results)
    )
  )
  class(batch) <- c("homer_motif_batch", "list")

  return(batch)
}


#' Rebuild batch combined/summary from individual results
#' @keywords internal
.homer_rebuild_batch <- function(batch) {
  rebuilt <- .homer_build_batch(
    results   = batch$results,
    genome    = batch$metadata$genome,
    bed_files = batch$metadata$bed_files
  )
  return(rebuilt)
}


# ==============================================================================
# Print methods
# ==============================================================================

#' @method print homer_motif
#' @export
print.homer_motif <- function(x, ...) {
  cat("HOMER Motif Enrichment Result\n")
  cat("------------------------------\n")
  cat("Output:", x$output_dir, "\n")
  if (!is.null(x$known)) {
    cat("Known motifs:", nrow(x$known), "\n")
    n_sig <- sum(x$known$q_value < 0.05, na.rm = TRUE)
    cat("Significant (q < 0.05):", n_sig, "\n")
    if (n_sig > 0) {
      top <- x$known[which.min(x$known$q_value), ]
      cat("Top motif:", top$motif_family, "(q =", formatC(top$q_value, format = "e", digits = 2), ")\n")
    }
  } else {
    cat("Known motifs: none loaded\n")
  }
  if (!is.null(x$denovo)) {
    cat("De novo motifs:", nrow(x$denovo), "\n")
  }
  if (!is.null(x$metadata$genome)) cat("Genome:", x$metadata$genome, "\n")
  if (!is.null(x$metadata$n_target_seqs)) cat("Target sequences:", x$metadata$n_target_seqs, "\n")
  if (!is.null(x$metadata$n_background_seqs)) cat("Background sequences:", x$metadata$n_background_seqs, "\n")
  invisible(x)
}


#' @method print homer_motif_batch
#' @export
print.homer_motif_batch <- function(x, ...) {
  cat("HOMER Motif Enrichment Batch Result\n")
  cat("------------------------------\n")
  cat("Peak sets:", x$metadata$n_sets, "\n")
  if (!is.null(x$combined)) {
    cat("Total known motif rows:", nrow(x$combined), "\n")
  }
  cat("\nPer-set summary:\n")
  cat(x$summary, row.names = FALSE)
  invisible(x)
}
