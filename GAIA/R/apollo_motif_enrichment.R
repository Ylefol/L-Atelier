# GAIA/Apollo/motif_enrichment.R
# HOMER motif enrichment analysis for ATAC-seq peaks
#
# Wraps HOMER's findMotifsGenome.pl for known motif enrichment and optional
# de novo motif discovery. Supports single and batch analysis of BED files.
#
# Requires HOMER to be installed and on PATH.
# See: http://homer.ucsd.edu/homer/


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
  if (verbose) cat("HOMER found:", homer_path, "\n")

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
    cat("Running HOMER findMotifsGenome.pl...\n")
    cat("  BED file:", bed_file, "\n")
    cat("  Genome:", genome, "\n")
    cat("  Output:", output_dir, "\n")
    cat("  De novo:", if (denovo) paste("yes (n =", denovo_n, ")") else "no", "\n")
    if (!is.null(bg)) cat("  Background:", bg, "\n")
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

  if (verbose) cat("HOMER completed in", round(elapsed, 1), "seconds.\n")

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
    cat("Running HOMER batch motif enrichment\n")
    cat("  Peak sets:", n_sets, "\n")
    cat("  Sets:", paste(set_names, collapse = ", "), "\n")
  }

  # Run each set
  results <- list()
  for (i in seq_along(bed_files)) {
    set_name <- set_names[i]
    set_dir <- file.path(output_dir, set_name)

    if (verbose) cat("\n--- [", i, "/", n_sets, "] ", set_name, " ---\n", sep = "")

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
    cat("\nBatch complete. Summary:\n")
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

  if (verbose) cat("Loading HOMER results from:", homer_dir, "\n")

  # Parse known motif results
  known_file <- file.path(homer_dir, "knownResults.txt")
  known <- .homer_parse_known(known_file)

  if (is.null(known)) {
    warning("No known motif results found in: ", homer_dir)
  } else if (verbose) {
    cat("  Known motifs:", nrow(known), "\n")
    n_sig <- sum(known$q_value < 0.05, na.rm = TRUE)
    cat("  Significant (q < 0.05):", n_sig, "\n")
  }

  # Parse de novo motif results (if available)
  denovo_file <- file.path(homer_dir, "homerMotifs.all.motifs")
  denovo <- .homer_parse_denovo(denovo_file)

  if (!is.null(denovo) && verbose) {
    cat("  De novo motifs:", nrow(denovo), "\n")
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
    cat("Loading HOMER batch results from:", output_dir, "\n")
    cat("  Sets found:", length(set_names), "\n")
  }

  # Load each set
  results <- list()
  for (set_name in set_names) {
    set_dir <- file.path(output_dir, set_name)
    if (verbose) cat("  Loading:", set_name, "\n")
    results[[set_name]] <- APOLLO_load_homer_results(set_dir, verbose = FALSE)
  }

  batch <- .homer_build_batch(results, genome = NA_character_, bed_files = NULL)

  if (verbose) {
    cat("Batch loaded. Summary:\n")
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
      if (verbose) cat("No known motifs to filter.\n")
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

    if (verbose) cat("Filtered:", n_before, "->", nrow(df), "known motifs\n")
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
      cat("Filtered combined:", n_before, "->", nrow(homer_result$combined), "rows\n")
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
  cat("  Output:", x$output_dir, "\n")
  if (!is.null(x$known)) {
    cat("  Known motifs:", nrow(x$known), "\n")
    n_sig <- sum(x$known$q_value < 0.05, na.rm = TRUE)
    cat("  Significant (q < 0.05):", n_sig, "\n")
    if (n_sig > 0) {
      top <- x$known[which.min(x$known$q_value), ]
      cat("  Top motif:", top$motif_family, "(q =", formatC(top$q_value, format = "e", digits = 2), ")\n")
    }
  } else {
    cat("  Known motifs: none loaded\n")
  }
  if (!is.null(x$denovo)) {
    cat("  De novo motifs:", nrow(x$denovo), "\n")
  }
  if (!is.null(x$metadata$genome)) cat("  Genome:", x$metadata$genome, "\n")
  if (!is.null(x$metadata$n_target_seqs)) cat("  Target sequences:", x$metadata$n_target_seqs, "\n")
  if (!is.null(x$metadata$n_background_seqs)) cat("  Background sequences:", x$metadata$n_background_seqs, "\n")
  invisible(x)
}


#' @method print homer_motif_batch
#' @export
print.homer_motif_batch <- function(x, ...) {
  cat("HOMER Motif Enrichment Batch Result\n")
  cat("  Peak sets:", x$metadata$n_sets, "\n")
  if (!is.null(x$combined)) {
    cat("  Total known motif rows:", nrow(x$combined), "\n")
  }
  cat("\nPer-set summary:\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
