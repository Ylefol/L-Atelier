# GAIA/Eleuthia/composition_export.R
# Export sequence composition results from APOLLO_sequence_composition()
#
# Saves the composition data.frame as CSV, optionally separates the
# sequence column, and writes a summary statistics text file.


#' Export Sequence Composition Results
#'
#' Saves the data.frame returned by \code{APOLLO_sequence_composition()} to
#' CSV, with an optional separate file for the raw DNA sequences and a
#' plain-text summary of composition statistics.
#'
#' @param composition A data.frame from \code{APOLLO_sequence_composition()},
#'   or a named list of such data.frames (one per peak set). When a list is
#'   provided, each element is exported to its own subdirectory under
#'   \code{output_dir} named after the list element.
#'   Expected columns: \code{peak_id}, \code{width}, \code{gc_percent},
#'   \code{at_percent}, \code{n_Nbases}. Optional: \code{longest_homopolymer},
#'   \code{homopolymer_base}, \code{longest_dinucleotide},
#'   \code{dinucleotide_motif}, \code{sequence}.
#' @param output_dir Character. Directory to save results. Created if needed.
#' @param prefix Character. Prefix for output filenames. Default: "composition".
#' @param save_sequences Logical. If \code{TRUE} and a \code{sequence} column
#'   exists, save sequences to a separate CSV and drop the column from the
#'   main export. Default: FALSE (sequences omitted from main CSV but not
#'   saved separately).
#' @param save_rds Logical. Save the full data.frame as RDS. Default: FALSE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Invisible character vector of file paths created.
#'
#' @details
#' Exports the following files:
#' \itemize{
#'   \item \strong{Composition CSV} — all numeric columns without the sequence
#'     string (\code{{prefix}_composition.csv})
#'   \item \strong{Sequences CSV} — two-column table (\code{peak_id},
#'     \code{sequence}) if \code{save_sequences = TRUE} and a \code{sequence}
#'     column is present (\code{{prefix}_sequences.csv})
#'   \item \strong{Summary TXT} — descriptive statistics for GC%, AT%,
#'     region width, and repeat metrics (\code{{prefix}_summary.txt})
#'   \item \strong{RDS} — full data.frame if \code{save_rds = TRUE}
#' }
#'
#' @examples
#' \dontrun{
#' comp <- APOLLO_sequence_composition(peaks, "genome.fa")
#' ELEUTHIA_export_sequence_composition(comp, "results/composition/")
#'
#' # Named list — each set gets its own subdirectory
#' comp_list <- APOLLO_sequence_composition(list(WT = wt, KO = ko), "genome.fa")
#' ELEUTHIA_export_sequence_composition(comp_list, "results/composition/")
#' # -> results/composition/WT/composition_composition.csv
#' # -> results/composition/KO/composition_composition.csv
#'
#' # Save sequences separately
#' ELEUTHIA_export_sequence_composition(comp, "results/composition/",
#'                                       save_sequences = TRUE)
#' }
#' @export
ELEUTHIA_export_sequence_composition <- function(composition,
                                                  output_dir,
                                                  prefix = "composition",
                                                  save_sequences = FALSE,
                                                  save_rds = FALSE,
                                                  verbose = TRUE) {

  # List dispatch: export each element to its own subdirectory
  if (is.list(composition) && !is.data.frame(composition)) {
    set_names <- names(composition)
    if (is.null(set_names)) {
      set_names <- paste0("set_", seq_along(composition))
    }
    if (verbose) {
      cat("=== Exporting Sequence Composition Results ===\n")
      cat("Sets:", paste(set_names, collapse = ", "), "\n\n")
    }
    all_files <- character(0)
    for (i in seq_along(composition)) {
      nm <- set_names[i]
      if (verbose) cat("--- ", nm, " ---\n", sep = "")
      set_files <- ELEUTHIA_export_sequence_composition(
        composition    = composition[[i]],
        output_dir     = file.path(output_dir, nm),
        prefix         = prefix,
        save_sequences = save_sequences,
        save_rds       = save_rds,
        verbose        = verbose
      )
      all_files <- c(all_files, set_files)
    }
    if (verbose) {
      cat("\n--- Export Summary ---\n")
      cat("Total files created:", length(all_files), "\n")
      cat("Output directory:", output_dir, "\n")
    }
    return(invisible(all_files))
  }

  if (!is.data.frame(composition)) {
    stop("'composition' must be a data.frame or named list of data.frames.",
         call. = FALSE)
  }

  required_cols <- c("peak_id", "gc_percent", "at_percent")
  missing_cols <- setdiff(required_cols, colnames(composition))
  if (length(missing_cols) > 0) {
    stop("'composition' is missing expected columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  n_regions <- nrow(composition)
  if (n_regions == 0) {
    stop("'composition' data.frame is empty.", call. = FALSE)
  }

  if (verbose) {
    cat("=== Exporting Sequence Composition Results ===\n")
    cat("Regions:", n_regions, "\n")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created directory:", output_dir, "\n")
  }

  files_created <- character(0)
  has_sequence <- "sequence" %in% colnames(composition)

  # --------------------------------------------------------------------------
  # Sequences CSV (optional)
  # --------------------------------------------------------------------------
  if (has_sequence && save_sequences) {
    seq_file <- file.path(output_dir, paste0(prefix, "_sequences.csv"))
    seq_df <- composition[, c("peak_id", "sequence"), drop = FALSE]
    write.csv(seq_df, seq_file, row.names = FALSE)
    files_created <- c(files_created, seq_file)
    if (verbose) cat("  Sequences:", seq_file, "\n")
  }

  # --------------------------------------------------------------------------
  # Composition CSV (without sequence column)
  # --------------------------------------------------------------------------
  comp_export <- composition[, setdiff(colnames(composition), "sequence"),
                               drop = FALSE]
  comp_file <- file.path(output_dir, paste0(prefix, "_composition.csv"))
  write.csv(comp_export, comp_file, row.names = FALSE)
  files_created <- c(files_created, comp_file)
  if (verbose) cat("  Composition:", comp_file, "\n")

  # --------------------------------------------------------------------------
  # Summary TXT
  # --------------------------------------------------------------------------
  .fmt_stat <- function(x, digits = 2) {
    x_valid <- x[!is.na(x)]
    if (length(x_valid) == 0) return("  (no data)")
    paste0(
      "  n = ", length(x_valid), "\n",
      "  mean = ", round(mean(x_valid), digits), "\n",
      "  median = ", round(median(x_valid), digits), "\n",
      "  range = [", round(min(x_valid), digits), ", ",
      round(max(x_valid), digits), "]"
    )
  }

  summary_lines <- c(
    "=== Sequence Composition Summary ===",
    paste("Regions:", n_regions),
    paste("Exported:", Sys.time()),
    "",
    "--- GC Content (%) ---",
    .fmt_stat(composition$gc_percent),
    "",
    "--- AT Content (%) ---",
    .fmt_stat(composition$at_percent)
  )

  if ("width" %in% colnames(composition)) {
    summary_lines <- c(summary_lines,
      "",
      "--- Region Width (bp) ---",
      .fmt_stat(composition$width, digits = 0)
    )
  }

  if ("n_Nbases" %in% colnames(composition)) {
    n_valid <- sum(!is.na(composition$n_Nbases))
    pct_with_N <- round(100 * sum(composition$n_Nbases > 0, na.rm = TRUE) /
                          n_valid, 1)
    summary_lines <- c(summary_lines,
      "",
      "--- N Bases ---",
      paste("  Regions with N bases:", sum(composition$n_Nbases > 0, na.rm = TRUE),
            paste0("(", pct_with_N, "%)"))
    )
  }

  if ("longest_homopolymer" %in% colnames(composition)) {
    summary_lines <- c(summary_lines,
      "",
      "--- Longest Homopolymer (bp) ---",
      .fmt_stat(composition$longest_homopolymer, digits = 1)
    )
    if ("homopolymer_base" %in% colnames(composition)) {
      base_tbl <- sort(table(composition$homopolymer_base), decreasing = TRUE)
      base_tbl <- base_tbl[base_tbl > 0]
      if (length(base_tbl) > 0) {
        summary_lines <- c(summary_lines, "  Dominant base frequency:")
        for (i in seq_along(base_tbl)) {
          summary_lines <- c(summary_lines,
            paste0("    ", names(base_tbl)[i], ": ", base_tbl[i])
          )
        }
      }
    }
  }

  if ("longest_dinucleotide" %in% colnames(composition)) {
    summary_lines <- c(summary_lines,
      "",
      "--- Longest Dinucleotide Repeat (bp) ---",
      .fmt_stat(composition$longest_dinucleotide, digits = 1)
    )
    if ("dinucleotide_motif" %in% colnames(composition)) {
      di_tbl <- sort(table(composition$dinucleotide_motif), decreasing = TRUE)
      di_tbl <- di_tbl[di_tbl > 0]
      top_di <- head(di_tbl, 5)
      if (length(top_di) > 0) {
        summary_lines <- c(summary_lines, "  Top dinucleotide motifs:")
        for (i in seq_along(top_di)) {
          summary_lines <- c(summary_lines,
            paste0("    ", names(top_di)[i], ": ", top_di[i])
          )
        }
      }
    }
  }

  summ_file <- file.path(output_dir, paste0(prefix, "_summary.txt"))
  writeLines(summary_lines, summ_file)
  files_created <- c(files_created, summ_file)
  if (verbose) cat("  Summary:", summ_file, "\n")

  # --------------------------------------------------------------------------
  # RDS
  # --------------------------------------------------------------------------
  if (save_rds) {
    rds_file <- file.path(output_dir, paste0(prefix, "_composition.rds"))
    saveRDS(composition, rds_file)
    files_created <- c(files_created, rds_file)
    if (verbose) cat("  RDS:", rds_file, "\n")
  }

  if (verbose) {
    cat("\n--- Export Summary ---\n")
    cat("Total files created:", length(files_created), "\n")
    cat("Output directory:", output_dir, "\n")
  }

  invisible(files_created)
}
