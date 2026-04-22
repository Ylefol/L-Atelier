#' Create a PARSE Biosciences sample sheet template
#'
#' Generates a CSV template for PARSE split-pipe runs. Each row represents one
#' sequencing run (library). Fill in paths and metadata, then pass to
#' \code{\link{HORIZON_validate_parse_sheet}} before running
#' \code{\link{HORIZON_run_splitpipe}}.
#'
#' Required columns:
#' \describe{
#'   \item{run_id}{Unique run identifier. Used to name the output sub-folder.}
#'   \item{fastq_r1}{Absolute path to R1 FASTQ file.}
#'   \item{fastq_r2}{Absolute path to R2 FASTQ file.}
#'   \item{output_dir}{Absolute path to root output directory. A sub-folder
#'     named \code{run_id} will be created inside it.}
#'   \item{genome_dir}{Absolute path to the split-pipe genome directory
#'     (built with \code{split-pipe --mode mkref}).}
#'   \item{chemistry}{Library chemistry version: \code{"v2"} or \code{"v3"}.}
#'   \item{kit}{Library kit: \code{"WT"}, \code{"WT_mini"}, or \code{"WT_mega"}.}
#' }
#'
#' Additional columns can be added freely; they are carried through as metadata.
#'
#' @param output_path Character or NULL. Path to write the template CSV.
#'   If NULL, the data.frame is returned without writing. Default NULL.
#' @param n_runs Integer. Number of placeholder rows to include. Default 1.
#'
#' @return A data.frame with the template, invisibly.
#' @export
HORIZON_create_parse_sheet <- function(output_path = NULL, n_runs = 1) {
  template <- data.frame(
    run_id     = paste0("run_", seq_len(n_runs)),
    fastq_r1   = rep("", n_runs),
    fastq_r2   = rep("", n_runs),
    output_dir = rep("", n_runs),
    genome_dir = rep("", n_runs),
    chemistry  = rep("v3", n_runs),
    kit        = rep("WT", n_runs),
    stringsAsFactors = FALSE
  )

  if (!is.null(output_path)) {
    write.csv(template, output_path, row.names = FALSE)
    cat("[PARSE] Sample sheet template written to:", output_path, "\n")
  }

  invisible(template)
}


#' Validate a PARSE Biosciences sample sheet
#'
#' Reads a PARSE sample sheet CSV (or accepts a data.frame) and validates it
#' before use. Checks required columns, unique run IDs, valid chemistry and kit
#' values, and optionally whether all paths exist on disk.
#'
#' @param path Character or data.frame. Path to a CSV, or an already-loaded
#'   data.frame.
#' @param check_files Logical. Whether to verify that FASTQ files and genome
#'   directories exist on disk. Default TRUE. Set FALSE if paths are on a drive
#'   not currently mounted.
#'
#' @return The validated sample sheet as a data.frame.
#' @export
HORIZON_validate_parse_sheet <- function(path, check_files = TRUE) {
  if (is.character(path)) {
    if (!file.exists(path))
      stop("Sample sheet file not found: ", path, call. = FALSE)
    ss <- read.csv(path, stringsAsFactors = FALSE)
  } else if (is.data.frame(path)) {
    ss <- path
  } else {
    stop("'path' must be a file path string or a data.frame.", call. = FALSE)
  }

  required_cols <- c("run_id", "fastq_r1", "fastq_r2",
                     "output_dir", "genome_dir", "chemistry", "kit")
  missing_cols  <- setdiff(required_cols, colnames(ss))
  if (length(missing_cols) > 0)
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)

  if (anyDuplicated(ss$run_id)) {
    dups <- unique(ss$run_id[duplicated(ss$run_id)])
    stop("Duplicate run_id values: ", paste(dups, collapse = ", "), call. = FALSE)
  }

  .validate_chemistry_kit(ss)

  if (check_files) {
    missing_r1 <- ss$fastq_r1[!file.exists(ss$fastq_r1)]
    if (length(missing_r1) > 0)
      stop("R1 FASTQ file(s) not found:\n",
           paste(" ", missing_r1, collapse = "\n"), call. = FALSE)

    missing_r2 <- ss$fastq_r2[!file.exists(ss$fastq_r2)]
    if (length(missing_r2) > 0)
      stop("R2 FASTQ file(s) not found:\n",
           paste(" ", missing_r2, collapse = "\n"), call. = FALSE)

    missing_genome <- ss$genome_dir[!dir.exists(ss$genome_dir)]
    if (length(missing_genome) > 0)
      stop("Genome directory/directories not found:\n",
           paste(" ", missing_genome, collapse = "\n"), call. = FALSE)
  }

  cat("[PARSE] Sample sheet OK —", nrow(ss), "run(s) validated.\n")
  ss
}


#' Run PARSE Biosciences split-pipe
#'
#' Wraps the \code{split-pipe} CLI to perform demultiplexing, alignment, and
#' cell/nucleus barcode counting for a single PARSE Biosciences run. Supports
#' both single-cell and single-nucleus RNA-seq data.
#'
#' \code{split-pipe} must be installed in a dedicated conda environment passed
#' via \code{conda_env}. It cannot share the standard \code{horizon_cli}
#' environment due to conflicting Python dependencies.
#'
#' Example installation:
#' \preformatted{
#'   conda create -n parse_env python=3.10
#'   conda activate parse_env
#'   pip install parsebiosciences
#' }
#'
#' \strong{Sample layout format:}
#' A named character vector mapping sample names to split-pipe well range
#' strings. Names become the split-pipe sample identifiers.
#' \preformatted{
#'   # All wells as a single sample
#'   sample_layout = c(all_cells = "A1-D12")
#'
#'   # Two conditions across separate well ranges
#'   sample_layout = c(control = "A1-B6", treatment = "C1-D6")
#' }
#'
#' The split-pipe output directory (\code{output_dir/run_id/}) contains a
#' per-sample DGE matrix compatible with
#' \code{TALARIA_load_parse()} in CAULDRON.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_parse_sheet}}. The output location for each
#'   run is controlled by the \code{output_dir} column in this sheet — a
#'   sub-folder named \code{run_id} is created inside it. There is no separate
#'   \code{output_dir} parameter on this function.
#' @param run_id Character. Run ID to process (must match a row in
#'   \code{sample_sheet}).
#' @param sample_layout Named character vector. Maps sample names to split-pipe
#'   well range strings (e.g. \code{c(control = "A1-B6", treatment = "C1-D6")}).
#' @param conda_env Character. Path to the conda environment containing
#'   \code{split-pipe}
#'   (e.g. \code{"/path/to/conda/envs/parse_env"}).
#' @param mode Character. split-pipe run mode. Default \code{"all"}
#'   (demultiplex + align + count in one pass). Use \code{"split"} for
#'   demultiplexing only, or \code{"count"} to recount from a prior split
#'   output.
#' @param threads Integer. Number of threads passed to split-pipe
#'   (\code{--nthreads}). Default \code{8L}.
#' @param force Logical. If \code{FALSE} (default), skip processing when the
#'   output directory already exists. Set \code{TRUE} to reprocess and
#'   overwrite.
#' @param ... Additional split-pipe flags passed verbatim to the command.
#'
#' @return Character. Path to the split-pipe output directory, invisibly.
#' @export
HORIZON_run_splitpipe <- function(sample_sheet,
                                   run_id,
                                   sample_layout,
                                   conda_env,
                                   mode    = "all",
                                   threads = 8L,
                                   force   = FALSE,
                                   ...) {

  if (!is.character(run_id) || length(run_id) != 1L)
    stop("run_id must be a single character string.", call. = FALSE)

  row <- .get_parse_row(sample_sheet, run_id)

  if (!is.character(sample_layout) || is.null(names(sample_layout)) ||
      any(!nzchar(names(sample_layout))))
    stop("sample_layout must be a named character vector ",
         "(e.g. c(control = 'A1-B6', treatment = 'C1-D6')).", call. = FALSE)

  if (!is.character(conda_env) || length(conda_env) != 1L)
    stop("conda_env must be a single path to the conda environment.",
         call. = FALSE)

  splitpipe_bin <- file.path(conda_env, "bin", "split-pipe")
  if (!file.exists(splitpipe_bin))
    stop("split-pipe not found at: ", splitpipe_bin,
         "\n  Install via: pip install parsebiosciences (inside parse_env)",
         call. = FALSE)

  valid_modes <- c("all", "split", "count")
  if (!mode %in% valid_modes)
    stop("mode must be one of: ", paste(valid_modes, collapse = ", "),
         call. = FALSE)

  out_dir <- file.path(row$output_dir, run_id)

  if (!isTRUE(force) && dir.exists(out_dir)) {
    cat("[PARSE] Output directory already exists for run:", run_id,
        "— skipping (use force=TRUE to reprocess)\n")
    return(invisible(out_dir))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Build one --sample <name> <wells> pair per entry
  sample_flags <- unlist(lapply(seq_along(sample_layout), function(i) {
    c("--sample", names(sample_layout)[i], sample_layout[i])
  }))

  args <- c(
    "--mode",       mode,
    "--chemistry",  row$chemistry,
    "--kit",        row$kit,
    "--nthreads",   as.integer(threads),
    "--fq1",        row$fastq_r1,
    "--fq2",        row$fastq_r2,
    "--output_dir", out_dir,
    "--genome_dir", row$genome_dir,
    sample_flags,
    c(...)
  )

  cat("[PARSE] Running split-pipe for run:", run_id, "\n")
  cat("    Mode:      ", mode, "\n")
  cat("    Chemistry: ", row$chemistry, "\n")
  cat("    Kit:       ", row$kit, "\n")
  cat("    Genome:    ", row$genome_dir, "\n")
  cat("    Output:    ", out_dir, "\n")
  cat("    Samples:   ",
      paste(names(sample_layout), sample_layout, sep = "=", collapse = ", "),
      "\n\n")

  ret <- system2(splitpipe_bin, args = args)

  if (ret != 0L)
    stop("split-pipe failed for run: ", run_id,
         " (exit code ", ret, ")", call. = FALSE)

  cat("[PARSE] split-pipe complete for run:", run_id, "\n")
  cat("    Output:", out_dir, "\n")
  invisible(out_dir)
}


#' Generate spliced/unspliced velocity matrices from PARSE split-pipe output
#'
#' Calls \code{inst/python/parse_velocity.py} to build per-run spliced and
#' unspliced sparse matrices from the \code{tscp_assignment.csv} files produced
#' by \code{\link{HORIZON_run_splitpipe}}, concatenates them across all runs,
#' filters to cells present in the combined split-pipe cell metadata, and writes
#' a single \code{adata_vel.h5ad} ready for RNA velocity analysis in CAULDRON.
#'
#' \strong{Prerequisites:}
#' \itemize{
#'   \item All runs in \code{sample_sheet} must have been processed by
#'     \code{\link{HORIZON_run_splitpipe}}.
#'   \item \code{tscp_assignment.csv} files must be uncompressed. If split-pipe
#'     wrote \code{.csv.gz} files, run \code{gunzip} on them first. This
#'     function will stop with an informative error if compressed files are
#'     found in place of the expected \code{.csv}.
#'   \item Required Python packages in \code{conda_env}: \code{scanpy},
#'     \code{scvelo}, \code{anndata}, \code{dask}, \code{pandas},
#'     \code{scipy}, \code{numpy}.
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_parse_sheet}}. All rows must share the same
#'   \code{output_dir}.
#' @param conda_env Character. Path to the conda environment containing
#'   \code{split-pipe} and the required Python packages — the same environment
#'   used for \code{\link{HORIZON_run_splitpipe}}.
#' @param cell_metadata_path Character. Path to the \code{cell_metadata.csv}
#'   produced by split-pipe's combined output (e.g.
#'   \code{"<output_dir>/combined/<sample>/DGE_unfiltered/cell_metadata.csv"}).
#' @param output_file Character or NULL. Path for the final
#'   \code{adata_vel.h5ad}. Defaults to \code{<output_dir>/adata_vel.h5ad}.
#' @param force Logical. If \code{FALSE} (default), skip processing when
#'   \code{output_file} already exists. Set \code{TRUE} to reprocess.
#'
#' @return Character. Path to the written \code{adata_vel.h5ad}, invisibly.
#' @export
HORIZON_run_parse_velocity <- function(sample_sheet,
                                        conda_env,
                                        cell_metadata_path,
                                        output_file = NULL,
                                        force       = FALSE) {

  # Resolve working_dir from sample sheet (output_dir must be consistent)
  output_dirs <- unique(sample_sheet$output_dir)
  if (length(output_dirs) > 1)
    stop("All rows in sample_sheet must share the same output_dir for velocity ",
         "processing. Found: ", paste(output_dirs, collapse = ", "), call. = FALSE)
  working_dir <- output_dirs

  run_ids <- sample_sheet$run_id

  if (!is.character(conda_env) || length(conda_env) != 1L)
    stop("conda_env must be a single path to the conda environment.", call. = FALSE)

  python_bin <- file.path(conda_env, "bin", "python")
  if (!file.exists(python_bin))
    stop("Python not found at: ", python_bin,
         "\n  Check that conda_env points to a valid environment.", call. = FALSE)

  if (is.null(output_file))
    output_file <- file.path(working_dir, "adata_vel.h5ad")

  if (!isTRUE(force) && file.exists(output_file)) {
    cat("[PARSE] Velocity AnnData already exists — skipping",
        "(use force=TRUE to reprocess)\n")
    cat("    Output:", output_file, "\n")
    return(invisible(output_file))
  }

  if (!file.exists(cell_metadata_path))
    stop("cell_metadata_path not found: ", cell_metadata_path, call. = FALSE)

  # Check for uncompressed tscp files; warn clearly if .gz found
  tscp_csv  <- file.path(working_dir, run_ids, "process", "tscp_assignment.csv")
  tscp_gz   <- paste0(tscp_csv, ".gz")
  missing   <- !file.exists(tscp_csv)
  has_gz    <- file.exists(tscp_gz)

  if (any(missing)) {
    gz_msg <- ifelse(has_gz[missing],
                     " (.gz found — run gunzip first)",
                     " (file not found)")
    stop("tscp_assignment.csv missing for run(s):\n",
         paste0("  ", run_ids[missing], gz_msg[missing], collapse = "\n"),
         call. = FALSE)
  }

  script <- system.file("python", "parse_velocity.py", package = "HORIZON")
  if (!nzchar(script))
    stop("parse_velocity.py not found in HORIZON installation.", call. = FALSE)

  args <- c(
    script,
    "--working_dir",   working_dir,
    "--run_ids",       paste(run_ids, collapse = ","),
    "--cell_metadata", cell_metadata_path,
    "--output_file",   output_file
  )

  cat("[PARSE] Generating velocity matrices for", length(run_ids), "run(s)\n")
  cat("    Working dir:     ", working_dir, "\n")
  cat("    Cell metadata:   ", cell_metadata_path, "\n")
  cat("    Output file:     ", output_file, "\n\n")

  ret <- system2(python_bin, args = args)

  if (ret != 0L)
    stop("parse_velocity.py failed (exit code ", ret, ")", call. = FALSE)

  cat("[PARSE] Velocity processing complete\n")
  cat("    Output:", output_file, "\n")
  invisible(output_file)
}


# Internal: extract a single run row from a PARSE sample sheet
.get_parse_row <- function(sample_sheet, run_id) {
  row <- sample_sheet[sample_sheet$run_id == run_id, ]
  if (nrow(row) == 0)
    stop("run_id not found in sample sheet: '", run_id, "'", call. = FALSE)
  if (nrow(row) > 1)
    stop("Duplicate run_id in sample sheet: '", run_id, "'", call. = FALSE)
  row
}
