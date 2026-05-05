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
#'   \item{kit}{Library kit: \code{"WT"}, \code{"WT_mini"}, \code{"WT_mega"},
#'     \code{"WT_mega_384"}, \code{"WT_penta"}, or \code{"WT_penta_384"}.}
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
#' @param kit_score_skip Logical. If \code{TRUE}, passes \code{--kit_score_skip}
#'   to split-pipe to bypass the automatic kit barcode check. Useful when the
#'   detected kit does not exactly match the specified kit but processing should
#'   proceed anyway. Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints the full split-pipe command
#'   before executing it. Useful for verifying parameter passing. Default
#'   \code{TRUE}.
#' @param ... Additional split-pipe flags passed verbatim to the command.
#'
#' @return Character. Path to the split-pipe output directory, invisibly.
#' @export
HORIZON_run_splitpipe <- function(sample_sheet,
                                   run_id,
                                   sample_layout,
                                   conda_env,
                                   mode           = "all",
                                   threads        = 8L,
                                   force          = FALSE,
                                   kit_score_skip = FALSE,
                                   verbose        = TRUE,
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
  conda_env <- path.expand(conda_env)

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
    if (isTRUE(kit_score_skip)) "--kit_score_skip",
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

  if (isTRUE(verbose))
    cat("[PARSE] Command: ", splitpipe_bin, paste(args, collapse = " "), "\n\n")

  ret <- .horizon_run_with_log(splitpipe_bin, args, conda_env)

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
#' and writes a single raw \code{adata_vel.h5ad}. No cell/gene filtering or
#' metadata addition is performed — use \code{\link{HORIZON_parse_DGE_filter}}
#' for that step.
#'
#' \strong{Prerequisites:}
#' \itemize{
#'   \item All runs in \code{sample_sheet} must have been processed by
#'     \code{\link{HORIZON_run_splitpipe}} and combined with
#'     \code{\link{HORIZON_combine_splitpipe}}.
#'   \item \code{tscp_assignment.csv} files may be gzip-compressed
#'     (\code{.csv.gz}); this function will decompress them automatically
#'     in-place before processing.
#'   \item Required Python packages in \code{conda_env}: \code{scanpy},
#'     \code{scvelo}, \code{anndata}, \code{dask}, \code{pandas},
#'     \code{scipy}, \code{numpy}.
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_parse_sheet}}. All rows must share the same
#'   \code{output_dir}.
#' @param conda_env Character. Path to the conda environment containing
#'   the required Python packages.
#' @param output_file Character or NULL. Path for the raw
#'   \code{adata_vel.h5ad}. Defaults to \code{<output_dir>/adata_vel.h5ad}.
#' @param force Logical. If \code{FALSE} (default), skip processing when
#'   \code{output_file} already exists. Set \code{TRUE} to reprocess.
#'
#' @return Character. Path to the written \code{adata_vel.h5ad}, invisibly.
#' @export
HORIZON_run_parse_velocity <- function(sample_sheet,
                                        conda_env,
                                        output_file = NULL,
                                        force       = FALSE) {

  output_dirs <- unique(sample_sheet$output_dir)
  if (length(output_dirs) > 1)
    stop("All rows in sample_sheet must share the same output_dir for velocity ",
         "processing. Found: ", paste(output_dirs, collapse = ", "), call. = FALSE)
  working_dir <- output_dirs

  run_ids <- sample_sheet$run_id

  if (!is.character(conda_env) || length(conda_env) != 1L)
    stop("conda_env must be a single path to the conda environment.", call. = FALSE)
  conda_env <- path.expand(conda_env)

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

  # Locate tscp files; auto-gunzip if only compressed versions are present
  tscp_csv <- file.path(working_dir, run_ids, "process", "tscp_assignment.csv")
  tscp_gz  <- paste0(tscp_csv, ".gz")
  missing  <- !file.exists(tscp_csv)
  has_gz   <- file.exists(tscp_gz)

  if (any(missing & has_gz)) {
    cat("[PARSE] Decompressing tscp_assignment.csv.gz file(s)...\n")
    for (gz in tscp_gz[missing & has_gz]) {
      cat("    gunzip:", gz, "\n")
      ret <- system2("gunzip", args = gz)
      if (ret != 0L)
        stop("gunzip failed for: ", gz, call. = FALSE)
    }
    missing <- !file.exists(tscp_csv)
  }

  if (any(missing)) {
    stop("tscp_assignment.csv missing for run(s):\n",
         paste0("  ", run_ids[missing], collapse = "\n"),
         call. = FALSE)
  }

  script <- system.file("python", "parse_velocity.py", package = "HORIZON")
  if (!nzchar(script))
    stop("parse_velocity.py not found in HORIZON installation.", call. = FALSE)

  args <- c(
    script,
    "--working_dir", working_dir,
    "--run_ids",     paste(run_ids, collapse = ","),
    "--output_file", output_file
  )

  cat("[PARSE] Generating velocity matrices for", length(run_ids), "run(s)\n")
  cat("    Working dir: ", working_dir, "\n")
  cat("    Output file: ", output_file, "\n\n")

  ret <- .horizon_run_with_log(python_bin, args, conda_env)

  if (ret != 0L)
    stop("parse_velocity.py failed (exit code ", ret, ")", call. = FALSE)

  cat("[PARSE] Velocity processing complete\n")
  cat("    Output:", output_file, "\n")
  cat("    Run HORIZON_parse_DGE_filter() to filter cells/genes and add metadata.\n")
  invisible(output_file)
}


#' Filter velocity AnnData using PARSE DGE output and add metadata
#'
#' Takes the raw \code{adata_vel.h5ad} from \code{\link{HORIZON_run_parse_velocity}}
#' and filters it to the cells and genes present in the split-pipe DGE output
#' (typically \code{DGE_filtered}), then adds cell and gene metadata from
#' \code{cell_metadata.csv} and \code{all_genes.csv}.
#'
#' @param velocity_h5ad Character. Path to the raw \code{adata_vel.h5ad}
#'   produced by \code{\link{HORIZON_run_parse_velocity}}.
#' @param dge_dir Character. Path to the split-pipe DGE directory (e.g.
#'   \code{"<output_dir>/combined/all-sample/DGE_filtered"}).
#'   Must contain \code{cell_metadata.csv} and \code{all_genes.csv}.
#' @param conda_env Character. Path to the conda environment used for
#'   \code{\link{HORIZON_run_parse_velocity}}.
#' @param output_file Character or NULL. Output path for the filtered h5ad.
#'   Defaults to \code{adata_vel_filtered.h5ad} in the same directory as
#'   \code{velocity_h5ad}.
#' @param force Logical. Reprocess even if \code{output_file} already exists.
#'
#' @return Character. Path to the filtered h5ad, invisibly.
#' @export
HORIZON_parse_DGE_filter <- function(velocity_h5ad,
                                      dge_dir,
                                      conda_env,
                                      output_file = NULL,
                                      force       = FALSE) {

  velocity_h5ad <- path.expand(velocity_h5ad)
  dge_dir       <- path.expand(dge_dir)
  conda_env     <- path.expand(conda_env)

  if (!file.exists(velocity_h5ad))
    stop("velocity_h5ad not found: ", velocity_h5ad, call. = FALSE)
  if (!dir.exists(dge_dir))
    stop("dge_dir not found: ", dge_dir, call. = FALSE)

  cell_meta <- file.path(dge_dir, "cell_metadata.csv")
  all_genes <- file.path(dge_dir, "all_genes.csv")

  if (!file.exists(cell_meta))
    stop("cell_metadata.csv not found in dge_dir: ", dge_dir, call. = FALSE)
  if (!file.exists(all_genes))
    stop("all_genes.csv not found in dge_dir: ", dge_dir, call. = FALSE)

  python_bin <- file.path(conda_env, "bin", "python")
  if (!file.exists(python_bin))
    stop("Python not found at: ", python_bin, call. = FALSE)

  if (is.null(output_file))
    output_file <- file.path(dirname(velocity_h5ad), "adata_vel_filtered.h5ad")

  if (!isTRUE(force) && file.exists(output_file)) {
    cat("[PARSE] Filtered AnnData already exists — skipping",
        "(use force=TRUE to reprocess)\n")
    cat("    Output:", output_file, "\n")
    return(invisible(output_file))
  }

  script <- system.file("python", "parse_dge_filter.py", package = "HORIZON")
  if (!nzchar(script))
    stop("parse_dge_filter.py not found in HORIZON installation.", call. = FALSE)

  args <- c(
    script,
    "--velocity_h5ad", velocity_h5ad,
    "--cell_metadata", cell_meta,
    "--all_genes",     all_genes,
    "--output_file",   output_file
  )

  cat("[PARSE] Filtering velocity AnnData\n")
  cat("    Input:         ", velocity_h5ad, "\n")
  cat("    DGE directory: ", dge_dir, "\n")
  cat("    Output:        ", output_file, "\n\n")

  ret <- .horizon_run_with_log(python_bin, args, conda_env)

  if (ret != 0L)
    stop("parse_dge_filter.py failed (exit code ", ret, ")", call. = FALSE)

  cat("[PARSE] DGE filtering complete\n")
  cat("    Output:", output_file, "\n")
  invisible(output_file)
}


#' Combine multiple PARSE split-pipe runs
#'
#' Wraps the \code{split-pipe --mode comb} step to merge all per-run outputs
#' produced by \code{\link{HORIZON_run_splitpipe}} into a single combined
#' dataset. This step is required before running
#' \code{\link{HORIZON_run_parse_velocity}} and before loading data into
#' CAULDRON via \code{TALARIA_load_parse()}.
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_parse_sheet}}. All rows must share the same
#'   \code{output_dir}. Each run's split-pipe output directory
#'   (\code{output_dir/run_id/}) is passed as a sublibrary to split-pipe.
#' @param conda_env Character. Path to the conda environment containing
#'   \code{split-pipe} — the same environment used for
#'   \code{\link{HORIZON_run_splitpipe}}.
#' @param output_subdir Character. Name of the subdirectory within
#'   \code{output_dir} where the combined output will be written. Default
#'   \code{"combined"}. The resulting path
#'   (\code{output_dir/combined/all-sample/}) is what
#'   \code{\link{HORIZON_run_parse_velocity}} expects for
#'   \code{cell_metadata_path}.
#' @param force Logical. If \code{FALSE} (default), skip if the combined
#'   output directory already exists. Set \code{TRUE} to reprocess.
#' @param verbose Logical. If \code{TRUE} (default), prints the full
#'   split-pipe command before executing.
#'
#' @return Character. Path to the combined output directory, invisibly.
#' @export
HORIZON_combine_splitpipe <- function(sample_sheet,
                                       conda_env,
                                       output_subdir = "combined",
                                       force         = FALSE,
                                       verbose       = TRUE) {

  output_dirs <- unique(sample_sheet$output_dir)
  if (length(output_dirs) > 1)
    stop("All rows in sample_sheet must share the same output_dir. Found: ",
         paste(output_dirs, collapse = ", "), call. = FALSE)

  combined_dir <- file.path(output_dirs, output_subdir)

  if (!is.character(conda_env) || length(conda_env) != 1L)
    stop("conda_env must be a single path to the conda environment.", call. = FALSE)
  conda_env <- path.expand(conda_env)

  splitpipe_bin <- file.path(conda_env, "bin", "split-pipe")
  if (!file.exists(splitpipe_bin))
    stop("split-pipe not found at: ", splitpipe_bin,
         "\n  Install via: pip install parsebiosciences (inside parse_env)",
         call. = FALSE)

  if (!isTRUE(force) && dir.exists(combined_dir)) {
    cat("[PARSE] Combined output already exists — skipping",
        "(use force=TRUE to reprocess)\n")
    cat("    Output:", combined_dir, "\n")
    return(invisible(combined_dir))
  }

  sublib_dirs <- file.path(output_dirs, sample_sheet$run_id)
  missing_dirs <- sublib_dirs[!dir.exists(sublib_dirs)]
  if (length(missing_dirs) > 0)
    stop("Per-run split-pipe output missing for:\n",
         paste(" ", missing_dirs, collapse = "\n"),
         "\n  Run HORIZON_run_splitpipe() for all runs first.", call. = FALSE)

  args <- c(
    "--mode",         "comb",
    "--output_dir",   combined_dir,
    "--sublibraries", sublib_dirs
  )

  cat("[PARSE] Combining", length(sublib_dirs), "split-pipe run(s)\n")
  cat("    Runs:  ", paste(sample_sheet$run_id, collapse = ", "), "\n")
  cat("    Output:", combined_dir, "\n\n")

  if (isTRUE(verbose))
    cat("[PARSE] Command:", splitpipe_bin, paste(args, collapse = " "), "\n\n")

  ret <- .horizon_run_with_log(splitpipe_bin, args, conda_env)

  if (ret != 0L)
    stop("split-pipe combine failed (exit code ", ret, ")", call. = FALSE)

  cat("[PARSE] Combine complete\n")
  cat("    Output:", combined_dir, "\n")
  invisible(combined_dir)
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
