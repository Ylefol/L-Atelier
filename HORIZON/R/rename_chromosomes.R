#' Rename chromosome identifiers in BAM, BED, or BigWig files
#'
#' Scans the sample output directory for files whose names end with
#' \code{file_suffix}, then renames chromosome identifiers according to
#' \code{chr_map}.  Useful when aligning to a T2T or NCBI-style reference that
#' uses accession-style names (e.g. \code{NC_060925.1}) rather than conventional
#' UCSC-style names (e.g. \code{chr1}).
#'
#' \strong{BAM}: Requires \code{samtools} on \code{PATH} (only the header
#' \code{@SQ SN:} fields are rewritten; the binary read data is unchanged).
#' The output BAM is automatically re-indexed.
#'
#' \strong{BigWig}: Requires the Bioconductor package \code{rtracklayer}.
#' Install via \code{BiocManager::install("rtracklayer")}.
#'
#' \strong{BED}: Handled in pure R (no extra dependencies).
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param sample_id Character. Sample ID to process.
#' @param chr_map Named character vector (names = old chromosome identifiers,
#'   values = new names), a path to a two-column tab-separated file (no header;
#'   column 1 = old, column 2 = new), or \code{"T2TCHM13"} for the built-in
#'   T2T-CHM13v2.0 NCBI-to-UCSC map
#'   (\code{NC_060925.1} \eqn{\to} \code{chr1}, \ldots).
#' @param file_suffix Character. Files whose basenames end with this string
#'   (case-sensitive) are processed.  Use a specific trailing fragment to target
#'   one file (e.g. \code{"RPKM.bw"}, \code{"_sorted.bam"}) or a broader
#'   ending to match several files of the same type.  The file type
#'   (\code{bam}, \code{bed}, or \code{bw}) is derived from the extension of
#'   \code{file_suffix} and determines the renaming strategy.
#' @param overwrite Logical. If \code{TRUE}, the renamed file replaces the
#'   original (via a safe temp-file swap).  Default \code{FALSE}: a new file is
#'   written with \code{rename_suffix} inserted before the extension.
#' @param rename_suffix Character. String appended to the base name before the
#'   extension when \code{overwrite = FALSE}.  Default \code{"_chrrenamed"}.
#'
#' @return Character vector of output file paths, invisibly.
#' @export
HORIZON_rename_chromosomes <- function(sample_sheet,
                                       sample_id,
                                       chr_map,
                                       file_suffix,
                                       overwrite     = FALSE,
                                       rename_suffix = "_chrrenamed") {

  # ── Built-in T2T-CHM13v2.0 map ─────────────────────────────────────────────
  .t2t_map <- c(
    "NC_060925.1" = "chr1",  "NC_060926.1" = "chr2",  "NC_060927.1" = "chr3",
    "NC_060928.1" = "chr4",  "NC_060929.1" = "chr5",  "NC_060930.1" = "chr6",
    "NC_060931.1" = "chr7",  "NC_060932.1" = "chr8",  "NC_060933.1" = "chr9",
    "NC_060934.1" = "chr10", "NC_060935.1" = "chr11", "NC_060936.1" = "chr12",
    "NC_060937.1" = "chr13", "NC_060938.1" = "chr14", "NC_060939.1" = "chr15",
    "NC_060940.1" = "chr16", "NC_060941.1" = "chr17", "NC_060942.1" = "chr18",
    "NC_060943.1" = "chr19", "NC_060944.1" = "chr20", "NC_060945.1" = "chr21",
    "NC_060946.1" = "chr22", "NC_060947.1" = "chrX",  "NC_060948.1" = "chrY",
    "NC_012920.1" = "chrM"
  )

  # ── Resolve chr_map ─────────────────────────────────────────────────────────
  if (identical(chr_map, "T2TCHM13")) {
    chr_map <- .t2t_map
  } else if (is.character(chr_map) && length(chr_map) == 1 &&
             file.exists(chr_map)) {
    tbl     <- read.table(chr_map, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE)
    chr_map <- setNames(as.character(tbl[[2]]), as.character(tbl[[1]]))
  } else if (!is.character(chr_map) || is.null(names(chr_map))) {
    stop("chr_map must be a named character vector, a path to a 2-column TSV, ",
         "or \"T2TCHM13\".", call. = FALSE)
  }

  # ── Derive file_type from file_suffix extension ──────────────────────────────
  ext       <- tolower(tools::file_ext(file_suffix))
  file_type <- switch(ext, bw = "bw", bigwig = "bw", bam = "bam", bed = "bed",
                      stop("Cannot determine file type from suffix '", file_suffix,
                           "'. Extension must be one of: bam, bed, bw.",
                           call. = FALSE))

  # ── External dependency checks ───────────────────────────────────────────────
  if (file_type == "bam" && !nzchar(Sys.which("samtools")))
    stop("'samtools' not found on PATH (required for BAM chromosome renaming). ",
         "Install via: conda install -c bioconda samtools", call. = FALSE)

  if (file_type == "bw" && !requireNamespace("rtracklayer", quietly = TRUE))
    stop("Package 'rtracklayer' is required for BigWig renaming. ",
         "Install via: BiocManager::install(\"rtracklayer\")", call. = FALSE)

  # ── Locate matching files ────────────────────────────────────────────────────
  row      <- .get_sample_row(sample_sheet, sample_id)
  root_dir <- file.path(row$output_dir, sample_id)

  if (!dir.exists(root_dir))
    stop("Sample output directory not found: ", root_dir, call. = FALSE)

  all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  targets   <- all_files[endsWith(basename(all_files), file_suffix)]

  if (length(targets) == 0)
    stop("No files ending with '", file_suffix, "' found under: ", root_dir,
         call. = FALSE)

  cat("Found ", length(targets), " file(s) matching '", file_suffix, "'.")

  # ── Process each file ────────────────────────────────────────────────────────
  out_paths <- character(length(targets))

  for (i in seq_along(targets)) {
    in_path  <- targets[[i]]
    out_path <- if (isTRUE(overwrite)) {
      in_path
    } else {
      ext  <- tools::file_ext(in_path)
      base <- tools::file_path_sans_ext(in_path)
      if (file.exists(paste0(base, rename_suffix, ".", ext)))
        warning("Output already exists and will be overwritten: ",
                paste0(base, rename_suffix, ".", ext))
      paste0(base, rename_suffix, ".", ext)
    }

    cat("Processing: ", basename(in_path))

    switch(file_type,
      bam = .rename_chr_bam(in_path, out_path, chr_map),
      bed = .rename_chr_bed(in_path, out_path, chr_map),
      bw  = .rename_chr_bw(in_path,  out_path, chr_map)
    )

    out_paths[[i]] <- out_path
  }

  invisible(out_paths)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.rename_chr_bam <- function(in_path, out_path, chr_map) {
  # Extract current header
  hdr <- system2("samtools", args = c("view", "-H", in_path), stdout = TRUE)

  # Rewrite SN: fields in @SQ lines only
  for (old in names(chr_map)) {
    escaped <- gsub(".", "\\.", old, fixed = TRUE)
    pattern <- paste0("(SN:)", escaped, "(?=\\t|$)")
    hdr     <- gsub(pattern, paste0("\\1", chr_map[[old]]), hdr, perl = TRUE)
  }

  n_renamed <- sum(grepl("^@SQ", hdr))  # approximate; real count from mapping
  if (n_renamed == 0)
    warning("No @SQ lines found in BAM header: ", basename(in_path))

  tmp_hdr <- tempfile(fileext = ".hdr")
  on.exit(unlink(tmp_hdr), add = TRUE)
  writeLines(hdr, tmp_hdr)

  # samtools reheader writes to stdout; use temp if overwriting in place
  write_to <- if (identical(in_path, out_path)) tempfile(fileext = ".bam") else out_path

  exit_code <- system2("samtools",
                       args   = c("reheader", tmp_hdr, in_path),
                       stdout = write_to)

  if (exit_code != 0)
    stop("samtools reheader failed (exit ", exit_code, "): ",
         basename(in_path), call. = FALSE)

  if (identical(in_path, out_path)) {
    file.rename(write_to, out_path)
    old_bai <- paste0(in_path, ".bai")
    if (file.exists(old_bai)) unlink(old_bai)
  }

  cat("  Indexing: ", basename(out_path))
  Rsamtools::indexBam(out_path)

  invisible(out_path)
}


.rename_chr_bed <- function(in_path, out_path, chr_map) {
  lines <- readLines(in_path)

  renamed <- vapply(lines, function(line) {
    # Pass through empty lines and BED/browser/track header lines unchanged
    if (!nzchar(line) || grepl("^#|^browser |^track ", line))
      return(line)
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (fields[1] %in% names(chr_map))
      fields[1] <- chr_map[[fields[1]]]
    paste(fields, collapse = "\t")
  }, character(1L), USE.NAMES = FALSE)

  write_to <- if (identical(in_path, out_path)) tempfile() else out_path
  writeLines(renamed, write_to)
  if (identical(in_path, out_path)) file.rename(write_to, out_path)

  invisible(out_path)
}


.rename_chr_bw <- function(in_path, out_path, chr_map) {
  gr      <- rtracklayer::import(in_path, format = "BigWig")
  present <- GenomeInfoDb::seqlevels(gr)
  relevant <- chr_map[names(chr_map) %in% present]

  if (length(relevant) == 0)
    warning("No chromosome names in chr_map matched those in: ", basename(in_path))

  if (length(relevant) > 0)
    gr <- GenomeInfoDb::renameSeqlevels(gr, relevant)

  write_to <- if (identical(in_path, out_path)) tempfile(fileext = ".bw") else out_path
  rtracklayer::export(gr, write_to, format = "BigWig")
  if (identical(in_path, out_path)) file.rename(write_to, out_path)

  invisible(out_path)
}
