# ==============================================================================
# ELEUTHIA - BigWig Generation
# ==============================================================================
# Functions for converting fragment-level BED files to BigWig coverage tracks.
# Not specific to any single assay type — applicable to any sequencing data
# where fragment-level BED files are available.
# ==============================================================================


#' Convert Fragment-Level BED to BigWig
#'
#' @description Converts a fragment-level BED file (one row per fragment) into
#' a BigWig coverage file suitable for genome browser visualization or use with
#' [AETHER_plot_coverage_tracks()]. Coverage is computed at a configurable bin
#' resolution and optionally normalized to CPM.
#'
#' @param bed_file Character. Path to the fragment-level BED file (0-based
#'   half-open coordinates: chr, start, end).
#' @param output_path Character. Output file path. Extension should be `.bw`
#'   or `.bigwig`; format is inferred by `rtracklayer` from the extension.
#' @param chrom_sizes Exact chromosome lengths. Accepts either:
#'   - A named numeric vector: `c(chr1 = 248956422L, chr2 = 242193529L)`
#'   - A data.frame with columns `chr` and `size`, as returned by
#'     `APOLLO_get_chromosome_sizes()`.
#'   Names/values must match chromosome names in the BED file exactly.
#'   No inferring from data.
#' @param normalize Character. Normalization method: `"CPM"` (default) scales
#'   signal by counts per million total fragments; `"raw"` leaves coverage as
#'   absolute fragment counts.
#' @param bin_size Integer. Bin width in base pairs for coverage aggregation.
#'   Larger values produce smoother signal and smaller files. Default `10`.
#' @param chromosomes Character vector. Restrict processing to these
#'   chromosomes. Default `NULL` processes all chromosomes present in both
#'   the BED file and `chrom_sizes`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return Invisibly returns `output_path`.
#'
#' @details
#' BED files use 0-based half-open coordinates. These are converted to
#' 1-based closed intervals for GRanges internally.
#'
#' Requires Bioconductor packages: `GenomicRanges`, `IRanges`,
#' `GenomeInfoDb`, `rtracklayer`. Install via `BiocManager::install()`.
#'
#' @export
ELEUTHIA_bed_to_bigwig <- function(bed_file,
                                   output_path,
                                   chrom_sizes,
                                   normalize  = "CPM",
                                   bin_size   = 10L,
                                   chromosomes = NULL,
                                   verbose    = TRUE) {

  # --- Input validation -------------------------------------------------------
  if (!file.exists(bed_file))
    stop("BED file not found: ", bed_file, call. = FALSE)

  # Accept APOLLO_get_chromosome_sizes data.frame (chr, size) or named vector
  if (is.data.frame(chrom_sizes)) {
    if (!all(c("chr", "size") %in% colnames(chrom_sizes)))
      stop("chrom_sizes data.frame must have columns 'chr' and 'size'",
           call. = FALSE)
    chrom_sizes <- setNames(as.numeric(chrom_sizes$size), chrom_sizes$chr)
  }

  if (!is.numeric(chrom_sizes) || is.null(names(chrom_sizes)))
    stop("chrom_sizes must be a named numeric vector or a data.frame with ",
         "'chr' and 'size' columns (e.g. from APOLLO_get_chromosome_sizes())",
         call. = FALSE)
  normalize <- match.arg(normalize, c("CPM", "raw"))
  bin_size  <- as.integer(bin_size)
  if (bin_size < 1L)
    stop("bin_size must be >= 1", call. = FALSE)

  # --- Package availability ---------------------------------------------------
  for (pkg in c("GenomicRanges", "IRanges", "GenomeInfoDb", "rtracklayer")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. ",
           "Install with: BiocManager::install('", pkg, "')",
           call. = FALSE)
  }

  # --- Read BED ---------------------------------------------------------------
  if (verbose) cat("[ELEUTHIA] Reading: ", basename(bed_file))

  if (requireNamespace("data.table", quietly = TRUE)) {
    raw <- data.table::fread(bed_file, header = FALSE, select = 1:3,
                             showProgress = FALSE, data.table = FALSE)
  } else {
    raw <- utils::read.table(bed_file, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)[, 1:3]
  }

  if (ncol(raw) < 3)
    stop("BED file must have at least 3 columns (chr, start, end)", call. = FALSE)

  colnames(raw) <- c("chr", "start", "end")
  raw$start <- as.integer(raw$start)
  raw$end   <- as.integer(raw$end)

  # --- Chromosome filtering ---------------------------------------------------
  valid_chrs <- intersect(names(chrom_sizes), unique(raw$chr))
  if (!is.null(chromosomes))
    valid_chrs <- intersect(valid_chrs, chromosomes)

  if (length(valid_chrs) == 0)
    stop("No chromosomes remain after filtering. ",
         "Check chromosome name compatibility between BED and chrom_sizes.",
         call. = FALSE)

  bed <- raw[raw$chr %in% valid_chrs, , drop = FALSE]
  total_fragments <- nrow(bed)

  if (verbose)
    cat("[ELEUTHIA] Summary\n")
    cat("    Fragments : ", format(total_fragments, big.mark = ","),
            "\n    Chromosomes : ", length(valid_chrs),
            "\n    Normalize : ", normalize,
            "\n    Bin size  : ", bin_size, " bp")

  # --- Build Seqinfo ----------------------------------------------------------
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames   = valid_chrs,
    seqlengths = as.integer(chrom_sizes[valid_chrs])
  )

  # --- GRanges (BED 0-based half-open → GRanges 1-based closed) --------------
  fragments_gr <- GenomicRanges::GRanges(
    seqnames = bed$chr,
    ranges   = IRanges::IRanges(start = bed$start + 1L, end = bed$end),
    seqinfo  = seqinfo
  )

  # --- Coverage ---------------------------------------------------------------
  if (verbose) cat("[ELEUTHIA]   Computing coverage...")
  cov_rle <- GenomicRanges::coverage(fragments_gr)

  # Tile genome into bins and compute per-bin average
  bins   <- GenomicRanges::tileGenome(
    seqlengths             = GenomeInfoDb::seqlengths(seqinfo),
    tilewidth              = bin_size,
    cut.last.tile.in.chrom = TRUE
  )
  binned <- GenomicRanges::binnedAverage(bins, cov_rle, "score")

  # --- Normalize --------------------------------------------------------------
  if (normalize == "CPM") {
    cpm_factor    <- 1e6 / total_fragments
    binned$score  <- binned$score * cpm_factor
    if (verbose)
      cat("[ELEUTHIA]   CPM factor: ", round(cpm_factor, 6))
  }

  # Drop zero-score bins (reduces file size meaningfully for sparse data)
  binned <- binned[binned$score > 0]

  # --- Export BigWig ----------------------------------------------------------
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  if (verbose) cat("[ELEUTHIA]   Writing: ", basename(output_path))
  rtracklayer::export.bw(binned, output_path)

  if (verbose) cat("    Done.")
  invisible(output_path)
}


#' Compute Coverage for a Genomic Region Directly from a BED File
#'
#' @description Extracts fragments overlapping a specific genomic region from
#' a fragment-level BED file and computes binned coverage — without generating
#' a full genome-wide BigWig file. Intended as the efficient backend for
#' [AETHER_plot_coverage_tracks()] when BigWig files are not pre-computed.
#'
#' @param bed_file Character. Path to the fragment-level BED file (0-based
#'   half-open coordinates: chr, start, end).
#' @param region Genomic region. One of:
#'   - Character string: `"chr1:1000000-2000000"`
#'   - Named character vector: `c(chr="chr1", start="1000000", end="2000000")`
#'   - A `GRanges` object (single range)
#' @param chrom_sizes Exact chromosome lengths. Accepts either:
#'   - A named numeric vector: `c(chr1 = 248956422L, chr2 = 242193529L)`
#'   - A data.frame with columns `chr` and `size`, as returned by
#'     `APOLLO_get_chromosome_sizes()`.
#' @param normalize Character. `"CPM"` (default) scales by counts per million
#'   total fragments in the file; `"raw"` returns absolute counts.
#' @param bin_size Integer. Bin width in base pairs. Default `10`.
#' @param verbose Logical. Print progress messages. Default `FALSE`.
#'
#' @return A data.frame with columns `start`, `end`, `score`, `midpoint`
#'   covering the requested region. Same format as `rtracklayer::import.bw()`
#'   output, making it interchangeable in downstream plotting code.
#'
#' @details
#' The entire BED file is read once (using `data.table`) to obtain both the
#' total fragment count (for CPM) and the region-filtered subset. Only bins
#' within the requested region are computed — no genome-wide tiling.
#'
#' Requires: `data.table`, `GenomicRanges`, `IRanges`, `GenomeInfoDb`.
#'
#' @export
ELEUTHIA_bed_region_coverage <- function(bed_file,
                                         region,
                                         chrom_sizes,
                                         normalize = "CPM",
                                         bin_size  = 10L,
                                         verbose   = FALSE) {

  # --- Input validation -------------------------------------------------------
  if (!file.exists(bed_file))
    stop("BED file not found: ", bed_file, call. = FALSE)

  if (is.data.frame(chrom_sizes)) {
    if (!all(c("chr", "size") %in% colnames(chrom_sizes)))
      stop("chrom_sizes data.frame must have columns 'chr' and 'size'", call. = FALSE)
    chrom_sizes <- setNames(as.numeric(chrom_sizes$size), chrom_sizes$chr)
  }
  if (!is.numeric(chrom_sizes) || is.null(names(chrom_sizes)))
    stop("chrom_sizes must be a named numeric vector or APOLLO_get_chromosome_sizes() data.frame",
         call. = FALSE)

  normalize <- match.arg(normalize, c("CPM", "raw"))
  bin_size  <- as.integer(bin_size)
  if (bin_size < 1L) stop("bin_size must be >= 1", call. = FALSE)

  # --- Package checks ---------------------------------------------------------
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' required. Install with: install.packages('data.table')",
         call. = FALSE)
  for (pkg in c("GenomicRanges", "IRanges", "GenomeInfoDb")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' required. Install with: BiocManager::install('", pkg, "')",
           call. = FALSE)
  }

  # --- Parse region -----------------------------------------------------------
  reg       <- .ELEUTHIA_parse_region(region)
  chr       <- reg$chr
  reg_start <- reg$start
  reg_end   <- reg$end

  if (!chr %in% names(chrom_sizes))
    stop("Chromosome '", chr, "' not found in chrom_sizes", call. = FALSE)

  if (verbose)
    cat("[ELEUTHIA] Region coverage: ", chr, ":", reg_start, "-", reg_end,
            " from ", basename(bed_file))

  # --- Read full BED (one pass for total count + region filter) ---------------
  all_frags <- data.table::fread(bed_file, header = FALSE, select = 1:3,
                                 showProgress = FALSE, data.table = FALSE)
  if (ncol(all_frags) < 3)
    stop("BED file must have at least 3 columns (chr, start, end)", call. = FALSE)

  colnames(all_frags) <- c("chr", "start", "end")
  total_fragments     <- nrow(all_frags)

  region_frags <- all_frags[
    all_frags$chr   == chr       &
      all_frags$end   >= reg_start &
      all_frags$start <= reg_end,
    , drop = FALSE
  ]

  if (verbose)
    cat("[ELEUTHIA]   Total: ", format(total_fragments, big.mark = ","),
            " | In region: ", format(nrow(region_frags), big.mark = ","))

  # Return empty frame if no fragments in region
  if (nrow(region_frags) == 0L) {
    warning("No fragments found in region ", chr, ":", reg_start, "-", reg_end,
            call. = FALSE)
    return(data.frame(start    = integer(0), end      = integer(0),
                      score    = numeric(0), midpoint = numeric(0)))
  }

  # --- Build coverage for the region only -------------------------------------
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames   = chr,
    seqlengths = as.integer(chrom_sizes[[chr]])
  )

  fragments_gr <- GenomicRanges::GRanges(
    seqnames = region_frags$chr,
    ranges   = IRanges::IRanges(start = region_frags$start + 1L,
                                end   = region_frags$end),
    seqinfo  = seqinfo
  )

  cov_rle <- GenomicRanges::coverage(fragments_gr)

  # Tile only the requested region (not the full chromosome)
  bin_starts <- seq(reg_start, reg_end, by = bin_size)
  bin_ends   <- pmin(bin_starts + bin_size - 1L, reg_end)
  bins <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = bin_starts, end = bin_ends),
    seqinfo  = seqinfo
  )

  binned <- GenomicRanges::binnedAverage(bins, cov_rle, "score")

  # --- Normalize --------------------------------------------------------------
  if (normalize == "CPM")
    binned$score <- binned$score * (1e6 / total_fragments)

  # --- Return data.frame matching rtracklayer::import.bw() format -------------
  df          <- as.data.frame(binned)[, c("start", "end", "score"), drop = FALSE]
  df$midpoint <- (df$start + df$end) / 2
  df
}


# Internal: parse region to a simple list(chr, start, end) — no GRanges needed
.ELEUTHIA_parse_region <- function(region) {

  if (inherits(region, "GRanges")) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("GenomicRanges required to use a GRanges region", call. = FALSE)
    if (length(region) != 1L)
      stop("region GRanges must contain exactly one range", call. = FALSE)
    return(list(
      chr   = as.character(GenomicRanges::seqnames(region)),
      start = GenomicRanges::start(region),
      end   = GenomicRanges::end(region)
    ))
  }

  if (is.character(region) && length(region) == 1L && is.null(names(region))) {
    parts <- regmatches(region,
                        regexec("^([^:]+):(\\d+)-(\\d+)$", region))[[1]]
    if (length(parts) != 4L)
      stop("Cannot parse region '", region, "'. Expected 'chr:start-end'",
           call. = FALSE)
    return(list(chr   = parts[2],
                start = as.integer(parts[3]),
                end   = as.integer(parts[4])))
  }

  if (is.character(region) && !is.null(names(region))) {
    if (!all(c("chr", "start", "end") %in% names(region)))
      stop("Named region vector must have 'chr', 'start', 'end' elements",
           call. = FALSE)
    return(list(chr   = region[["chr"]],
                start = as.integer(region[["start"]]),
                end   = as.integer(region[["end"]])))
  }

  stop("region must be a GRanges, 'chr:start-end' string, ",
       "or named character vector with chr/start/end", call. = FALSE)
}
