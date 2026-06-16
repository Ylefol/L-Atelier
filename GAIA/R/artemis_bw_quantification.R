#' Prepare BigWig Signal Data for 2D Density Scatter Comparison
#'
#' Extracts mean BigWig signal over genomic regions for two groups, handles
#' replicate averaging, optional log2 ratio computation, and low-signal
#' filtering. Designed to feed directly into
#' \code{\link{AETHER_plot_density_scatter}}.
#'
#' @param sample_sheet data.frame. Must include columns identified by
#'   \code{sample_col} and \code{bw_col}.
#' @param group_x Named list defining the x-axis. See Details.
#' @param group_y Named list defining the y-axis. See Details.
#' @param regions One of: \code{"bins"} for genome-wide fixed-width bins,
#'   a path to a BED file (3-column, 0-based half-open), or a
#'   \code{GRanges} object.
#' @param bin_size Integer. Bin width in bp. Only used when
#'   \code{regions = "bins"}. Default 1000.
#' @param genome Character or \code{NULL}. Genome identifier (e.g.
#'   \code{"hg38"}). Used for logging when \code{regions = "bins"}.
#'   Chromosome sizes are derived from the BigWig seqinfo; this parameter
#'   does not load a BSgenome package. Default \code{NULL}.
#' @param ratio Logical. If \code{TRUE}, compute log2 ratio for any axis
#'   whose list has exactly 2 names. Default \code{FALSE}.
#' @param pseudocount Numeric. Added to group means before log2 ratio:
#'   \code{log2((numerator + pseudocount) / (denominator + pseudocount))}.
#'   Prevents \code{log(0)} and dampens instability at near-zero signal.
#'   Symmetric — does not introduce directional bias. Default 1.
#' @param min_signal Numeric or \code{NULL}. Regions where the representative
#'   signal on \emph{both} axes falls below this threshold are removed before
#'   ratio computation. Representative = single group mean (1-name list) or
#'   \code{max(group1, group2)} (2-name list). Regions enriched in only one
#'   group are retained. Default \code{NULL} (no filtering).
#' @param sample_col Character. Column in \code{sample_sheet} holding sample
#'   IDs. Default \code{"sample_id"}.
#' @param bw_col Character. Column in \code{sample_sheet} holding BigWig file
#'   paths. Default \code{"bw_path"}.
#' @param coverage \code{NULL} or an \code{artemis_bw_coverage} object returned
#'   by \code{\link{ARTEMIS_load_bw_coverage}}. When provided, BigWig files are
#'   not read from disk — signal is extracted from the pre-loaded RLE coverage
#'   objects. Speeds up multi-region workflows by loading each file only once.
#'   Default \code{NULL}.
#'
#' @details
#' \strong{Input list logic:}\cr
#' Each of \code{group_x} and \code{group_y} is a named list mapping group
#' labels to character vectors of sample IDs:
#' \itemize{
#'   \item \strong{1 name} — direct mean signal across all samples.
#'         E.g. \code{list(CTRL = c("s1", "s2"))}.
#'   \item \strong{2 names + \code{ratio = TRUE}} — log2 ratio. The
#'         \emph{first name is the denominator}, the \emph{second is the
#'         numerator}: \code{log2((2nd + pc) / (1st + pc))}.
#'         E.g. \code{list(CTRL = c("s1","s2"), TREAT = c("s3","s4"))}
#'         yields \code{log2((TREAT + pc) / (CTRL + pc))}.
#' }
#'
#' Mixed axes are valid (one direct, one ratio).
#'
#' \strong{Order of operations:}
#' validate inputs → resolve paths → define regions → extract mean signal
#' per BigWig → average replicates per group → filter (min_signal) →
#' compute log2 ratio → return data frame.
#'
#' @return data.frame with columns \code{x}, \code{y}, \code{chr},
#'   \code{start}, \code{end} (one row per region). Attributes
#'   \code{x_label} and \code{y_label} store descriptive axis labels
#'   that \code{AETHER_plot_density_scatter} uses as defaults.
#'
#' @export
ARTEMIS_prepare_bw_comparison <- function(
  sample_sheet,
  group_x,
  group_y,
  regions,
  bin_size    = 1000L,
  genome      = NULL,
  ratio       = FALSE,
  pseudocount = 1,
  min_signal  = NULL,
  sample_col  = "sample_id",
  bw_col      = "bw_path",
  coverage    = NULL
) {
  # ---- Input validation -------------------------------------------------------
  if (!is.data.frame(sample_sheet))
    stop("[ARTEMIS] sample_sheet must be a data.frame", call. = FALSE)
  if (!sample_col %in% colnames(sample_sheet))
    stop("[ARTEMIS] sample_col '", sample_col, "' not found in sample_sheet",
         call. = FALSE)
  if (!bw_col %in% colnames(sample_sheet))
    stop("[ARTEMIS] bw_col '", bw_col, "' not found in sample_sheet",
         call. = FALSE)

  .check_group <- function(g, gname) {
    if (!is.list(g) || is.null(names(g)) || any(nchar(names(g)) == 0L))
      stop("[ARTEMIS] ", gname, " must be a named list ",
           "(e.g. list(CTRL = c('s1', 's2')))", call. = FALSE)
    n <- length(g)
    if (!n %in% c(1L, 2L))
      stop("[ARTEMIS] ", gname, " must have 1 or 2 names (got ", n, ")",
           call. = FALSE)
    if (n == 2L && !ratio)
      stop("[ARTEMIS] ", gname, " has 2 names but ratio = FALSE. ",
           "Set ratio = TRUE to compute a log2 ratio on this axis.", call. = FALSE)
  }
  .check_group(group_x, "group_x")
  .check_group(group_y, "group_y")

  all_samples <- unique(unlist(c(group_x, group_y)))
  sheet_ids   <- as.character(sample_sheet[[sample_col]])
  missing_ids <- setdiff(all_samples, sheet_ids)
  if (length(missing_ids))
    stop("[ARTEMIS] Sample ID(s) not found in sample_sheet: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)

  # ---- Validate or resolve coverage/paths ------------------------------------
  if (!is.null(coverage)) {
    if (!inherits(coverage, "artemis_bw_coverage"))
      stop("[ARTEMIS] coverage must be an object returned by ARTEMIS_load_bw_coverage()",
           call. = FALSE)
    missing_cov <- setdiff(all_samples, names(coverage$coverage))
    if (length(missing_cov))
      stop("[ARTEMIS] Sample(s) not found in coverage object: ",
           paste(missing_cov, collapse = ", "), call. = FALSE)
  } else {
    idx       <- match(all_samples, sheet_ids)
    bw_paths  <- setNames(as.character(sample_sheet[[bw_col]])[idx], all_samples)
    bad_files <- names(bw_paths)[!file.exists(bw_paths)]
    if (length(bad_files))
      stop("[ARTEMIS] BigWig file(s) not found on disk: ",
           paste(bad_files, collapse = ", "), call. = FALSE)
  }

  # ---- Define regions ---------------------------------------------------------
  cat("[ARTEMIS] Defining regions...\n")
  if (inherits(regions, "GRanges")) {
    query_regions <- regions
    cat("    GRanges object:", length(query_regions), "regions\n")

  } else if (is.character(regions) && length(regions) == 1L && regions == "bins") {
    genome_tag <- if (!is.null(genome)) paste0(" (", genome, ")") else ""
    cat("    Mode: genome-wide", bin_size, "bp bins", genome_tag, "\n")
    if (!is.null(coverage)) {
      si_len <- lengths(coverage$coverage[[all_samples[1]]])
    } else {
      ref_bw <- rtracklayer::import.bw(bw_paths[[1]])
      si_len <- GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(ref_bw))
    }
    keep_chrs <- !grepl("_", names(si_len))  # drop alt/random/patch contigs
    si_len    <- si_len[keep_chrs]
    if (!length(si_len))
      stop("[ARTEMIS] No canonical chromosomes found in BigWig seqinfo",
           call. = FALSE)
    cat("    Chromosomes:", length(si_len),
        "(", paste(utils::head(names(si_len), 5), collapse = ", "),
        if (length(si_len) > 5) "..." else "", ")\n")
    tiles         <- GenomicRanges::tileGenome(
      si_len,
      tilewidth              = bin_size,
      cut.last.tile.in.chrom = TRUE
    )
    query_regions <- as(tiles, "GRanges")
    cat("    Total bins:", length(query_regions), "\n")

  } else if (is.character(regions) && length(regions) == 1L) {
    if (!file.exists(regions))
      stop("[ARTEMIS] BED file not found: ", regions, call. = FALSE)
    cat("    Reading BED file:", basename(regions), "\n")
    bed <- read.table(regions, header = FALSE, stringsAsFactors = FALSE,
                      comment.char = "#")
    if (ncol(bed) < 3L)
      stop("[ARTEMIS] BED file must have at least 3 columns (chr, start, end)",
           call. = FALSE)
    # Drop non-data rows: track/browser headers or column-name lines
    # (rows where column 2 cannot be parsed as integer)
    keep_rows <- !is.na(suppressWarnings(as.integer(bed[[2]])))
    if (!all(keep_rows)) {
      cat("    Skipping", sum(!keep_rows), "non-data row(s) (header/track lines)\n")
      bed <- bed[keep_rows, ]
    }
    if (nrow(bed) == 0L)
      stop("[ARTEMIS] BED file contains no valid data rows", call. = FALSE)
    # Detect column layout: standard BED (chr|start|end) vs name-first BED (name|chr|start|end)
    # Chromosome names are short: start with "chr", are plain digits, or are X/Y/MT/M.
    # Peak names (e.g. "K27me3_C1_peak_1") fail all of those criteria.
    sample_col1 <- as.character(utils::head(bed[[1]], min(20L, nrow(bed))))
    is_chr_col1 <- all(grepl("^(chr|[0-9]{1,3}$|X$|Y$|MT?$|Un)", sample_col1))
    if (!is_chr_col1 && ncol(bed) >= 4L) {
      cat("    Name-first BED format detected: using columns 2-4 for chr/start/end\n")
      chr_col   <- 2L
      start_col <- 3L
      end_col   <- 4L
    } else {
      chr_col   <- 1L
      start_col <- 2L
      end_col   <- 3L
    }
    # Re-validate: start column must be numeric after format detection
    keep_rows <- !is.na(suppressWarnings(as.integer(bed[[start_col]])))
    if (!all(keep_rows)) {
      bed <- bed[keep_rows, ]
    }
    query_regions <- GenomicRanges::GRanges(
      seqnames = bed[[chr_col]],
      ranges   = IRanges::IRanges(
        start = as.integer(bed[[start_col]]) + 1L,  # BED 0-based -> 1-based
        end   = as.integer(bed[[end_col]])
      )
    )
    cat("    Regions loaded:", length(query_regions), "\n")

  } else {
    stop('[ARTEMIS] "regions" must be "bins", a BED file path, or a GRanges object',
         call. = FALSE)
  }

  # ---- Extract mean signal per BigWig -----------------------------------------
  n_regions  <- length(query_regions)
  signal_mat <- matrix(0, nrow = n_regions, ncol = length(all_samples),
                       dimnames = list(NULL, all_samples))
  cat("[ARTEMIS] Extracting signal (", length(all_samples), " file(s), ",
      n_regions, " regions)\n", sep = "")

  rchrs   <- as.character(GenomicRanges::seqnames(query_regions))
  rstarts <- GenomicRanges::start(query_regions)
  rends   <- GenomicRanges::end(query_regions)

  # ---- Chromosome name harmonization ------------------------------------------
  # For 'bins' mode, regions are tiled directly from BigWig seqinfo — no mismatch.
  # For BED/GRanges modes, detect and auto-fix chr1 vs 1 naming discrepancies.
  extract_rchrs <- rchrs
  if (!(is.character(regions) && length(regions) == 1L && regions == "bins")) {
    bw_chrs <- if (!is.null(coverage)) {
      coverage$bw_chrnames
    } else {
      names(GenomeInfoDb::seqlengths(
        GenomicRanges::seqinfo(rtracklayer::BigWigFile(bw_paths[[1]]))
      ))
    }
    n_overlap <- sum(unique(rchrs) %in% bw_chrs)
    if (n_overlap == 0L) {
      reg_has_chr <- any(grepl("^chr", rchrs))
      bw_has_chr  <- any(grepl("^chr", bw_chrs))
      if (reg_has_chr && !bw_has_chr) {
        extract_rchrs <- sub("^chr", "", rchrs)
        cat("[ARTEMIS] Chr name harmonization: stripped 'chr' prefix to match BigWig (NCBI style)\n")
      } else if (!reg_has_chr && bw_has_chr) {
        extract_rchrs <- paste0("chr", rchrs)
        cat("[ARTEMIS] Chr name harmonization: added 'chr' prefix to match BigWig (UCSC style)\n")
      } else {
        warning(
          "[ARTEMIS] No chromosome overlap between regions and BigWig files.",
          " Signal will be all zeros.\n",
          "  Region chrs (first 5): ", paste(utils::head(unique(rchrs), 5), collapse = ", "), "\n",
          "  BigWig chrs (first 5): ", paste(utils::head(bw_chrs, 5), collapse = ", ")
        )
      }
    }
  }

  for (samp in all_samples) {
    cat("    ", samp, "...\n", sep = "")
    cov_rle <- if (!is.null(coverage)) {
      coverage$coverage[[samp]]
    } else {
      bw_gr <- rtracklayer::import.bw(bw_paths[[samp]])
      GenomicRanges::coverage(bw_gr, weight = "score")
    }

    for (chr in names(cov_rle)) {
      chr_mask <- extract_rchrs == chr
      if (!any(chr_mask)) next
      chr_len  <- length(cov_rle[[chr]])
      s        <- pmin(pmax(rstarts[chr_mask], 1L), chr_len)
      e        <- pmin(rends[chr_mask],         chr_len)
      valid    <- s <= e
      if (!any(valid)) next
      views    <- IRanges::Views(cov_rle[[chr]], start = s[valid], end = e[valid])
      means    <- IRanges::viewMeans(views)
      row_idx  <- which(chr_mask)[valid]
      signal_mat[row_idx, samp] <- means
    }
  }

  # ---- Compute per-group means ------------------------------------------------
  .grp_mean <- function(samples) {
    if (length(samples) == 1L) signal_mat[, samples, drop = TRUE]
    else rowMeans(signal_mat[, samples, drop = FALSE])
  }

  x_single <- length(group_x) == 1L
  y_single <- length(group_y) == 1L

  if (x_single) {
    x_vals  <- .grp_mean(group_x[[1]])
    x_label <- paste0("Mean signal: ", names(group_x)[1])
  } else {
    x_den   <- .grp_mean(group_x[[1]])
    x_num   <- .grp_mean(group_x[[2]])
    x_label <- paste0("log2(", names(group_x)[2], " / ", names(group_x)[1], ")")
  }

  if (y_single) {
    y_vals  <- .grp_mean(group_y[[1]])
    y_label <- paste0("Mean signal: ", names(group_y)[1])
  } else {
    y_den   <- .grp_mean(group_y[[1]])
    y_num   <- .grp_mean(group_y[[2]])
    y_label <- paste0("log2(", names(group_y)[2], " / ", names(group_y)[1], ")")
  }

  # ---- min_signal filter ------------------------------------------------------
  if (!is.null(min_signal)) {
    cat("[ARTEMIS] Applying min_signal filter (threshold =", min_signal, ")\n")
    # Representative signal per axis: single mean, or max of two groups
    x_rep <- if (x_single) x_vals else pmax(x_den, x_num)
    y_rep <- if (y_single) y_vals else pmax(y_den, y_num)
    keep  <- !(x_rep < min_signal & y_rep < min_signal)
    cat("    Before:", sum(rep(TRUE, length(keep))), "| kept:", sum(keep),
        "| removed:", sum(!keep), "\n")
    query_regions <- query_regions[keep]
    if (x_single) x_vals <- x_vals[keep]
    else          { x_den <- x_den[keep]; x_num <- x_num[keep] }
    if (y_single) y_vals <- y_vals[keep]
    else          { y_den <- y_den[keep]; y_num <- y_num[keep] }
  }

  # ---- Log2 ratio (computed after filtering) ----------------------------------
  if (!x_single)
    x_vals <- log2((x_num + pseudocount) / (x_den + pseudocount))
  if (!y_single)
    y_vals <- log2((y_num + pseudocount) / (y_den + pseudocount))

  # ---- Assemble output --------------------------------------------------------
  out <- data.frame(
    x     = x_vals,
    y     = y_vals,
    chr   = as.character(GenomicRanges::seqnames(query_regions)),
    start = GenomicRanges::start(query_regions),
    end   = GenomicRanges::end(query_regions),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  attr(out, "x_label") <- x_label
  attr(out, "y_label") <- y_label

  cat("[ARTEMIS] Done.", nrow(out), "regions returned.\n")
  out
}


#' Pre-load BigWig Coverage for Fast Multi-Region Extraction
#'
#' Reads all BigWig files in a sample sheet once and stores their coverage as
#' RLE objects. The result can be passed to \code{\link{ARTEMIS_prepare_bw_comparison}}
#' via the \code{coverage} argument to skip redundant file I/O when extracting
#' signal over many region sets from the same samples.
#'
#' @param sample_sheet data.frame. Must include columns \code{sample_col} and
#'   \code{bw_col}. Duplicate sample IDs are deduplicated automatically.
#' @param sample_col Character. Column holding sample IDs. Default
#'   \code{"sample_id"}.
#' @param bw_col Character. Column holding BigWig file paths. Default
#'   \code{"bw_path"}.
#'
#' @return An \code{artemis_bw_coverage} object (a list) with:
#'   \describe{
#'     \item{\code{coverage}}{Named list of RleList objects, one per sample.}
#'     \item{\code{bw_chrnames}}{Character vector of chromosome names from the
#'       first BigWig, used for chromosome name harmonization in
#'       \code{ARTEMIS_prepare_bw_comparison}.}
#'   }
#'
#' @export
ARTEMIS_load_bw_coverage <- function(
  sample_sheet,
  sample_col = "sample_id",
  bw_col     = "bw_path"
) {
  if (!is.data.frame(sample_sheet))
    stop("[ARTEMIS] sample_sheet must be a data.frame", call. = FALSE)
  if (!sample_col %in% colnames(sample_sheet))
    stop("[ARTEMIS] sample_col '", sample_col, "' not found in sample_sheet",
         call. = FALSE)
  if (!bw_col %in% colnames(sample_sheet))
    stop("[ARTEMIS] bw_col '", bw_col, "' not found in sample_sheet",
         call. = FALSE)

  sample_ids <- as.character(sample_sheet[[sample_col]])
  bw_paths   <- setNames(as.character(sample_sheet[[bw_col]]), sample_ids)
  bw_paths   <- bw_paths[!duplicated(names(bw_paths))]

  bad_files <- names(bw_paths)[!file.exists(bw_paths)]
  if (length(bad_files))
    stop("[ARTEMIS] BigWig file(s) not found: ",
         paste(bad_files, collapse = ", "), call. = FALSE)

  cat("[ARTEMIS] Loading BigWig coverage (", length(bw_paths), " file(s))\n",
      sep = "")

  cov_list    <- vector("list", length(bw_paths))
  names(cov_list) <- names(bw_paths)
  bw_chrnames <- NULL

  for (samp in names(bw_paths)) {
    cat("    ", samp, "...\n", sep = "")
    bw_gr            <- rtracklayer::import.bw(bw_paths[[samp]])
    cov_list[[samp]] <- GenomicRanges::coverage(bw_gr, weight = "score")
    if (is.null(bw_chrnames))
      bw_chrnames <- names(cov_list[[samp]])
  }

  cat("[ARTEMIS] Coverage loaded for", length(cov_list), "sample(s).\n")
  structure(
    list(coverage = cov_list, bw_chrnames = bw_chrnames),
    class = "artemis_bw_coverage"
  )
}
