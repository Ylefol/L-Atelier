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
  bw_col      = "bw_path"
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

  # ---- Resolve BigWig paths ---------------------------------------------------
  idx       <- match(all_samples, sheet_ids)
  bw_paths  <- setNames(as.character(sample_sheet[[bw_col]])[idx], all_samples)
  bad_files <- names(bw_paths)[!file.exists(bw_paths)]
  if (length(bad_files))
    stop("[ARTEMIS] BigWig file(s) not found on disk: ",
         paste(bad_files, collapse = ", "), call. = FALSE)

  # ---- Define regions ---------------------------------------------------------
  cat("[ARTEMIS] Defining regions...\n")
  if (inherits(regions, "GRanges")) {
    query_regions <- regions
    cat("    GRanges object:", length(query_regions), "regions\n")

  } else if (is.character(regions) && length(regions) == 1L && regions == "bins") {
    genome_tag <- if (!is.null(genome)) paste0(" (", genome, ")") else ""
    cat("    Mode: genome-wide", bin_size, "bp bins", genome_tag, "\n")
    ref_bw    <- rtracklayer::import.bw(bw_paths[[1]])
    si_len    <- GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(ref_bw))
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
    query_regions <- GenomicRanges::GRanges(
      seqnames = bed[[1]],
      ranges   = IRanges::IRanges(
        start = as.integer(bed[[2]]) + 1L,  # BED 0-based -> 1-based
        end   = as.integer(bed[[3]])
      )
    )
    cat("    Regions loaded:", length(query_regions), "\n")

  } else {
    stop('[ARTEMIS] "regions" must be "bins", a BED file path, or a GRanges object',
         call. = FALSE)
  }

  # ---- Extract mean signal per BigWig -----------------------------------------
  n_regions  <- length(query_regions)
  n_samples  <- length(bw_paths)
  signal_mat <- matrix(0, nrow = n_regions, ncol = n_samples,
                       dimnames = list(NULL, names(bw_paths)))
  cat("[ARTEMIS] Extracting signal (", n_samples, " file(s), ",
      n_regions, " regions)\n", sep = "")

  rchrs   <- as.character(GenomicRanges::seqnames(query_regions))
  rstarts <- GenomicRanges::start(query_regions)
  rends   <- GenomicRanges::end(query_regions)

  for (samp in names(bw_paths)) {
    cat("    ", samp, "...\n", sep = "")
    bw_gr   <- rtracklayer::import.bw(bw_paths[[samp]])
    cov_rle <- GenomicRanges::coverage(bw_gr, weight = "score")

    for (chr in names(cov_rle)) {
      chr_mask <- rchrs == chr
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
