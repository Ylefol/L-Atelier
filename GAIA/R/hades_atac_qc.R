# ==============================================================================
# HADES - ATAC-seq QC functions
# ==============================================================================


#' Compute TSS Enrichment Score from BigWig files
#'
#' @description
#' Computes per-sample TSS enrichment scores from BigWig files. The TSS
#' enrichment score is a standard ATAC-seq QC metric (ENCODE guidelines) that
#' quantifies how specifically Tn5 inserted at open chromatin near transcription
#' start sites relative to flanking background signal.
#'
#' The score is computed by:
#' \enumerate{
#'   \item Aggregating normalised BigWig signal across all TSS windows into a
#'     single mean profile.
#'   \item Normalising by the mean signal in the outermost \code{flank_bins}
#'     bins on each side of the window (background estimate).
#'   \item Reporting the maximum of the normalised profile as the enrichment
#'     score (typically the centre bin over the TSS).
#' }
#'
#' ENCODE recommends a score >= 6 for human ATAC-seq data as a passing
#' threshold. This value was derived empirically and may not be universally
#' appropriate — verify for your cell type and coverage depth (Corces et al.,
#' 2018, Nature Methods, doi:10.1038/s41592-018-0105-7).
#'
#' @param sample_sheet data.frame. Must contain at minimum columns
#'   \code{sample_id} and \code{bw_loc} (path to the BigWig file for each
#'   sample). An optional \code{condition} column is carried through to the
#'   result and used for colouring in \code{\link{AETHER_plot_tss_enrichment}}.
#'   Exactly one of \code{sample_sheet} or \code{bw_file} must be provided.
#' @param bw_file Character. Path to a single BigWig file. Use when processing
#'   one sample without a sample sheet. Exactly one of \code{sample_sheet} or
#'   \code{bw_file} must be provided.
#' @param txdb TxDb object. Used to extract TSS positions. Create with
#'   \code{APOLLO_make_txdb()} or load a pre-built TxDb annotation package.
#' @param window Integer. Base pairs upstream and downstream of each TSS to
#'   include in the profile window. Default \code{3000L}.
#' @param bin_size Integer. Approximate width of each bin in bp. The window is
#'   divided into \code{2 * floor(window / bin_size) + 1} equal bins.
#'   Default \code{10L}.
#' @param flank_bins Integer. Number of bins at each edge of the window used to
#'   estimate background signal. Background is the mean of the outermost
#'   \code{flank_bins} on each side. At the default settings (bin_size = 10,
#'   flank_bins = 10) this represents 100 bp of background on each flank.
#'   Default \code{10L}.
#' @param threshold Numeric. Score threshold used to flag samples as passing
#'   QC. The ENCODE recommended value of 6 is a guideline for human data; no
#'   universally established threshold exists — treat this as a reference point,
#'   not an absolute standard. Default \code{6}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{hades_tss_enrichment} S3 object containing:
#'   \describe{
#'     \item{\code{$scores}}{data.frame with columns \code{sample_id},
#'       \code{tss_enrichment_score}, \code{pass} (score >= threshold), and
#'       \code{condition} if present in the sample sheet.}
#'     \item{\code{$profiles}}{Named list of normalised mean signal profiles
#'       (one numeric vector per sample, length = \code{n_bins}).}
#'     \item{\code{$params}}{List of run parameters: \code{window},
#'       \code{bin_size}, \code{flank_bins}, \code{n_bins}, \code{n_tss},
#'       \code{threshold}.}
#'   }
#'
#' @seealso \code{\link{AETHER_plot_tss_enrichment}}
#' @export
HADES_tss_enrichment <- function(
  sample_sheet = NULL,
  bw_file      = NULL,
  txdb,
  window       = 3000L,
  bin_size     = 10L,
  flank_bins   = 10L,
  threshold    = 6,
  verbose      = TRUE
) {

  # --- Package check -----------------------------------------------------------
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("Package 'GenomicFeatures' is required. ",
         "Install with: BiocManager::install('GenomicFeatures')", call. = FALSE)

  # --- Input validation --------------------------------------------------------
  if (!is.null(sample_sheet) && !is.null(bw_file))
    stop("Provide either 'sample_sheet' or 'bw_file', not both.", call. = FALSE)
  if (is.null(sample_sheet) && is.null(bw_file))
    stop("One of 'sample_sheet' or 'bw_file' must be provided.", call. = FALSE)

  window     <- as.integer(window)
  bin_size   <- as.integer(bin_size)
  flank_bins <- as.integer(flank_bins)

  # Number of bins: symmetric around TSS, centre bin sits on the TSS
  n_bins_half <- as.integer(window / bin_size)
  n_bins      <- 2L * n_bins_half + 1L

  if (flank_bins * 2L >= n_bins)
    stop("flank_bins (", flank_bins, ") must be less than half of n_bins (",
         n_bins, "). Reduce flank_bins or increase window / bin_size.",
         call. = FALSE)

  # --- Build sample table ------------------------------------------------------
  if (!is.null(sample_sheet)) {
    if (!is.data.frame(sample_sheet))
      stop("'sample_sheet' must be a data.frame.", call. = FALSE)
    if (!all(c("sample_id", "bw_loc") %in% colnames(sample_sheet)))
      stop("'sample_sheet' must contain columns 'sample_id' and 'bw_loc'.",
           call. = FALSE)
    # Carry all sample sheet columns through (except bw_loc, renamed to bw_path)
    meta_cols <- setdiff(colnames(sample_sheet), c("sample_id", "bw_loc"))
    samples   <- cbind(
      data.frame(sample_id = sample_sheet$sample_id,
                 bw_path   = sample_sheet$bw_loc,
                 stringsAsFactors = FALSE),
      sample_sheet[, meta_cols, drop = FALSE]
    )
  } else {
    if (!file.exists(bw_file))
      stop("BigWig file not found: ", bw_file, call. = FALSE)
    samples <- data.frame(
      sample_id = tools::file_path_sans_ext(basename(bw_file)),
      bw_path   = bw_file,
      stringsAsFactors = FALSE
    )
  }

  missing_bw <- samples$bw_path[!file.exists(samples$bw_path)]
  if (length(missing_bw) > 0L)
    stop("BigWig file(s) not found:\n  ",
         paste(missing_bw, collapse = "\n  "), call. = FALSE)

  # --- Extract TSS positions ---------------------------------------------------
  if (verbose) cat("[HADES] Extracting TSS positions from TxDb...\n")

  all_genes <- tryCatch(
    GenomicFeatures::genes(txdb),
    error = function(e)
      stop("Failed to retrieve genes from TxDb: ", conditionMessage(e),
           call. = FALSE)
  )

  si           <- GenomicRanges::seqinfo(all_genes)
  chr_lens     <- seqlengths(si)
  strands      <- as.character(GenomicRanges::strand(all_genes))
  starts       <- GenomicRanges::start(all_genes)
  ends         <- GenomicRanges::end(all_genes)
  chrs         <- as.character(GenomicRanges::seqnames(all_genes))
  tss_pos      <- ifelse(strands == "-", ends, starts)

  # Keep only TSSs with complete windows (not near chromosome boundaries)
  chr_len_per_gene <- chr_lens[chrs]
  complete <- !is.na(chr_len_per_gene) &
              (tss_pos - window) >= 1L &
              (tss_pos + window) <= chr_len_per_gene

  pos_df <- data.frame(
    gene_id = names(all_genes)[complete],
    chr     = chrs[complete],
    tss     = tss_pos[complete],
    stringsAsFactors = FALSE
  )

  n_tss <- nrow(pos_df)
  if (n_tss == 0L)
    stop("No TSS positions with complete windows found. ",
         "Check that chromosome names in the TxDb match those in the BigWig files.",
         call. = FALSE)

  if (verbose)
    cat("     TSSs in TxDb: ", length(all_genes),
        " | with complete windows: ", n_tss, "\n", sep = "")

  # --- Import regions (reduce overlapping windows to minimise BigWig queries) --
  tss_windows <- GenomicRanges::GRanges(
    seqnames = pos_df$chr,
    ranges   = IRanges::IRanges(
      start = pos_df$tss - window,
      end   = pos_df$tss + window
    )
  )
  import_gr <- GenomicRanges::reduce(tss_windows)

  # --- Per-sample processing ---------------------------------------------------
  score_vec    <- numeric(nrow(samples))
  profile_list <- vector("list", nrow(samples))
  names(profile_list) <- samples$sample_id

  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    bw  <- samples$bw_path[i]

    if (verbose) cat("[HADES] Processing: ", sid, "\n", sep = "")

    sig <- tryCatch(
      rtracklayer::import.bw(bw, which = import_gr),
      error = function(e)
        stop("Failed to read BigWig for '", sid, "': ",
             conditionMessage(e), call. = FALSE)
    )

    if (length(sig) == 0L) {
      warning("No signal found in BigWig for '", sid,
              "' — score set to NA.", call. = FALSE)
      score_vec[i]        <- NA_real_
      profile_list[[sid]] <- rep(NA_real_, n_bins)
      next
    }

    # Per-base weighted coverage from BigWig score column (RleList)
    sig_cov   <- GenomicRanges::coverage(sig, weight = "score")
    bw_chroms <- names(sig_cov)

    # Build signal matrix: rows = TSS, cols = bins
    # Strategy: for each chromosome, bin the Rle at bin_size resolution using
    # IRanges::viewMeans (compiled C — fast), then extract each TSS window from
    # the binned vector via vectorised matrix indexing.
    # This avoids creating millions of GRanges objects (tile + binnedAverage).
    score_mat <- matrix(0.0, nrow = n_tss, ncol = n_bins)

    for (chr in bw_chroms) {
      chr_tss_idx <- which(pos_df$chr == chr)
      if (length(chr_tss_idx) == 0L) next

      chr_rle <- sig_cov[[chr]]
      chr_len <- length(chr_rle)

      # Bin the chromosome Rle into bin_size-wide blocks
      bin_starts_chr <- seq(1L, chr_len, by = bin_size)
      bin_ends_chr   <- pmin(bin_starts_chr + bin_size - 1L, chr_len)
      chr_views      <- IRanges::Views(chr_rle,
                                        start = bin_starts_chr,
                                        end   = bin_ends_chr)
      bin_means      <- IRanges::viewMeans(chr_views)
      n_chr_bins     <- length(bin_means)

      # Map each TSS window start to its bin index in bin_means
      win_starts     <- pos_df$tss[chr_tss_idx] - window
      bin_start_idx  <- floor((win_starts - 1L) / bin_size) + 1L

      # Discard TSSs whose window falls outside the binned range
      in_range       <- bin_start_idx >= 1L &
                        (bin_start_idx + n_bins - 1L) <= n_chr_bins
      if (!any(in_range)) next

      valid_starts  <- bin_start_idx[in_range]
      valid_row_idx <- chr_tss_idx[in_range]

      # Vectorised extraction: outer() builds an index matrix where row i
      # contains the bin_means indices for TSS i's window
      col_indices <- outer(valid_starts, seq(0L, n_bins - 1L), "+")
      score_mat[valid_row_idx, ] <- matrix(bin_means[col_indices],
                                            nrow = length(valid_starts),
                                            ncol = n_bins)
    }

    mean_profile <- colMeans(score_mat, na.rm = TRUE)

    # Normalise by mean of outermost flank_bins on each side
    background <- mean(
      c(mean_profile[seq_len(flank_bins)],
        mean_profile[seq.int(n_bins - flank_bins + 1L, n_bins)]),
      na.rm = TRUE
    )

    if (is.na(background) || background <= 0) {
      warning("Background signal is zero or NA for '", sid,
              "'. Check that the BigWig has signal in flanking regions. ",
              "Score set to NA.", call. = FALSE)
      score_vec[i]        <- NA_real_
      profile_list[[sid]] <- mean_profile
      next
    }

    norm_profile        <- mean_profile / background
    score_vec[i]        <- max(norm_profile, na.rm = TRUE)
    profile_list[[sid]] <- norm_profile

    if (verbose)
      cat("    TSS enrichment score: ", round(score_vec[i], 3L), "\n", sep = "")
  }

  # --- Assemble result ---------------------------------------------------------
  scores_df <- data.frame(
    sample_id            = samples$sample_id,
    tss_enrichment_score = score_vec,
    pass                 = !is.na(score_vec) & score_vec >= threshold,
    stringsAsFactors     = FALSE
  )
  # Attach all sample metadata columns (everything except bw_path)
  meta_cols <- setdiff(colnames(samples), c("sample_id", "bw_path"))
  if (length(meta_cols) > 0L)
    scores_df <- cbind(scores_df, samples[, meta_cols, drop = FALSE])

  n_pass <- sum(scores_df$pass, na.rm = TRUE)
  if (verbose) {
    cat("[HADES] TSS enrichment complete\n")
    cat("    Threshold: ", threshold, " | Passing: ", n_pass, " / ",
        nrow(scores_df), "\n", sep = "")
  }

  structure(
    list(
      scores   = scores_df,
      profiles = profile_list,
      params   = list(
        window     = window,
        bin_size   = bin_size,
        flank_bins = flank_bins,
        n_bins     = n_bins,
        n_tss      = n_tss,
        threshold  = threshold
      )
    ),
    class = "hades_tss_enrichment"
  )
}


#' @export
print.hades_tss_enrichment <- function(x, ...) {
  p      <- x$params
  scores <- x$scores
  cat("hades_tss_enrichment\n")
  cat("  Window:      ±", p$window, " bp  |  Bin size: ~", p$bin_size,
      " bp  |  n_bins: ", p$n_bins, "\n", sep = "")
  cat("  TSSs used:   ", p$n_tss, "\n", sep = "")
  cat("  Threshold:   ", p$threshold, "\n", sep = "")
  cat("  Samples (", nrow(scores), "):\n", sep = "")
  for (i in seq_len(nrow(scores))) {
    sc   <- scores$tss_enrichment_score[i]
    flag <- if (isTRUE(scores$pass[i])) "PASS" else if (is.na(sc)) " NA " else "FAIL"
    cond <- if ("condition" %in% colnames(scores))
              paste0("  [", scores$condition[i], "]") else ""
    cat("    ", flag, "  ", scores$sample_id[i], cond,
        "  score = ", if (is.na(sc)) "NA" else round(sc, 3L), "\n", sep = "")
  }
  invisible(x)
}
