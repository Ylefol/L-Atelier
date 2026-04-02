# ==============================================================================
# AETHER - TSS/TES Profile Heatmap (scale-regions layout)
# ==============================================================================
# deepTools computeMatrix scale-regions style figure.
# X-axis layout per panel:
#   [TSS - window] ... TSS ... (gene body, scaled to body_bins) ... TES ... [TES + window]
# Rows = genes, sorted by descending total signal (highest at top).
# Each sample produces one panel (profile curve + heatmap).
# ==============================================================================


#' Plot Signal Profile Heatmap (scale-regions layout)
#'
#' @description Produces a deepTools \code{computeMatrix scale-regions}-style
#' figure. Each panel shows a mean signal profile curve above a row-sorted
#' heatmap. The x-axis spans from \code{window} bp upstream of the TSS through
#' the scaled gene body to \code{window} bp downstream of the TES. The gene
#' body is scaled to \code{body_bins} equal bins regardless of gene length, so
#' genes of different sizes are directly comparable. Rows are sorted by
#' descending total signal (highest at top). Multiple samples are displayed
#' side-by-side in one composite figure.
#'
#' @param bigwig_files Named list of BigWig file paths. Each element is a
#'   character vector. If multiple files are supplied they are averaged into
#'   one track. Element names become the panel labels.
#'   Example: `list(WT = c("wt1.bw","wt2.bw"), KO = "ko.bw")`
#' @param txdb TxDb object (from [APOLLO_make_txdb()]). Used to extract TSS
#'   and TES positions for all (or selected) genes.
#' @param genes Optional character vector of gene identifiers to include. If
#'   `NULL` (default) all genes in `txdb` are used. Matched against raw TxDb
#'   gene IDs and NCBI "gene-" stripped forms.
#' @param peaks Optional `GRanges` object or path to a BED/narrowPeak file of
#'   called peaks (e.g. from ATAC or ChIP). When supplied, only genes whose
#'   extended region (body ± `window`) overlaps at least one peak are retained.
#'   This is more principled than the `max_genes` heuristic: rather than
#'   guessing the top-N genes by TSS signal, you provide the exact loci of
#'   interest. `max_genes` still acts as a hard cap on the result.
#' @param window Integer. Base pairs upstream of TSS and downstream of TES to
#'   include as flanking regions. Default `3000`.
#' @param bin_size Integer. Bin width in bp for the flanking regions. The
#'   number of flank bins on each side = `ceiling(window / bin_size)`.
#'   Default `50`.
#' @param body_bins Integer. Number of bins the gene body is scaled to,
#'   regardless of gene length. Default `100`.
#' @param color_low Character. Fill colour for zero/low signal.
#'   Default `"#FFF7EC"` (near-white cream).
#' @param color_high Character. Fill colour for maximum signal.
#'   Default `"#7F0000"` (dark red).
#' @param profile_color Character. Line and fill colour for the profile curve.
#'   Default `"#2166AC"` (blue).
#' @param cap_quantile Numeric in (0,1). Signal values above this quantile are
#'   capped before display (shared cap across all panels). Default `0.99`.
#' @param max_genes Integer or `NULL`. Hard cap on the number of genes
#'   displayed. When `peaks` is supplied this limits the peak-filtered set;
#'   when `peaks` is `NULL` the TSS-signal pre-filter keeps the top
#'   `max_genes` genes. Set `NULL` to disable. Default `5000`.
#' @param profile_height Numeric. Relative height of the profile panel.
#'   Default `1`.
#' @param heatmap_height Numeric. Relative height of the heatmap panel.
#'   Default `5`.
#' @param title Optional character string. Overall figure title.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A `patchwork` composite plot.
#'
#' @details
#' The gene body is scaled independently per gene: a 500 bp gene and a 50 kb
#' gene both contribute \code{body_bins} columns to the matrix. The flanking
#' regions use a fixed \code{bin_size}. Minus-strand genes are flipped so that
#' upstream is always on the left. TSS and TES positions are marked with dashed
#' vertical lines on both the profile and heatmap panels.
#'
#' Requires: `GenomicRanges`, `IRanges`, `S4Vectors`, `rtracklayer`,
#' `GenomicFeatures`, `EnrichedHeatmap`, `patchwork`.
#'
#' @export
AETHER_plot_profile_heatmap <- function(bigwig_files,
                                         txdb,
                                         genes          = NULL,
                                         peaks          = NULL,
                                         window         = 3000L,
                                         bin_size       = 50L,
                                         body_bins      = 100L,
                                         color_low      = "#FFF7EC",
                                         color_high     = "#7F0000",
                                         profile_color  = "#2166AC",
                                         cap_quantile   = 0.99,
                                         max_genes      = 5000L,
                                         profile_height = 1,
                                         heatmap_height = 5,
                                         title          = NULL,
                                         verbose        = TRUE) {

  # --- Package checks ---------------------------------------------------------
  for (pkg in c("GenomicRanges", "IRanges", "S4Vectors",
                "rtracklayer", "GenomicFeatures", "EnrichedHeatmap",
                "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. ",
           "Install with: BiocManager::install('", pkg, "')", call. = FALSE)
  }

  # --- Input validation -------------------------------------------------------
  window    <- as.integer(window)
  bin_size  <- as.integer(bin_size)
  body_bins <- as.integer(body_bins)
  n_up      <- as.integer(ceiling(window / bin_size))
  n_down    <- n_up
  total_bins <- n_up + body_bins + n_down

  if (is.character(bigwig_files)) {
    if (is.null(names(bigwig_files)) || any(names(bigwig_files) == ""))
      stop("bigwig_files must be a fully named list or named character vector",
           call. = FALSE)
    bigwig_files <- as.list(bigwig_files)
  }
  if (!is.list(bigwig_files) || is.null(names(bigwig_files)) ||
      any(names(bigwig_files) == ""))
    stop("bigwig_files must be a named list", call. = FALSE)

  # --- Gene positions ---------------------------------------------------------
  if (verbose) cat("[AETHER] Extracting gene positions from TxDb...")
  pos_df <- .aether_phm_get_positions(txdb, genes, verbose)
  if (nrow(pos_df) == 0L)
    stop("No genes found — check txdb and the genes argument", call. = FALSE)
  if (verbose) cat("[AETHER]   Genes: ", nrow(pos_df))

  # --- Gene filtering ---------------------------------------------------------
  if (!is.null(peaks)) {
    # --- Peak-based filter (preferred) --------------------------------------
    # Resolve peaks to GRanges from various input types
    if (is.character(peaks)) {
      if (!file.exists(peaks))
        stop("peaks file not found: ", peaks, call. = FALSE)
      peaks <- rtracklayer::import(peaks)
    } else if (is.data.frame(peaks)) {
      # Accept data.frame from APOLLO_annotate_peaks (has chr/start/end columns)
      if (!all(c("chr", "start", "end") %in% colnames(peaks)))
        stop("peaks data.frame must have columns: chr, start, end", call. = FALSE)
      peaks <- GenomicRanges::GRanges(
        seqnames = peaks$chr,
        ranges   = IRanges::IRanges(start = peaks$start + 1L, end = peaks$end)
      )
    }
    if (!inherits(peaks, "GRanges"))
      stop("peaks must be a GRanges object, a data.frame with chr/start/end ",
           "columns (from APOLLO_annotate_peaks), or a path to a BED/narrowPeak file.",
           call. = FALSE)

    if (verbose)
      cat("[AETHER] Filtering genes by peak overlap (", length(peaks),
              " peaks)...")

    # Extended gene regions: body ± window
    gene_ext <- GenomicRanges::GRanges(
      seqnames = pos_df$chr,
      ranges   = IRanges::IRanges(
        start = pmax(1L, pmin(pos_df$tss, pos_df$tes) - window),
        end   = pmax(pos_df$tss, pos_df$tes) + window
      )
    )
    has_peak <- GenomicRanges::countOverlaps(gene_ext, peaks) > 0L
    pos_df   <- pos_df[has_peak, ]

    if (nrow(pos_df) == 0L)
      stop("No genes overlap the supplied peaks. Check that peaks and TxDb ",
           "use the same chromosome naming style.", call. = FALSE)
    if (verbose)
      cat("[AETHER]   Genes with overlapping peaks: ", nrow(pos_df))

    # Still apply max_genes as a hard cap if needed
    if (!is.null(max_genes) && nrow(pos_df) > max_genes) {
      if (verbose)
        cat("[AETHER]   Capping to ", max_genes,
                " genes by TSS signal pre-filter...")
      pos_df <- .aether_phm_prefilter(
        pos_df    = pos_df,
        bw_path   = as.list(bigwig_files[[1L]])[[1L]],
        max_genes = max_genes,
        score_win = min(500L, window),
        verbose   = verbose
      )
    }

  } else if (!is.null(max_genes) && nrow(pos_df) > max_genes) {
    # --- TSS-signal pre-filter (fallback when no peaks supplied) ------------
    if (verbose)
      cat("[AETHER] No peaks supplied — pre-filtering ", nrow(pos_df),
              " genes to top ", max_genes, " by TSS signal...")
    pos_df <- .aether_phm_prefilter(
      pos_df    = pos_df,
      bw_path   = as.list(bigwig_files[[1L]])[[1L]],
      max_genes = max_genes,
      score_win = min(500L, window),
      verbose   = verbose
    )
    if (verbose) cat("[AETHER]   Retained: ", nrow(pos_df), " genes")
  }

  # --- Extract signal matrices ------------------------------------------------
  if (verbose) cat("[AETHER] Extracting signal matrices...")

  raw_results <- setNames(lapply(names(bigwig_files), function(samp) {
    bw_files <- as.list(bigwig_files[[samp]])
    if (verbose) cat("[AETHER]   ", samp)
    .aether_phm_bw_to_matrix_sr(bw_files, pos_df, window, n_up, body_bins, n_down)
  }), names(bigwig_files))

  # Unpack matrices and sync layout with what normalizeToMatrix actually produced
  all_mats <- lapply(raw_results, `[[`, "mat")
  actual   <- raw_results[[1L]]$layout
  n_up      <- actual$n_up
  body_bins <- actual$body_bins
  n_down    <- actual$n_down
  total_bins <- n_up + body_bins + n_down

  # --- Row order (sort by descending total signal, first sample) --------------
  row_order <- order(rowSums(all_mats[[1L]], na.rm = TRUE), decreasing = TRUE)
  n_display <- length(row_order)

  # --- Global color cap -------------------------------------------------------
  all_pos_vals <- unlist(lapply(all_mats, function(m) {
    v <- m[row_order, ]; v[v > 0]
  }))
  cap_val <- if (length(all_pos_vals) > 0L)
    quantile(all_pos_vals, cap_quantile, na.rm = TRUE)
  else 1.0

  # --- Shared profile y-limit (same scale across all panels) -----------------
  profile_ylim <- max(
    sapply(all_mats, function(m)
      max(colMeans(m[row_order, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    ), na.rm = TRUE
  ) * 1.15   # 15% headroom, consistent with previous expansion(mult=c(0,0.15))

  # --- X-axis ticks -----------------------------------------------------------
  # Bin coordinates: 1 = first upstream bin, n_up = last upstream bin (TSS),
  # n_up + body_bins = TES, total_bins = last downstream bin
  tss_bin  <- n_up + 0.5          # boundary between upstream and body
  tes_bin  <- n_up + body_bins + 0.5  # boundary between body and downstream

  x_breaks <- c(1, tss_bin, (tss_bin + tes_bin) / 2, tes_bin, total_bins)
  x_labels <- c(paste0("-", window / 1000L, "kb"), "TSS", "", "TES",
                paste0("+", window / 1000L, "kb"))

  # --- Build panels -----------------------------------------------------------
  n_panels <- length(bigwig_files)
  panels   <- vector("list", n_panels)

  for (i in seq_along(names(bigwig_files))) {
    samp <- names(bigwig_files)[i]
    mat  <- all_mats[[samp]][row_order, , drop = FALSE]
    mat  <- pmin(mat, cap_val)
    panels[[i]] <- .aether_phm_build_panel_sr(
      mat           = mat,
      panel_title   = samp,
      n_up          = n_up,
      body_bins     = body_bins,
      n_down        = n_down,
      tss_bin       = tss_bin,
      tes_bin       = tes_bin,
      x_breaks      = x_breaks,
      x_labels      = x_labels,
      color_low     = color_low,
      color_high    = color_high,
      profile_color = profile_color,
      cap_val       = cap_val,
      profile_ylim  = profile_ylim,
      show_y_axis   = (i == 1L),
      show_legend   = (i == n_panels)
    )
  }

  # --- Assemble ---------------------------------------------------------------
  combined <- lapply(panels, function(p)
    patchwork::wrap_plots(list(p$profile, p$heatmap),
                          ncol    = 1,
                          heights = c(profile_height, heatmap_height))
  )

  result <- patchwork::wrap_plots(combined, nrow = 1)
  plot_title <- if (!is.null(title))
    paste0(title, " (n = ", n_display, ")")
  else
    paste0("n = ", n_display)

  result <- result +
    patchwork::plot_annotation(
      title = plot_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold",
                                              hjust = 0.5))
    )
  result
}


# ------------------------------------------------------------------------------
# Internal: extract TSS/TES positions from TxDb (shared with previous version)
# ------------------------------------------------------------------------------
.aether_phm_get_positions <- function(txdb, genes, verbose) {

  all_genes <- tryCatch(
    GenomicFeatures::genes(txdb),
    error = function(e)
      stop("Failed to retrieve genes from TxDb: ", conditionMessage(e),
           call. = FALSE)
  )

  if (!is.null(genes)) {
    raw_ids     <- names(all_genes)
    cleaned_ids <- sub("^gene-", "", raw_ids)
    keep        <- raw_ids %in% genes | cleaned_ids %in% genes
    n_miss      <- length(setdiff(genes, c(raw_ids[keep], cleaned_ids[keep])))
    if (n_miss > 0L && verbose)
      cat("[AETHER]   ", n_miss, " supplied gene(s) not found in TxDb")
    all_genes <- all_genes[keep]
  }

  if (length(all_genes) == 0L)
    return(data.frame(gene_id = character(), chr = character(),
                      tss = integer(), tes = integer(), strand = character(),
                      stringsAsFactors = FALSE))

  strands <- as.character(GenomicRanges::strand(all_genes))
  starts  <- GenomicRanges::start(all_genes)
  ends    <- GenomicRanges::end(all_genes)

  data.frame(
    gene_id = names(all_genes),
    chr     = as.character(GenomicRanges::seqnames(all_genes)),
    tss     = ifelse(strands == "-", ends,   starts),
    tes     = ifelse(strands == "-", starts, ends),
    strand  = strands,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# Internal: fast TSS-signal pre-filter
#
# Imports only ±score_win bp around each TSS (a tiny fraction of the full
# gene-body regions), sums signal per gene, and returns the top max_genes rows
# of pos_df sorted by descending TSS signal.
# Using reduce() on TSS windows collapses overlapping regions into fewer import
# queries, making the BigWig import much faster than querying per gene.
# ------------------------------------------------------------------------------
.aether_phm_prefilter <- function(pos_df, bw_path, max_genes,
                                    score_win = 500L, verbose = TRUE) {
  n <- nrow(pos_df)

  tss_gr <- GenomicRanges::GRanges(
    seqnames = pos_df$chr,
    ranges   = IRanges::IRanges(
      start = pmax(1L, pos_df$tss - score_win),
      end   = pos_df$tss + score_win
    )
  )

  # Merge overlapping windows → fewer import queries
  query_gr <- GenomicRanges::reduce(tss_gr)

  sig <- tryCatch(
    rtracklayer::import.bw(bw_path, which = query_gr),
    error = function(e) NULL
  )

  if (is.null(sig) || length(sig) == 0L) {
    if (verbose) cat("[AETHER]   No TSS signal found — keeping first ",
                         max_genes, " genes")
    return(pos_df[seq_len(min(max_genes, n)), ])
  }

  # Sum BigWig signal within each per-gene TSS window
  hits   <- GenomicRanges::findOverlaps(tss_gr, sig)
  scores <- numeric(n)
  if (length(hits) > 0L) {
    agg <- tapply(sig$score[S4Vectors::subjectHits(hits)],
                  S4Vectors::queryHits(hits), sum)
    scores[as.integer(names(agg))] <- as.numeric(agg)
  }

  top_idx <- order(scores, decreasing = TRUE)[seq_len(min(max_genes, n))]
  pos_df[top_idx, ]
}


# ------------------------------------------------------------------------------
# Internal: scale-regions signal matrix extraction via EnrichedHeatmap
#
# Delegates to EnrichedHeatmap::normalizeToMatrix() which handles:
#   - Gene body scaling to body_bins columns
#   - Fixed flank bins on each side (window / bin_size each)
#   - Strand-aware extraction (minus-strand genes flipped automatically)
#   - Weighted mean per bin (mean_mode="w0": 0 for uncovered positions)
#
# Returns matrix(n_genes x (n_up + body_bins + n_down)), averaged over files.
# ------------------------------------------------------------------------------
.aether_phm_bw_to_matrix_sr <- function(bw_files, pos_df, window,
                                          n_up, body_bins, n_down) {
  if (!requireNamespace("EnrichedHeatmap", quietly = TRUE))
    stop("Package 'EnrichedHeatmap' is required. ",
         "Install with: BiocManager::install('EnrichedHeatmap')", call. = FALSE)

  # Gene body GRanges — normalizeToMatrix handles strand flipping internally
  target_gr <- GenomicRanges::GRanges(
    seqnames = pos_df$chr,
    ranges   = IRanges::IRanges(
      start = pmin(pos_df$tss, pos_df$tes),
      end   = pmax(pos_df$tss, pos_df$tes)
    ),
    strand = pos_df$strand
  )
  names(target_gr) <- pos_df$gene_id

  # Flank bin width; derived from window / n_up to stay consistent with
  # however n_up was computed (ceiling(window / bin_size))
  bin_w <- as.integer(round(window / n_up))

  # Import regions: gene bodies + flanks on both sides (restrict BigWig import)
  import_gr <- GenomicRanges::GRanges(
    seqnames = pos_df$chr,
    ranges   = IRanges::IRanges(
      start = pmax(1L, pmin(pos_df$tss, pos_df$tes) - window),
      end   = pmax(pos_df$tss, pos_df$tes) + window
    )
  )

  mat_list <- lapply(bw_files, function(f) {
    if (!file.exists(f))
      stop("BigWig not found: ", f, call. = FALSE)

    # Import as GRanges with score column — normalizeToMatrix requires GRanges,
    # not RleList, when include_target=TRUE (gene body mode)
    sig_gr <- tryCatch(
      rtracklayer::import.bw(f, which = import_gr),
      error = function(e)
        stop("Failed to read '", basename(f), "': ", conditionMessage(e),
             call. = FALSE)
    )

    nm <- EnrichedHeatmap::normalizeToMatrix(
      signal         = sig_gr,
      target         = target_gr,
      value_column   = "score",
      extend         = window,
      w              = bin_w,
      k              = body_bins,
      include_target = TRUE,
      mean_mode      = "w0",
      background     = 0,
      verbose        = FALSE
    )

    m <- unclass(nm)         # strip NormalizedMatrix S4 class → plain matrix
    m[is.na(m)] <- 0.0

    # Capture the actual layout from NormalizedMatrix attributes on first file
    if (!exists("layout_out")) {
      layout_out <<- list(
        n_up      = length(attr(nm, "upstream_index")),
        body_bins = length(attr(nm, "target_index")),
        n_down    = length(attr(nm, "downstream_index"))
      )
    }
    m
  })

  mat <- Reduce("+", mat_list) / length(mat_list)
  list(mat = mat, layout = layout_out)
}


# ------------------------------------------------------------------------------
# Internal: build one profile + heatmap panel (scale-regions layout)
# Returns list(profile = ggplot, heatmap = ggplot)
# ------------------------------------------------------------------------------
.aether_phm_build_panel_sr <- function(mat, panel_title, n_up, body_bins, n_down,
                                        tss_bin, tes_bin, x_breaks, x_labels,
                                        color_low, color_high, profile_color,
                                        cap_val, profile_ylim, show_y_axis,
                                        show_legend) {

  n_genes    <- nrow(mat)
  total_bins <- n_up + body_bins + n_down
  x_ctrs     <- seq_len(total_bins)   # bin indices as x coordinates

  # --- Profile curve ----------------------------------------------------------
  mean_sig <- colMeans(mat, na.rm = TRUE)

  p_profile <- ggplot(data.frame(x = x_ctrs, y = mean_sig),
                      aes(x = .data$x, y = .data$y)) +
    geom_area(fill = profile_color, alpha = 0.25) +
    geom_line(color = profile_color, linewidth = 0.6) +
    geom_vline(xintercept = tss_bin, color = "grey40", linetype = "dashed",
               linewidth = 0.35) +
    geom_vline(xintercept = tes_bin, color = "grey40", linetype = "dashed",
               linewidth = 0.35) +
    scale_x_continuous(limits = c(0.5, total_bins + 0.5), expand = c(0, 0),
                       breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, profile_ylim), expand = c(0, 0)) +
    labs(title = panel_title, x = NULL, y = if (show_y_axis) "Mean" else NULL) +
    theme_bw() +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.title.x       = element_blank(),
      axis.title.y       = if (show_y_axis) element_text(size = 7) else element_blank(),
      axis.text.y        = if (show_y_axis) element_text(size = 6) else element_blank(),
      axis.ticks.y       = if (show_y_axis) element_line() else element_blank(),
      plot.title         = element_text(size = 8, hjust = 0.5, face = "bold"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4, 4, 0, 4)
    )

  # --- Heatmap ----------------------------------------------------------------
  # y is reversed: row 1 (highest signal) at the top (y = n_genes in reversed axis)
  heat_df <- data.frame(
    x    = rep(x_ctrs,          each  = n_genes),
    y    = rep(seq_len(n_genes), times = total_bins),
    fill = pmax(0.0, as.vector(mat))
  )

  p_heat <- ggplot(heat_df, aes(x = .data$x, y = .data$y, fill = .data$fill)) +
    geom_raster(interpolate = FALSE) +
    geom_vline(xintercept = tss_bin, color = "white", linetype = "dashed",
               linewidth = 0.4) +
    geom_vline(xintercept = tes_bin, color = "white", linetype = "dashed",
               linewidth = 0.4) +
    scale_fill_gradient(
      low    = color_low,
      high   = color_high,
      limits = c(0, cap_val),
      oob    = function(x, range) pmin(pmax(x, range[1]), range[2]),
      name   = "Signal"
    ) +
    scale_x_continuous(limits = c(0.5, total_bins + 0.5), expand = c(0, 0),
                       breaks = x_breaks, labels = x_labels) +
    # Reversed y: row 1 (top of sorted matrix = highest signal) appears at top
    scale_y_reverse(expand = c(0, 0),
                    breaks = c(1L, as.integer(n_genes / 2L), n_genes),
                    labels = c("1", as.character(as.integer(n_genes / 2L)),
                               as.character(n_genes))) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.x       = element_text(size = 7),
      axis.title.x      = element_blank(),
      legend.position   = if (show_legend) "right" else "none",
      legend.key.width  = unit(0.3, "cm"),
      legend.key.height = unit(0.8, "cm"),
      legend.text       = element_text(size = 6),
      legend.title      = element_text(size = 7),
      panel.border      = element_rect(colour = "grey40"),
      plot.margin       = margin(0, 4, 4, 4)
    )

  list(profile = p_profile, heatmap = p_heat)
}
