# ==============================================================================
# AETHER - Coverage Track Plots
# ==============================================================================
# Genome browser-style coverage track plots from BigWig files.
# Multiple tracks stacked vertically over a constrained genomic region,
# with optional peak annotation overlay and replicate variance display.
# ==============================================================================


#' Plot Genome Coverage Tracks
#'
#' @description Produces a genome browser-style stacked coverage track plot
#' over a specified genomic region. Each named track shows the read pileup
#' signal as a filled area curve. Multiple files per track are averaged into
#' a single signal; an optional variance ribbon can be displayed when
#' replicates are provided.
#'
#' Two input modes are supported — provide exactly one:
#' - **BigWig mode** (`bigwig_files`): reads signal from pre-computed BigWig
#'   files. Use when BigWigs already exist (e.g. from [ELEUTHIA_bed_to_bigwig()]).
#' - **BED mode** (`bed_files`): computes coverage on-the-fly for the
#'   requested region directly from fragment-level BED files via
#'   `ELEUTHIA_bed_region_coverage()`. More efficient when only a small region
#'   is needed and BigWig files have not been pre-generated.
#'
#' An optional peak annotation track can be added at the bottom.
#'
#' @param bigwig_files Named list of character vectors, one element per track.
#'   Each element is a character vector of BigWig file paths to merge into
#'   that track. A named character vector is accepted as shorthand where each
#'   track contains a single file. Mutually exclusive with `bed_files`.
#'
#' @param bed_files Named list of character vectors, one element per track.
#'   Each element is a character vector of fragment-level BED file paths.
#'   Coverage is computed on-the-fly for the requested region only.
#'   A named character vector is accepted as shorthand. Mutually exclusive
#'   with `bigwig_files`. Requires `chrom_sizes`.
#'
#' @param region Genomic region to display. One of:
#'   - Character string: `"chr1:1000000-2000000"`
#'   - Named character vector: `c(chr="chr1", start="1000000", end="2000000")`
#'   - A `GRanges` object (single range)
#' @param chrom_sizes Required when using `bed_files`. Exact chromosome
#'   lengths as a named numeric vector or `APOLLO_get_chromosome_sizes()`
#'   data.frame. Ignored in BigWig mode.
#' @param normalize Character. Normalization for BED mode only: `"CPM"`
#'   (default) or `"raw"`. Ignored in BigWig mode (signal is already
#'   normalized in the file).
#' @param bin_size Integer. Bin width in bp for BED mode only. Default `10`.
#'   Ignored in BigWig mode.
#' @param peak_bed Optional. Peak annotation track displayed at the bottom.
#'   Either a path to a BED file or a data.frame with at least three columns
#'   (chr, start, end).
#' @param group_colors Optional named character vector mapping track names to
#'   colors. If `NULL`, colors are assigned automatically — known group names
#'   (WT, KO, etc.) use the AETHER default palette; others are generated.
#' @param scale_y Character. Y-axis scaling: `"free"` (default, each track
#'   autoscales independently) or `"fixed"` (all tracks share the same
#'   y-axis maximum).
#' @param show_variance Logical. When a track contains multiple files, display
#'   a shaded ribbon showing the spread across replicates. Default `FALSE`.
#'   Has no effect on single-file tracks.
#' @param variance_type Character. Type of variance ribbon when
#'   `show_variance = TRUE`: `"sd"` (default, mean ± 1 SD) or `"range"`
#'   (min to max across replicates).
#' @param track_height Numeric. Relative height unit for each coverage track.
#'   Default `1`.
#' @param annotation_height Numeric. Relative height unit for the peak
#'   annotation track. Default `0.3`.
#' @param peak_color Character. Fill color for peak rectangles in the
#'   annotation track. Default `"steelblue4"`.
#' @param title Optional character string. Plot title placed above the first
#'   track.
#' @param verbose Logical. Print progress messages. Default `FALSE`.
#'
#' @return A `patchwork` plot object (combination of `ggplot` panels).
#'
#' @details
#' When multiple files are provided for a track, per-bin mean is computed
#' across all files. Bins absent in a file are treated as zero before
#' averaging.
#'
#' BigWig mode requires: `GenomicRanges`, `IRanges`, `rtracklayer`.
#' BED mode requires: `data.table`, `GenomicRanges`, `IRanges`, `GenomeInfoDb`.
#' Both modes require: `patchwork`.
#'
#' @export
AETHER_plot_coverage_tracks <- function(bigwig_files      = NULL,
                                        bed_files         = NULL,
                                        region,
                                        chrom_sizes       = NULL,
                                        normalize         = "CPM",
                                        bin_size          = 10L,
                                        peak_bed          = NULL,
                                        group_colors      = NULL,
                                        scale_y           = "free",
                                        show_variance     = FALSE,
                                        variance_type     = "sd",
                                        track_height      = 1,
                                        annotation_height = 0.3,
                                        peak_color        = "steelblue4",
                                        title             = NULL,
                                        verbose           = FALSE) {

  # --- Mode validation --------------------------------------------------------
  if (is.null(bigwig_files) && is.null(bed_files))
    stop("Provide either bigwig_files (BigWig mode) or bed_files (BED mode)",
         call. = FALSE)
  if (!is.null(bigwig_files) && !is.null(bed_files))
    stop("bigwig_files and bed_files are mutually exclusive — provide only one",
         call. = FALSE)

  use_bed <- !is.null(bed_files)

  if (use_bed && is.null(chrom_sizes))
    stop("chrom_sizes is required when using bed_files mode", call. = FALSE)

  # --- Package checks ---------------------------------------------------------
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required. Install with: install.packages('patchwork')",
         call. = FALSE)

  if (!use_bed) {
    for (pkg in c("GenomicRanges", "IRanges", "rtracklayer")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' required. ",
             "Install with: BiocManager::install('", pkg, "')",
             call. = FALSE)
    }

  }

  # --- Normalise input --------------------------------------------------------
  # Accept named character vector as shorthand for single-file list
  .normalise_files <- function(x, arg) {
    if (is.character(x)) {
      if (is.null(names(x)) || any(names(x) == ""))
        stop(arg, " must be a fully named character vector or named list",
             call. = FALSE)
      x <- as.list(x)
    }
    if (!is.list(x) || is.null(names(x)) || any(names(x) == ""))
      stop(arg, " must be a named list or named character vector", call. = FALSE)
    x
  }

  source_files  <- if (use_bed) {
    .normalise_files(bed_files,    "bed_files")
  } else {
    .normalise_files(bigwig_files, "bigwig_files")
  }

  scale_y       <- match.arg(scale_y,       c("free", "fixed"))
  variance_type <- match.arg(variance_type, c("sd", "range"))
  normalize     <- match.arg(normalize,     c("CPM", "raw"))

  # --- Parse region -----------------------------------------------------------
  region_gr  <- .AETHER_parse_region(region)
  chr        <- as.character(GenomicRanges::seqnames(region_gr))
  reg_start  <- GenomicRanges::start(region_gr)
  reg_end    <- GenomicRanges::end(region_gr)
  region_lab <- paste0(chr, ":",
                       format(reg_start, big.mark = ",", scientific = FALSE),
                       "-",
                       format(reg_end,   big.mark = ",", scientific = FALSE))

  if (verbose) message("[AETHER] Mode: ", if (use_bed) "BED" else "BigWig",
                       " | Region: ", region_lab)

  # --- Read and merge signal per track ----------------------------------------
  cov_list <- lapply(names(source_files), function(nm) {
    files <- source_files[[nm]]

    dfs <- lapply(files, function(f) {
      if (use_bed) {
        if (verbose) message("[AETHER]   Reading [", nm, "]: ", basename(f))
        ELEUTHIA_bed_region_coverage(f, region_gr, chrom_sizes,
                                     normalize = normalize,
                                     bin_size  = bin_size,
                                     verbose   = FALSE)
      } else {
        if (!file.exists(f))
          stop("BigWig not found: ", f, call. = FALSE)
        if (verbose) message("[AETHER]   Reading [", nm, "]: ", basename(f))
        sig <- tryCatch(
          rtracklayer::import.bw(f, which = region_gr),
          error = function(e) {
            stop("rtracklayer failed to read BigWig file: ", basename(f), "\n",
                 "  Original error: ", conditionMessage(e), "\n",
                 "  If this is an S4 dispatch or seqinfo error, your Bioconductor\n",
                 "  packages may be out of sync. Run BiocManager::valid() to\n",
                 "  identify and reinstall any outdated packages.",
                 call. = FALSE)
          }
        )
        df  <- as.data.frame(sig)[, c("start", "end", "score"), drop = FALSE]
        df$midpoint <- (df$start + df$end) / 2
        df
      }
    })

    if (length(dfs) == 1L) {
      # Single file — no averaging needed
      df <- dfs[[1]]
      data.frame(
        midpoint   = df$midpoint,
        sample     = nm,
        mean_score = df$score,
        sd_score   = NA_real_,
        min_score  = df$score,
        max_score  = df$score,
        stringsAsFactors = FALSE
      )
    } else {
      # Multiple files — merge by midpoint, fill missing bins with 0
      all_mid <- sort(unique(unlist(lapply(dfs, `[[`, "midpoint"))))
      score_mat <- sapply(dfs, function(df) {
        s <- df$score[match(all_mid, df$midpoint)]
        s[is.na(s)] <- 0
        s
      })
      data.frame(
        midpoint   = all_mid,
        sample     = nm,
        mean_score = rowMeans(score_mat),
        sd_score   = apply(score_mat, 1, sd),
        min_score  = apply(score_mat, 1, min),
        max_score  = apply(score_mat, 1, max),
        stringsAsFactors = FALSE
      )
    }
  })
  names(cov_list) <- names(source_files)

  # --- Assign track colors ----------------------------------------------------
  n_tracks <- length(source_files)

  if (is.null(group_colors)) {
    known_pal <- AETHER_default_palette("groups")
    if (all(names(source_files) %in% names(known_pal))) {
      group_colors <- known_pal[names(source_files)]
    } else {
      base_colors <- c("#2166AC", "#D6604D", "#1A9641", "#762A83",
                       "#F4A582", "#4DAC26", "#CA0020", "#5E4FA2")
      if (n_tracks <= length(base_colors)) {
        group_colors <- setNames(base_colors[seq_len(n_tracks)], names(source_files))
      } else {
        group_colors <- setNames(
          grDevices::colorRampPalette(base_colors)(n_tracks),
          names(source_files)
        )
      }
    }
  }

  # --- Shared y-axis limit if fixed -------------------------------------------
  y_lim <- NULL
  if (scale_y == "fixed") {
    y_max <- max(sapply(cov_list, function(df) {
      if (show_variance && !all(is.na(df$sd_score))) {
        if (variance_type == "sd")    max(df$mean_score + df$sd_score, na.rm = TRUE)
        else                          max(df$max_score,                 na.rm = TRUE)
      } else {
        max(df$mean_score, na.rm = TRUE)
      }
    }), na.rm = TRUE)
    y_lim <- c(0, y_max * 1.05)
  }

  # --- Build one ggplot per coverage track ------------------------------------
  has_annotation <- !is.null(peak_bed)

  track_plots <- lapply(seq_along(names(source_files)), function(i) {
    nm    <- names(source_files)[i]
    df    <- cov_list[[nm]]
    color <- group_colors[[nm]]

    n_files        <- length(source_files[[nm]])
    can_variance   <- show_variance && n_files > 1L && !all(is.na(df$sd_score))
    show_xaxis     <- (i == n_tracks) && !has_annotation

    # Build geom layers
    layers <- list(
      geom_line(data = df,
                aes(x = .data$midpoint, y = .data$mean_score),
                color = color, linewidth = 0.6)
    )

    if (can_variance) {
      ribbon <- if (variance_type == "sd") {
        geom_ribbon(
          data = df,
          aes(x    = .data$midpoint,
              ymin = pmax(0, .data$mean_score - .data$sd_score),
              ymax = .data$mean_score + .data$sd_score),
          fill = color, alpha = 0.25
        )
      } else {
        geom_ribbon(
          data = df,
          aes(x    = .data$midpoint,
              ymin = .data$min_score,
              ymax = .data$max_score),
          fill = color, alpha = 0.25
        )
      }
      # Insert ribbon first so it renders behind the area and line
      layers <- c(list(ribbon), layers)
    }

    p <- ggplot() +
      layers +
      scale_x_continuous(
        limits = c(reg_start, reg_end),
        expand = c(0, 0),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)
      ) +
      labs(
        y = nm,
        x = if (show_xaxis) region_lab else NULL
      ) +
      theme_bw() +
      theme(
        axis.text.x      = if (show_xaxis) element_text(size = 8) else element_blank(),
        axis.ticks.x     = if (show_xaxis) element_line()         else element_blank(),
        axis.title.x     = if (show_xaxis) element_text(size = 8) else element_blank(),
        axis.title.y     = element_text(size = 8, angle = 0, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(1, 5, 1, 5)
      )

    if (!is.null(y_lim)) {
      p <- p + scale_y_continuous(limits = y_lim, expand = c(0, 0))
    } else {
      p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    }

    p
  })

  # Add title to first track
  if (!is.null(title)) {
    track_plots[[1]] <- track_plots[[1]] +
      labs(title = title) +
      theme(plot.title = element_text(size = 10, face = "bold"))
  }

  # --- Build peak annotation track (optional) ---------------------------------
  heights <- rep(track_height, n_tracks)

  if (has_annotation) {
    if (is.character(peak_bed) && length(peak_bed) == 1L) {
      if (!file.exists(peak_bed))
        stop("peak_bed file not found: ", peak_bed, call. = FALSE)
      peaks_raw <- utils::read.table(peak_bed, header = FALSE, sep = "\t",
                                     stringsAsFactors = FALSE)
    } else {
      peaks_raw <- as.data.frame(peak_bed)
    }

    if (ncol(peaks_raw) < 3)
      stop("peak_bed must have at least 3 columns (chr, start, end)", call. = FALSE)

    peaks_raw         <- peaks_raw[, 1:3, drop = FALSE]
    colnames(peaks_raw) <- c("chr", "start", "end")

    peaks <- peaks_raw[
      peaks_raw$chr   == chr       &
        peaks_raw$end   >= reg_start &
        peaks_raw$start <= reg_end,
      , drop = FALSE
    ]
    peaks$start <- pmax(peaks$start, reg_start)
    peaks$end   <- pmin(peaks$end,   reg_end)

    p_anno <- ggplot() +
      scale_x_continuous(
        limits = c(reg_start, reg_end),
        expand = c(0, 0),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = region_lab, y = "Peaks") +
      theme_bw() +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0.5),
        axis.text.x  = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid   = element_blank(),
        plot.margin  = margin(1, 5, 1, 5)
      )

    if (nrow(peaks) > 0) {
      p_anno <- p_anno +
        geom_rect(
          data = peaks,
          aes(xmin = .data$start, xmax = .data$end, ymin = 0.1, ymax = 0.9),
          fill = peak_color, color = NA
        )
    }

    # Remove x-axis from last coverage track — annotation takes it over
    track_plots[[n_tracks]] <- track_plots[[n_tracks]] +
      labs(x = NULL) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank())

    all_plots <- c(track_plots, list(p_anno))
    heights   <- c(heights, annotation_height)
  } else {
    all_plots <- track_plots
  }

  # --- Assemble with patchwork ------------------------------------------------
  patchwork::wrap_plots(all_plots, ncol = 1, heights = heights)
}


# Internal: parse region specification to a single GRanges
.AETHER_parse_region <- function(region) {

  if (inherits(region, "GRanges")) {
    if (length(region) != 1L)
      stop("region GRanges must contain exactly one range", call. = FALSE)
    return(region)
  }

  if (is.character(region) && length(region) == 1L && is.null(names(region))) {
    # "chr1:1000000-2000000"
    parts <- regmatches(region,
                        regexec("^([^:]+):(\\d+)-(\\d+)$", region))[[1]]
    if (length(parts) != 4L)
      stop("Cannot parse region '", region,
           "'. Expected format: 'chr:start-end'", call. = FALSE)
    return(GenomicRanges::GRanges(
      parts[2],
      IRanges::IRanges(as.integer(parts[3]), as.integer(parts[4]))
    ))
  }

  if (is.character(region) && !is.null(names(region))) {
    # c(chr="chr1", start="1000000", end="2000000")
    required <- c("chr", "start", "end")
    if (!all(required %in% names(region)))
      stop("Named region vector must contain 'chr', 'start', and 'end' elements",
           call. = FALSE)
    return(GenomicRanges::GRanges(
      region[["chr"]],
      IRanges::IRanges(as.integer(region[["start"]]), as.integer(region[["end"]]))
    ))
  }

  stop("region must be a GRanges object, a 'chr:start-end' string, ",
       "or a named character vector with chr/start/end elements",
       call. = FALSE)
}
