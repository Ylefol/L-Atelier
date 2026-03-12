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
#' An optional peak annotation track can be added at the bottom. Additional
#' annotation tracks (e.g. CpG islands, methylation sites) can be added via
#' `annotation_tracks`.
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
#' @param peak_bed Optional. Peak annotation track displayed below coverage
#'   tracks. Either a path to a BED file or a data.frame with at least three
#'   columns (chr, start, end). Equivalent to adding a `"region"` type entry
#'   to `annotation_tracks` named `"Peaks"`.
#' @param annotation_tracks Optional named list of additional annotation tracks
#'   to display below coverage tracks (and below `peak_bed` if provided). Each
#'   element is a named list describing one track. Two types are supported:
#'   \describe{
#'     \item{`type = "region"`}{Genomic intervals shown as filled rectangles
#'       (e.g. CpG islands, repeats). Required field: `data` (BED file path or
#'       data.frame with chr/start/end columns). Optional fields: `color`
#'       (fill colour, default `"steelblue4"`), `height` (relative patchwork
#'       height, default = `annotation_height`).}
#'     \item{`type = "score"`}{Per-site scores shown as a lollipop plot — a
#'       vertical segment from 0 to the score, capped with a point, both
#'       coloured by score value (e.g. methylation fractions). Required field:
#'       `data` (BigWig file path). Optional fields: `color_low` (colour at
#'       score=0, default `"steelblue"`), `color_high` (colour at score=1,
#'       default `"firebrick"`), `score_limits` (numeric(2), y-axis and colour
#'       scale limits; auto-detected from data if `NULL`), `height` (relative
#'       patchwork height, default = `annotation_height * 2`).}
#'     \item{`type = "genes"`}{Collapsed IGV-style gene model track. Shows
#'       the intron backbone, exon blocks (UTR height), CDS blocks (taller),
#'       strand direction arrows, and italic gene name labels. Multiple genes
#'       are automatically stacked into rows to avoid overlap. Required field:
#'       either `txdb` (TxDb object from `APOLLO_make_txdb()` — fast, preferred)
#'       or `gtf_file` (path to GTF — slower, reads genome-wide). Optional
#'       fields: `chr_mapping` (e.g. `"T2T"`, applied when using `gtf_file`),
#'       `show_cds` (logical, default `TRUE`), `label_genes` (logical, default
#'       `TRUE`), `label_size` (default `2.5`), `color` (default `"black"`),
#'       `height` (default = `annotation_height * 3`).}
#'   }
#'   Example:
#'   ```r
#'   annotation_tracks = list(
#'     "CpG Islands" = list(type="region", data="cpg.bed", color="darkgreen"),
#'     "GpC Meth"    = list(type="score",  data="meth.bw",
#'                          color_low="steelblue", color_high="firebrick")
#'   )
#'   ```
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
#' @param annotation_height Numeric. Default relative height unit for
#'   annotation tracks (`peak_bed` and `annotation_tracks` entries that do not
#'   specify their own `height`). Default `0.2`.
#' @param peak_color Character. Fill color for peak rectangles in the
#'   `peak_bed` annotation track. Default `"steelblue4"`.
#' @param vline Optional numeric vector of genomic positions (bp) at which to
#'   draw a vertical reference line through every track and annotation panel.
#'   Each position is drawn as an independent line within its panel, so panels
#'   are not visually connected. Accepts either plain numerics or
#'   `"chr:position"` strings (e.g. `"chr10:74073596"`) — the chromosome is
#'   ignored since all panels share the same region. Default `NULL` (no lines).
#' @param vline_color Character. Colour of the vertical reference line(s).
#'   Default `"black"`.
#' @param vline_type Character. Line type of the vertical reference line(s)
#'   (any value accepted by `linetype`). Default `"dotted"`.
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
                                        annotation_tracks = NULL,
                                        group_colors      = NULL,
                                        scale_y           = "free",
                                        show_variance     = FALSE,
                                        variance_type     = "sd",
                                        track_height      = 1,
                                        annotation_height = 0.2,
                                        peak_color        = "steelblue4",
                                        vline             = NULL,
                                        vline_color       = "black",
                                        vline_type        = "dotted",
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

  # --- Parse vline positions --------------------------------------------------
  vline_pos <- NULL
  if (!is.null(vline)) {
    vline_pos <- vapply(vline, function(v) {
      if (is.character(v)) {
        # Accept "chr:pos" — strip chromosome, keep position
        parts <- strsplit(v, ":")[[1]]
        as.numeric(parts[length(parts)])
      } else {
        as.numeric(v)
      }
    }, numeric(1))
    # Silently drop positions outside the plotted region
    vline_pos <- vline_pos[vline_pos >= reg_start & vline_pos <= reg_end]
    if (length(vline_pos) == 0) vline_pos <- NULL
  }

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

  # --- Build unified annotation track list ------------------------------------
  # peak_bed is kept for backwards compatibility; internally it becomes a
  # "region" type annotation track prepended before annotation_tracks.
  anno_list <- list()

  if (!is.null(peak_bed)) {
    anno_list[["Peaks"]] <- list(
      data   = peak_bed,
      type   = "region",
      color  = peak_color,
      height = annotation_height
    )
  }

  if (!is.null(annotation_tracks)) {
    if (!is.list(annotation_tracks) || is.null(names(annotation_tracks)) ||
        any(names(annotation_tracks) == ""))
      stop("annotation_tracks must be a fully named list", call. = FALSE)
    for (nm in names(annotation_tracks)) {
      at <- annotation_tracks[[nm]]
      if (!is.list(at))
        stop("annotation_tracks[['", nm, "']] must be a list", call. = FALSE)
      if (is.null(at$type))   at$type   <- "region"
      if (is.null(at$height)) at$height <- if (at$type == "score")
        annotation_height * 2 else if (at$type == "genes")
        annotation_height * 3 else annotation_height
      anno_list[[nm]] <- at
    }
  }

  has_annotation <- length(anno_list) > 0

  # --- Build one ggplot per coverage track ------------------------------------
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

  # Remove x-axis from last coverage track when annotation panels follow
  if (has_annotation) {
    track_plots[[n_tracks]] <- track_plots[[n_tracks]] +
      labs(x = NULL) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank())
  }

  # --- Build annotation panels ------------------------------------------------
  n_anno       <- length(anno_list)
  anno_plots   <- vector("list", n_anno)
  anno_heights <- numeric(n_anno)

  for (i in seq_len(n_anno)) {
    nm <- names(anno_list)[i]
    at <- anno_list[[nm]]
    is_last <- i == n_anno

    if (at$type == "region") {
      anno_plots[[i]] <- .AETHER_build_region_track(
        nm, at, chr, reg_start, reg_end, region_lab, is_last
      )
    } else if (at$type == "score") {
      anno_plots[[i]] <- .AETHER_build_score_track(
        nm, at, region_gr, chr, reg_start, reg_end, region_lab, is_last, verbose
      )
    } else if (at$type == "genes") {
      anno_plots[[i]] <- .AETHER_build_genes_track(
        nm, at, region_gr, chr, reg_start, reg_end, region_lab, is_last, verbose
      )
    } else {
      stop("annotation_tracks[['", nm, "']]: unknown type '", at$type,
           "'. Use \"region\", \"score\", or \"genes\".", call. = FALSE)
    }

    anno_heights[i] <- at$height
  }

  # --- Apply vlines to every panel --------------------------------------------
  if (!is.null(vline_pos)) {
    vline_layer <- lapply(vline_pos, function(x)
      geom_vline(xintercept = x, color = vline_color,
                 linetype = vline_type, linewidth = 0.4)
    )
    all_plots_pre <- c(track_plots, anno_plots)
    all_plots <- lapply(all_plots_pre, function(p) {
      for (l in vline_layer) p <- p + l
      p
    })
  } else {
    all_plots <- c(track_plots, anno_plots)
  }

  # --- Assemble with patchwork ------------------------------------------------
  heights <- c(rep(track_height, n_tracks), anno_heights)

  patchwork::wrap_plots(all_plots, ncol = 1, heights = heights)
}


# ------------------------------------------------------------------------------
# Internal: region annotation track (filled rectangles — peaks, CpG islands…)
# ------------------------------------------------------------------------------
.AETHER_build_region_track <- function(nm, at, chr, reg_start, reg_end,
                                       region_lab, is_last) {

  # --- Load interval data -----------------------------------------------------
  if (is.character(at$data) && length(at$data) == 1L) {
    if (!file.exists(at$data))
      stop("annotation_tracks[['", nm, "']]: file not found: ", at$data,
           call. = FALSE)
    raw <- utils::read.table(at$data, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)
  } else {
    raw <- as.data.frame(at$data)
  }

  if (ncol(raw) < 3)
    stop("annotation_tracks[['", nm, "']]: data must have at least 3 columns",
         " (chr, start, end)", call. = FALSE)

  raw <- raw[, 1:3, drop = FALSE]
  colnames(raw) <- c("chr", "start", "end")

  regions <- raw[raw$chr == chr & raw$end >= reg_start & raw$start <= reg_end,
                 , drop = FALSE]
  regions$start <- pmax(regions$start, reg_start)
  regions$end   <- pmin(regions$end,   reg_end)

  color <- if (is.null(at$color)) "steelblue4" else at$color

  # --- Build plot -------------------------------------------------------------
  p <- ggplot() +
    scale_x_continuous(
      limits = c(reg_start, reg_end),
      expand = c(0, 0),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = if (is_last) region_lab else NULL, y = nm) +
    theme_bw() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0.5),
      axis.text.x  = if (is_last) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (is_last) element_line()         else element_blank(),
      axis.title.x = if (is_last) element_text(size = 8) else element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(1, 5, 1, 5)
    )

  if (nrow(regions) > 0) {
    p <- p +
      geom_rect(
        data = regions,
        aes(xmin = .data$start, xmax = .data$end, ymin = 0.2, ymax = 0.8),
        fill = color, color = NA
      )
  }

  p
}


# ------------------------------------------------------------------------------
# Internal: score annotation track (lollipop — per-site methylation etc.)
# ------------------------------------------------------------------------------
.AETHER_build_score_track <- function(nm, at, region_gr, chr, reg_start, reg_end,
                                      region_lab, is_last, verbose) {

  if (is.null(at$data))
    stop("annotation_tracks[['", nm, "']]: 'data' field is required",
         call. = FALSE)

  color_low  <- if (is.null(at$color_low))  "steelblue" else at$color_low
  color_high <- if (is.null(at$color_high)) "firebrick" else at$color_high

  # --- Read BigWig ------------------------------------------------------------
  if (!file.exists(at$data))
    stop("annotation_tracks[['", nm, "']]: file not found: ", at$data,
         call. = FALSE)
  if (verbose) message("[AETHER]   Reading score track [", nm, "]: ",
                       basename(at$data))

  sig <- tryCatch(
    rtracklayer::import.bw(at$data, which = region_gr),
    error = function(e)
      stop("Failed to read BigWig for annotation_tracks[['", nm, "']]: ",
           conditionMessage(e), call. = FALSE)
  )

  df      <- as.data.frame(sig)[, c("start", "end", "score"), drop = FALSE]
  df$pos  <- (df$start + df$end) / 2L

  # --- Score limits (colour scale + y-axis) -----------------------------------
  if (!is.null(at$score_limits)) {
    s_lim <- at$score_limits
  } else {
    s_range <- range(df$score, na.rm = TRUE)
    s_lim   <- c(0, max(s_range[2] * 1.05, s_range[2] + 0.01))
  }

  # --- Build plot -------------------------------------------------------------
  p <- ggplot(df) +
    geom_segment(
      aes(x    = .data$pos,
          xend = .data$pos,
          y    = 0,
          yend = .data$score,
          color = .data$score),
      linewidth = 0.3
    ) +
    geom_point(
      aes(x = .data$pos, y = .data$score, color = .data$score),
      size = 0.7
    ) +
    scale_color_gradient(
      low    = color_low,
      high   = color_high,
      limits = s_lim,
      guide  = "none"
    ) +
    scale_x_continuous(
      limits = c(reg_start, reg_end),
      expand = c(0, 0),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_y_continuous(
      limits = s_lim,
      expand = expansion(mult = c(0, 0.02)),
      breaks = pretty(s_lim, n = 3)
    ) +
    labs(x = if (is_last) region_lab else NULL, y = nm) +
    theme_bw() +
    theme(
      axis.text.y        = element_text(size = 6),
      axis.ticks.y       = element_line(),
      axis.title.y       = element_text(size = 8, angle = 0, hjust = 1, vjust = 0.5),
      axis.text.x        = if (is_last) element_text(size = 8) else element_blank(),
      axis.ticks.x       = if (is_last) element_line()         else element_blank(),
      axis.title.x       = if (is_last) element_text(size = 8) else element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(1, 5, 1, 5)
    )

  p
}


# ------------------------------------------------------------------------------
# Internal: parse region specification to a single GRanges
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# Internal: gene model annotation track (IGV-style collapsed gene models)
# ------------------------------------------------------------------------------
#
# at fields:
#   txdb        — TxDb object (preferred; fast, indexed)
#   gtf_file    — path to GTF/GFF (alternative; reads exons genome-wide, slower)
#   chr_mapping — passed to APOLLO_get_chr_mapping() when using gtf_file
#   show_cds    — distinguish CDS (tall) from UTR (short) blocks; default TRUE
#   label_genes — show gene name labels; default TRUE
#   label_size  — text size for gene labels; default 2.5
#   color       — gene model color; default "black"
# ------------------------------------------------------------------------------
.AETHER_build_genes_track <- function(nm, at, region_gr, chr, reg_start, reg_end,
                                       region_lab, is_last, verbose) {

  if (is.null(at$txdb) && is.null(at$gtf_file))
    stop("annotation_tracks[['", nm, "']]: 'genes' type requires either a ",
         "'txdb' field (TxDb object) or a 'gtf_file' field (path).",
         call. = FALSE)

  color       <- if (is.null(at$color))       "black" else at$color
  label_size  <- if (is.null(at$label_size))  2.5     else at$label_size
  show_cds    <- if (is.null(at$show_cds))    TRUE    else at$show_cds
  label_genes <- if (is.null(at$label_genes)) TRUE    else at$label_genes

  # --- Extract gene models ----------------------------------------------------
  if (!is.null(at$txdb)) {
    gm <- .AETHER_gene_models_txdb(at$txdb, region_gr, show_cds, verbose)
  } else {
    gm <- .AETHER_gene_models_gtf(at$gtf_file, region_gr, chr,
                                   at$chr_mapping, show_cds, verbose)
  }

  genes_df <- gm$genes   # gene_id, gene_name, gene_start, gene_end, strand
  exons_df <- gm$exons   # gene_id, start, end
  cds_df   <- gm$cds     # gene_id, start, end (or NULL)

  # Common x scale used by all layer builds below
  x_sc <- scale_x_continuous(
    limits = c(reg_start, reg_end), expand = c(0, 0),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE)
  )
  base_theme <- theme_bw() + theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0.5),
    axis.text.x  = if (is_last) element_text(size = 8) else element_blank(),
    axis.ticks.x = if (is_last) element_line()         else element_blank(),
    axis.title.x = if (is_last) element_text(size = 8) else element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(1, 5, 1, 5)
  )

  # Empty region — return blank track
  if (is.null(genes_df) || nrow(genes_df) == 0) {
    p <- ggplot() + x_sc +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = if (is_last) region_lab else NULL, y = nm) +
      base_theme
    return(p)
  }

  # --- Clip to region ---------------------------------------------------------
  genes_df$gene_start <- pmax(genes_df$gene_start, reg_start)
  genes_df$gene_end   <- pmin(genes_df$gene_end,   reg_end)

  exons_df$start <- pmax(exons_df$start, reg_start)
  exons_df$end   <- pmin(exons_df$end,   reg_end)
  exons_df <- exons_df[exons_df$end > exons_df$start, , drop = FALSE]

  if (!is.null(cds_df) && nrow(cds_df) > 0) {
    cds_df$start <- pmax(cds_df$start, reg_start)
    cds_df$end   <- pmin(cds_df$end,   reg_end)
    cds_df <- cds_df[cds_df$end > cds_df$start, , drop = FALSE]
  }

  # --- Assign stacking rows (greedy interval scheduling) ----------------------
  buffer           <- (reg_end - reg_start) * 0.015
  genes_df$row     <- .AETHER_assign_gene_rows(genes_df$gene_start,
                                                genes_df$gene_end, buffer)
  row_map          <- setNames(genes_df$row, genes_df$gene_id)
  exons_df$row     <- row_map[exons_df$gene_id]

  if (!is.null(cds_df) && nrow(cds_df) > 0)
    cds_df$row <- row_map[cds_df$gene_id]

  n_rows <- max(genes_df$row)

  # --- Direction arrows (placed only on intron / non-exon backbone) -----------
  arrows_df <- .AETHER_gene_arrows(genes_df, exons_df, reg_start, reg_end)

  # --- Visual constants -------------------------------------------------------
  utr_h <- 0.28   # half-height of UTR exon blocks
  cds_h <- 0.42   # half-height of CDS blocks (taller)
  exon_h <- if (show_cds && !is.null(cds_df) && nrow(cds_df) > 0) utr_h else cds_h

  # --- Build plot -------------------------------------------------------------
  p <- ggplot() +
    # Intron backbone
    geom_segment(
      data = genes_df,
      aes(x = .data$gene_start, xend = .data$gene_end,
          y = .data$row, yend = .data$row),
      color = color, linewidth = 0.5
    ) +
    # Exon blocks
    geom_rect(
      data = exons_df,
      aes(xmin = .data$start, xmax = .data$end,
          ymin = .data$row - exon_h, ymax = .data$row + exon_h),
      fill = color, color = NA
    ) +
    x_sc +
    scale_y_continuous(limits = c(0.2, n_rows + 0.8), expand = c(0, 0)) +
    labs(x = if (is_last) region_lab else NULL, y = nm) +
    base_theme

  # CDS blocks on top (taller than UTR)
  if (show_cds && !is.null(cds_df) && nrow(cds_df) > 0) {
    p <- p + geom_rect(
      data = cds_df,
      aes(xmin = .data$start, xmax = .data$end,
          ymin = .data$row - cds_h, ymax = .data$row + cds_h),
      fill = color, color = NA
    )
  }

  # Direction arrows along intron backbone
  if (!is.null(arrows_df) && nrow(arrows_df) > 0) {
    p <- p + geom_text(
      data    = arrows_df,
      mapping = aes(x = .data$x, y = .data$row, label = .data$arrow),
      size    = 2.2, color = color, alpha = 0.55
    )
  }

  # Gene name labels (italic, just above the backbone)
  if (label_genes) {
    genes_df$label_x <- (genes_df$gene_start + genes_df$gene_end) / 2
    p <- p + geom_text(
      data    = genes_df,
      mapping = aes(x = .data$label_x, y = .data$row + 0.52,
                    label = .data$gene_name),
      size     = label_size, color = color,
      vjust    = 0, fontface = "italic"
    )
  }

  p
}


# Greedy interval scheduling: assign genes to the lowest available row
.AETHER_assign_gene_rows <- function(starts, ends, buffer = 0) {
  n <- length(starts)
  if (n == 0) return(integer(0))
  ord      <- order(starts)
  rows     <- integer(n)
  row_ends <- numeric(0)   # rightmost end in each row

  for (idx in ord) {
    placed <- FALSE
    for (r in seq_along(row_ends)) {
      if (starts[idx] > row_ends[r] + buffer) {
        rows[idx]   <- r
        row_ends[r] <- ends[idx]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      row_ends  <- c(row_ends, ends[idx])
      rows[idx] <- length(row_ends)
    }
  }
  rows
}


# Generate direction arrow positions along intron backbone regions
.AETHER_gene_arrows <- function(genes_df, exons_df, reg_start, reg_end) {
  spacing <- (reg_end - reg_start) / 30   # ~30 candidate positions across window

  do.call(rbind, lapply(seq_len(nrow(genes_df)), function(i) {
    g       <- genes_df[i, ]
    g_exons <- exons_df[exons_df$gene_id == g$gene_id, , drop = FALSE]

    vis_start  <- max(g$gene_start, reg_start)
    vis_end    <- min(g$gene_end,   reg_end)
    if (vis_start >= vis_end) return(NULL)
    seq_from   <- vis_start + spacing / 2
    if (seq_from > vis_end) return(NULL)
    candidates <- seq(seq_from, vis_end, by = spacing)
    if (length(candidates) == 0) return(NULL)

    # Keep only positions not inside an exon block
    if (nrow(g_exons) > 0) {
      in_exon    <- vapply(candidates, function(pos)
        any(pos >= g_exons$start & pos <= g_exons$end), logical(1))
      candidates <- candidates[!in_exon]
    }
    if (length(candidates) == 0) return(NULL)

    data.frame(
      x     = candidates,
      row   = g$row,
      arrow = if (g$strand == "-") "\u25C2" else "\u25B8",  # ◂ or ▸
      stringsAsFactors = FALSE
    )
  }))
}


# Extract gene models from a TxDb object (fast — uses indexed SQL queries)
.AETHER_gene_models_txdb <- function(txdb, region_gr, show_cds, verbose) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("GenomicFeatures required: BiocManager::install('GenomicFeatures')",
         call. = FALSE)

  if (verbose) message("[AETHER]   Loading gene models from TxDb...")

  all_genes <- tryCatch(
    GenomicFeatures::genes(txdb),
    error = function(e) stop("Failed to get genes from TxDb: ",
                             conditionMessage(e), call. = FALSE)
  )
  genes_in_reg <- IRanges::subsetByOverlaps(all_genes, region_gr,
                                             ignore.strand = TRUE)
  if (length(genes_in_reg) == 0)
    return(list(genes = data.frame(), exons = data.frame(), cds = NULL))

  gene_ids   <- names(genes_in_reg)
  # Strip common NCBI "gene-" prefix for cleaner display labels
  gene_names <- sub("^gene-", "", gene_ids)

  genes_df <- data.frame(
    gene_id    = gene_ids,
    gene_name  = gene_names,
    gene_start = GenomicRanges::start(genes_in_reg),
    gene_end   = GenomicRanges::end(genes_in_reg),
    strand     = as.character(GenomicRanges::strand(genes_in_reg)),
    stringsAsFactors = FALSE
  )

  # Exons grouped by gene — unique positions to collapse isoforms
  exons_by_gene <- tryCatch(
    GenomicFeatures::exonsBy(txdb, by = "gene")[gene_ids],
    error = function(e) { warning("Could not retrieve exons: ", e$message); NULL }
  )

  exons_df <- if (!is.null(exons_by_gene) && length(exons_by_gene) > 0) {
    exons_unl <- unlist(exons_by_gene)
    df <- data.frame(
      gene_id = rep(names(exons_by_gene), lengths(exons_by_gene)),
      start   = GenomicRanges::start(exons_unl),
      end     = GenomicRanges::end(exons_unl),
      stringsAsFactors = FALSE
    )
    unique(df)
  } else {
    data.frame(gene_id = character(), start = integer(), end = integer())
  }

  # CDS grouped by gene
  cds_df <- NULL
  if (show_cds) {
    cds_by_gene <- tryCatch(
      GenomicFeatures::cdsBy(txdb, by = "gene")[gene_ids],
      error = function(e) NULL
    )
    if (!is.null(cds_by_gene) && length(cds_by_gene) > 0) {
      cds_unl <- unlist(cds_by_gene)
      cds_df  <- unique(data.frame(
        gene_id = rep(names(cds_by_gene), lengths(cds_by_gene)),
        start   = GenomicRanges::start(cds_unl),
        end     = GenomicRanges::end(cds_unl),
        stringsAsFactors = FALSE
      ))
    }
  }

  list(genes = genes_df, exons = exons_df, cds = cds_df)
}


# Extract gene models from a GTF file (reads genome-wide exons — slower)
.AETHER_gene_models_gtf <- function(gtf_file, region_gr, chr, chr_mapping,
                                     show_cds, verbose) {
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("rtracklayer required: BiocManager::install('rtracklayer')", call. = FALSE)

  if (!file.exists(gtf_file))
    stop("gtf_file not found: ", gtf_file, call. = FALSE)

  if (verbose) message("[AETHER]   Reading gene models from GTF (may be slow)...")

  # Read gene- and exon-level records
  genes_gr <- rtracklayer::import(gtf_file, feature.type = "gene")
  exons_gr <- rtracklayer::import(gtf_file, feature.type = "exon")
  cds_gr   <- if (show_cds)
    rtracklayer::import(gtf_file, feature.type = "CDS") else NULL

  # Apply chromosome renaming
  .rename_gr <- function(gr) {
    if (is.null(chr_mapping) || is.null(gr)) return(gr)
    map_vec <- if (is.character(chr_mapping) && length(chr_mapping) == 1L &&
                   is.null(names(chr_mapping)))
      APOLLO_get_chr_mapping(chr_mapping) else chr_mapping
    cur  <- GenomeInfoDb::seqlevels(gr)
    hits <- intersect(cur, names(map_vec))
    if (length(hits) > 0) {
      rv <- map_vec[hits]; names(rv) <- hits
      gr <- GenomeInfoDb::renameSeqlevels(gr, rv)
    }
    gr
  }
  genes_gr <- .rename_gr(genes_gr)
  exons_gr <- .rename_gr(exons_gr)
  cds_gr   <- .rename_gr(cds_gr)

  # Filter to region chromosome + overlap
  .chr_filter <- function(gr) {
    if (is.null(gr)) return(gr)
    gr[as.character(GenomicRanges::seqnames(gr)) == chr]
  }
  genes_gr <- IRanges::subsetByOverlaps(.chr_filter(genes_gr), region_gr,
                                         ignore.strand = TRUE)
  exons_gr <- .chr_filter(exons_gr)
  cds_gr   <- .chr_filter(cds_gr)

  if (length(genes_gr) == 0)
    return(list(genes = data.frame(), exons = data.frame(), cds = NULL))

  mc      <- as.data.frame(GenomicRanges::mcols(genes_gr))
  sym_col <- intersect(c("gene_name", "gene"), colnames(mc))[1]
  id_col  <- intersect(c("gene_id", "ID"),     colnames(mc))[1]

  gene_names <- if (!is.na(sym_col)) as.character(mc[[sym_col]]) else
    as.character(mc[[if (!is.na(id_col)) id_col else sym_col]])
  gene_ids <- if (!is.na(id_col)) as.character(mc[[id_col]]) else gene_names

  genes_df <- data.frame(
    gene_id    = gene_ids,
    gene_name  = gene_names,
    gene_start = GenomicRanges::start(genes_gr),
    gene_end   = GenomicRanges::end(genes_gr),
    strand     = as.character(GenomicRanges::strand(genes_gr)),
    stringsAsFactors = FALSE
  )
  genes_df <- unique(genes_df)

  # Match exons to genes via gene_id / gene_name attribute
  .extract_feature_df <- function(gr, gene_ids_vec, gene_names_vec) {
    if (is.null(gr) || length(gr) == 0)
      return(data.frame(gene_id=character(), start=integer(), end=integer()))
    fmc     <- as.data.frame(GenomicRanges::mcols(gr))
    id_col  <- intersect(c("gene_id", "ID"),            colnames(fmc))[1]
    sym_col <- intersect(c("gene_name", "gene"),         colnames(fmc))[1]
    # Prefer gene_id for matching; fall back to gene_name
    if (!is.na(id_col) && any(as.character(fmc[[id_col]]) %in% gene_ids_vec)) {
      keep   <- as.character(fmc[[id_col]]) %in% gene_ids_vec
      gid    <- as.character(fmc[[id_col]])[keep]
    } else if (!is.na(sym_col)) {
      keep   <- as.character(fmc[[sym_col]]) %in% gene_names_vec
      # Map gene_name back to gene_id
      nm_map <- setNames(gene_ids_vec, gene_names_vec)
      gid    <- nm_map[as.character(fmc[[sym_col]])[keep]]
    } else {
      return(data.frame(gene_id=character(), start=integer(), end=integer()))
    }
    gr_filt <- gr[keep]
    unique(data.frame(
      gene_id = gid,
      start   = GenomicRanges::start(gr_filt),
      end     = GenomicRanges::end(gr_filt),
      stringsAsFactors = FALSE
    ))
  }

  exons_df <- .extract_feature_df(exons_gr, gene_ids, gene_names)
  cds_df   <- if (show_cds)
    .extract_feature_df(cds_gr, gene_ids, gene_names) else NULL
  if (!is.null(cds_df) && nrow(cds_df) == 0) cds_df <- NULL

  list(genes = genes_df, exons = exons_df, cds = cds_df)
}
