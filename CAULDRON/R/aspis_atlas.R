# ==============================================================================
# ASPIS - Atlas Visualisation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Extends v5 with satellite mini-embeddings in the four device corners.
# Each satellite shows the same UMAP/tSNE as the central scatter, coloured by
# one meta_col entry, using the same colour scale as the corresponding circos
# track.  Up to four satellites (one per corner) can be shown simultaneously.
#
# Track layout (outside → inside):
#   1. Label track      (optional, outermost)
#   2. Metadata tracks  (optional) — one ring per meta_col element.
#   3. Colour ring      — one sector per category, coloured by palette.
#   4. Dendrogram track (optional) — per-sector hclust of cells / bins.
#   5. Central scatter  — UMAP or tSNE via grid viewport.
#   6. Satellite panels — mini scatter in device corners (TL, TR, BL, BR).
# ==============================================================================


# ------------------------------------------------------------------------------
# Helper: density-aware label positions for the central scatter
# ------------------------------------------------------------------------------
# Tries k-means splits (k = 2..k_max) within each colour_by category.
# A split is accepted when the best mean silhouette score >= min_sil.
# Returns a data.frame: x, y, label, size, angle — one row per placed label.
# ------------------------------------------------------------------------------
.aspis_split_labels <- function(emb,            # n×2 matrix (all cells)
                                col_vals,        # character vector, one per cell
                                cats,            # ordered category names
                                pal,             # named colour vector
                                label_size,      # base text size (used when size_range=NULL)
                                min_sil,         # silhouette threshold for split
                                k_max,           # max sub-clusters per category
                                size_scale,      # scale text by sub-cluster size?
                                angle_labels,    # tilt text toward centroid?
                                size_range,      # c(min_sz, max_sz) absolute ggplot units
                                min_cells,       # min cells to place any label
                                nudge) {         # PC2 shift in units of PC2 SD

  if (!requireNamespace("cluster", quietly = TRUE))
    stop("Package 'cluster' is required for label_split. ",
         "Install via: install.packages('cluster')", call. = FALSE)

  # Resolve size range: NULL → derive from label_size
  sz_min <- size_range[1L]
  sz_max <- size_range[2L]

  max_sil_n  <- 500L   # cap for silhouette distance matrix (O(n^2) cost)
  k_max_use  <- max(2L, as.integer(k_max))
  min_cells_use <- max(1L, as.integer(min_cells))

  # Returns list(angle, dx, dy) for a sub-cluster's label placement.
  # angle — PC1 direction (text baseline), constrained to -90..90°
  # dx/dy — nudge offset along PC2 scaled by PC2 SD; PC2 is always
  #          oriented to have a positive-y component so that positive
  #          nudge consistently shifts the label "upward" on screen.
  .cluster_orientation <- function(sub_pts) {
    zero <- list(angle = 0, dx = 0, dy = 0)
    if (nrow(sub_pts) < 3L) return(zero)
    cov_mat <- stats::cov(sub_pts)
    if (anyNA(cov_mat) || abs(det(cov_mat)) < .Machine$double.eps) return(zero)
    eig  <- eigen(cov_mat, symmetric = TRUE)
    pc1  <- eig$vectors[, 1L]
    pc2  <- eig$vectors[, 2L]
    sd2  <- sqrt(max(0, eig$values[2L]))   # spread along PC2

    # Orient PC2 so positive nudge always goes "upward" on screen
    if (pc2[2L] < 0) pc2 <- -pc2

    a <- if (angle_labels) {
      ang <- atan2(pc1[2L], pc1[1L]) * 180 / pi
      if (ang >  90) ang <- ang - 180
      if (ang < -90) ang <- ang + 180
      ang
    } else 0

    list(angle = a,
         dx    = nudge * sd2 * pc2[1L],
         dy    = nudge * sd2 * pc2[2L])
  }

  # Map sub-cluster fraction to a size within [sz_min, sz_max]
  .scaled_size <- function(n_sub, n_cat) {
    if (!size_scale) return(sz_max)
    frac <- sqrt(n_sub / n_cat)   # 0..1; largest sub-cluster → sz_max
    sz_min + (sz_max - sz_min) * frac
  }

  rows <- lapply(cats, function(cat) {
    idx <- which(col_vals == cat)
    pts <- emb[idx, , drop = FALSE]
    n   <- nrow(pts)

    # Whole category below threshold — no label at all
    if (n < min_cells_use) return(NULL)

    # Single centroid fallback
    .single <- function() {
      ori <- .cluster_orientation(pts)
      data.frame(x     = mean(pts[, 1L]) + ori$dx,
                 y     = mean(pts[, 2L]) + ori$dy,
                 label = cat,
                 size  = sz_max,
                 angle = ori$angle,
                 stringsAsFactors = FALSE)
    }

    if (n < 4L) return(.single())

    # Try k = 2..min(k_max, n-1), evaluate mean silhouette on subsampled pts
    sil_idx  <- if (n > max_sil_n) sample(seq_len(n), max_sil_n) else seq_len(n)
    sil_pts  <- pts[sil_idx, , drop = FALSE]
    sil_dist <- stats::dist(sil_pts)

    best_k   <- 1L
    best_sil <- -Inf
    best_km  <- NULL

    for (k in seq(2L, min(k_max_use, n - 1L))) {
      km <- tryCatch(
        stats::kmeans(pts, centers = k, nstart = 10L, iter.max = 50L),
        error = function(e) NULL
      )
      if (is.null(km)) next
      cl_sub <- km$cluster[sil_idx]
      if (length(unique(cl_sub)) < k) next   # degenerate: empty cluster

      sil      <- cluster::silhouette(cl_sub, sil_dist)
      mean_sil <- mean(sil[, "sil_width"])

      if (mean_sil > best_sil) {
        best_sil <- mean_sil
        best_k   <- k
        best_km  <- km
      }
    }

    # Reject split if silhouette below threshold
    if (best_k == 1L || best_sil < min_sil || is.null(best_km))
      return(.single())

    # One row per accepted sub-cluster; drop any below min_cells threshold
    sub_rows <- lapply(seq_len(best_k), function(kk) {
      sub_idx <- which(best_km$cluster == kk)
      n_sub   <- length(sub_idx)
      if (n_sub < min_cells_use) return(NULL)
      sub_pts <- pts[sub_idx, , drop = FALSE]
      ori <- .cluster_orientation(sub_pts)
      data.frame(x     = mean(sub_pts[, 1L]) + ori$dx,
                 y     = mean(sub_pts[, 2L]) + ori$dy,
                 label = cat,
                 size  = .scaled_size(n_sub, n),
                 angle = ori$angle,
                 stringsAsFactors = FALSE)
    })
    sub_rows <- Filter(Negate(is.null), sub_rows)

    # If every sub-cluster was filtered out, fall back to single centroid
    if (length(sub_rows) == 0L) return(.single())
    do.call(rbind, sub_rows)
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame(x=numeric(0), y=numeric(0),
                                            label=character(0), size=numeric(0),
                                            angle=numeric(0), stringsAsFactors=FALSE))
  do.call(rbind, rows)
}


# ------------------------------------------------------------------------------
# Helper: extract a legend grob for a satellite panel
# ------------------------------------------------------------------------------
# Builds a minimal legend-only ggplot for the given meta_track entry and
# extracts the guide-box grob via ggplotGrob().  Returns NULL on failure.
# ------------------------------------------------------------------------------
.aspis_sat_legend_grob <- function(tr, display_name, title_cex, text_cex,
                                    horizontal = FALSE, size = 1) {
  leg_pos   <- if (horizontal) "bottom" else "right"
  leg_theme <- theme_void() +
    theme(
      legend.position    = leg_pos,
      legend.title       = element_text(size  = title_cex * 11,
                                        face  = "bold",
                                        hjust = 0),
      legend.text        = element_text(size  = text_cex * 11),
      legend.background  = element_rect(fill  = NA, colour = NA),
      legend.key         = element_rect(fill  = NA, colour = NA)
    )

  if (!tr$is_continuous) {
    cats     <- names(tr$pal)
    dummy_df <- data.frame(x   = rep(1L, length(cats)),
                           y   = seq_along(cats),
                           cat = factor(cats, levels = cats),
                           stringsAsFactors = FALSE)
    p_leg <- ggplot(dummy_df, aes(x = x, y = y, colour = cat)) +
      geom_point(size = 2) +
      scale_colour_manual(values = tr$pal, name = display_name) +
      guides(colour = guide_legend(override.aes  = list(size = 3),
                                   direction      = if (horizontal) "horizontal" else "vertical",
                                   nrow           = if (horizontal) 1L else NULL)) +
      leg_theme +
      theme(legend.key.size = unit(0.3 * size, "cm"))
  } else {
    dummy_df <- data.frame(x = 0, y = 0, val = mean(tr$expr_range))
    p_leg <- ggplot(dummy_df, aes(x = x, y = y, colour = val)) +
      geom_point(size = 0) +
      scale_colour_gradientn(
        colours = tr$col_ramp(256L),
        limits  = tr$expr_range,
        name    = display_name
      ) +
      guides(colour = guide_colourbar(
        barwidth  = if (horizontal) unit(1.5 * size, "cm") else unit(0.3 * size, "cm"),
        barheight = if (horizontal) unit(0.3 * size, "cm") else unit(1.5 * size, "cm"),
        direction = if (horizontal) "horizontal" else "vertical",
        title.position = "top",
        ticks     = FALSE
      )) +
      leg_theme
  }

  g         <- ggplotGrob(p_leg)
  # ggplot2 < 3.5:  single "guide-box" entry
  # ggplot2 >= 3.5: positional entries "guide-box-right", "guide-box-left", etc.
  # In both cases the active legend is the only non-zeroGrob among them.
  guide_idx <- which(grepl("guide-box", g$layout$name, fixed = TRUE))
  guide_idx <- guide_idx[!vapply(guide_idx,
                                  function(i) inherits(g$grobs[[i]], "zeroGrob"),
                                  logical(1L))]
  if (length(guide_idx) == 0L) return(NULL)
  g$grobs[[guide_idx[1L]]]
}


#' Atlas visualisation: metadata tracks + colour ring + dendrogram + scatter
#'
#' Draws a circos-style figure with a central UMAP/tSNE scatter surrounded by
#' up to three concentric track groups (outside → inside):
#' \enumerate{
#'   \item \strong{Metadata tracks} (optional) — one ring per \code{meta_col}
#'     entry.  Entries can be discrete \code{colData} columns, numeric
#'     \code{colData} columns, or gene names (expression from \code{meta_assay}).
#'     Track style is controlled by \code{meta_style}.
#'   \item \strong{Colour ring} — one sector per category in \code{colour_by},
#'     width proportional to cell/bin count.
#'   \item \strong{Dendrogram track} (optional) — per-sector hierarchical
#'     clustering of cells or bin centroids.
#' }
#'
#' @param sce A \code{SingleCellExperiment} with a UMAP or tSNE embedding and,
#'   ideally, a PCA.
#' @param colour_by Character.  Discrete \code{colData} column for the scatter
#'   colouring and ring sectors.
#' @param dimred Character.  \code{"UMAP"} (default) or \code{"tSNE"}.
#'
#' @param style Character.  \code{"points"} (default), \code{"contour"}, or
#'   \code{"both"}.
#' @param point_size Numeric.  Point size in the central scatter.  Default
#'   \code{0.5}.
#' @param point_alpha Numeric.  Point transparency (0–1).  Default \code{0.5}.
#' @param ring_palette Character vector or \code{NULL}.  Colours for
#'   \code{colour_by} categories — applied consistently to both the central
#'   scatter points and the colour ring segments so they always match.
#'   \code{NULL} uses the built-in 20-colour palette.
#' @param label_clusters Logical.  Overlay centroid labels on the scatter.
#'   One label per \code{colour_by} category, placed at the overall centroid.
#'   Default \code{FALSE}.
#' @param label_size Numeric.  Base centroid label size (ggplot size units).
#'   Default \code{3}.
#' @param label_split Logical.  Enable density-aware sub-cluster labelling.
#'   When \code{TRUE}, each \code{colour_by} category is analysed independently
#'   in 2D embedding space: k-means splits (k = 2 to \code{label_split_k_max})
#'   are tried and the best split is accepted if its mean silhouette score
#'   exceeds \code{label_split_min_sil}.  Accepted splits place one label per
#'   sub-cluster, all showing the same category name.  Overrides
#'   \code{label_clusters} when active.  Default \code{FALSE}.
#' @param label_split_min_sil Numeric (0–1).  Minimum mean silhouette score
#'   required to accept a k > 1 split.  Higher values are stricter (require
#'   more clearly separated sub-clouds); lower values allow splits of weakly
#'   separated distributions.  A value around \code{0.25} is permissive;
#'   \code{0.5} is strict.  Default \code{0.35}.
#' @param label_split_k_max Integer.  Maximum number of sub-clusters to try
#'   per category.  Default \code{5L}.
#' @param label_split_size_scale Logical.  Scale label text size proportionally
#'   to the square-root of the sub-cluster's fraction of the category's cells.
#'   The largest sub-cluster gets the upper bound of \code{label_split_size_range};
#'   smaller sub-clusters are mapped linearly within the range.  Default
#'   \code{TRUE}.
#' @param label_split_size_range Numeric vector of length 2: \code{c(min, max)}
#'   in ggplot text-size units.  Controls the absolute size bounds when
#'   \code{label_split_size_scale = TRUE}.  When \code{FALSE} all labels get
#'   the maximum value.  \code{NULL} (default) derives the range automatically
#'   from \code{label_size}: \code{c(label_size * 0.5, label_size)}.
#' @param label_split_min_cells Integer.  Minimum number of cells a category
#'   (or an individual sub-cluster after splitting) must contain for a label to
#'   be placed.  Categories below this threshold are skipped entirely;
#'   sub-clusters below it are dropped and the remaining ones retain their
#'   labels.  If all sub-clusters of a category are dropped the category falls
#'   back to a single centroid label (provided the category itself meets the
#'   threshold).  Default \code{50L}.
#' @param label_split_nudge Numeric.  Shift each label away from its centroid
#'   along the axis perpendicular to the text (PC2), in units of that
#'   sub-cluster's PC2 standard deviation.  Positive values push the label
#'   "upward" on screen (PC2 is always oriented toward positive y); negative
#'   values push it "downward".  Useful for moving labels off the dense core of
#'   a cluster while keeping the angle unchanged.  Default \code{0} (no shift).
#' @param label_split_angle Logical.  Tilt each label along the principal axis
#'   (PC1) of the sub-cluster's point cloud.  Text is constrained to the
#'   −90°–90° range so it remains readable.  Default \code{TRUE}.
#' @param n_levels Integer or \code{"dynamic"}.  Contour levels.  Default
#'   \code{8}.
#' @param contour_alpha Numeric.  Contour opacity.  Default \code{0.7}.
#' @param contour_adjust Numeric.  KDE bandwidth multiplier for contour
#'   silhouettes.  Values above \code{1} produce smoother, more spread-out
#'   cluster outlines; values below \code{1} tighten the silhouettes to the
#'   data.  Only applies when \code{style} is \code{"contour"} or
#'   \code{"both"}.  Default \code{1} (automatic bandwidth).
#' @param boundary_pad Numeric.  Fractional axis expansion for KDE.  Default
#'   \code{0.10}.
#' @param inner_pad Numeric.  Additional scatter viewport expansion beyond the
#'   data range.  Default \code{0.25}.
#'
#' @param gap_degrees Numeric.  White-space degrees between sectors.  Default
#'   \code{2}.
#' @param inner_scale Numeric (0–1).  Fraction of the innermost ring boundary
#'   used for the scatter viewport.  Default \code{0.95}.
#'
#' @param show_ring Logical.  Draw the colour ring track.  Set to \code{FALSE}
#'   to hide it (useful when only the metadata tracks and/or dendrogram are
#'   needed).  Default \code{TRUE}.
#' @param ring_height Numeric (0–1).  Radial thickness of the colour ring.
#'   Default \code{0.15}.
#' @param ring_labels Logical.  Draw category labels inside the colour ring.
#'   Default \code{TRUE}.
#' @param ring_label_cex Numeric.  Ring label size.  Default \code{0.65}.
#' @param ring_label_col Character.  Ring label text colour.  Default
#'   \code{"white"}.
#' @param ring_gap Numeric (0–1).  Radial thickness of the empty spacer
#'   between the metadata tracks and the colour ring.  Default \code{0.01}.
#' @param ring_border_col Character or \code{NA}.  Border colour for the
#'   colour ring rectangles.  When \code{NA} (default) and
#'   \code{meta_style = "heatmap"}, automatically inherits
#'   \code{meta_border_col}.
#' @param min_seg_frac Numeric (0–1) or \code{NULL}.  Minimum arc fraction
#'   for any sector.  When set, sectors smaller than this floor are expanded
#'   and larger sectors are proportionally compressed so the total arc remains
#'   360°.  The internal x-axis of each sector is unaffected, so dendrograms
#'   and metadata dots fill their sector correctly regardless.  \code{NULL}
#'   (default) keeps sectors strictly proportional to cell count.
#' @param bottom_gap_degrees Numeric.  Angular width (degrees) of the fixed
#'   gap placed at the visual bottom (6 o'clock) of the circos plot.  This gap
#'   is always present and slightly wider than \code{gap_degrees} so that
#'   metadata track labels can be drawn without overlap.  Default \code{15}.
#'   Note: fixing this gap at the bottom overrides the embedding-alignment
#'   rotation (sectors retain their relative order but not absolute orientation).
#'
#' @param show_meta_labels Logical.  Draw a short label for each metadata track
#'   in the bottom gap, at the radial level of that track.  A short tick
#'   connects each label to its track.  Default \code{FALSE}.
#' @param meta_label_names Named character vector or \code{NULL}.  Override
#'   display names for metadata tracks.  Names must match entries in
#'   \code{meta_col}; values are the text shown.  Useful for shortening long
#'   column or gene names.  \code{NULL} uses the original \code{meta_col}
#'   entries verbatim.
#' @param meta_label_cex Numeric.  Text size for metadata track labels.
#'   Default \code{0.65}.
#' @param meta_label_col Character.  Colour for metadata track label text and
#'   connector ticks.  Default \code{"black"}.
#' @param meta_label_lwd Numeric.  Line width of the connector ticks.
#'   Default \code{0.5}.
#'
#' @param show_label_track Logical.  Draw an outermost text label track with
#'   one wrapped label per sector.  Labels are the \code{colour_by} category
#'   names.  Text is wrapped to avoid spilling beyond sector boundaries.
#'   Default \code{FALSE}.
#' @param label_track_height Numeric (0–1).  Radial thickness of the label
#'   track.  Increase to accommodate more lines of wrapped text.  Default
#'   \code{0.12}.
#' @param label_track_cex Numeric.  Label text size.  Default \code{0.55}.
#' @param label_track_col Character.  Label text colour.  Default
#'   \code{"black"}.
#' @param label_track_font Integer.  Font face: \code{1} plain (default),
#'   \code{2} bold, \code{3} italic.
#'
#' @param show_dendrogram Logical.  Draw the per-sector dendrogram track.
#'   Default \code{TRUE}.
#' @param dend_height Numeric (0–1).  Radial thickness of the dendrogram
#'   track.  Default \code{0.15}.
#' @param dend_clip Numeric (0–1].  Fraction of the global maximum dendrogram
#'   height to display.  Values below 1 zoom into the lower portion of the
#'   height range: the track shows only heights 0–\code{dend_clip × max_h},
#'   stretching that region to fill the full track while clipping the root at
#'   the boundary.  Useful for inspecting fine leaf structure.  Default
#'   \code{1} (no clipping).
#' @param dend_power Numeric (> 0).  Power exponent for a non-linear height
#'   transformation applied before drawing.  Compresses the visual prominence
#'   of the root arch without clipping it.  The transformation is
#'   \code{h_new = (h / max_h)^power × max_h}, which is concave for
#'   \code{power < 1}: small heights (leaf forks) expand radially while the
#'   gap between the last cluster-level fork and the root shrinks.  A value
#'   of \code{0.5} (square-root) gives moderate root compression;
#'   \code{0.3} gives stronger compression.  Default \code{1} (no
#'   transformation).
#' @param dend_n_pcs Integer.  PCA dimensions used for cell/bin distances.
#'   Default \code{20}.
#' @param dend_max_cells Integer or \code{NULL}.  Maximum cells per sector
#'   used for dendrogram construction when \strong{not} binning
#'   (\code{n_bin = NULL}).  Sectors exceeding this count are randomly
#'   subsampled before \code{hclust()} to keep memory and runtime manageable
#'   (full \eqn{n \times n} distance matrix).  Has no effect when \code{n_bin}
#'   is set, because in that case clustering operates on bin centroids rather
#'   than individual cells.  \code{NULL} uses all cells.  Default \code{500}.
#' @param dend_col Character.  Dendrogram line colour.  Default
#'   \code{"grey30"}.
#' @param dend_lwd Numeric.  Dendrogram line width.  Default \code{0.5}.
#'
#' @param n_bin Integer or \code{NULL}.  Cells per bin.  \code{NULL} or
#'   \code{1} disables binning (one leaf per cell).  When set, cells are
#'   sorted by PC1 and grouped into consecutive bins; dendrograms and metadata
#'   tracks operate on bin centroids / aggregated values.  Default \code{NULL}.
#'
#' @param meta_col Character vector or \code{NULL}.  One or more track
#'   specifications drawn from outermost to innermost.  Each entry can be:
#'   \itemize{
#'     \item A \strong{discrete} \code{colData} column (character/factor) —
#'       per-bin majority vote, categorical palette.
#'     \item A \strong{numeric} \code{colData} column (integer/double) —
#'       per-bin mean, colour gradient.
#'     \item A \strong{gene name} present in \code{rownames(sce)} — per-bin
#'       mean expression from \code{meta_assay}, colour gradient.
#'   }
#'   \code{NULL} (default) skips all metadata tracks.
#' @param meta_assay Character.  Assay used for gene expression tracks.
#'   Default \code{"logcounts"}.
#' @param meta_palette Named list or \code{NULL}.  Colour specification per
#'   track, keyed by column/gene name.  For discrete tracks: a character
#'   vector.  For continuous/gene tracks: a two-element vector
#'   \code{c(low, high)} for \code{colorRampPalette}.  Missing entries use
#'   defaults (built-in categorical palette or \code{c("grey90", "#2171B5")}).
#' @param meta_height Numeric (0–1).  Radial thickness of each metadata track.
#'   Default \code{0.04}.
#' @param meta_style Character.  Visual style for all metadata tracks.
#'   \code{"dot"} (default) — filled circles; \code{"dot_border"} — filled
#'   circles with an outline; \code{"heatmap"} — filled colour tiles spanning
#'   the full bin width and track height.
#' @param meta_cex Numeric.  Point size for \code{"dot"} and
#'   \code{"dot_border"} styles.  Default \code{0.3}.
#' @param meta_border_col Character, \code{NA}, or \code{NULL}.  Border colour
#'   applied to metadata track elements: outline of dots
#'   (\code{"dot_border"} style) or tile borders (\code{"heatmap"} style).
#'   \code{NULL} or \code{NA} draws no border.  Also inherited by
#'   \code{ring_border_col} when that is \code{NA} and
#'   \code{meta_style = "heatmap"}.  Default \code{"white"}.
#' @param meta_border_lwd Numeric.  Border line width for \code{"dot_border"}
#'   style.  Default \code{0.3}.
#' @param meta_scale_q Numeric vector of length 2.  Quantiles used to define
#'   the low and high anchors of the colour gradient for continuous tracks
#'   (numeric \code{colData} columns or genes).  Values outside this range are
#'   clamped to the nearest anchor colour.  Default \code{c(0, 1)} (full range).
#'   For sparse genes expressed in only a subset of cells, try \code{c(0, 0.95)}
#'   or \code{c(0, 0.99)} to prevent a few zero-expressing cells from washing
#'   out the colour scale.
#'
#' @param show_satellites Logical.  Draw satellite mini-embedding panels in the
#'   device corners.  Set to \code{FALSE} to suppress all satellites.  Default
#'   \code{TRUE}.
#' @param sat_col Character vector or \code{NULL}.  Which \code{meta_col}
#'   entries to display as satellite mini-embeddings in the four device corners
#'   (top-left, top-right, bottom-left, bottom-right).  Must be a subset of
#'   \code{meta_col}; at most four entries are used.  \code{NULL} (default)
#'   uses the first four entries of \code{meta_col}.  If \code{meta_col} is
#'   \code{NULL} no satellites are drawn regardless of this parameter.
#' @param sat_custom Named list or \code{NULL}.  Escape hatch for satellite
#'   panels that aren't a simple metadata/gene scatter -- e.g. the output of
#'   \code{\link{ASPIS_plot_velocity_stream}}.  Names must be one or more of
#'   \code{"TL"}, \code{"TR"}, \code{"BL"}, \code{"BR"}; each value must be a
#'   \code{ggplot} object.  A named corner is drawn exactly as supplied
#'   (structural chrome -- panel border, plot margin, background -- is
#'   stripped to match the other satellites, but the plot's own title and
#'   legend are left untouched, since there is no metadata track to derive a
#'   title from or an automatic legend for) instead of a \code{sat_col}-driven
#'   scatter for that corner.  Any corner not named here still falls back to
#'   the usual \code{sat_col}/\code{meta_col} behaviour.  Default
#'   \code{list()} (no custom corners).
#' @param sat_size Numeric (0–1).  Fraction of the device canvas occupied by
#'   each satellite panel.  Default \code{0.22}.
#' @param sat_legend Logical or character vector.  Controls which satellite
#'   panels display a legend.  \code{FALSE} (default) shows no legends.
#'   \code{TRUE} shows a legend for every satellite.  A character vector of
#'   \code{meta_col} names shows legends only for the named satellites.
#'   Discrete legends are automatically suppressed when the number of
#'   categories exceeds \code{sat_legend_max_cats}.
#' @param sat_legend_max_cats Integer.  Maximum number of categories a discrete
#'   satellite track may have before its legend is automatically hidden.
#'   Continuous/gradient tracks are always shown regardless of this limit.
#'   Default \code{10L}.
#' @param sat_legend_title_cex Numeric.  Legend title size multiplier (×11 pt).
#'   Default \code{0.6}.
#' @param sat_legend_text_cex Numeric.  Legend item text size multiplier
#'   (×11 pt).  Default \code{0.55}.
#' @param sat_legend_pad Numeric (0–1).  Padding between each legend and its
#'   nearest device edge, in npc units.  TL and TR legends are shifted inward
#'   from the left and right edges respectively; BL and BR legends are shifted
#'   upward from the bottom edge.  Default \code{0.01}.
#' @param sat_legend_size Numeric.  Size multiplier applied uniformly to all
#'   legend geometry: colourbar bar dimensions and discrete key size.  Values
#'   above 1 enlarge; below 1 shrink.  Text sizes are unaffected (use
#'   \code{sat_legend_title_cex} / \code{sat_legend_text_cex} for those).
#'   Default \code{1}.
#' @param sat_point_size Numeric.  Point size in satellite scatters.  Default
#'   \code{0.3}.
#' @param sat_point_alpha Numeric.  Point transparency in satellite scatters.
#'   Default \code{0.5}.
#' @param sat_title_cex Numeric.  Title text size for satellite panels,
#'   expressed as a multiplier on the base ggplot size (11 pt).  Default
#'   \code{0.7} (≈ 7.7 pt).
#' @param sat_title_col Character.  Satellite title colour.  Default
#'   \code{"black"}.
#'
#' @param output_file Character or \code{NULL}.  Path for automatic PNG device
#'   management.  \code{NULL} (default) uses the current device.
#' @param width Numeric.  Output width in inches.  Default \code{8}.
#' @param height Numeric.  Output height in inches.  Default \code{8}.
#' @param res Integer.  Output resolution in dpi.  Default \code{150}.
#'
#' @return \code{invisible(NULL)}.  Output is drawn to the current device.
#' @export
ASPIS_plot_atlas <- function(sce,
                                 colour_by,
                                 dimred            = c("UMAP", "tSNE"),

                                 # Central scatter
                                 style             = c("points", "contour", "both"),
                                 point_size        = 0.5,
                                 point_alpha       = 0.5,
                                 ring_palette      = NULL,
                                 label_clusters    = FALSE,
                                 label_size        = 3,
                                 label_split       = FALSE,
                                 label_split_min_sil = 0.35,
                                 label_split_k_max = 5L,
                                 label_split_size_scale = TRUE,
                                 label_split_size_range = NULL,
                                 label_split_min_cells = 50L,
                                 label_split_nudge = 0,
                                 label_split_angle = TRUE,
                                 n_levels          = 8L,
                                 contour_alpha     = 0.7,
                                 contour_adjust    = 1,
                                 boundary_pad      = 0.10,
                                 inner_pad         = 0.25,

                                 # Circos layout
                                 gap_degrees       = 2,
                                 inner_scale       = 0.95,

                                 # Label track (outermost)
                                 show_label_track  = FALSE,
                                 label_track_height = 0.12,
                                 label_track_cex   = 0.55,
                                 label_track_col   = "black",
                                 label_track_font  = 1L,

                                 # Colour ring
                                 show_ring         = TRUE,
                                 ring_height       = 0.15,
                                 ring_labels       = TRUE,
                                 ring_label_cex    = 0.65,
                                 ring_label_col    = "white",
                                 ring_gap          = 0.01,
                                 ring_border_col   = NA,
                                 min_seg_frac      = NULL,
                                 bottom_gap_degrees = 15,

                                 # Metadata track labels (drawn in bottom gap)
                                 show_meta_labels  = FALSE,
                                 meta_label_names  = NULL,
                                 meta_label_cex    = 0.65,
                                 meta_label_col    = "black",
                                 meta_label_lwd    = 0.5,

                                 # Dendrogram
                                 show_dendrogram   = TRUE,
                                 dend_height       = 0.15,
                                 dend_clip         = 1,
                                 dend_power        = 1,
                                 dend_n_pcs        = 20L,
                                 dend_max_cells    = 500L,
                                 dend_col          = "grey30",
                                 dend_lwd          = 0.5,

                                 # Binning
                                 n_bin             = NULL,

                                 # Metadata tracks
                                 meta_col          = NULL,
                                 meta_assay        = "logcounts",
                                 meta_palette      = NULL,
                                 meta_height       = 0.04,
                                 meta_style        = c("dot", "dot_border", "heatmap"),
                                 meta_cex          = 0.3,
                                 meta_border_col   = "white",
                                 meta_border_lwd   = 0.3,
                                 meta_scale_q      = c(0, 1),

                                 # Satellite panels
                                 show_satellites   = TRUE,
                                 sat_col           = NULL,
                                 sat_custom        = list(),
                                 sat_size          = 0.22,
                                 sat_point_size    = 0.3,
                                 sat_point_alpha   = 0.5,
                                 show_sat_title    = TRUE,
                                 sat_title_cex     = 0.7,
                                 sat_title_col     = "black",
                                 sat_legend        = FALSE,
                                 sat_legend_max_cats = 10L,
                                 sat_legend_title_cex = 0.6,
                                 sat_legend_text_cex  = 0.55,
                                 sat_legend_pad    = 0.01,
                                 sat_legend_size   = 1,

                                 # Output
                                 output_file       = NULL,
                                 width             = 8,
                                 height            = 8,
                                 res               = 150) {

  # ── Dependencies ──────────────────────────────────────────────────────────────
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required. Install via: install.packages('circlize')",
         call. = FALSE)

  # ── Input validation ──────────────────────────────────────────────────────────
  if (!colour_by %in% names(colData(sce)))
    stop("colour_by '", colour_by, "' not found in colData(sce).", call. = FALSE)

  if (!is.null(meta_col)) {
    meta_col <- as.character(meta_col)
    bad      <- setdiff(meta_col, c(names(colData(sce)), rownames(sce)))
    if (length(bad))
      stop("meta_col entry/entries not found in colData(sce) or rownames(sce): ",
           paste(bad, collapse = ", "), call. = FALSE)
    if (!meta_assay %in% assayNames(sce))
      stop("meta_assay '", meta_assay, "' not found in assayNames(sce).", call. = FALSE)
  }

  col_vals <- colData(sce)[[colour_by]]
  if (!is.factor(col_vals) && !is.character(col_vals))
    stop("colour_by must reference a discrete (character or factor) column.",
         call. = FALSE)

  dimred     <- match.arg(dimred)
  dimred_key <- if (dimred == "UMAP") "UMAP" else "tSNE"
  if (!dimred_key %in% reducedDimNames(sce))
    stop(dimred, " not found. Run TALOS_run_", tolower(dimred), "() first.",
         call. = FALSE)

  style      <- match.arg(style)
  meta_style <- match.arg(meta_style)

  col_vals_chr <- as.character(col_vals)
  n_bin_use    <- if (is.null(n_bin) || as.integer(n_bin) <= 1L) 1L else as.integer(n_bin)

  # NULL meta_border_col means no border (NA is what circlize expects)
  meta_border_eff <- if (is.null(meta_border_col)) NA else meta_border_col

  # Auto-match ring border to meta border when using heatmap style
  eff_ring_border <- if (is.na(ring_border_col) && meta_style == "heatmap") {
    meta_border_eff
  } else {
    ring_border_col
  }

  # ── Palettes ──────────────────────────────────────────────────────────────────
  cats_alpha <- levels(as.factor(col_vals_chr))
  n_cat      <- length(cats_alpha)
  base_pal   <- if (!is.null(ring_palette)) ring_palette else .aspis_discrete_palette()
  pal        <- setNames(rep(base_pal, length.out = n_cat), cats_alpha)

  # ── Metadata track definitions ────────────────────────────────────────────────
  # Each entry: list(col, is_continuous, vals, pal, col_ramp, expr_range)
  meta_tracks <- NULL
  if (!is.null(meta_col)) {
    meta_tracks <- lapply(meta_col, function(col) {
      cd_val    <- if (col %in% names(colData(sce))) colData(sce)[[col]] else NULL
      is_num_cd <- !is.null(cd_val) && (is.numeric(cd_val) || is.integer(cd_val))

      if (!is.null(cd_val) && !is_num_cd) {
        # ── Discrete colData column ────────────────────────────────────────
        vals    <- as.character(cd_val)
        cats    <- levels(as.factor(vals))
        base_p  <- if (!is.null(meta_palette) && !is.null(meta_palette[[col]])) {
          meta_palette[[col]]
        } else {
          .aspis_discrete_palette()
        }
        pal_vec <- setNames(rep(base_p, length.out = length(cats)), cats)
        list(col = col, is_continuous = FALSE, vals = vals,
             pal = pal_vec, col_ramp = NULL, expr_range = NULL)

      } else {
        # ── Continuous: numeric colData or gene expression ─────────────────
        expr       <- if (is_num_cd) as.numeric(cd_val) else
                      as.numeric(assay(sce, meta_assay)[col, ])
        expr_range <- as.numeric(quantile(expr, probs = meta_scale_q, na.rm = TRUE))
        grad       <- if (!is.null(meta_palette) && !is.null(meta_palette[[col]])) {
          meta_palette[[col]]
        } else {
          c("grey90", "#2171B5")
        }
        list(col = col, is_continuous = TRUE, vals = expr,
             pal = NULL, col_ramp = grDevices::colorRampPalette(grad),
             grad_colors = grad, expr_range = expr_range)
      }
    })
  }

  # ── Satellite column resolution ───────────────────────────────────────────────
  sat_col_use <- NULL
  if (!is.null(meta_col)) {
    if (!is.null(sat_col)) {
      sat_col <- as.character(sat_col)
      bad_sat <- setdiff(sat_col, meta_col)
      if (length(bad_sat))
        stop("sat_col entries not found in meta_col: ",
             paste(bad_sat, collapse = ", "), call. = FALSE)
      sat_col_use <- sat_col[seq_len(min(4L, length(sat_col)))]
    } else {
      sat_col_use <- meta_col[seq_len(min(4L, length(meta_col)))]
    }
  }

  # ── Custom satellite corner resolution ───────────────────────────────────────
  # sat_custom claims specific corners outright; sat_col_use fills whatever
  # corners are left over, in TL, TR, BL, BR order (setdiff() preserves the
  # order of its first argument, so this stays deterministic).
  corner_labels <- c("TL", "TR", "BL", "BR")
  if (length(sat_custom) > 0L) {
    if (is.null(names(sat_custom)) || any(names(sat_custom) == ""))
      stop("sat_custom must be a named list keyed by \"TL\", \"TR\", \"BL\", or \"BR\".",
           call. = FALSE)
    bad_corners <- setdiff(names(sat_custom), corner_labels)
    if (length(bad_corners) > 0L)
      stop("sat_custom names must be one of \"TL\", \"TR\", \"BL\", \"BR\": bad entries: ",
           paste(bad_corners, collapse = ", "), call. = FALSE)
    if (anyDuplicated(names(sat_custom)) > 0L)
      stop("sat_custom has duplicate corner names.", call. = FALSE)
    not_gg <- names(sat_custom)[!vapply(sat_custom, inherits, logical(1L), "ggplot")]
    if (length(not_gg) > 0L)
      stop("sat_custom entries must be ggplot objects: ",
           paste(not_gg, collapse = ", "), call. = FALSE)
  }
  open_corners    <- setdiff(corner_labels, names(sat_custom))
  corner_meta_col <- setNames(vector("list", 4L), corner_labels)
  if (!is.null(sat_col_use) && length(open_corners) > 0L) {
    n_fill <- min(length(sat_col_use), length(open_corners))
    if (n_fill > 0L)
      for (k in seq_len(n_fill)) corner_meta_col[[open_corners[k]]] <- sat_col_use[k]
  }

  # ── Satellite legend column resolution ───────────────────────────────────────
  sat_legend_use <- character(0L)
  if (!isFALSE(sat_legend) && !is.null(sat_col_use) && !is.null(meta_tracks)) {
    candidates <- if (isTRUE(sat_legend)) sat_col_use else
                    intersect(as.character(sat_legend), sat_col_use)
    # Suppress discrete tracks with too many categories
    sat_legend_use <- Filter(function(col) {
      tr <- meta_tracks[[which(meta_col == col)[1L]]]
      tr$is_continuous || length(tr$pal) <= as.integer(sat_legend_max_cats)
    }, candidates)
  }

  # ── Embedding + sector order (CW from top by centroid angle) ─────────────────
  emb    <- reducedDim(sce, dimred_key)
  emb_cx <- mean(emb[, 1L])
  emb_cy <- mean(emb[, 2L])

  centroid_angles <- vapply(cats_alpha, function(cat) {
    idx <- col_vals_chr == cat
    dx  <- mean(emb[idx, 1L]) - emb_cx
    dy  <- mean(emb[idx, 2L]) - emb_cy
    (pi / 2 - atan2(dy, dx)) %% (2 * pi)
  }, numeric(1L))

  ord       <- order(centroid_angles)
  cat_order <- cats_alpha[ord]

  # ── Rotate cat_order to align gap with the embedding ─────────────────────────
  # With start_deg = 270 - G/2 fixed, the first sector starts at 262.5° in
  # circlize convention.  Converting: circlize_angle = 90 - cwt_deg, so the
  # sector starting at 262.5° should have CWT angle = 90 - 262.5 = -172.5 ≡
  # 187.5° ≡ ~3.27 rad.  We rotate cat_order so the sector with CWT nearest
  # that target is first, keeping the relative clockwise order intact.
  target_cwt <- ((90 - (270 - bottom_gap_degrees / 2)) %% 360) * pi / 180
  cwt_sorted <- centroid_angles[ord]
  ang_dist   <- abs((cwt_sorted - target_cwt + pi) %% (2 * pi) - pi)
  best_i     <- which.min(ang_dist)
  if (best_i > 1L)
    cat_order <- c(cat_order[seq(best_i, n_cat)],
                   cat_order[seq_len(best_i - 1L)])

  # ── PCA matrix ───────────────────────────────────────────────────────────────
  has_pca   <- "PCA" %in% reducedDimNames(sce)
  n_pcs_use <- if (has_pca) min(as.integer(dend_n_pcs),
                                ncol(reducedDim(sce, "PCA"))) else ncol(emb)
  if ((show_dendrogram || !is.null(meta_col)) && !has_pca)
    message("PCA not found; using embedding coordinates for cell distances.")
  pca_or_emb <- if (has_pca) reducedDim(sce, "PCA")[, seq_len(n_pcs_use), drop = FALSE
                ] else emb

  # ── Per-sector info ───────────────────────────────────────────────────────────
  sector_info <- lapply(cat_order, function(cat) {
    idx   <- which(col_vals_chr == cat)
    n_sub <- length(idx)
    mat   <- pca_or_emb[idx, , drop = FALSE]

    if (n_sub < 2L)
      return(list(dend = NULL, n_leaves = n_sub, bin_mode = FALSE,
                  leaf_cell_idx = idx, leaf_bins = NULL))

    if (n_bin_use <= 1L) {
      max_n <- if (is.null(dend_max_cells)) Inf else as.integer(dend_max_cells)
      if (n_sub > max_n) {
        keep  <- sample(seq_len(n_sub), max_n)
        idx   <- idx[keep]; mat <- mat[keep, , drop = FALSE]; n_sub <- length(idx)
      }
      if (!show_dendrogram)
        return(list(dend = NULL, n_leaves = n_sub, bin_mode = FALSE,
                    leaf_cell_idx = idx, leaf_bins = NULL))
      dend        <- as.dendrogram(hclust(dist(mat), method = "average"))
      idx_ordered <- idx[order.dendrogram(dend)]
      list(dend = dend, n_leaves = length(idx_ordered), bin_mode = FALSE,
           leaf_cell_idx = idx_ordered, leaf_bins = NULL)

    } else {
      if (n_sub < n_bin_use) {
        min_cells_needed <- n_bin_use * 2L   # need at least 2 bins
        stop("Sector '", cat, "' has only ", n_sub, " cell(s) but n_bin = ", n_bin_use,
             " requires at least ", min_cells_needed, " cells per sector to form 2 or more bins.\n",
             "  Fix: reduce n_bin (current: ", n_bin_use, ") to at most floor(", n_sub, " / 2) = ",
             floor(n_sub / 2L), ", or set n_bin = NULL to disable binning.",
             call. = FALSE)
      }
      pc1_ord    <- order(mat[, 1L])
      idx_sorted <- idx[pc1_ord]; mat_sorted <- mat[pc1_ord, , drop = FALSE]
      n_bins     <- max(2L, ceiling(n_sub / n_bin_use))
      bin_ids    <- pmin(as.integer(ceiling(seq_len(n_sub) / n_bin_use)), n_bins)
      bin_cents  <- t(vapply(seq_len(n_bins), function(b)
        colMeans(mat_sorted[bin_ids == b, , drop = FALSE]), numeric(ncol(mat_sorted))))
      dend_bins    <- as.dendrogram(hclust(dist(bin_cents), method = "average"))
      bin_leaf_ord <- order.dendrogram(dend_bins)
      leaf_bins    <- lapply(bin_leaf_ord, function(b) idx_sorted[bin_ids == b])
      list(dend = dend_bins, n_leaves = n_bins, bin_mode = TRUE,
           leaf_cell_idx = NULL, leaf_bins = leaf_bins)
    }
  })
  names(sector_info) <- cat_order

  n_leaves_vec <- as.integer(vapply(sector_info, function(x) x$n_leaves, numeric(1L)))
  xlim_mat     <- matrix(c(rep(0, n_cat), n_leaves_vec), ncol = 2)

  dend_heights <- vapply(sector_info,
                         function(x) if (is.null(x$dend)) 0 else attr(x$dend, "height"),
                         numeric(1L))
  global_max_h <- if (any(dend_heights > 0)) max(dend_heights) else 1

  # ── Dendrogram height transformation (dend_power) ─────────────────────────────
  # Applies h_new = (h / max_h)^power * max_h before drawing.
  # power < 1 is concave: expands leaf-level forks radially while compressing
  # the gap between the last cluster fork and the root (tall arch).
  # power = 1 (default) is identity; root stays at max_h in both cases.
  if (show_dendrogram && as.numeric(dend_power) != 1 && global_max_h > 0) {
    .pow_t   <- as.numeric(dend_power)
    .max_h_t <- global_max_h
    sector_info <- setNames(lapply(names(sector_info), function(nm) {
      si <- sector_info[[nm]]
      if (!is.null(si$dend))
        si$dend <- dendrapply(si$dend, function(node) {
          h <- attr(node, "height")
          if (!is.null(h) && h > 0)
            attr(node, "height") <- (h / .max_h_t)^.pow_t * .max_h_t
          node
        })
      si
    }), names(sector_info))
    # Recompute in case float rounding shifted the root slightly
    dend_heights <- vapply(sector_info,
                           function(x) if (is.null(x$dend)) 0 else attr(x$dend, "height"),
                           numeric(1L))
    global_max_h <- if (any(dend_heights > 0)) max(dend_heights) else 1
  }

  # ── Sector arc widths (apply min_seg_frac floor if requested) ────────────────
  raw_fracs <- n_leaves_vec / sum(n_leaves_vec)
  if (!is.null(min_seg_frac) && min_seg_frac > 0) {
    adj_fracs <- pmax(raw_fracs, as.numeric(min_seg_frac))
    adj_fracs <- adj_fracs / sum(adj_fracs)   # renormalise to sum = 1
  } else {
    adj_fracs <- raw_fracs
  }
  # sector.width controls arc; xlim stays as actual leaf counts (0..n_leaves)
  # so dendrograms and dots fill their sector correctly regardless of arc scaling

  # ── Ring start angle (fixed so the bottom gap always sits at 6 o'clock) ───────
  # With clock.wise=TRUE, the gap between the last and first sector is placed
  # at start_deg going counter-clockwise, so its centre lands at:
  #   start_deg - (360 - bottom_gap_degrees) = start_deg - 360 + G
  # Setting that centre to 270° (visual bottom):  start_deg = 270 - G/2
  start_deg <- 270 - bottom_gap_degrees / 2

  # ── Gap vector: uniform small gaps, except large gap after last sector ────────
  gap_vec              <- setNames(rep(gap_degrees, n_cat), cat_order)
  gap_vec[cat_order[n_cat]] <- bottom_gap_degrees

  # ── Per-sector arc widths in degrees (for label wrapping) ────────────────────
  # Total arc = 360 - (n-1) small gaps - 1 bottom gap
  total_arc_deg    <- 360 - gap_degrees * (n_cat - 1L) - bottom_gap_degrees
  arc_degrees_vec  <- setNames(adj_fracs * total_arc_deg, cat_order)

  # ── Build central scatter ─────────────────────────────────────────────────────
  # When label_split is active we add our own label layer below, so suppress
  # the simple centroid labels inside .aspis_plot_dimred.
  p_scatter <- .aspis_plot_dimred(
    sce, dimred_key,
    colour_by      = colour_by,
    point_size     = point_size,
    point_alpha    = point_alpha,
    palette        = pal,
    title          = NULL,
    label_clusters = if (label_split) FALSE else label_clusters,
    label_size     = label_size,
    assay_name     = "logcounts",
    ncol           = NULL,
    style          = style,
    n_levels       = n_levels,
    contour_alpha  = contour_alpha,
    contour_adjust = contour_adjust,
    show_legend    = FALSE,
    boundary_pad   = boundary_pad
  )

  x_rng  <- range(emb[, 1L])
  y_rng  <- range(emb[, 2L])
  tpad_x <- diff(x_rng) * (boundary_pad + inner_pad)
  tpad_y <- diff(y_rng) * (boundary_pad + inner_pad)

  p_scatter <- p_scatter +
    coord_cartesian(xlim = c(x_rng[1L] - tpad_x, x_rng[2L] + tpad_x),
                    ylim = c(y_rng[1L] - tpad_y, y_rng[2L] + tpad_y)) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(plot.margin      = unit(c(0, 0, 0, 0), "pt"),
          plot.background  = element_rect(fill = NA, colour = NA),
          panel.background = element_rect(fill = NA, colour = NA),
          panel.border     = element_blank(),
          axis.line        = element_blank(),
          axis.ticks       = element_blank(),
          axis.text        = element_blank())

  # ── Density-aware split labels ────────────────────────────────────────────────
  # When label_clusters is explicitly set to FALSE by the caller, treat it as a
  # master "no labels" switch and suppress split labels too.
  eff_label_split <- label_split && (missing(label_clusters) || isTRUE(label_clusters))
  if (eff_label_split) {
    # Resolve size_range: NULL → derive from label_size
    eff_size_range <- if (!is.null(label_split_size_range)) {
      as.numeric(label_split_size_range)[1:2]
    } else {
      c(label_size * 0.5, label_size)
    }

    lbl_df <- .aspis_split_labels(
      emb          = emb,
      col_vals     = col_vals_chr,
      cats         = cats_alpha,
      pal          = pal,
      label_size   = label_size,
      min_sil      = as.numeric(label_split_min_sil),
      k_max        = as.integer(label_split_k_max),
      size_scale   = isTRUE(label_split_size_scale),
      angle_labels = isTRUE(label_split_angle),
      size_range   = eff_size_range,
      min_cells    = as.integer(label_split_min_cells),
      nudge        = as.numeric(label_split_nudge)
    )
    p_scatter <- p_scatter +
      geom_text(data        = lbl_df,
                aes(x = x, y = y, label = label, size = size, angle = angle),
                colour      = "black",
                fontface    = "bold",
                inherit.aes = FALSE) +
      scale_size_identity(guide = "none")
  }

  # ── Satellite mini-scatter plots ──────────────────────────────────────────────
  # Built here (before device opens) so ggplot construction errors surface early.
  # Returns one entry per corner_labels position (TL, TR, BL, BR) -- NULL for a
  # corner that has neither a sat_custom entry nor a sat_col_use column -- so
  # the printing loop below can index it positionally against corner_x/corner_y.
  sat_plots <- NULL
  if (show_satellites && (!is.null(sat_col_use) || length(sat_custom) > 0L)) {
    sat_plots <- lapply(corner_labels, function(corner) {

      if (!is.null(sat_custom[[corner]])) {
        # ── User-supplied plot (e.g. ASPIS_plot_velocity_stream()) ───────────
        # Only structural chrome is stripped so it sits flush in its corner
        # like the other satellites; the caller's own title/legend (if any)
        # are left untouched -- there's no metadata track here to derive a
        # display name or an automatic legend from.
        return(
          sat_custom[[corner]] +
            theme(
              plot.margin      = unit(c(0, 0, 0, 0), "pt"),
              plot.background  = element_rect(fill = NA, colour = NA),
              panel.background = element_rect(fill = NA, colour = NA),
              panel.border     = element_blank()
            )
        )
      }

      col <- corner_meta_col[[corner]]
      if (is.null(col)) return(NULL)

      tr_idx <- which(meta_col == col)[1L]
      tr     <- meta_tracks[[tr_idx]]

      display_name <- if (!is.null(meta_label_names) && !is.null(meta_label_names[[col]]))
        meta_label_names[[col]] else col

      sat_pal <- if (!tr$is_continuous) tr$pal else tr$grad_colors

      p <- .aspis_plot_dimred(
        sce, dimred_key,
        colour_by      = col,
        point_size     = sat_point_size,
        point_alpha    = sat_point_alpha,
        palette        = sat_pal,
        title          = NULL,
        label_clusters = FALSE,
        label_size     = label_size,
        assay_name     = meta_assay,
        ncol           = NULL,
        style          = "points",
        n_levels       = n_levels,
        contour_alpha  = contour_alpha,
        show_legend    = FALSE,
        boundary_pad   = boundary_pad
      )

      p +
        coord_cartesian(xlim = c(x_rng[1L] - tpad_x, x_rng[2L] + tpad_x),
                        ylim = c(y_rng[1L] - tpad_y, y_rng[2L] + tpad_y)) +
        labs(title = if (show_sat_title) display_name else NULL, x = NULL, y = NULL) +
        theme(
          plot.title       = element_text(size   = sat_title_cex * 11,
                                          hjust  = 0.5,
                                          colour = sat_title_col,
                                          margin = margin(b = 1, t = 2)),
          plot.margin      = unit(c(0, 0, 0, 0), "pt"),
          plot.background  = element_rect(fill = NA, colour = NA),
          panel.background = element_rect(fill = NA, colour = NA),
          panel.border     = element_blank(),
          axis.line        = element_blank(),
          axis.ticks       = element_blank(),
          axis.text        = element_blank(),
          legend.position  = "none"
        )
    })
  }

  # ── Satellite legend grobs ────────────────────────────────────────────────────
  # Built before device opens; one entry per corner_labels position (NULL if no
  # legend for that corner). Indexed by corner rather than sat_col_use position,
  # since sat_custom may have shifted metadata columns out of their "natural"
  # TL/TR/BL/BR slot. Custom-plot corners never get an auto legend here -- there
  # is no metadata track to build one from -- so any legend for those must be
  # baked into the ggplot object passed via sat_custom.
  sat_legend_grobs <- vector("list", 4L)
  if (length(sat_legend_use) > 0L) {
    for (.i in seq_along(corner_labels)) {
      .col <- corner_meta_col[[corner_labels[.i]]]
      if (is.null(.col) || !.col %in% sat_legend_use) next
      .tr   <- meta_tracks[[which(meta_col == .col)[1L]]]
      .dn   <- if (!is.null(meta_label_names) && !is.null(meta_label_names[[.col]]))
                 meta_label_names[[.col]] else .col
      # Bottom satellites (slots 3 and 4 = BL, BR) get horizontal colourbar/legend
      .horiz <- .i >= 3L
      sat_legend_grobs[[.i]] <- tryCatch(
        .aspis_sat_legend_grob(.tr, .dn, sat_legend_title_cex, sat_legend_text_cex,
                               horizontal = .horiz,
                               size       = as.numeric(sat_legend_size)),
        error = function(e) { warning("Could not build legend for '", .col, "': ", e$message); NULL }
      )
    }
  }

  # ── Device management ─────────────────────────────────────────────────────────
  if (!is.null(output_file)) {
    png(output_file, width = width, height = height, units = "in", res = res,
        bg = "white")
    on.exit({ try(circlize::circos.clear(), silent = TRUE); dev.off() })
  } else {
    on.exit(try(circlize::circos.clear(), silent = TRUE))
  }
  old_par <- par(bty = "n")
  on.exit(par(old_par), add = TRUE)

  # ── Circos initialisation ─────────────────────────────────────────────────────
  circlize::circos.par(
    canvas.xlim             = c(-1, 1),
    canvas.ylim             = c(-1, 1),
    start.degree            = start_deg,
    clock.wise              = TRUE,
    gap.after               = gap_vec,
    track.margin            = c(0, 0),
    cell.padding            = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )
  circlize::circos.initialize(factors      = cat_order,
                              xlim         = xlim_mat,
                              sector.width = adj_fracs)

  # ── Label track (outermost) ───────────────────────────────────────────────────
  if (show_label_track) {
    arc_deg_local      <- arc_degrees_vec
    lt_cex_local       <- label_track_cex
    lt_col_local       <- label_track_col
    lt_font_local      <- label_track_font

    circlize::circos.track(
      ylim         = c(0, 1),
      track.height = label_track_height,
      bg.border    = NA,
      bg.col       = NA,
      panel.fun    = function(x, y) {
        si  <- circlize::get.cell.meta.data("sector.index")
        xl  <- circlize::get.cell.meta.data("cell.xlim")
        yl  <- circlize::get.cell.meta.data("cell.ylim")

        # Estimate chars per line: ~2.2 chars per degree at cex=1 (empirical)
        arc_deg    <- arc_deg_local[[si]]
        chars_line <- max(1L, floor(arc_deg * 2.2 / lt_cex_local))
        lines      <- strwrap(si, width = chars_line)
        n_lines    <- length(lines)

        # Place lines evenly within the track (outer → inner)
        y_top    <- yl[2L] - 0.05 * diff(yl)
        y_bot    <- yl[1L] + 0.05 * diff(yl)
        y_pos    <- if (n_lines == 1L) mean(yl) else
                    seq(y_top, y_bot, length.out = n_lines)

        for (i in seq_len(n_lines)) {
          circlize::circos.text(
            mean(xl), y_pos[i],
            labels     = lines[i],
            facing     = "bending.inside",
            niceFacing = TRUE,
            cex        = lt_cex_local,
            col        = lt_col_local,
            font       = lt_font_local
          )
        }
      }
    )
  }

  # ── Metadata tracks ───────────────────────────────────────────────────────────
  if (!is.null(meta_tracks)) {
    si_info_local    <- sector_info
    meta_style_local <- meta_style
    meta_cex_local   <- meta_cex
    meta_bcol_local  <- meta_border_eff
    meta_blwd_local  <- meta_border_lwd

    for (.tr in meta_tracks) {
      local({
        tr <- .tr
        circlize::circos.track(
          ylim         = c(0, 1),
          track.height = meta_height,
          bg.border    = NA,
          bg.col       = NA,
          panel.fun    = function(x, y) {
            si   <- circlize::get.cell.meta.data("sector.index")
            info <- si_info_local[[si]]
            yl   <- circlize::get.cell.meta.data("cell.ylim")
            n_i  <- info$n_leaves
            x_pos <- seq_len(n_i) - 0.5

            # Resolve per-leaf fill colours
            if (!tr$is_continuous) {
              meta_i <- if (!info$bin_mode) {
                tr$vals[info$leaf_cell_idx]
              } else {
                vapply(info$leaf_bins, function(ci) {
                  tbl <- table(tr$vals[ci]); names(tbl)[which.max(tbl)]
                }, character(1L))
              }
              fills <- tr$pal[meta_i]
            } else {
              leaf_expr <- if (!info$bin_mode) {
                tr$vals[info$leaf_cell_idx]
              } else {
                vapply(info$leaf_bins, function(ci) mean(tr$vals[ci], na.rm = TRUE),
                       numeric(1L))
              }
              span  <- tr$expr_range[2L] - tr$expr_range[1L]
              norm  <- if (span > 0) (leaf_expr - tr$expr_range[1L]) / span else
                       rep(0.5, n_i)
              norm  <- pmin(pmax(norm, 0), 1)
              ramp  <- tr$col_ramp(256L)
              fills <- ramp[pmax(1L, pmin(256L, as.integer(norm * 255L) + 1L))]
            }

            if (meta_style_local == "heatmap") {
              for (i in seq_len(n_i))
                circlize::circos.rect(i - 1L, yl[1L], i, yl[2L],
                                      col = fills[i], border = meta_bcol_local)
            } else if (meta_style_local == "dot_border") {
              old_lwd <- graphics::par(lwd = meta_blwd_local)
              circlize::circos.points(x_pos, rep(0.5, n_i),
                                      col = meta_bcol_local, bg = fills,
                                      pch = 21L, cex = meta_cex_local)
              graphics::par(old_lwd)
            } else {
              circlize::circos.points(x_pos, rep(0.5, n_i),
                                      col = fills, pch = 16L, cex = meta_cex_local)
            }
          }
        )
      })
    }

    # Spacer between metadata tracks and colour ring
    if (show_ring && ring_gap > 0)
      circlize::circos.track(ylim = c(0, 1), track.height = ring_gap,
                             bg.border = NA, bg.col = NA,
                             panel.fun = function(x, y) {})
  }

  # ── Colour ring ───────────────────────────────────────────────────────────────
  if (show_ring) {
    pal_local        <- pal
    ring_border_eff  <- eff_ring_border
    show_labels      <- ring_labels
    label_cex        <- ring_label_cex
    label_col        <- ring_label_col

    circlize::circos.track(
      ylim         = c(0, 1),
      track.height = ring_height,
      bg.border    = NA,
      bg.col       = NA,
      panel.fun    = function(x, y) {
        si <- circlize::get.cell.meta.data("sector.index")
        xl <- circlize::get.cell.meta.data("cell.xlim")
        yl <- circlize::get.cell.meta.data("cell.ylim")
        circlize::circos.rect(xl[1L], yl[1L], xl[2L], yl[2L],
                              col = pal_local[si], border = ring_border_eff)
        if (show_labels)
          circlize::circos.text(mean(xl), mean(yl), labels = si,
                                facing = "bending.inside", niceFacing = TRUE,
                                cex = label_cex, col = label_col, font = 2L)
      }
    )
  }

  # ── Dendrogram track ──────────────────────────────────────────────────────────
  if (show_dendrogram) {
    si_info_dend <- sector_info
    # dend_clip clips the root: only heights 0..clip_max_h are shown.
    # Any branch above clip_max_h is drawn beyond the track ylim and circlize
    # clips it at the cell boundary, shrinking the tall root visually.
    clip_max_h   <- global_max_h * max(min(as.numeric(dend_clip), 1), 1e-6)

    circlize::circos.track(
      ylim         = c(0, clip_max_h),
      track.height = dend_height,
      bg.border    = NA,
      bg.col       = NA,
      panel.fun    = function(x, y) {
        si     <- circlize::get.cell.meta.data("sector.index")
        dend_i <- si_info_dend[[si]]$dend
        if (!is.null(dend_i))
          circlize::circos.dendrogram(dend_i, facing = "outside",
                                      max_height = clip_max_h)
      }
    )
  }

  circlize::circos.clear()
  graphics::box(col = "white", lwd = 3)

  # ── Metadata track labels (drawn in bottom gap using base graphics) ───────────
  # Labels sit at the left side of the bottom gap, at each track's radial level.
  # Radial midpoint of meta track i (1-indexed, outermost first):
  #   r_i = 1 - label_space - (i - 0.5) * meta_height
  # At the bottom (270°), the left edge of the gap (where the first sector
  # starts) is at canvas coords: x = -r * sin(G/2),  y = -r * cos(G/2).
  if (show_meta_labels && !is.null(meta_tracks)) {
    label_space_l <- if (show_label_track) label_track_height else 0
    g_half_rad    <- (bottom_gap_degrees / 2) * pi / 180

    display_names <- vapply(seq_along(meta_tracks), function(i) {
      col_i <- meta_tracks[[i]]$col
      if (!is.null(meta_label_names) && !is.null(meta_label_names[[col_i]]))
        meta_label_names[[col_i]]
      else
        col_i
    }, character(1L))

    for (i in seq_along(meta_tracks)) {
      r_i    <- 1 - label_space_l - (i - 0.5) * meta_height
      # Left boundary of the gap (where the first sector starts).
      # At 270° - G/2: x = -r*sin(G/2),  y = -r*cos(G/2).
      # The tick goes rightward INTO the gap; text is left-aligned so it
      # also reads rightward, staying within the empty gap space.
      x_edge     <- -r_i * sin(g_half_rad)
      y_edge     <- -r_i * cos(g_half_rad)
      x_tick_end <- x_edge + 0.01        # one step into the gap (rightward)
      x_text     <- x_tick_end + 0.01    # label starts just past the tick
      graphics::segments(x0  = x_edge,     y0 = y_edge,
                         x1  = x_tick_end,  y1 = y_edge,
                         col = meta_label_col, lwd = meta_label_lwd)
      graphics::text(x      = x_text, y      = y_edge,
                     labels = display_names[i],
                     adj    = c(0, 0.5),
                     cex    = meta_label_cex,
                     col    = meta_label_col)
    }
  }

  # ── Central scatter via grid viewport ─────────────────────────────────────────
  label_space   <- if (show_label_track) label_track_height else 0
  ring_gap_used <- if (!is.null(meta_tracks) && show_ring) ring_gap else 0
  meta_space    <- if (!is.null(meta_tracks)) length(meta_tracks) * meta_height +
                                              ring_gap_used else 0
  ring_space    <- if (show_ring) ring_height else 0
  dend_space    <- if (show_dendrogram) dend_height else 0
  vp_size       <- (1 - label_space - meta_space - ring_space - dend_space) * inner_scale

  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
    width = grid::unit(vp_size, "npc"), height = grid::unit(vp_size, "npc")
  ))
  print(p_scatter, newpage = FALSE)
  grid::popViewport()

  # ── Satellite viewports (device corners: TL, TR, BL, BR) ──────────────────────
  # Legend placement (semi-asymmetric):
  #   TL(1): legend below satellite, left-pinned  → hangs from bottom-left
  #   TR(2): legend below satellite, right-pinned → hangs from bottom-right
  #   BL(3): legend right of satellite (inner edge), bottom-pinned → toward centre
  #   BR(4): legend left of satellite (inner edge), bottom-pinned → toward centre
  if (!is.null(sat_plots)) {
    corner_x     <- c(0, 1, 0, 1)
    corner_y     <- c(1, 1, 0, 0)
    corner_justx <- c("left", "right", "left", "right")
    corner_justy <- c("top",  "top",   "bottom", "bottom")

    # Legend anchor positions in npc (x, y, justx, justy)
    # sat_legend_pad shifts each legend away from its nearest hard edge:
    #   TL → rightward from left edge  |  TR → leftward from right edge
    #   BL/BR → upward from bottom edge
    .lp   <- as.numeric(sat_legend_pad)
    leg_x <- c(.lp,             1 - .lp,       sat_size,       1 - sat_size)
    leg_y <- c(1 - sat_size,    1 - sat_size,  .lp,            .lp)
    leg_jx   <- c("left",      "right",       "left",         "right")
    leg_jy   <- c("top",       "top",         "bottom",       "bottom")

    for (i in seq_along(sat_plots)) {
      if (is.null(sat_plots[[i]])) next

      # Scatter panel
      grid::pushViewport(grid::viewport(
        x      = grid::unit(corner_x[i],    "npc"),
        y      = grid::unit(corner_y[i],    "npc"),
        width  = grid::unit(sat_size, "npc"),
        height = grid::unit(sat_size, "npc"),
        just   = c(corner_justx[i], corner_justy[i])
      ))
      print(sat_plots[[i]], newpage = FALSE)
      grid::popViewport()

      # Legend panel (if present for this slot)
      leg_grob <- if (i <= length(sat_legend_grobs)) sat_legend_grobs[[i]] else NULL
      if (!is.null(leg_grob)) {
        grid::pushViewport(grid::viewport(
          x      = grid::unit(leg_x[i], "npc"),
          y      = grid::unit(leg_y[i], "npc"),
          width  = grid::grobWidth(leg_grob),
          height = grid::grobHeight(leg_grob),
          just   = c(leg_jx[i], leg_jy[i])
        ))
        grid::grid.draw(leg_grob)
        grid::popViewport()
      }
    }
  }

  invisible(NULL)
}
