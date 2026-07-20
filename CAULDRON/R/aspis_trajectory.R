# ==============================================================================
# ASPIS - Trajectory & Velocity Visualisation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Plots for RNA velocity, pseudotime, and trajectory results.
#   - ASPIS_plot_velocity        per-cell coloured arrows + origin circles
#   - ASPIS_plot_velocity_stream KDE blob background + Euler streamlines
#   - ASPIS_plot_pseudotime      pseudotime gradient on embedding
#   - ASPIS_plot_trajectory      principal graph overlay on embedding
#
# All public functions prefixed: ASPIS_
# Internal helpers prefixed: .aspis_
# ==============================================================================


# ==============================================================================
# Internal shared helpers
# ==============================================================================

# Validate SCE, resolve colour_by, build base data frame
.aspis_vel_prep <- function(sce, dimred, colour_by, assay_name) {
  if (!isTRUE(metadata(sce)$velocity_run))
    stop("[ASPIS] TRIPODES_run_velocity() has not been run on this SCE.",
         call. = FALSE)

  vel_name <- paste0("velocity_", dimred)

  if (!dimred %in% reducedDimNames(sce))
    stop("[ASPIS] Embedding '", dimred, "' not found. Available: ",
         paste(reducedDimNames(sce), collapse = ", "), ".", call. = FALSE)

  if (!vel_name %in% reducedDimNames(sce))
    stop("[ASPIS] Velocity embedding '", vel_name, "' not found. ",
         "Re-run TRIPODES_run_velocity() with '", dimred,
         "' present at that time.", call. = FALSE)

  is_gene <- colour_by %in% rownames(sce)
  is_meta <- colour_by %in% names(colData(sce))
  if (!is_gene && !is_meta)
    stop("[ASPIS] '", colour_by,
         "' not found in colData(sce) or rownames(sce).", call. = FALSE)

  colour_vals <- if (is_gene) {
    as.numeric(assay(sce, assay_name)[colour_by, ])
  } else {
    colData(sce)[[colour_by]]
  }

  coords <- as.data.frame(reducedDim(sce, dimred))
  vel    <- as.data.frame(reducedDim(sce, vel_name))
  colnames(coords) <- c("dim1", "dim2")
  colnames(vel)    <- c("vel1", "vel2")

  data.frame(
    dim1   = coords$dim1,
    dim2   = coords$dim2,
    vel1   = vel$vel1,
    vel2   = vel$vel2,
    colour = colour_vals,
    stringsAsFactors = FALSE
  )
}

# Build colour scale and named palette from a prepared data frame
.aspis_build_cs <- function(df, palette, colour_by) {
  is_discrete <- is.factor(df$colour) || is.character(df$colour)
  if (is_discrete) {
    lvls <- if (is.factor(df$colour)) levels(df$colour) else sort(unique(df$colour))
    pal  <- if (!is.null(palette)) palette else .aspis_discrete_palette()
    pal  <- stats::setNames(rep(pal, length.out = length(lvls)), lvls)
    list(
      is_discrete = TRUE,
      pal         = pal,
      scale       = scale_colour_manual(values = pal, name = colour_by)
    )
  } else {
    list(
      is_discrete = FALSE,
      pal         = NULL,
      scale       = if (!is.null(palette))
        scale_colour_gradientn(colours = palette, name = colour_by)
      else
        scale_colour_viridis_c(name = colour_by)
    )
  }
}

# Shared theme for all velocity plots
.aspis_vel_theme <- function() {
  theme_classic() +
    theme(
      legend.position = "right",
      plot.title      = element_text(size = 11, face = "bold"),
      axis.title      = element_text(size = 9),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.line       = element_blank(),
      panel.border    = element_rect(colour = "grey70", fill = NA,
                                     linewidth = 0.5)
    )
}


# ==============================================================================
# Velocity field helpers (for stream plot)
# ==============================================================================

# Compute Gaussian-weighted velocity field on a regular grid.
# Returns list(grid_x, grid_y, dx, dy, has_vel, xlim, ylim).
#
# `ws` (the raw, unnormalised sum of kernel weights reaching a grid point) is
# also a local density estimate. A Gaussian kernel has infinite support, so
# with thousands of cells `ws` clears a bare numerical-epsilon threshold at
# essentially every grid point -- including corners of the bounding box with
# no nearby cells at all, on a non-convex embedding shape. That let both
# has_vel and dx/dy be "valid" (but noise-driven, dominated by whichever
# handful of cells are least-far-away) in genuinely empty space, so seeds
# landed there and .aspis_interp_vel() never saw an NA to stop a streamline
# at. min_density masks grid points whose local density falls below a
# fraction of the peak grid density (mirrors scVelo's own min_mass masking
# in velocity_embedding_stream) -- both has_vel AND dx/dy are masked so
# streamlines actually terminate at the edge of real data support instead of
# wandering through empty margins.
.aspis_velocity_field <- function(x, y, vx, vy, grid_res, bandwidth,
                                   min_density = 0.05) {
  pad_x <- diff(range(x)) * 0.05
  pad_y <- diff(range(y)) * 0.05
  xlim  <- range(x) + c(-pad_x,  pad_x)
  ylim  <- range(y) + c(-pad_y,  pad_y)
  gx    <- seq(xlim[1], xlim[2], length.out = grid_res)
  gy    <- seq(ylim[1], ylim[2], length.out = grid_res)

  if (is.null(bandwidth))
    bandwidth <- (diff(xlim) + diff(ylim)) / 2 / grid_res

  h2      <- 2 * bandwidth^2
  dx_mat  <- matrix(NA_real_, grid_res, grid_res)
  dy_mat  <- matrix(NA_real_, grid_res, grid_res)
  ws_mat  <- matrix(0,        grid_res, grid_res)
  has_vel <- matrix(FALSE,    grid_res, grid_res)

  for (i in seq_len(grid_res)) {
    for (j in seq_len(grid_res)) {
      d2 <- (x - gx[i])^2 + (y - gy[j])^2
      w  <- exp(-d2 / h2)
      ws <- sum(w)
      ws_mat[i, j] <- ws
      if (ws > .Machine$double.eps * 100) {
        dx_mat[i, j]  <- sum(w * vx) / ws
        dy_mat[i, j]  <- sum(w * vy) / ws
        has_vel[i, j] <- TRUE
      }
    }
  }

  # Mask out low-density grid points -- NA the vectors (not just has_vel) so
  # streamline interpolation actually stops there instead of extrapolating.
  density_thresh <- min_density * max(ws_mat)
  low_density    <- ws_mat < density_thresh
  dx_mat[low_density]  <- NA_real_
  dy_mat[low_density]  <- NA_real_
  has_vel[low_density] <- FALSE

  list(grid_x  = gx, grid_y  = gy,
       dx      = dx_mat, dy      = dy_mat,
       has_vel = has_vel,
       xlim    = xlim,   ylim    = ylim)
}

# Bilinear interpolation of velocity at point (px, py).
.aspis_interp_vel <- function(px, py, field) {
  gx <- field$grid_x;  gy <- field$grid_y
  ix <- findInterval(px, gx, rightmost.closed = TRUE)
  iy <- findInterval(py, gy, rightmost.closed = TRUE)
  if (ix < 1L || ix >= length(gx) || iy < 1L || iy >= length(gy))
    return(c(NA_real_, NA_real_))
  tx <- (px - gx[ix]) / (gx[ix + 1L] - gx[ix])
  ty <- (py - gy[iy]) / (gy[iy + 1L] - gy[iy])
  blerp <- function(mat) {
    (1 - tx) * (1 - ty) * mat[ix,      iy     ] +
    tx        * (1 - ty) * mat[ix + 1L, iy     ] +
    (1 - tx)  * ty       * mat[ix,      iy + 1L] +
    tx        * ty       * mat[ix + 1L, iy + 1L]
  }
  c(blerp(field$dx), blerp(field$dy))
}

# Trace a single streamline by Euler integration (unit-normalised direction).
# Returns data.frame(x, y) or NULL if < 3 steps.
.aspis_trace_streamline <- function(x0, y0, field, dt, max_steps, min_vel_sq) {
  px <- numeric(max_steps + 1L)
  py <- numeric(max_steps + 1L)
  px[1L] <- x0;  py[1L] <- y0
  n <- 1L
  for (i in seq_len(max_steps)) {
    v   <- .aspis_interp_vel(px[n], py[n], field)
    if (any(is.na(v))) break
    vsq <- v[1L]^2 + v[2L]^2
    if (vsq <= min_vel_sq) break
    vn   <- v / sqrt(vsq)           # unit direction
    nx   <- px[n] + dt * vn[1L]
    ny   <- py[n] + dt * vn[2L]
    if (nx < field$xlim[1L] || nx > field$xlim[2L] ||
        ny < field$ylim[1L] || ny > field$ylim[2L]) break
    n      <- n + 1L
    px[n]  <- nx;  py[n] <- ny
  }
  if (n < 3L) return(NULL)
  data.frame(x = px[seq_len(n)], y = py[seq_len(n)])
}

# Select n_seeds seed points from non-empty grid cells (with light jitter).
.aspis_stream_seeds <- function(field, n_seeds) {
  idx <- which(field$has_vel, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(matrix(NA_real_, 0L, 2L))
  sel     <- idx[sample(nrow(idx), min(n_seeds, nrow(idx))), , drop = FALSE]
  dx_step <- diff(field$grid_x[1:2]) * 0.4
  dy_step <- diff(field$grid_y[1:2]) * 0.4
  cbind(
    field$grid_x[sel[, 1L]] + stats::runif(nrow(sel), -dx_step, dx_step),
    field$grid_y[sel[, 2L]] + stats::runif(nrow(sel), -dy_step, dy_step)
  )
}

# Trace all streamlines; optionally attach per-stream colour from seed_colours.
# seed_colours: character vector length nrow(seeds), or NULL for fixed colour.
.aspis_trace_all_streams <- function(seeds, field, max_steps, dt,
                                      min_vel_frac, seed_colours) {
  vel_sq_vals <- c(field$dx^2 + field$dy^2)[field$has_vel]
  min_vel_sq  <- if (length(vel_sq_vals) > 0L)
    (min_vel_frac^2) * mean(vel_sq_vals) else 0

  stream_list <- vector("list", nrow(seeds))
  id <- 0L
  for (s in seq_len(nrow(seeds))) {
    tr <- .aspis_trace_streamline(seeds[s, 1L], seeds[s, 2L],
                                   field, dt, max_steps, min_vel_sq)
    if (!is.null(tr)) {
      id         <- id + 1L
      tr$stream_id <- id
      if (!is.null(seed_colours)) tr$stream_colour <- seed_colours[s]
      stream_list[[s]] <- tr
    }
  }
  do.call(rbind, Filter(Negate(is.null), stream_list))
}


# ==============================================================================
# ASPIS_plot_velocity — per-cell coloured arrows + origin circle
# ==============================================================================

#' Plot RNA velocity with coloured per-cell arrows and origin circles
#'
#' Overlays per-cell velocity arrows on a UMAP or tSNE embedding.  Each cell's
#' arrow and its semi-transparent origin circle share the cell's
#' cluster/feature colour, replicating the \code{scv.pl.velocity_embedding}
#' style.  For smooth streamlines see \code{\link{ASPIS_plot_velocity_stream}}.
#'
#' @param sce A \code{SingleCellExperiment} with velocity results from
#'   \code{\link{TRIPODES_run_velocity}}.
#' @param dimred Character. Embedding to use. Default \code{"UMAP"}.
#' @param colour_by Character. \code{colData} column or gene name for cell
#'   colour. Default \code{"cluster"}.
#' @param arrow_scale Numeric. Multiplier on arrow length. Default \code{1}.
#' @param arrow_density Numeric in (0,1]. Fraction of cells that get arrows.
#'   Reduce for large datasets. Default \code{0.3}.
#' @param arrow_alpha Numeric. Arrow opacity. Default \code{0.9}.
#' @param arrow_size Numeric. Arrow line width. Default \code{0.4}.
#' @param arrow_length Numeric. Arrowhead length in cm. Default \code{0.08}.
#' @param circle_size Numeric. Size of the semi-transparent origin circle.
#'   Default \code{3}.
#' @param circle_alpha Numeric. Opacity of the origin circle. Default \code{0.2}.
#' @param point_size Numeric. Cell point size. Default \code{0.5}.
#' @param point_alpha Numeric. Cell point opacity. Default \code{0.8}.
#' @param palette Character vector or \code{NULL} for defaults.
#' @param assay_name Character. Assay for gene expression. Default
#'   \code{"logcounts"}.
#' @param title Character or \code{NULL}. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_velocity <- function(sce,
                                 dimred        = "UMAP",
                                 colour_by     = "cluster",
                                 arrow_scale   = 1,
                                 arrow_density = 0.3,
                                 arrow_alpha   = 0.9,
                                 arrow_size    = 0.4,
                                 arrow_length  = 0.08,
                                 circle_size   = 3,
                                 circle_alpha  = 0.2,
                                 point_size    = 0.5,
                                 point_alpha   = 0.8,
                                 palette       = NULL,
                                 assay_name    = "logcounts",
                                 title         = NULL) {

  df <- .aspis_vel_prep(sce, dimred, colour_by, assay_name)
  cs <- .aspis_build_cs(df, palette, colour_by)

  df$xend <- df$dim1 + arrow_scale * df$vel1
  df$yend <- df$dim2 + arrow_scale * df$vel2

  arrow_density <- min(max(arrow_density, 0), 1)
  arrow_idx <- if (arrow_density < 1)
    sample(nrow(df), max(1L, round(nrow(df) * arrow_density)))
  else
    seq_len(nrow(df))
  df_arrows <- df[arrow_idx, ]

  auto_title <- paste0("RNA Velocity \u2014 ", dimred,
                       " (", metadata(sce)$velocity_mode, ")")

  ggplot(df, aes(x = dim1, y = dim2)) +
    geom_point(aes(colour = colour),
               size = point_size, alpha = point_alpha) +
    geom_point(data  = df_arrows,
               aes(colour = colour),
               size  = circle_size,
               alpha = circle_alpha,
               shape = 16) +
    geom_segment(data      = df_arrows,
                 aes(x = dim1, y = dim2, xend = xend, yend = yend,
                     colour = colour),
                 alpha     = arrow_alpha,
                 linewidth = arrow_size,
                 arrow     = arrow(length = unit(arrow_length, "cm"),
                                   type   = "closed")) +
    cs$scale +
    labs(title  = if (is.null(title)) auto_title else title,
         x      = paste0(dimred, " 1"),
         y      = paste0(dimred, " 2"),
         colour = colour_by) +
    .aspis_vel_theme()
}


# ==============================================================================
# ASPIS_plot_velocity_stream — KDE blob background + streamlines
# ==============================================================================

#' Plot RNA velocity as smooth streamlines over a KDE cluster background
#'
#' Replicates the \code{scv.pl.velocity_embedding_stream} style.  Per-cluster
#' KDE density blobs provide the colour context; streamlines show the inferred
#' directional flow of RNA velocity across the embedding.
#'
#' Streamlines are computed by interpolating per-cell velocity vectors onto a
#' regular grid (Gaussian kernel weighting) and then integrating forward along
#' unit-normalised direction vectors using Euler steps.
#'
#' @param sce A \code{SingleCellExperiment} with velocity results from
#'   \code{\link{TRIPODES_run_velocity}}.
#' @param dimred Character. Embedding to use. Default \code{"UMAP"}.
#' @param colour_by Character. \code{colData} column used for blob and
#'   (optionally) streamline colour.  Must be discrete. Default \code{"cluster"}.
#' @param grid_res Integer. Resolution of the velocity grid (grid_res × grid_res
#'   cells). Higher values give smoother fields but are slower to compute.
#'   Default \code{30}.
#' @param n_streams Integer. Number of streamline seed points. Default
#'   \code{200}.
#' @param stream_length Integer. Maximum integration steps per streamline.
#'   Default \code{50}.
#' @param stream_colour Character colour string, or \code{NULL}.  A fixed string
#'   (e.g. \code{"black"}) renders all streamlines in that colour.  \code{NULL}
#'   colours each streamline by the nearest cell's \code{colour_by} value,
#'   matching the blob palette. Default \code{"black"}.
#' @param stream_size Numeric. Streamline line width. Default \code{0.4}.
#' @param stream_alpha Numeric. Streamline opacity. Default \code{0.8}.
#' @param bandwidth Numeric or \code{NULL}. Gaussian kernel bandwidth for
#'   velocity field interpolation.  \code{NULL} auto-scales to grid spacing.
#' @param min_density Numeric in [0, 1]. Grid points whose local cell density
#'   falls below this fraction of the grid's peak density are masked out of
#'   the velocity field entirely (no seeding, streamlines cannot cross them).
#'   Without this mask, the Gaussian kernel's infinite support means distant
#'   grid points (e.g. in empty corners of a non-convex embedding) still get
#'   a spuriously "valid" direction dominated by whichever cells are
#'   least-far-away, producing streamlines that wander through empty space
#'   with no relation to the actual data. Raise this value if streamlines
#'   still appear disconnected from the point cloud; lower it if streamlines
#'   are being cut off too early near the edges of dense clusters.
#'   Default \code{0.05}.
#' @param blob_alpha Numeric. Maximum opacity of KDE polygon fills. Default
#'   \code{0.25}.
#' @param blob_bins Integer. Number of density contour levels per cluster.
#'   Fewer bins produce softer blobs. Default \code{4}.
#' @param show_points Logical. Whether to overlay individual cell points.
#'   Default \code{FALSE}.
#' @param point_size Numeric. Cell point size (if \code{show_points=TRUE}).
#'   Default \code{0.3}.
#' @param point_alpha Numeric. Cell point opacity. Default \code{0.3}.
#' @param min_vel_frac Numeric. Fraction of mean velocity magnitude below which
#'   a streamline step is considered stalled and tracing stops. Default
#'   \code{0.05}.
#' @param palette Character vector or \code{NULL} for defaults.
#' @param assay_name Character. Assay for gene expression. Default
#'   \code{"logcounts"}.
#' @param title Character or \code{NULL}. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_velocity_stream <- function(sce,
                                        dimred        = "UMAP",
                                        colour_by     = "cluster",
                                        grid_res      = 30,
                                        n_streams     = 200,
                                        stream_length = 50,
                                        stream_colour = "black",
                                        stream_size   = 0.4,
                                        stream_alpha  = 0.8,
                                        bandwidth     = NULL,
                                        min_density   = 0.05,
                                        blob_alpha    = 0.25,
                                        blob_bins     = 4,
                                        show_points   = FALSE,
                                        point_size    = 0.3,
                                        point_alpha   = 0.3,
                                        min_vel_frac  = 0.05,
                                        palette       = NULL,
                                        assay_name    = "logcounts",
                                        title         = NULL) {

  df <- .aspis_vel_prep(sce, dimred, colour_by, assay_name)
  cs <- .aspis_build_cs(df, palette, colour_by)

  if (!cs$is_discrete)
    stop("[ASPIS] ASPIS_plot_velocity_stream requires a discrete colour_by. ",
         "For continuous variables use ASPIS_plot_velocity_v2.", call. = FALSE)

  # ── Velocity field & streamlines ──────────────────────────────────────────
  field <- .aspis_velocity_field(df$dim1, df$dim2, df$vel1, df$vel2,
                                   grid_res    = grid_res,
                                   bandwidth   = bandwidth,
                                   min_density = min_density)

  # dt: half a grid-cell width in data units (unit-normalised steps)
  dt <- (diff(field$xlim) + diff(field$ylim)) / 2 / grid_res * 0.5

  seeds <- .aspis_stream_seeds(field, n_streams)

  # Assign per-seed colours when stream_colour = NULL
  seed_colours <- NULL
  if (is.null(stream_colour) && nrow(seeds) > 0L) {
    seed_colours <- vapply(seq_len(nrow(seeds)), function(s) {
      d2   <- (df$dim1 - seeds[s, 1L])^2 + (df$dim2 - seeds[s, 2L])^2
      near <- which.min(d2)
      as.character(df$colour[near])
    }, character(1L))
  }

  stream_df <- .aspis_trace_all_streams(seeds, field, stream_length, dt,
                                         min_vel_frac, seed_colours)

  # ── Build plot ─────────────────────────────────────────────────────────────
  auto_title <- paste0("RNA Velocity (stream) \u2014 ", dimred,
                       " (", metadata(sce)$velocity_mode, ")")

  p <- ggplot(df, aes(x = dim1, y = dim2))

  # KDE blob background — one stat_density_2d layer per cluster
  for (cl in names(cs$pal)) {
    df_cl <- df[as.character(df$colour) == cl, ]
    if (nrow(df_cl) < 10L) next
    p <- p + stat_density_2d(
      data            = df_cl,
      aes(x = dim1, y = dim2),
      geom            = "polygon",
      fill            = cs$pal[[cl]],
      colour          = NA,
      alpha           = blob_alpha,
      contour_var     = "ndensity",
      bins            = blob_bins
    )
  }

  # Optional cell points
  if (show_points)
    p <- p + geom_point(aes(colour = colour),
                        size = point_size, alpha = point_alpha)

  # Invisible point layer to force colour scale + legend
  p <- p + geom_point(aes(colour = colour), alpha = 0, size = 0) +
    cs$scale

  # Streamlines
  if (!is.null(stream_df) && nrow(stream_df) > 0L) {
    if (is.null(stream_colour)) {
      # Coloured by nearest-cell cluster
      p <- p + geom_path(
        data      = stream_df,
        aes(x = x, y = y, group = stream_id,
            colour = stream_colour),
        linewidth = stream_size,
        alpha     = stream_alpha,
        arrow     = arrow(length = unit(0.1, "cm"), type = "closed",
                          ends = "last")
      )
    } else {
      # Fixed colour
      p <- p + geom_path(
        data      = stream_df,
        aes(x = x, y = y, group = stream_id),
        colour    = stream_colour,
        linewidth = stream_size,
        alpha     = stream_alpha,
        arrow     = arrow(length = unit(0.1, "cm"), type = "closed",
                          ends = "last")
      )
    }
  }

  p +
    labs(title  = if (is.null(title)) auto_title else title,
         x      = paste0(dimred, " 1"),
         y      = paste0(dimred, " 2"),
         colour = colour_by) +
    .aspis_vel_theme()
}


# ==============================================================================
# ASPIS_plot_pseudotime — pseudotime gradient on embedding
# ==============================================================================

#' Plot pseudotime as a colour gradient on a 2D embedding
#'
#' Colours cells by their Monocle3 pseudotime value on a UMAP or tSNE
#' embedding.  Cells with \code{NA} pseudotime (those in partitions without a
#' root) are drawn in \code{na_colour}.
#'
#' @param sce A \code{SingleCellExperiment} with pseudotime results from
#'   \code{\link{TRIPODES_run_monocle}}.
#' @param dimred Character. Embedding to plot on. \code{NULL} (default) uses
#'   the embedding recorded in \code{metadata(sce)$monocle_dimred}; otherwise
#'   specify any name present in \code{reducedDimNames(sce)}.
#' @param pseudotime_col Character. \code{colData} column containing pseudotime
#'   values. Default \code{"monocle_pseudotime"}.
#' @param na_colour Character. Colour for cells with \code{NA} pseudotime.
#'   Default \code{"grey80"}.
#' @param palette Character vector of colours for the gradient, or \code{NULL}
#'   for viridis (default).
#' @param point_size Numeric. Cell point size. Default \code{0.5}.
#' @param point_alpha Numeric. Cell point opacity. Default \code{0.7}.
#' @param title Character or \code{NULL}. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_pseudotime <- function(sce,
                                   dimred         = NULL,
                                   pseudotime_col = "monocle_pseudotime",
                                   na_colour      = "grey80",
                                   palette        = NULL,
                                   point_size     = 0.5,
                                   point_alpha    = 0.7,
                                   title          = NULL) {

  # ── Validate ─────────────────────────────────────────────────────────────────
  if (!isTRUE(metadata(sce)$monocle_run))
    stop("[ASPIS] TRIPODES_run_monocle() has not been run on this SCE.",
         call. = FALSE)

  if (!pseudotime_col %in% names(colData(sce)))
    stop("[ASPIS] '", pseudotime_col, "' not found in colData(sce).",
         call. = FALSE)

  dimred <- dimred %||% metadata(sce)$monocle_dimred %||% "UMAP"

  if (!dimred %in% reducedDimNames(sce))
    stop("[ASPIS] Embedding '", dimred, "' not found. Available: ",
         paste(reducedDimNames(sce), collapse = ", "), ".", call. = FALSE)

  # ── Build data frame ──────────────────────────────────────────────────────────
  coords <- as.data.frame(reducedDim(sce, dimred))
  colnames(coords) <- c("dim1", "dim2")
  coords$pseudotime <- colData(sce)[[pseudotime_col]]

  # Draw NA cells first (bottom layer), then ordered cells on top
  df_na  <- coords[ is.na(coords$pseudotime), ]
  df_ord <- coords[!is.na(coords$pseudotime), ]

  colour_scale <- if (!is.null(palette)) {
    scale_colour_gradientn(colours = palette, name = pseudotime_col,
                           na.value = na_colour)
  } else {
    scale_colour_viridis_c(name = pseudotime_col, na.value = na_colour,
                           option = "C")
  }

  auto_title <- paste0("Pseudotime \u2014 ", dimred)

  p <- ggplot(coords, aes(x = dim1, y = dim2))

  if (nrow(df_na) > 0)
    p <- p + geom_point(data  = df_na,
                        aes(x = dim1, y = dim2),
                        colour = na_colour,
                        size   = point_size,
                        alpha  = point_alpha)

  p +
    geom_point(data  = df_ord,
               aes(x = dim1, y = dim2, colour = pseudotime),
               size  = point_size,
               alpha = point_alpha) +
    colour_scale +
    labs(title  = if (is.null(title)) auto_title else title,
         x      = paste0(dimred, " 1"),
         y      = paste0(dimred, " 2"),
         colour = pseudotime_col) +
    .aspis_vel_theme()
}


# ==============================================================================
# ASPIS_plot_trajectory — principal graph overlay on embedding
# ==============================================================================

#' Plot Monocle3 principal graph trajectory on a 2D embedding
#'
#' Overlays the Monocle3 principal graph edges on a cell embedding, with cells
#' optionally coloured by pseudotime, cluster, or any other feature.  The
#' principal graph node coordinates stored by \code{\link{TRIPODES_run_monocle}}
#' are used directly, so the graph is always aligned to the correct embedding.
#'
#' @param sce A \code{SingleCellExperiment} with trajectory results from
#'   \code{\link{TRIPODES_run_monocle}}.
#' @param dimred Character. Embedding to plot on. \code{NULL} (default) uses
#'   the embedding recorded in \code{metadata(sce)$monocle_dimred}.  If a
#'   different embedding is specified, a warning is issued because the principal
#'   graph node coordinates were computed for the original embedding.
#' @param colour_by Character. \code{colData} column or gene name for cell
#'   colour.  Default \code{"monocle_pseudotime"}.
#' @param na_colour Character. Colour for cells with \code{NA} values (e.g.
#'   unordered cells when \code{colour_by = "monocle_pseudotime"}).
#'   Default \code{"grey80"}.
#' @param edge_colour Character. Colour of principal graph edges.
#'   Default \code{"black"}.
#' @param edge_size Numeric. Line width of graph edges. Default \code{0.8}.
#' @param edge_alpha Numeric. Opacity of graph edges. Default \code{0.8}.
#' @param point_size Numeric. Cell point size. Default \code{0.5}.
#' @param point_alpha Numeric. Cell point opacity. Default \code{0.7}.
#' @param palette Character vector or \code{NULL} for defaults.
#' @param assay_name Character. Assay for gene expression when
#'   \code{colour_by} is a gene. Default \code{"logcounts"}.
#' @param title Character or \code{NULL}. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_trajectory <- function(sce,
                                   dimred       = NULL,
                                   colour_by    = "monocle_pseudotime",
                                   na_colour    = "grey80",
                                   edge_colour  = "black",
                                   edge_size    = 0.8,
                                   edge_alpha   = 0.8,
                                   point_size   = 0.5,
                                   point_alpha  = 0.7,
                                   palette      = NULL,
                                   assay_name   = "logcounts",
                                   title        = NULL) {

  # ── Validate ──────────────────────────────────────────────────────────────────
  if (!isTRUE(metadata(sce)$monocle_run))
    stop("[ASPIS] TRIPODES_run_monocle() has not been run on this SCE.",
         call. = FALSE)

  if (is.null(metadata(sce)$monocle_graph) ||
      is.null(metadata(sce)$monocle_graph_nodes))
    stop("[ASPIS] Principal graph not found in metadata. ",
         "Re-run TRIPODES_run_monocle().", call. = FALSE)

  monocle_dimred <- metadata(sce)$monocle_dimred %||% "UMAP"
  dimred <- dimred %||% monocle_dimred

  if (!dimred %in% reducedDimNames(sce))
    stop("[ASPIS] Embedding '", dimred, "' not found. Available: ",
         paste(reducedDimNames(sce), collapse = ", "), ".", call. = FALSE)

  if (!identical(dimred, monocle_dimred))
    warning("[ASPIS] Plotting on '", dimred, "' but principal graph node ",
            "coordinates were computed for '", monocle_dimred, "'. ",
            "The graph overlay will not align correctly.", call. = FALSE)

  # ── Resolve colour_by ────────────────────────────────────────────────────────
  is_gene <- colour_by %in% rownames(sce)
  is_meta <- colour_by %in% names(colData(sce))
  if (!is_gene && !is_meta)
    stop("[ASPIS] '", colour_by,
         "' not found in colData(sce) or rownames(sce).", call. = FALSE)

  colour_vals <- if (is_gene) {
    as.numeric(assay(sce, assay_name)[colour_by, ])
  } else {
    colData(sce)[[colour_by]]
  }

  # ── Build cell data frame ────────────────────────────────────────────────────
  coords <- as.data.frame(reducedDim(sce, dimred))
  colnames(coords) <- c("dim1", "dim2")
  coords$colour <- colour_vals

  # ── Build principal graph edge data frame ────────────────────────────────────
  pg_nodes <- metadata(sce)$monocle_graph_nodes   # nodes × 2
  pg_graph <- metadata(sce)$monocle_graph

  edge_list <- igraph::as_edgelist(pg_graph, names = FALSE)
  edge_df   <- data.frame(
    x    = pg_nodes[edge_list[, 1L], 1L],
    y    = pg_nodes[edge_list[, 1L], 2L],
    xend = pg_nodes[edge_list[, 2L], 1L],
    yend = pg_nodes[edge_list[, 2L], 2L]
  )

  # ── Colour scale ──────────────────────────────────────────────────────────────
  is_discrete <- is.factor(coords$colour) || is.character(coords$colour)

  colour_scale <- if (is_discrete) {
    pal <- if (!is.null(palette)) palette else .aspis_discrete_palette()
    lvls <- if (is.factor(coords$colour)) levels(coords$colour)
            else sort(unique(coords$colour))
    pal <- stats::setNames(rep(pal, length.out = length(lvls)), lvls)
    scale_colour_manual(values = pal, name = colour_by, na.value = na_colour)
  } else {
    if (!is.null(palette))
      scale_colour_gradientn(colours = palette, name = colour_by,
                             na.value = na_colour)
    else
      scale_colour_viridis_c(name = colour_by, na.value = na_colour,
                             option = "C")
  }

  auto_title <- paste0("Trajectory \u2014 ", dimred)

  # ── Build plot ───────────────────────────────────────────────────────────────
  ggplot(coords, aes(x = dim1, y = dim2)) +
    geom_point(aes(colour = colour),
               size  = point_size,
               alpha = point_alpha) +
    geom_segment(data      = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour    = edge_colour,
                 linewidth = edge_size,
                 alpha     = edge_alpha) +
    colour_scale +
    labs(title  = if (is.null(title)) auto_title else title,
         x      = paste0(dimred, " 1"),
         y      = paste0(dimred, " 2"),
         colour = colour_by) +
    .aspis_vel_theme()
}
