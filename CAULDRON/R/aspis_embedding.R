# ==============================================================================
# ASPIS - Embedding & QC Visualisation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The Shield of Achilles, forged by Hephaestus and described in extraordinary
# detail in the Iliad (Book 18) — it depicted all of human civilisation on
# its surface: cities, fields, weddings, battles. The shield that shows
# everything. ASPIS makes the invisible structure of data visible.
#
# This file covers dimensionality-reduction and QC plots:
#   - Elbow plot (PCA variance)
#   - UMAP and tSNE scatter plots
#   - Embedding parameter-sweep grids
#   - QC metric violin plots
#   - Shared internal helpers (.aspis_plot_dimred, .aspis_discrete_palette, …)
#
# All public functions prefixed: ASPIS_
# ==============================================================================


#' Elbow plot of PCA variance explained
#'
#' Plots the percentage of variance explained by each principal component,
#' with an optional cumulative variance overlay, a user-specified cut-off line,
#' and algorithmic suggestions to help choose the number of components.
#'
#' Two suggestion methods are available and can be overlaid simultaneously:
#'
#' \describe{
#'   \item{\code{"elbow"}}{Finds the PC where the second derivative of the
#'     variance curve is maximised — i.e., where the rate of decline itself
#'     decelerates most sharply.  Tends to suggest fewer PCs for data with a
#'     quick early drop-off.}
#'   \item{\code{"cumulative"}}{Finds the minimum number of PCs needed to
#'     explain at least \code{cum_threshold}\% of total variance.  Scales
#'     naturally with how spread-out the variance is across components.}
#' }
#'
#' Suggestions are printed to the console and drawn as distinct vertical lines
#' (green = elbow, purple = cumulative) so they can be compared against one
#' another and against the user's chosen \code{n_pcs}.
#'
#' @param sce A \code{SingleCellExperiment} with \code{reducedDims(sce)[["PCA"]]}
#'   populated by \code{\link{TALOS_run_pca}}.
#' @param n_show Integer. Number of components to display. Default \code{50}.
#' @param n_pcs Integer. If provided, draws a grey dashed line marking the
#'   user's chosen cut-off. Default \code{NULL}.
#' @param show_cumulative Logical. Overlay cumulative variance as a red dashed
#'   line. Default \code{TRUE}.
#' @param suggest Character or \code{NULL}.  Algorithmic cut-off suggestions to
#'   display.  One of \code{"elbow"}, \code{"cumulative"}, \code{"all"}, or
#'   \code{NULL} (default, no suggestion).
#' @param cum_threshold Numeric. Cumulative variance target (%) for the
#'   \code{"cumulative"} method. Default \code{80}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_elbow <- function(sce,
                              n_show          = 50L,
                              n_pcs           = NULL,
                              show_cumulative = TRUE,
                              suggest         = NULL,
                              cum_threshold   = 80) {

  if (!"PCA" %in% reducedDimNames(sce))
    stop("PCA not found. Run TALOS_run_pca() first.", call. = FALSE)

  pct_var <- attr(reducedDim(sce, "PCA"), "percentVar")

  if (is.null(pct_var))
    stop("No percentVar attribute found on the PCA result. ",
         "Re-run TALOS_run_pca() to ensure it is stored.", call. = FALSE)

  if (!is.null(suggest))
    suggest <- match.arg(suggest, c("elbow", "cumulative", "all"))

  n_show <- min(as.integer(n_show), length(pct_var))

  df <- data.frame(
    pc      = seq_len(n_show),
    var_exp = pct_var[seq_len(n_show)],
    cum_var = cumsum(pct_var)[seq_len(n_show)]
  )

  # ── Base plot ────────────────────────────────────────────────────────────────
  p <- ggplot(df, aes(x = pc)) +
    geom_col(aes(y = var_exp), fill = "#4E79A7", alpha = 0.7, width = 0.7) +
    geom_point(aes(y = var_exp), colour = "#4E79A7", size = 1.5) +
    labs(x = "Principal component", y = "Variance explained (%)",
         title = "PCA elbow plot") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (isTRUE(show_cumulative))
    p <- p +
      geom_line(aes(y = cum_var), colour = "#E15759", linewidth = 0.7,
                linetype = "dashed") +
      geom_point(aes(y = cum_var), colour = "#E15759", size = 1.5, shape = 17)

  # ── User-chosen cut-off (grey) ───────────────────────────────────────────────
  if (!is.null(n_pcs)) {
    n_pcs <- min(as.integer(n_pcs), n_show)
    p <- p +
      geom_vline(xintercept = n_pcs, linetype = "dashed",
                 colour = "grey30", linewidth = 0.7) +
      annotate("text", x = n_pcs + 0.4, y = max(df$var_exp) * 0.97,
               label = paste0("chosen\nPC", n_pcs),
               hjust = 0, size = 3, colour = "grey30")
  }

  # ── Algorithmic suggestions ──────────────────────────────────────────────────
  if (!is.null(suggest)) {
    sugg <- .aspis_suggest_pcs(pct_var, suggest, cum_threshold)

    if (!is.null(sugg$elbow)) {
      e  <- min(sugg$elbow, n_show)
      ce <- df$cum_var[e]
      p  <- p +
        geom_vline(xintercept = e, linetype = "dotdash",
                   colour = "#59A14F", linewidth = 0.8) +
        annotate("text", x = e + 0.4, y = max(df$var_exp) * 0.80,
                 label = paste0("elbow\nPC", e),
                 hjust = 0, size = 3, colour = "#59A14F")
      cat(sprintf("Elbow suggestion      : PC%d  (%.1f%% cum. variance)\n", e, ce))
    }

    if (!is.null(sugg$cumulative)) {
      cv  <- min(sugg$cumulative, n_show)
      p   <- p +
        geom_vline(xintercept = cv, linetype = "dotdash",
                   colour = "#B07AA1", linewidth = 0.8) +
        annotate("text", x = cv + 0.4, y = max(df$var_exp) * 0.63,
                 label = paste0(">=", cum_threshold, "%\nPC", cv),
                 hjust = 0, size = 3, colour = "#B07AA1")
      cat(sprintf("Cumulative suggestion : PC%d  (>= %g%% variance explained)\n",
                      cv, cum_threshold))
    }
  }

  # ── Caption ──────────────────────────────────────────────────────────────────
  caption_parts <- "Bars/blue: individual variance"
  if (isTRUE(show_cumulative))
    caption_parts <- paste0(caption_parts, "  |  Red dashed: cumulative variance")
  if (!is.null(suggest) && suggest %in% c("elbow", "all"))
    caption_parts <- paste0(caption_parts, "  |  Green dot-dash: elbow suggestion")
  if (!is.null(suggest) && suggest %in% c("cumulative", "all"))
    caption_parts <- paste0(caption_parts,
                            "  |  Purple dot-dash: cumulative suggestion")

  p + labs(caption = caption_parts)
}


#' UMAP embedding plot
#'
#' Plots the UMAP embedding stored in \code{reducedDims(sce)[["UMAP"]]},
#' coloured by a \code{colData} column or gene expression value.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"UMAP"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_umap}} first).
#' @param colour_by Character scalar or vector.  One or more \code{colData}
#'   column names or gene names present in \code{rownames(sce)}.  A single
#'   value returns a \code{ggplot}; a vector produces one panel per element
#'   arranged in a grid (requires \code{gridExtra}).  Default \code{"cluster"}.
#' @param point_size Numeric. Point size. Default \code{0.8}.
#' @param point_alpha Numeric. Point transparency (0–1). Default \code{0.6}.
#' @param palette Character vector or \code{NULL}.  For discrete variables:
#'   a vector of colours (recycled as needed).  For continuous variables or
#'   gene expression: a 2-element vector \code{c(low, high)}.
#'   \code{NULL} (default) uses built-in palettes.
#' @param title Character or \code{NULL}.  Plot title.  \code{NULL} auto-generates
#'   \code{"UMAP — <colour_by>"}.  Ignored when \code{colour_by} is a vector.
#' @param label_clusters Logical.  Overlay cluster centroid labels.  Only
#'   applied when \code{colour_by} resolves to a discrete variable.
#'   Default \code{FALSE}.
#' @param label_size Numeric.  Size of centroid labels.  Default \code{4}.
#' @param assay_name Character.  Assay used when \code{colour_by} is a gene.
#'   Default \code{"logcounts"}.
#' @param ncol Integer or \code{NULL}.  Number of columns in the panel grid
#'   when \code{colour_by} is a vector.  \code{NULL} (default) uses
#'   \code{min(length(colour_by), 3)}.
#' @param style Character.  Visual style: \code{"points"} (default, classic
#'   scatter), \code{"contour"} (density ring lines only, no points), or
#'   \code{"both"} (points with density ring lines overlaid).  Topographic
#'   styles require a discrete \code{colour_by}.
#' @param n_levels Integer or \code{"dynamic"}.  Number of contour levels when
#'   using a fixed spacing.  Pass \code{"dynamic"} to compute per-cluster
#'   normalised KDE contours — each cluster is contoured relative to its own
#'   density range, avoiding smearing on dense clusters while still showing
#'   structure in dispersed ones.  Default \code{8}.
#' @param contour_alpha Numeric.  Opacity of the contour lines.  Default \code{0.7}.
#' @param show_legend Logical.  Whether to display the colour legend.
#'   Default \code{TRUE}.
#' @param boundary_pad Numeric.  Fractional expansion of the plot boundaries
#'   beyond the data range (e.g. \code{0.10} = 10\%).  Increase if contour
#'   lines are clipped at the edges; decrease if there is too much whitespace.
#'   Default \code{0.10}.
#'
#' @return A \code{ggplot} (single \code{colour_by}) or a \code{gtable} grid
#'   (multiple \code{colour_by}).
#' @export
ASPIS_plot_umap <- function(sce,
                             colour_by      = "cluster",
                             point_size     = 0.8,
                             point_alpha    = 0.6,
                             palette        = NULL,
                             title          = NULL,
                             label_clusters = FALSE,
                             label_size     = 4,
                             assay_name     = "logcounts",
                             ncol           = NULL,
                             style          = c("points", "contour", "both"),
                             n_levels       = 8L,
                             contour_alpha  = 0.7,
                             show_legend    = TRUE,
                             boundary_pad   = 0.10) {

  style <- match.arg(style)
  if (!"UMAP" %in% reducedDimNames(sce))
    stop("UMAP not found. Run TALOS_run_umap() first.", call. = FALSE)

  .aspis_plot_dimred(sce, "UMAP",
                     colour_by      = colour_by,
                     point_size     = point_size,
                     point_alpha    = point_alpha,
                     palette        = palette,
                     title          = title,
                     label_clusters = label_clusters,
                     label_size     = label_size,
                     assay_name     = assay_name,
                     ncol           = ncol,
                     style          = style,
                     n_levels       = n_levels,
                     contour_alpha  = contour_alpha,
                     show_legend    = show_legend,
                     boundary_pad   = boundary_pad)
}


#' tSNE embedding plot
#'
#' Plots the tSNE embedding stored in \code{reducedDims(sce)[["tSNE"]]},
#' coloured by a \code{colData} column or gene expression value.
#'
#' @param sce A \code{SingleCellExperiment} with \code{"tSNE"} in
#'   \code{reducedDims} (run \code{\link{TALOS_run_tsne}} first).
#' @param colour_by Character scalar or vector.  One or more \code{colData}
#'   column names or gene names present in \code{rownames(sce)}.  A single
#'   value returns a \code{ggplot}; a vector produces one panel per element
#'   arranged in a grid (requires \code{gridExtra}).  Default \code{"cluster"}.
#' @param point_size Numeric. Point size. Default \code{0.8}.
#' @param point_alpha Numeric. Point transparency (0–1). Default \code{0.6}.
#' @param palette Character vector or \code{NULL}.  For discrete variables:
#'   a vector of colours (recycled as needed).  For continuous variables or
#'   gene expression: a 2-element vector \code{c(low, high)}.
#'   \code{NULL} (default) uses built-in palettes.
#' @param title Character or \code{NULL}.  Plot title.  \code{NULL} auto-generates
#'   \code{"tSNE — <colour_by>"}.  Ignored when \code{colour_by} is a vector.
#' @param label_clusters Logical.  Overlay cluster centroid labels.  Only
#'   applied when \code{colour_by} resolves to a discrete variable.
#'   Default \code{FALSE}.
#' @param label_size Numeric.  Size of centroid labels.  Default \code{4}.
#' @param assay_name Character.  Assay used when \code{colour_by} is a gene.
#'   Default \code{"logcounts"}.
#' @param ncol Integer or \code{NULL}.  Number of columns in the panel grid
#'   when \code{colour_by} is a vector.  \code{NULL} (default) uses
#'   \code{min(length(colour_by), 3)}.
#' @param style Character.  Visual style: \code{"points"} (default, classic
#'   scatter), \code{"contour"} (density ring lines only, no points), or
#'   \code{"both"} (points with density ring lines overlaid).  Topographic
#'   styles require a discrete \code{colour_by}.
#' @param n_levels Integer or \code{"dynamic"}.  Number of contour levels when
#'   using a fixed spacing.  Pass \code{"dynamic"} to compute per-cluster
#'   normalised KDE contours.  Default \code{8}.
#' @param contour_alpha Numeric.  Opacity of the contour lines.  Default \code{0.7}.
#' @param show_legend Logical.  Whether to display the colour legend.
#'   Default \code{TRUE}.
#' @param boundary_pad Numeric.  Fractional expansion of the plot boundaries
#'   beyond the data range (e.g. \code{0.10} = 10\%).  Increase if contour
#'   lines are clipped at the edges; decrease if there is too much whitespace.
#'   Default \code{0.10}.
#'
#' @return A \code{ggplot} (single \code{colour_by}) or a \code{gtable} grid
#'   (multiple \code{colour_by}).
#' @export
ASPIS_plot_tsne <- function(sce,
                             colour_by      = "cluster",
                             point_size     = 0.8,
                             point_alpha    = 0.6,
                             palette        = NULL,
                             title          = NULL,
                             label_clusters = FALSE,
                             label_size     = 4,
                             assay_name     = "logcounts",
                             ncol           = NULL,
                             style          = c("points", "contour", "both"),
                             n_levels       = 8L,
                             contour_alpha  = 0.7,
                             show_legend    = TRUE,
                             boundary_pad   = 0.10) {

  style <- match.arg(style)
  if (!"tSNE" %in% reducedDimNames(sce))
    stop("tSNE not found. Run TALOS_run_tsne() first.", call. = FALSE)

  .aspis_plot_dimred(sce, "tSNE",
                     colour_by      = colour_by,
                     point_size     = point_size,
                     point_alpha    = point_alpha,
                     palette        = palette,
                     title          = title,
                     label_clusters = label_clusters,
                     label_size     = label_size,
                     assay_name     = assay_name,
                     ncol           = ncol,
                     style          = style,
                     n_levels       = n_levels,
                     contour_alpha  = contour_alpha,
                     show_legend    = show_legend,
                     boundary_pad   = boundary_pad)
}


#' Grid of embedding plots from a parameter sweep
#'
#' Takes the output of \code{\link{TALOS_tune_umap}} or
#' \code{\link{TALOS_tune_tsne}} and plots every parameter combination as a
#' small embedding, arranged in a grid.  The best combination (by composite
#' score) is marked with a \code{★} in its panel title.
#'
#' Embeddings are read directly from the sweep object (pre-computed during the
#' tuning run) so no re-computation is needed.
#'
#' @param sweep A \code{talos_embedding_sweep} object produced by
#'   \code{\link{TALOS_tune_umap}} or \code{\link{TALOS_tune_tsne}}.
#' @param sce A \code{SingleCellExperiment} used only for \code{colData} and
#'   assay access when resolving \code{colour_by}.
#' @param colour_by Character.  A \code{colData} column name (e.g.
#'   \code{"cluster"}, \code{"region"}, \code{"condition"}) or a gene name.
#'   Default \code{"cluster"}.
#' @param point_size Numeric.  Point size.  Smaller values work better in a
#'   dense grid.  Default \code{0.3}.
#' @param point_alpha Numeric.  Point transparency.  Default \code{0.5}.
#' @param palette Character vector or \code{NULL}.  Passed to the colour
#'   scale — see \code{\link{ASPIS_plot_umap}} for details.  Default
#'   \code{NULL} uses built-in palettes.
#' @param assay_name Character.  Assay used when \code{colour_by} is a gene.
#'   Default \code{"logcounts"}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_embedding_grid <- function(sweep,
                                       sce,
                                       colour_by   = "cluster",
                                       point_size  = 0.3,
                                       point_alpha = 0.5,
                                       palette     = NULL,
                                       assay_name  = "logcounts") {

  if (!inherits(sweep, "talos_embedding_sweep"))
    stop("'sweep' must be a talos_embedding_sweep object from ",
         "TALOS_tune_umap() or TALOS_tune_tsne().", call. = FALSE)

  if (is.null(sweep$embeddings))
    stop("No embeddings found in sweep object. ",
         "Re-run TALOS_tune_umap() or TALOS_tune_tsne().", call. = FALSE)

  is_umap  <- sweep$params$type == "umap"
  res_df   <- sweep$results
  n_runs   <- nrow(res_df)
  dim_name <- if (is_umap) "UMAP" else "tSNE"

  # ── Resolve colour_by ────────────────────────────────────────────────────────
  is_gene <- colour_by %in% rownames(sce)
  is_meta <- colour_by %in% names(colData(sce))

  if (!is_gene && !is_meta)
    stop("'", colour_by, "' not found in colData(sce) or rownames(sce).",
         call. = FALSE)

  if (is_gene) {
    if (!assay_name %in% assayNames(sce))
      stop("Assay '", assay_name, "' not found.", call. = FALSE)
    colour_vals <- as.numeric(assay(sce, assay_name)[colour_by, ])
    is_discrete <- FALSE
  } else {
    colour_vals <- colData(sce)[[colour_by]]
    is_discrete <- is.factor(colour_vals) || is.character(colour_vals)
    if (is_discrete) colour_vals <- as.factor(colour_vals)
  }

  # ── Build panel label order (natural sort, best marked with ★) ───────────────
  if (is_umap) {
    ord  <- order(res_df$n_neighbors, res_df$min_dist)
    labs <- vapply(ord, function(i) {
      nn <- res_df$n_neighbors[i]; md <- res_df$min_dist[i]
      is_best <- nn == sweep$best_params$n_neighbors &&
                 md == sweep$best_params$min_dist
      sprintf("%snn=%d | md=%.2f", if (is_best) "\u2605 " else "", nn, md)
    }, character(1L))
  } else {
    ord  <- order(res_df$perplexity)
    labs <- vapply(ord, function(i) {
      pp      <- res_df$perplexity[i]
      is_best <- pp == sweep$best_params$perplexity
      sprintf("%sperplexity=%.0f", if (is_best) "\u2605 " else "", pp)
    }, character(1L))
  }

  # ── Assemble data from stored embeddings ─────────────────────────────────────
  df_list <- vector("list", n_runs)

  for (idx in seq_along(ord)) {
    i   <- ord[idx]
    emb <- sweep$embeddings[[i]]
    df_list[[idx]] <- data.frame(
      dim1        = emb[, 1L],
      dim2        = emb[, 2L],
      colour_val  = colour_vals,
      panel_label = factor(labs[idx], levels = labs)
    )
  }

  df_all <- do.call(rbind, df_list)

  # ── ncol: for UMAP lay out as n_neighbors × min_dist grid ───────────────────
  ncols <- if (is_umap) length(unique(res_df$min_dist)) else
              min(3L, n_runs)

  # ── Plot ─────────────────────────────────────────────────────────────────────
  p <- ggplot(df_all, aes(x = dim1, y = dim2, colour = colour_val)) +
    geom_point(size = point_size, alpha = point_alpha) +
    facet_wrap(~ panel_label, ncol = ncols) +
    labs(x      = paste(dim_name, "1"),
         y      = paste(dim_name, "2"),
         title  = sprintf("%s parameter grid  |  coloured by: %s  |  \u2605 = best",
                          dim_name, colour_by),
         colour = colour_by) +
    theme_bw(base_size = 9) +
    theme(panel.grid  = element_blank(),
          axis.ticks  = element_blank(),
          axis.text   = element_blank(),
          strip.text  = element_text(size = 7.5))

  # ── Colour scale ─────────────────────────────────────────────────────────────
  if (is_discrete) {
    n_lev <- nlevels(df_all$colour_val)
    pal   <- if (!is.null(palette)) rep(palette, length.out = n_lev) else
               rep(.aspis_discrete_palette(), length.out = n_lev)
    p <- p + scale_color_manual(values = pal)
  } else {
    low  <- if (!is.null(palette) && length(palette) >= 1L) palette[1L] else "grey90"
    high <- if (!is.null(palette) && length(palette) >= 2L) palette[2L] else
              if (is_gene) "#2166AC" else "#4E79A7"
    p <- p + scale_color_gradient(low = low, high = high, name = colour_by)
  }

  p
}


#' QC metric violin plots
#'
#' Visualises the distribution of per-cell quality control metrics computed by
#' \code{\link{AEGIS_compute_qc_metrics}}.  Each metric is shown as a violin
#' with an optional jitter overlay, arranged in a faceted grid.  Threshold
#' lines can be drawn to inspect prospective filtering cut-offs before calling
#' \code{\link{AEGIS_filter_cells}}.
#'
#' The following colData columns are produced by
#' \code{\link{AEGIS_compute_qc_metrics}} and can be passed to
#' \code{metrics}:
#'
#' \describe{
#'   \item{\code{sum}}{Total UMI count per cell (library size).  Low values
#'     indicate empty droplets or dead cells; extremely high values may indicate
#'     doublets.  Typically inspected on a log scale.}
#'   \item{\code{detected}}{Number of genes with at least one count.  Follows
#'     a similar distribution to \code{sum} but is less sensitive to a few
#'     highly expressed genes.  Typically inspected on a log scale.}
#'   \item{\code{subsets_mt_percent}}{Percentage of counts from mitochondrial
#'     genes.  High values suggest compromised cell membranes (cytoplasmic RNA
#'     lost, mitochondrial RNA retained).  Inspected on a linear scale.}
#'   \item{\code{subsets_ribo_percent}}{Percentage of counts from ribosomal
#'     protein genes.  Unusually high values may indicate stressed or
#'     proliferating cells.  Only present if ribosomal genes were detected.}
#' }
#'
#' @param sce A \code{SingleCellExperiment} with QC columns in
#'   \code{colData} (run \code{\link{AEGIS_compute_qc_metrics}} first).
#' @param metrics Character vector of \code{colData} column names to plot.
#'   Default \code{c("sum", "detected", "subsets_mt_percent")}.  Any numeric
#'   \code{colData} column is accepted.
#' @param group_by Character.  A \code{colData} column used to split violins
#'   by group (e.g. \code{"sample_id"}, \code{"condition"}).  \code{NULL}
#'   (default) shows a single violin per metric.
#' @param thresholds Named list of threshold values to draw as red dashed
#'   horizontal lines.  Names must match entries in \code{metrics}
#'   (on the original, untransformed scale).  Example:
#'   \code{list(subsets_mt_percent = 20, sum = 500)}.  \code{NULL} (default)
#'   draws no lines.
#' @param log_scale Character vector of metric names to display on a log10
#'   scale.  Values are transformed as \code{log10(x + 1)} and axis labels
#'   updated accordingly.  Default \code{c("sum", "detected")}.
#' @param show_points Logical.  Overlay individual cell points as jitter.
#'   Automatically suppressed when the number of cells exceeds
#'   \code{max_points} to avoid overplotting.  Default \code{TRUE}.
#' @param max_points Integer.  Maximum number of cells for which jitter points
#'   are drawn.  If \code{ncol(sce) > max_points}, a random subsample of
#'   \code{max_points} cells is shown.  Default \code{5000}.
#' @param point_size Numeric.  Jitter point size.  Default \code{0.3}.
#' @param point_alpha Numeric.  Jitter point transparency.  Default \code{0.3}.
#' @param ncol Integer or \code{NULL}.  Number of columns in the facet grid.
#'   \code{NULL} (default) uses \code{min(length(metrics), 3)}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_qc <- function(sce,
                           metrics     = c("sum", "detected", "subsets_mt_percent"),
                           group_by    = NULL,
                           thresholds  = NULL,
                           log_scale   = c("sum", "detected"),
                           show_points = TRUE,
                           max_points  = 5000L,
                           point_size  = 0.3,
                           point_alpha = 0.3,
                           ncol        = NULL) {

  available <- names(colData(sce))
  missing_m <- setdiff(metrics, available)
  if (length(missing_m) > 0L)
    stop("Metrics not found in colData: ", paste(missing_m, collapse = ", "),
         "\nAvailable colData columns: ", paste(available, collapse = ", "),
         call. = FALSE)

  if (!is.null(group_by) && !group_by %in% available)
    stop("group_by '", group_by, "' not found in colData.", call. = FALSE)

  log_scale <- intersect(log_scale, metrics)
  n_cells   <- ncol(sce)

  # ── Build group vector ───────────────────────────────────────────────────────
  group_vec <- if (!is.null(group_by)) {
    as.factor(colData(sce)[[group_by]])
  } else {
    factor(rep("cells", n_cells))
  }

  # ── Build long-format data frame ─────────────────────────────────────────────
  metric_labels <- setNames(vapply(metrics, function(m) {
    if (m %in% log_scale) paste0(m, " (log\u2081\u2080 + 1)") else m
  }, character(1L)), metrics)

  df_list <- lapply(metrics, function(m) {
    vals <- as.numeric(colData(sce)[[m]])
    if (m %in% log_scale) vals <- log10(vals + 1)
    data.frame(group  = group_vec,
               metric = unname(metric_labels[m]),
               value  = vals,
               stringsAsFactors = FALSE)
  })

  df <- do.call(rbind, df_list)
  df$metric <- factor(df$metric, levels = metric_labels)

  # ── Jitter subsampling ───────────────────────────────────────────────────────
  do_jitter <- isTRUE(show_points)
  if (do_jitter && n_cells > as.integer(max_points)) {
    jitter_idx <- sample.int(n_cells, as.integer(max_points))
    df_jitter  <- df[rep(jitter_idx, length(metrics)) +
                       rep(seq(0, (length(metrics) - 1L) * n_cells, by = n_cells),
                           each = length(jitter_idx)), ]
  } else {
    df_jitter <- df
  }

  has_groups <- !is.null(group_by)
  n_groups   <- nlevels(group_vec)
  ncol_use   <- if (!is.null(ncol)) as.integer(ncol) else min(length(metrics), 3L)

  # ── Colour palette ───────────────────────────────────────────────────────────
  pal <- if (has_groups) {
    rep(.aspis_discrete_palette(), length.out = n_groups)
  } else {
    "#4E79A7"
  }

  # ── Base plot ────────────────────────────────────────────────────────────────
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.75, linewidth = 0.4) +
    scale_fill_manual(values = pal) +
    facet_wrap(~ metric, scales = "free_y", ncol = ncol_use) +
    labs(x     = if (has_groups) group_by else NULL,
         y     = NULL,
         title = "QC metrics",
         fill  = if (has_groups) group_by else NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = if (has_groups && n_groups > 1L) "right" else "none")

  if (!has_groups)
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())

  # ── Jitter overlay ───────────────────────────────────────────────────────────
  if (do_jitter)
    p <- p + geom_jitter(data    = df_jitter,
                         width   = 0.15,
                         size    = point_size,
                         alpha   = point_alpha,
                         colour  = "grey20",
                         inherit.aes = TRUE)

  # ── Threshold lines ──────────────────────────────────────────────────────────
  if (!is.null(thresholds)) {
    thresh_rows <- lapply(names(thresholds), function(m) {
      if (!m %in% metrics) return(NULL)
      lab <- metric_labels[m]
      val <- if (m %in% log_scale) log10(thresholds[[m]] + 1) else thresholds[[m]]
      data.frame(metric = factor(lab, levels = levels(df$metric)), yintercept = val)
    })
    thresh_df <- do.call(rbind, Filter(Negate(is.null), thresh_rows))
    if (nrow(thresh_df) > 0L)
      p <- p + geom_hline(data     = thresh_df,
                          aes(yintercept = yintercept),
                          linetype = "dashed",
                          colour   = "#E15759",
                          linewidth = 0.7,
                          inherit.aes = FALSE)
  }

  p
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Compute algorithmic PC cut-off suggestions.
# - "elbow":      PC where the second derivative of variance is maximised
#                 (sharpest deceleration in the drop-off)
# - "cumulative": minimum PCs needed to reach cum_threshold% total variance
.aspis_suggest_pcs <- function(pct_var, method, cum_threshold = 80) {

  elbow_pc <- NULL
  cum_pc   <- NULL

  if (method %in% c("elbow", "all")) {
    n <- length(pct_var)
    if (n > 2) {
      d2       <- diff(diff(pct_var))   # length n-2; d2[i] <-> PC i+1
      elbow_pc <- as.integer(which.max(d2) + 1L)
    } else {
      elbow_pc <- 1L
    }
  }

  if (method %in% c("cumulative", "all")) {
    hits   <- which(cumsum(pct_var) >= cum_threshold)
    cum_pc <- if (length(hits)) as.integer(hits[1L]) else as.integer(length(pct_var))
  }

  list(elbow = elbow_pc, cumulative = cum_pc)
}


# Shared embedding plot engine (UMAP and tSNE).
# Also called by ASPIS_plot_atlas() in aspis_atlas.R.
.aspis_plot_dimred <- function(sce, dimred, colour_by, point_size, point_alpha,
                                palette, title, label_clusters, label_size,
                                assay_name, ncol = NULL,
                                style = "points", n_levels = 8L,
                                contour_alpha = 0.7, show_legend = TRUE,
                                boundary_pad = 0.10) {

  # ── Multi-panel: one plot per colour_by element ──────────────────────────────
  if (length(colour_by) > 1L) {
    if (!requireNamespace("gridExtra", quietly = TRUE))
      stop("Package 'gridExtra' is required for multi-panel plots. ",
           "Install via: install.packages(\"gridExtra\")", call. = FALSE)
    plots    <- lapply(colour_by, function(cb)
      .aspis_plot_dimred(sce, dimred, colour_by = cb, point_size = point_size,
                         point_alpha = point_alpha, palette = palette,
                         title = NULL, label_clusters = label_clusters,
                         label_size = label_size, assay_name = assay_name,
                         ncol = NULL, style = style, n_levels = n_levels,
                         contour_alpha = contour_alpha, show_legend = show_legend,
                         boundary_pad = boundary_pad))
    ncol_use <- if (!is.null(ncol)) as.integer(ncol) else min(length(colour_by), 3L)
    return(gridExtra::grid.arrange(grobs = plots, ncol = ncol_use))
  }

  # ── Extract embedding ────────────────────────────────────────────────────────
  emb <- reducedDim(sce, dimred)
  df  <- data.frame(dim1 = emb[, 1L], dim2 = emb[, 2L])

  # ── Resolve colour_by source ─────────────────────────────────────────────────
  is_gene <- colour_by %in% rownames(sce)
  is_meta <- colour_by %in% names(colData(sce))

  if (!is_gene && !is_meta)
    stop("'", colour_by, "' not found in colData(sce) or rownames(sce).",
         call. = FALSE)

  if (is_gene) {
    if (!assay_name %in% assayNames(sce))
      stop("Assay '", assay_name, "' not found.", call. = FALSE)
    df$colour_val <- as.numeric(assay(sce, assay_name)[colour_by, ])
    is_discrete   <- FALSE
    is_gene_val   <- TRUE
  } else {
    df$colour_val <- colData(sce)[[colour_by]]
    is_discrete   <- is.factor(df$colour_val) || is.character(df$colour_val)
    is_gene_val   <- FALSE
    if (is_discrete)
      df$colour_val <- as.factor(df$colour_val)
  }

  # ── Warn and fall back if topographic style requested for continuous data ─────
  if (style %in% c("contour", "both") && !is_discrete) {
    warning("Topographic style requires a discrete colour_by. ",
            "Falling back to style = 'points'.", call. = FALSE)
    style <- "points"
  }

  # ── Colour palette ───────────────────────────────────────────────────────────
  if (is_discrete) {
    n_lev <- nlevels(df$colour_val)
    pal   <- rep(if (!is.null(palette)) palette else .aspis_discrete_palette(),
                 length.out = n_lev)
  }

  # ── Axis / title labels ──────────────────────────────────────────────────────
  display_name <- if (dimred == "UMAP") "UMAP" else "tSNE"
  dim_labels   <- paste(display_name, 1:2)
  auto_title   <- if (!is.null(title)) title else
                    sprintf("%s \u2014 %s", display_name, colour_by)

  # ── Base canvas (no geoms yet — layer order matters) ─────────────────────────
  p <- ggplot(df, aes(x = dim1, y = dim2, colour = colour_val)) +
    labs(x      = dim_labels[1],
         y      = dim_labels[2],
         title  = auto_title,
         colour = colour_by) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank())

  # ── Axis limits — pre-expand by boundary_pad to ensure contours aren't clipped
  x_pad <- diff(range(df$dim1)) * boundary_pad
  y_pad <- diff(range(df$dim2)) * boundary_pad
  xlim  <- c(range(df$dim1)[1] - x_pad, range(df$dim1)[2] + x_pad)
  ylim  <- c(range(df$dim2)[1] - y_pad, range(df$dim2)[2] + y_pad)

  # ── Contour rings (drawn first so points render on top) ──────────────────────
  if (style %in% c("contour", "both")) {
    if (identical(n_levels, "dynamic")) {
      contour_df <- .aspis_dynamic_contours(df, n_levels = 10L, pad = boundary_pad)
      # Further extend to include actual contour coordinates if they exceed 15%
      xlim <- range(c(xlim, contour_df$dim1))
      ylim <- range(c(ylim, contour_df$dim2))
      p <- p + geom_path(
        data        = contour_df,
        aes(x       = dim1,
            y       = dim2,
            colour  = colour_val,
            group   = line_id),
        linewidth   = 0.25,
        alpha       = contour_alpha,
        inherit.aes = FALSE
      )
    } else {
      p <- p +
        stat_density_2d(
          aes(group = colour_val, colour = colour_val),
          bins      = as.integer(n_levels),
          linewidth = 0.25,
          alpha     = contour_alpha
        )
    }
  }

  # ── Points layer ─────────────────────────────────────────────────────────────
  if (style %in% c("points", "both")) {
    pt_size  <- if (style == "both") point_size * 0.6 else point_size
    pt_alpha <- if (style == "both") point_alpha * 0.7 else point_alpha
    p <- p + geom_point(size = pt_size, alpha = pt_alpha)
  }

  # ── Colour scale ─────────────────────────────────────────────────────────────
  if (is_discrete) {
    p <- p + scale_color_manual(values = pal)
  } else {
    low  <- if (!is.null(palette) && length(palette) >= 1L) palette[1L] else "grey90"
    high <- if (!is.null(palette) && length(palette) >= 2L) palette[2L] else
              if (is_gene_val) "#2166AC" else "#4E79A7"
    p <- p + scale_color_gradient(low = low, high = high, name = colour_by)
  }

  # ── Cluster centroid labels ──────────────────────────────────────────────────
  if (isTRUE(label_clusters) && is_discrete) {
    centroids <- stats::aggregate(cbind(dim1, dim2) ~ colour_val,
                                   data = df, FUN = mean)
    p <- p + geom_text(data        = centroids,
                       aes(x = dim1, y = dim2, label = colour_val),
                       colour      = "black",
                       size        = label_size,
                       fontface    = "bold",
                       inherit.aes = FALSE)
  }

  # ── Plot boundaries — add 3% breathing room around the computed limits ───────
  x_pad <- diff(xlim) * 0.03
  y_pad <- diff(ylim) * 0.03
  p <- p + coord_cartesian(
    xlim = c(xlim[1] - x_pad, xlim[2] + x_pad),
    ylim = c(ylim[1] - y_pad, ylim[2] + y_pad)
  )

  # ── Legend ───────────────────────────────────────────────────────────────────
  if (!isTRUE(show_legend))
    p <- p + theme(legend.position = "none")

  p
}


# Per-cluster normalised KDE contours.
# For each cluster: compute 2D KDE, normalise density to [0,1], extract
# contour lines at evenly spaced quantiles.  Returns a data.frame suitable
# for geom_path() with columns dim1, dim2, colour_val, line_id.
.aspis_dynamic_contours <- function(df, n_levels = 10L, pad = 0.10) {
  lvls    <- levels(df$colour_val)
  breaks  <- seq(0.05, 0.95, length.out = as.integer(n_levels))

  # Use global limits for all clusters so the KDE surface tapers to near-zero
  # at the grid edges — prevents contour lines being cut off by a cluster-local
  # rectangular boundary.  pad matches boundary_pad so the KDE grid and the
  # plot window are always consistent.
  x_pad  <- diff(range(df$dim1)) * pad
  y_pad  <- diff(range(df$dim2)) * pad
  g_lims <- c(range(df$dim1)[1] - x_pad, range(df$dim1)[2] + x_pad,
              range(df$dim2)[1] - y_pad, range(df$dim2)[2] + y_pad)

  out_list <- lapply(lvls, function(lv) {
    sub <- df[df$colour_val == lv, c("dim1", "dim2")]
    if (nrow(sub) < 10L) return(NULL)

    kde    <- MASS::kde2d(sub$dim1, sub$dim2, n = 200L, lims = g_lims)
    z      <- kde$z
    z_norm <- (z - min(z)) / (max(z) - min(z))

    clines <- grDevices::contourLines(kde$x, kde$y, z_norm, levels = breaks)
    if (length(clines) == 0L) return(NULL)

    do.call(rbind, lapply(seq_along(clines), function(i) {
      cl <- clines[[i]]
      data.frame(
        dim1       = cl$x,
        dim2       = cl$y,
        colour_val = lv,
        line_id    = paste0(lv, "_", i),
        stringsAsFactors = FALSE
      )
    }))
  })

  out            <- do.call(rbind, Filter(Negate(is.null), out_list))
  out$colour_val <- factor(out$colour_val, levels = lvls)
  out
}


# 20-colour categorical palette (Tableau-inspired).
.aspis_discrete_palette <- function() {
  c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#79706E", "#D4A6C8", "#86BCB6", "#FFBE7D", "#8CD17D",
    "#499894", "#E6D16A", "#D37295", "#FABFD2", "#B6992D")
}
