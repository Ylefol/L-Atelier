# Resolves plot colours with priority: an explicit call-time `palette`
# argument > colours carried through from KERAUNOS_propeller_proportions()
# (propeller_result$params$cluster_colors / $group_colors, if the caller
# supplied them there) > .aspis_discrete_palette() default. Always returns a
# vector named by `levels`, in `levels` order -- ready to hand straight to
# scale_*_manual(values = ...).
.aspis_resolve_palette <- function(explicit, stored, levels) {
  pal <- explicit %||% stored %||% .aspis_discrete_palette()

  if (is.null(names(pal))) {
    # Positional palette (the default set, or a caller-supplied unnamed
    # vector) -- recycle by position, then assign names to match `levels`.
    pal <- rep_len(pal, length(levels))
    names(pal) <- levels
  }
  # else: already a category->colour map (cluster_colors/group_colors from
  # KERAUNOS, or a named `palette` argument) -- go straight to a by-name
  # lookup below. Do NOT rep_len() a named vector first: rep_len() silently
  # drops names even when the length is already correct, which previously
  # caused pal[levels] to reassign the ORIGINAL (wrong) colour order onto
  # `levels`' names instead of doing a real by-name lookup.

  pal[levels]
}


#' Plot a propeller cell-type proportion test result
#'
#' Visualises the output of \code{\link{KERAUNOS_propeller_proportions}}: one
#' facet per cluster, showing each biological sample's proportion as a point
#' (coloured by group) with a crossbar at the group mean, plus the test's
#' FDR (or raw p-value) printed in the top-left corner of each facet.
#'
#' A boxplot is deliberately not used here -- with as few as 2 samples per
#' group (as in designs like Nico_midbrain), a boxplot's quartiles are
#' degenerate and visually imply a distribution the data doesn't have.
#' Individual points plus a mean crossbar show exactly what was tested
#' without overstating it.
#'
#' @param propeller_result A \code{keraunos_propeller} object from
#'   \code{\link{KERAUNOS_propeller_proportions}}.
#' @param palette Character vector of colours, one per group level.
#'   Overrides any \code{group_colors} carried on \code{propeller_result}.
#'   Default \code{NULL}: uses \code{propeller_result$params$group_colors}
#'   if it was supplied to \code{\link{KERAUNOS_propeller_proportions}},
#'   otherwise falls back to the built-in 20-colour Tableau-inspired set
#'   (same as \code{\link{ASPIS_plot_umap}}).
#' @param ncol Integer. Facets per row. Default \code{NULL} lets
#'   \code{facet_wrap} choose.
#' @param point_size Numeric. Sample point size. Default \code{2.5}.
#' @param jitter_width Numeric. Horizontal jitter applied to points so
#'   overlapping samples stay visible. Default \code{0.08}.
#' @param show_stats Logical. Annotate each facet with its test statistic.
#'   Default \code{TRUE}.
#' @param stat_label Character. \code{"fdr"} (default) or \code{"pvalue"} --
#'   which column from \code{propeller_result$results} to display.
#' @param order_by Character. \code{"significance"} (default) orders facets
#'   by ascending P.Value (most significant first, matching
#'   \code{propeller_result$results}); \code{"cluster"} orders alphabetically/
#'   numerically by cluster label.
#' @param title Character. Plot title. Default \code{NULL}.
#'
#' @return A ggplot object, faceted by cluster.
#' @export
ASPIS_plot_propeller <- function(propeller_result, palette = NULL, ncol = NULL,
                                  point_size = 2.5, jitter_width = 0.08,
                                  show_stats = TRUE,
                                  stat_label = c("fdr", "pvalue"),
                                  order_by = c("significance", "cluster"),
                                  title = NULL) {

  if (!inherits(propeller_result, "keraunos_propeller"))
    stop("'propeller_result' must be a keraunos_propeller object (see KERAUNOS_propeller_proportions()).",
         call. = FALSE)
  stat_label <- match.arg(stat_label)
  order_by   <- match.arg(order_by)

  df      <- propeller_result$proportions
  results <- propeller_result$results

  cluster_order <- if (order_by == "significance")
    as.character(results$cluster)
  else
    sort(unique(as.character(df$cluster)))

  df$cluster <- factor(as.character(df$cluster), levels = cluster_order)

  group_levels <- names(propeller_result$params$n_per_group)
  df$group <- factor(as.character(df$group), levels = group_levels)

  pal <- .aspis_resolve_palette(palette, propeller_result$params$group_colors, group_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = proportion, colour = group)) +
    ggplot2::geom_point(size = point_size,
                         position = ggplot2::position_jitter(width = jitter_width, height = 0)) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                           colour = "grey30", linewidth = 0.4, fatten = 0) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::facet_wrap(~cluster, scales = "free_y", ncol = ncol) +
    ggplot2::labs(title = title, x = NULL, y = "Proportion of cells",
                  colour = propeller_result$params$group_col) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
    )

  if (isTRUE(show_stats)) {
    stat_col <- if (stat_label == "fdr") "FDR" else "P.Value"
    stat_df <- data.frame(
      cluster = factor(as.character(results$cluster), levels = cluster_order),
      value   = results[[stat_col]]
    )
    stat_df$label <- ifelse(
      is.na(stat_df$value), "n/a",
      paste0(if (stat_label == "fdr") "FDR = " else "p = ",
             formatC(stat_df$value, format = "g", digits = 2))
    )
    p <- p + ggplot2::geom_text(
      data = stat_df,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 3.2, colour = "grey20"
    )
  }

  p
}


#' Plot a 100% stacked bar chart of per-sample cell-type composition
#'
#' Complements \code{\link{ASPIS_plot_propeller}}: that function shows one
#' facet per cluster (a statistical test result in isolation), which never
#' makes the compositional nature of the data visible. This shows the full
#' picture instead -- one bar per sample, cluster proportions stacked to the
#' full bar height, so it's immediately clear that every sample's cell-type
#' proportions sum to 100% (an increase in one cluster necessarily comes at
#' the expense of others).
#'
#' Each legend entry is two lines: the cluster label, then \code{FDR = ...}
#' from the same propeller test -- so the statistical context from
#' \code{\link{ASPIS_plot_propeller}} carries over without needing to
#' cross-reference a second plot. Legend key spacing is widened accordingly
#' so the two-line entries don't crowd each other.
#'
#' @param propeller_result A \code{keraunos_propeller} object from
#'   \code{\link{KERAUNOS_propeller_proportions}} (uses its
#'   \code{$proportions} and \code{$results$baseline_prop}/\code{$FDR}
#'   fields).
#' @param palette Character vector of colours, one per cluster. Overrides
#'   any \code{cluster_colors} carried on \code{propeller_result}. Default
#'   \code{NULL}: uses \code{propeller_result$params$cluster_colors} if it
#'   was supplied to \code{\link{KERAUNOS_propeller_proportions}}, otherwise
#'   falls back to the built-in 20-colour Tableau-inspired set (recycled if
#'   there are more clusters than colours).
#' @param facet_by_group Logical. Split samples into group-labelled facet
#'   panels so they're grouped visually by experimental group rather than in
#'   arbitrary order. Default \code{TRUE}.
#' @param title Character. Plot title. Default \code{NULL}.
#'
#' @return A ggplot object.
#' @export
ASPIS_plot_composition <- function(propeller_result, palette = NULL,
                                    facet_by_group = TRUE, title = NULL) {

  if (!inherits(propeller_result, "keraunos_propeller"))
    stop("'propeller_result' must be a keraunos_propeller object (see KERAUNOS_propeller_proportions()).",
         call. = FALSE)

  df      <- propeller_result$proportions
  results <- propeller_result$results

  # Largest overall cluster stacked at the bottom -- a stable reference
  # ordering (baseline proportion doesn't depend on the group comparison),
  # unlike ASPIS_plot_propeller's significance-based ordering.
  cluster_order <- results$cluster[order(-results$baseline_prop)]
  df$cluster <- factor(as.character(df$cluster), levels = cluster_order)

  group_levels <- names(propeller_result$params$n_per_group)
  df$group <- factor(as.character(df$group), levels = group_levels)

  pal <- .aspis_resolve_palette(palette, propeller_result$params$cluster_colors, cluster_order)

  # Legend entries double as the FDR readout for each cluster, so the
  # composition plot carries the same statistical context as
  # ASPIS_plot_propeller without needing a second plot to cross-reference.
  fdr_by_cluster <- results$FDR[match(cluster_order, results$cluster)]
  fdr_text <- ifelse(is.na(fdr_by_cluster), "n/a",
                      formatC(fdr_by_cluster, format = "g", digits = 2))
  legend_labels <- setNames(paste0(cluster_order, "\nFDR = ", fdr_text), cluster_order)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = proportion, fill = cluster)) +
    ggplot2::geom_col(position = "fill", width = 0.8) +
    ggplot2::scale_fill_manual(values = pal, breaks = cluster_order,
                                labels = legend_labels[cluster_order],
                                name = propeller_result$params$cluster_col) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                                 expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(title = title, x = propeller_result$params$sample_col,
                  y = "Proportion of cells") +
    ggplot2::guides(fill = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.ticks         = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "grey95", colour = NA),
      legend.text         = ggplot2::element_text(size = 9, lineheight = 0.9,
                                                    margin = ggplot2::margin(t = 2, b = 8)),
      legend.key.height    = grid::unit(1.1, "lines"),
      legend.spacing.y     = grid::unit(0.3, "cm")
    )

  if (isTRUE(facet_by_group))
    p <- p + ggplot2::facet_grid(~group, scales = "free_x", space = "free_x")

  p
}
