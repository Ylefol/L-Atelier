#' Aether - Statistical Comparison Plots
#'
#' @description Visualisation functions for statistical comparison results,
#' including group-level comparisons of distributions and proportions.


###############################################################################
########### Peak Annotation Distribution Comparison Plot ###########
###############################################################################

#' Plot Peak Annotation Distribution Comparison Between Groups
#'
#' @description
#' Visualises the output of \code{ARTEMIS_compare_annotation_distribution()}.
#' For each genomic feature, shows per-sample proportions as dots (coloured by
#' group) with a horizontal mean line, and overlays significance brackets
#' where the group comparison passes the chosen threshold.
#'
#' @param annotation_test List returned by
#'   \code{ARTEMIS_compare_annotation_distribution()}.
#' @param colors Named character vector mapping group labels to colours.
#'   If \code{NULL}, uses a default blue/red palette.
#' @param title Character. Plot title.
#' @param p_threshold Numeric. p.adj threshold for significance brackets.
#'   Default: \code{0.05}.
#' @param show_ns Logical. If \code{TRUE} (default), label non-significant
#'   features with "ns". If \code{FALSE}, omit brackets for non-significant
#'   features entirely.
#' @param nrow Integer. Number of rows in the facet grid. Default: \code{2}.
#'
#' @return A ggplot object.
#'
#' @details
#' Significance is encoded as: *** p.adj < 0.001, ** < 0.01, * < 0.05,
#' ns >= 0.05.
#'
#' The y-axis is free per facet so that features with very different overall
#' proportions (e.g., Intergenic ~40% vs 5' UTR ~1%) remain legible.
#'
#' @seealso \code{\link{ARTEMIS_compare_annotation_distribution}},
#'   \code{\link{AETHER_plot_annotation_bar}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test_result <- ARTEMIS_compare_annotation_distribution(
#'   annotated_list = annotated_list,
#'   metadata       = meta,
#'   group_col      = "group"
#' )
#' p <- AETHER_plot_annotation_comparison(test_result)
#' ggsave("annotation_comparison.pdf", p, width = 14, height = 8)
#' }
AETHER_plot_annotation_comparison <- function(annotation_test,
                                               colors      = NULL,
                                               title       = "Peak Annotation Distribution by Group",
                                               p_threshold = 0.05,
                                               show_ns     = TRUE,
                                               nrow        = 2L) {

  prop_df   <- annotation_test$counts
  results   <- annotation_test$results
  groups    <- annotation_test$groups
  group_col <- annotation_test$group_col

  if (is.null(colors)) {
    colors <- c("#377EB8", "#E41A1C")
    names(colors) <- groups
  } else if (is.null(names(colors))) {
    names(colors) <- groups
  }

  # Feature display order (consistent with AETHER_plot_annotation_bar)
  preferred_order <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR",
                       "Downstream", "Intergenic", "Other")
  all_features  <- unique(prop_df$feature)
  feature_order <- c(intersect(preferred_order, all_features),
                     setdiff(all_features, preferred_order))

  prop_df$feature      <- factor(prop_df$feature,      levels = feature_order)
  prop_df[[group_col]] <- factor(prop_df[[group_col]], levels = groups)
  results$feature      <- factor(results$feature,      levels = feature_order)

  # Significance labels
  .sig_label <- function(p_adj) {
    ifelse(is.na(p_adj),        "",
    ifelse(p_adj < 0.001,       "***",
    ifelse(p_adj < 0.01,        "**",
    ifelse(p_adj < p_threshold, "*",
           if (show_ns) "ns" else ""))))
  }
  results$sig_label <- .sig_label(results$p.adj)

  # ---------------------------------------------------------------------------
  # Base plot: dots + mean crossbar per group per feature
  # ---------------------------------------------------------------------------
  p <- ggplot(prop_df,
              aes(x     = .data[[group_col]],
                  y     = count,
                  color = .data[[group_col]])) +
    geom_jitter(width = 0.12, size = 2.4, alpha = 0.85) +
    stat_summary(fun      = mean,
                 geom     = "errorbar",
                 aes(ymin = after_stat(y), ymax = after_stat(y)),
                 width    = 0.45,
                 linewidth = 1.0,
                 color    = "grey20") +
    scale_color_manual(values = colors, name = group_col) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    facet_wrap(~ feature, scales = "free_y", nrow = nrow) +
    labs(title = title, x = NULL, y = "Peak count") +
    theme_bw(base_size = 12) +
    theme(
      strip.background   = element_rect(fill = "grey92"),
      strip.text         = element_text(face = "bold"),
      axis.text.x        = element_text(angle = 30, hjust = 1),
      legend.position    = "bottom",
      plot.title         = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )

  # ---------------------------------------------------------------------------
  # Significance brackets (per-facet, respecting free_y)
  # ---------------------------------------------------------------------------
  sig_ann <- results[results$sig_label != "", c("feature", "sig_label"), drop = FALSE]

  if (nrow(sig_ann) > 0L) {
    # Per-feature max to position brackets just above the data
    feat_max <- aggregate(count ~ feature, data = prop_df,
                          FUN = function(x) max(x, na.rm = TRUE))
    colnames(feat_max)[2] <- "feat_max"

    sig_ann <- merge(sig_ann, feat_max, by = "feature")
    sig_ann$feature   <- factor(sig_ann$feature, levels = feature_order)
    sig_ann$bracket_y <- sig_ann$feat_max * 1.10
    sig_ann$label_y   <- sig_ann$feat_max * 1.22
    # x in discrete-scale coordinates: group 1 = 1, group 2 = 2
    sig_ann$x1   <- 1L
    sig_ann$x2   <- 2L
    sig_ann$xmid <- 1.5

    p <- p +
      geom_segment(data        = sig_ann,
                   aes(x = x1, xend = x2, y = bracket_y, yend = bracket_y),
                   inherit.aes = FALSE,
                   color       = "grey40",
                   linewidth   = 0.4) +
      geom_text(data        = sig_ann,
                aes(x = xmid, y = label_y, label = sig_label),
                inherit.aes = FALSE,
                size        = 3.5,
                color       = "grey20")
  }

  return(p)
}
