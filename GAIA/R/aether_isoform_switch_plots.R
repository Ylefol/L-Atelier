###############################################################################
########### Isoform Switch Visualizations ###########
###############################################################################

# ==============================================================================
# Internal helpers
# ==============================================================================

#' Extract per-isoform results data.frame from an artemis_isoform_switch object
#' @keywords internal
.iso_extract_df <- function(switch_result, required) {
  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  df <- switch_result$isoform_results
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Required column(s) not found in isoform_results: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  df
}


#' Assign four-category significance labels to an isoform switch data.frame
#'
#' Adds a .cat column: "up", "down", "low_dIF", "non_sig". "up"/"down" refer to
#' increased/decreased isoform usage in the experiment group. NA q-values are
#' treated as non-significant.
#' @keywords internal
.iso_assign_categories <- function(df, q_col, alpha, dIF_cutoff) {
  qval <- df[[q_col]]
  dIF  <- df$dIF

  sig  <- !is.na(qval) & qval < alpha
  up   <- sig & !is.na(dIF) & dIF >  dIF_cutoff
  down <- sig & !is.na(dIF) & dIF < -dIF_cutoff
  low  <- sig & !is.na(dIF) & abs(dIF) <= dIF_cutoff

  df$.cat <- "non_sig"
  df$.cat[up]   <- "up"
  df$.cat[down] <- "down"
  df$.cat[low]  <- "low_dIF"
  df
}


# ==============================================================================
# Volcano plot
# ==============================================================================

#' Volcano Plot for Isoform Switch Results
#'
#' Creates a four-category volcano plot (-log10 q-value vs delta Isoform
#' Fraction) from an \code{artemis_isoform_switch} object, with embedded
#' counts in the legend and optional gene/isoform labeling. Categories:
#' increased usage, decreased usage, low-dIF (significant but below the dIF
#' cutoff), and non-significant.
#'
#' @param switch_result An \code{artemis_isoform_switch} object from
#'   \code{ARTEMIS_isoform_switch()}.
#' @param label_col Character. Column in \code{isoform_results} to use for
#'   point labels. Default: \code{"gene_id"}.
#' @param genes_of_interest Character vector. Values in \code{label_col} to
#'   always label. Default: NULL.
#' @param show_non_sig_interest Logical. If FALSE, genes of interest that do
#'   not meet both thresholds are not labeled. Default: TRUE.
#' @param label_top_n Integer. Label the top N isoforms by \code{q_col}
#'   regardless of direction. Default: 0 (disabled).
#' @param q_col Character. Q-value column to use for significance filtering
#'   and the y-axis: \code{"isoform_switch_q_value"} (default, isoform-level)
#'   or \code{"gene_switch_q_value"} (gene-level).
#' @param dIF_cutoff Numeric or NULL. Minimum absolute dIF to call a switch
#'   "increased"/"decreased" usage rather than "low-dIF". Default: NULL, which
#'   uses \code{switch_result$params$dIF_cutoff}.
#' @param alpha Numeric or NULL. Significance threshold on \code{q_col}.
#'   Default: NULL, which uses \code{switch_result$params$alpha}.
#' @param title Character or NULL. Plot title. Default: NULL, which builds
#'   "<experiment> vs <reference>" from \code{switch_result$params}.
#' @param colors Named character vector with colors for \code{"up"},
#'   \code{"down"}, \code{"low_dIF"}, \code{"non_sig"}. Default uses
#'   red/blue/green/gray.
#' @param point_size Numeric. Point size. Default: 0.8.
#' @param point_alpha Numeric. Point transparency. Default: 0.7.
#'
#' @return A ggplot object.
#'
#' @details
#' The four categories are:
#' \itemize{
#'   \item \strong{increased usage}: \code{dIF > dIF_cutoff} AND
#'     \code{q_col < alpha}
#'   \item \strong{decreased usage}: \code{dIF < -dIF_cutoff} AND
#'     \code{q_col < alpha}
#'   \item \strong{low-dIF}: \code{|dIF| <= dIF_cutoff} AND \code{q_col < alpha}
#'   \item \strong{non-significant}: \code{q_col >= alpha} or NA
#' }
#' A horizontal dashed line marks \code{alpha} on the \code{-log10(q_col)}
#' scale. Vertical dashed lines mark \code{±dIF_cutoff}.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_isoform_switch_volcano(switch_result)
#'
#' p <- AETHER_plot_isoform_switch_volcano(switch_result, label_top_n = 10,
#'                                         genes_of_interest = c("ENSG00000141510"))
#'
#' p <- AETHER_plot_isoform_switch_volcano(switch_result,
#'                                         q_col = "gene_switch_q_value")
#' }
#' @export
AETHER_plot_isoform_switch_volcano <- function(switch_result,
                                                label_col = "gene_id",
                                                genes_of_interest = NULL,
                                                show_non_sig_interest = TRUE,
                                                label_top_n = 0,
                                                q_col = "isoform_switch_q_value",
                                                dIF_cutoff = NULL,
                                                alpha = NULL,
                                                title = NULL,
                                                colors = NULL,
                                                point_size = 0.8,
                                                point_alpha = 0.7) {

  df <- .iso_extract_df(switch_result, required = c("dIF", q_col, label_col))

  if (is.null(dIF_cutoff)) dIF_cutoff <- switch_result$params$dIF_cutoff
  if (is.null(alpha))      alpha      <- switch_result$params$alpha
  if (is.null(title))
    title <- paste(switch_result$params$experiment, "vs", switch_result$params$reference)

  df <- .iso_assign_categories(df, q_col, alpha, dIF_cutoff)
  df <- df[!is.na(df[[q_col]]), , drop = FALSE]

  if (is.null(colors)) {
    colors <- c(up = "#B31B21", down = "#1465AC", low_dIF = "green", non_sig = "darkgray")
  }

  n_up      <- sum(df$.cat == "up")
  n_down    <- sum(df$.cat == "down")
  n_low_dIF <- sum(df$.cat == "low_dIF")
  n_non_sig <- sum(df$.cat == "non_sig")

  lbl <- c(
    up      = paste0("increased usage | ", q_col, "<", alpha, " (n=", n_up, ")"),
    down    = paste0("decreased usage | ", q_col, "<", alpha, " (n=", n_down, ")"),
    low_dIF = paste0("low-dIF | ", q_col, "<", alpha, " (n=", n_low_dIF, ")"),
    non_sig = paste0("non-significant | ", q_col, "≥", alpha, " (n=", n_non_sig, ")")
  )

  df$Significance <- factor(lbl[df$.cat],
                             levels = c(lbl["up"], lbl["down"], lbl["low_dIF"], lbl["non_sig"]))
  df$.cat <- NULL
  # Non-sig drawn first (background), significant on top
  df <- df[order(df$Significance, decreasing = TRUE), ]

  # Label data: genes of interest
  labs_interest <- df[0, ]
  if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
    labs_interest <- df[df[[label_col]] %in% genes_of_interest, , drop = FALSE]
    if (!show_non_sig_interest) {
      labs_interest <- labs_interest[labs_interest$Significance != lbl["non_sig"], , drop = FALSE]
    }
  }
  labs_interest$.label <- labs_interest[[label_col]]

  # Label data: top N
  labs_top <- df[0, ]
  if (label_top_n > 0) {
    labs_top <- head(df[order(df[[q_col]]), ], label_top_n)
  }
  labs_top$.label <- labs_top[[label_col]]

  p <- ggplot(df, aes(x = dIF, y = -log10(.data[[q_col]]), color = Significance)) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_vline(xintercept = c(-dIF_cutoff, dIF_cutoff),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(
      values = setNames(unname(colors[c("up", "down", "low_dIF", "non_sig")]),
                        unname(lbl[c("up", "down", "low_dIF", "non_sig")])),
      breaks = unname(lbl),
      name   = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab(expression(Delta*"Isoform Fraction (dIF)")) +
    ylab(bquote(-log[10]*"(" * .(q_col) * ")")) +
    ggtitle(title) +
    theme_light() +
    theme(
      text             = element_text(size = 10),
      plot.title       = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.position  = "bottom",
      legend.direction = "vertical"
    )

  if (nrow(labs_interest) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_interest, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 1, segment.colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  if (nrow(labs_top) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labs_top, mapping = aes(label = .label),
      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
      force = 0.5, colour = "black", show.legend = FALSE,
      label.size = 0.5, size = 3
    )
  }

  return(p)
}


# ==============================================================================
# Per-gene switch plot
# ==============================================================================

#' Per-Gene Isoform Switch Plot
#'
#' Thin convenience wrapper around
#' \code{IsoformSwitchAnalyzeR::switchPlot()} that fills in
#' \code{condition1}/\code{condition2}/\code{dIFcutoff}/\code{alphas} from an
#' \code{artemis_isoform_switch} object's \code{$params}, so only a gene (or
#' isoform) needs to be supplied.
#'
#' @param switch_result An \code{artemis_isoform_switch} object from
#'   \code{ARTEMIS_isoform_switch()}.
#' @param gene_id Character. A single \code{gene_id} to plot, matching the
#'   \code{gene_id} column of \code{switch_result$isoform_results} (not a gene
#'   name/symbol). Mutually exclusive with \code{isoform_id}.
#' @param isoform_id Character vector. One or more isoform IDs (from the same
#'   gene) to plot, matching \code{switch_result$isoform_results$isoform_id}.
#'   Alternative to \code{gene_id}.
#' @param IF_cutoff Numeric. Minimum isoform-fraction contribution to gene
#'   expression (in at least one condition) required for an isoform to be
#'   plotted. This is \code{IsoformSwitchAnalyzeR}'s own display filter, unrelated
#'   to \code{switch_result$params$dIF_cutoff}. Default: \code{0.05}.
#' @param dIF_cutoff Numeric or NULL. dIF cutoff used to annotate
#'   increased/decreased usage on the transcript plot. Default: NULL, which
#'   uses \code{switch_result$params$dIF_cutoff}.
#' @param alphas Numeric vector of length two, or NULL. Q-value thresholds for
#'   "*" and "***" significance annotation, respectively. Default: NULL, which
#'   uses \code{c(switch_result$params$alpha, 0.001)}.
#' @param ... Additional arguments passed to
#'   \code{IsoformSwitchAnalyzeR::switchPlot()} (e.g. \code{rescaleTranscripts},
#'   \code{logYaxis}, \code{localTheme}, \code{additionalArguments}).
#'
#' @return Invisibly returns the result of
#'   \code{IsoformSwitchAnalyzeR::switchPlot()}. As with the underlying
#'   function, the composite plot is drawn directly to the current graphics
#'   device as a side effect — it is not a single combinable ggplot object.
#'
#' @details
#' \code{condition1}/\code{condition2} are taken from
#' \code{switch_result$params$reference}/\code{$experiment} respectively —
#' this matches how \code{ARTEMIS_isoform_switch()} builds
#' \code{comparisonsToMake} internally (\code{condition_1 = reference},
#' \code{condition_2 = experiment}), so orientation is guaranteed consistent
#' with the rest of the \code{artemis_isoform_switch} object (e.g. the sign
#' of \code{dIF} in \code{AETHER_plot_isoform_switch_volcano()}).
#'
#' \code{gene_id}/\code{isoform_id} are validated against
#' \code{switch_result$isoform_results} before calling
#' \code{switchPlot()}, so a typo or a gene removed during \code{preFilter()}
#' fails with a clear message rather than an opaque error from
#' \pkg{IsoformSwitchAnalyzeR}.
#'
#' @examples
#' \dontrun{
#' AETHER_plot_isoform_switch_gene(switch_result, gene_id = "ENSG00000141510")
#'
#' AETHER_plot_isoform_switch_gene(switch_result,
#'                                 isoform_id = c("ENST00000269305", "ENST00000504290"))
#' }
#' @export
AETHER_plot_isoform_switch_gene <- function(switch_result,
                                             gene_id = NULL,
                                             isoform_id = NULL,
                                             IF_cutoff = 0.05,
                                             dIF_cutoff = NULL,
                                             alphas = NULL,
                                             ...) {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (is.null(gene_id) == is.null(isoform_id))
    stop("Supply exactly one of gene_id or isoform_id.", call. = FALSE)

  if (is.null(dIF_cutoff)) dIF_cutoff <- switch_result$params$dIF_cutoff
  if (is.null(alphas))     alphas     <- c(switch_result$params$alpha, 0.001)

  if (!is.null(gene_id)) {
    known_genes <- unique(switch_result$isoform_results$gene_id)
    if (!gene_id %in% known_genes)
      stop("gene_id '", gene_id, "' not found in switch_result$isoform_results$gene_id ",
           "(it may have been removed by preFilter, or this expects a gene_id ",
           "rather than a gene name/symbol).", call. = FALSE)
  } else {
    known_isoforms <- unique(switch_result$isoform_results$isoform_id)
    missing_iso <- setdiff(isoform_id, known_isoforms)
    if (length(missing_iso) > 0)
      stop("isoform_id(s) not found in switch_result$isoform_results$isoform_id: ",
           paste(missing_iso, collapse = ", "), call. = FALSE)
  }

  invisible(IsoformSwitchAnalyzeR::switchPlot(
    switchAnalyzeRlist = switch_result$switch_list,
    gene       = gene_id,
    isoform_id = isoform_id,
    condition1 = switch_result$params$reference,
    condition2 = switch_result$params$experiment,
    IFcutoff   = IF_cutoff,
    dIFcutoff  = dIF_cutoff,
    alphas     = alphas,
    ...
  ))
}
