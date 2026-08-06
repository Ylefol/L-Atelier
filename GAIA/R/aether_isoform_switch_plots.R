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
#'   "<reference> vs <experiment>" from \code{switch_result$params}.
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
    title <- paste(switch_result$params$reference, "vs", switch_result$params$experiment)

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


# ==============================================================================
# Consequence summary plot
# ==============================================================================

#' Functional Consequence Summary Bar Plot
#'
#' Bar plot of the number of switches with each functional consequence type
#' (domain gain/loss, coding potential, signal peptide, subcellular
#' location, IDR, etc.), faceted by feature. Built on
#' \code{IsoformSwitchAnalyzeR::extractConsequenceSummary()}, called live on
#' \code{switch_result$switch_list} -- this is a cheap tabulation over the
#' consequence annotations already computed by
#' \code{ARTEMIS_isoform_switch_consequences()}, so nothing needs to be
#' pre-stored under a separate field the way \code{$splicing_summary}/
#' \code{$splicing_enrichment} are.
#'
#' @param switch_result An \code{artemis_isoform_switch} object that has
#'   already been run through \code{ARTEMIS_isoform_switch_consequences()}.
#' @param count_by Character. \code{"genes"} (default) plots
#'   \code{nrGenesWithConsequences}; \code{"isoforms"} plots
#'   \code{nrIsoWithConsequences}.
#' @param consequences_to_analyze Character vector. Passed through to
#'   \code{extractConsequenceSummary()}'s \code{consequencesToAnalyze}.
#'   Default: \code{"all"} (every consequence type present in
#'   \code{switch_result$switch_list}).
#' @param title Character or NULL. Plot title. Default: NULL, which builds
#'   "<reference> vs <experiment>" from \code{switch_result$params}.
#' @param fill_color Character. Bar fill color. Default: \code{"#3A3A3A"}.
#'
#' @return A ggplot object.
#'
#' @details
#' Unlike \code{AETHER_plot_splicing_summary()}, consequence categories
#' don't reduce to a clean two-way "more/less" split -- a feature like
#' "Domains identified" has gain/loss/switch outcomes, "Coding potential"
#' has coding/non-coding, "Sub cell location" has gain/loss/switch, etc. --
#' so bars are a single color and faceted by \code{featureCompared} (one
#' panel per feature) with a free x scale, rather than dodged by direction.
#' Panels only appear for consequence types that were actually annotated
#' (i.e. the external tool files supplied to
#' \code{ARTEMIS_isoform_switch_consequences()}).
#'
#' \strong{Direction}: every \code{switchConsequence} label (e.g. "Domain
#' gain", "Domain loss", "Transcript is Noncoding") describes the isoform
#' \emph{upregulated in \code{switch_result$params$experiment}} relative to
#' the isoform used more in \code{switch_result$params$reference} --
#' \code{IsoformSwitchAnalyzeR}'s own convention
#' (\code{?analyzeSwitchConsequences}: "switchConsequence ... a short
#' description of the features of the upregulated isoform", where
#' "upregulated" always means higher usage in condition_2/experiment). So
#' e.g. a tall "Domain gain" bar means many genes have the
#' experiment-dominant isoform gaining a domain the reference-dominant
#' isoform didn't have.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_consequence_summary(switch_result)
#' p <- AETHER_plot_consequence_summary(switch_result, count_by = "isoforms")
#' }
#' @export
AETHER_plot_consequence_summary <- function(switch_result,
                                             count_by = c("genes", "isoforms"),
                                             consequences_to_analyze = "all",
                                             title = NULL,
                                             fill_color = "#3A3A3A") {

  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  if (is.null(switch_result$consequence_summary))
    stop("switch_result has no consequence_summary -- run ",
         "ARTEMIS_isoform_switch_consequences() first.", call. = FALSE)

  count_by <- match.arg(count_by)

  df <- IsoformSwitchAnalyzeR::extractConsequenceSummary(
    switch_result$switch_list,
    consequencesToAnalyze = consequences_to_analyze,
    plot         = FALSE,
    returnResult = TRUE
  )

  count_col <- if (count_by == "genes") "nrGenesWithConsequences" else "nrIsoWithConsequences"
  ylab_text <- if (count_by == "genes") "Number of genes" else "Number of isoforms"

  if (is.null(title))
    title <- paste(switch_result$params$reference, "vs", switch_result$params$experiment)
  subtitle <- paste0("Direction is relative to ", switch_result$params$experiment,
                      " -- e.g. \"gain\" = gained in the ", switch_result$params$experiment,
                      "-dominant isoform")

  p <- ggplot(df, aes(x = switchConsequence, y = .data[[count_col]])) +
    geom_col(fill = fill_color, width = 0.7) +
    xlab("Consequence of switch (features of the upregulated isoform)") +
    ylab(ylab_text) +
    ggtitle(title, subtitle = subtitle) +
    theme_light() +
    theme(
      text            = element_text(size = 10),
      plot.title      = element_text(size = 16, face = "bold"),
      plot.subtitle   = element_text(size = 9, face = "italic", color = "gray30"),
      axis.text.x     = element_text(angle = 40, hjust = 1),
      legend.position = "none"
    )

  if (length(unique(df$Comparison)) > 1) {
    p <- p + facet_grid(Comparison ~ featureCompared, scales = "free_x", space = "free_x")
  } else {
    p <- p + facet_wrap(~featureCompared, scales = "free_x", nrow = 1)
  }

  return(p)
}


# ==============================================================================
# Splicing summary / enrichment plots
# ==============================================================================

#' Extract a top-level data.frame field from an artemis_isoform_switch object
#' @keywords internal
.splicing_extract_df <- function(switch_result, field, from_fn) {
  if (!inherits(switch_result, "artemis_isoform_switch"))
    stop("switch_result must be an artemis_isoform_switch object from ",
         "ARTEMIS_isoform_switch()", call. = FALSE)

  df <- switch_result[[field]]
  if (is.null(df))
    stop("switch_result$", field, " not found -- run ", from_fn,
         "() first.", call. = FALSE)
  df
}


#' Splicing Event Summary Bar Plot
#'
#' Bar plot of the total number of alternative splicing events per type
#' (IR, A5, A3, ATSS, ATTS, ES, MES, MEE), split by whether the event
#' occurs in the isoform used \emph{more} or \emph{less} in the experiment
#' group, from \code{switch_result$splicing_summary}
#' (\code{IsoformSwitchAnalyzeR::extractSplicingSummary()} output).
#'
#' @param switch_result An \code{artemis_isoform_switch} object with
#'   splicing classification already run
#'   (\code{ARTEMIS_isoform_switch_splicing()}).
#' @param count_by Character. \code{"genes"} (default) plots
#'   \code{nrGenesWithConsequences}; \code{"isoforms"} plots
#'   \code{nrIsoWithConsequences}.
#' @param title Character or NULL. Plot title. Default: NULL, which builds
#'   "<reference> vs <experiment>" from \code{switch_result$params}.
#' @param colors Named character vector with colors for \code{"more"} and
#'   \code{"less"}. Default uses red/blue.
#'
#' @return A ggplot object.
#'
#' @details
#' This is the total count of events per type -- it does not indicate
#' whether the events belong to significant switches specifically (see
#' \code{AETHER_plot_splicing_enrichment()} for the significance-tested
#' gain/loss comparison). Bars are grouped by \code{AStype} and dodged by
#' direction (isoform used more vs. less in the experiment group).
#'
#' \strong{Direction}: "isoform used more/less" follows
#' \code{IsoformSwitchAnalyzeR}'s \code{isoformUpregulated}/
#' \code{isoformDownregulated} convention from
#' \code{analyzeSwitchConsequences()} -- "more" is the isoform with higher
#' usage in \code{switch_result$params$experiment} (positive dIF, relative
#' to \code{switch_result$params$reference}), not an arbitrary pairwise
#' label. So e.g. a tall "IR in isoform used more" bar means many genes gain
#' intron retention specifically in the experiment-dominant isoform.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_splicing_summary(switch_result)
#' p <- AETHER_plot_splicing_summary(switch_result, count_by = "isoforms")
#' }
#' @export
AETHER_plot_splicing_summary <- function(switch_result,
                                          count_by = c("genes", "isoforms"),
                                          title = NULL,
                                          colors = NULL) {

  count_by <- match.arg(count_by)
  df <- .splicing_extract_df(switch_result, "splicing_summary",
                              "ARTEMIS_isoform_switch_splicing")

  count_col <- if (count_by == "genes") "nrGenesWithConsequences" else "nrIsoWithConsequences"
  ylab_text <- if (count_by == "genes") "Number of genes" else "Number of isoforms"

  if (is.null(title))
    title <- paste(switch_result$params$reference, "vs", switch_result$params$experiment)
  subtitle <- paste0("\"more\"/\"less\" = isoform usage in ", switch_result$params$experiment,
                      " relative to ", switch_result$params$reference)

  if (is.null(colors)) colors <- c(more = "#B31B21", less = "#1465AC")

  df$Direction <- ifelse(grepl("more$", df$splicingResult), "more", "less")
  df$Direction <- factor(df$Direction, levels = c("more", "less"))

  p <- ggplot(df, aes(x = AStype, y = .data[[count_col]], fill = Direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = colors,
      labels = c(more = "isoform used more", less = "isoform used less"),
      name   = NULL
    ) +
    xlab("Alternative splicing type") +
    ylab(ylab_text) +
    ggtitle(title, subtitle = subtitle) +
    theme_light() +
    theme(
      text            = element_text(size = 10),
      plot.title      = element_text(size = 16, face = "bold"),
      plot.subtitle   = element_text(size = 9, face = "italic", color = "gray30"),
      legend.position = "bottom"
    )

  if (length(unique(df$Comparison)) > 1) p <- p + facet_wrap(~Comparison)

  return(p)
}


#' Splicing Event Enrichment Plot
#'
#' Point-range plot of the proportion of gain-vs-loss events per splicing
#' type (\code{propUp}, with confidence interval), from
#' \code{switch_result$splicing_enrichment}
#' (\code{IsoformSwitchAnalyzeR::extractSplicingEnrichment()} output). Unlike
#' \code{AETHER_plot_splicing_summary()}, this compares the two isoforms
#' within each switch directly and is significance-tested per type
#' (\code{prop.test()}, FDR-corrected).
#'
#' @param switch_result An \code{artemis_isoform_switch} object with
#'   splicing classification already run
#'   (\code{ARTEMIS_isoform_switch_splicing()}).
#' @param min_events Integer. Minimum total events (\code{nUp + nDown})
#'   required for a splicing type to be plotted. Default: 10, matching
#'   \code{extractSplicingEnrichment()}'s own \code{minEventsForPlotting}
#'   default.
#' @param title Character or NULL. Plot title. Default: NULL, which builds
#'   "<reference> vs <experiment>" from \code{switch_result$params}.
#' @param colors Named character vector with colors for \code{"TRUE"} and
#'   \code{"FALSE"} significance. Default uses red/gray.
#'
#' @return A ggplot object.
#'
#' @details
#' \code{propUp} is the fraction of gain/loss events (per type) that are a
#' \emph{gain} -- e.g. for \code{IR}, the fraction of IR gain-or-loss events
#' that are a gain of intron retention in the isoform used more (higher
#' usage, positive dIF) in \code{switch_result$params$experiment} relative
#' to \code{switch_result$params$reference} -- \code{IsoformSwitchAnalyzeR}'s
#' \code{isoformUpregulated} convention, same as
#' \code{AETHER_plot_splicing_summary()}. A dashed reference line at 0.5
#' marks "no directional preference"; points whose confidence interval
#' excludes 0.5 are the ones \code{prop.test()} calls significant
#' (FDR-corrected, \code{Significant} column). Point size is scaled by the
#' total number of events (\code{nUp + nDown}) so sparsely-supported types
#' are visually de-emphasized.
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_splicing_enrichment(switch_result)
#' p <- AETHER_plot_splicing_enrichment(switch_result, min_events = 5)
#' }
#' @export
AETHER_plot_splicing_enrichment <- function(switch_result,
                                             min_events = 10,
                                             title = NULL,
                                             colors = NULL) {

  df <- .splicing_extract_df(switch_result, "splicing_enrichment",
                              "ARTEMIS_isoform_switch_splicing")

  df$n_total <- df$nUp + df$nDown
  df <- df[df$n_total >= min_events, , drop = FALSE]
  if (nrow(df) == 0)
    stop("No splicing types have nUp + nDown >= min_events (", min_events,
         "). Lower min_events to plot anything.", call. = FALSE)

  if (is.null(title))
    title <- paste(switch_result$params$reference, "vs", switch_result$params$experiment)
  subtitle <- paste0("Proportion is of the isoform used more in ", switch_result$params$experiment,
                      " (vs ", switch_result$params$reference, ")")

  if (is.null(colors)) colors <- c(`TRUE` = "#B31B21", `FALSE` = "darkgray")

  df$Significant <- as.character(df$Significant)

  p <- ggplot(df, aes(x = AStype, y = propUp, color = Significant)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = propUpCiLo, ymax = propUpCiHi), width = 0.2) +
    geom_point(aes(size = n_total)) +
    scale_color_manual(values = colors, name = "Significant (FDR)") +
    scale_size_continuous(name = "Total events\n(nUp + nDown)") +
    coord_flip() +
    xlab("Alternative splicing type") +
    ylab("Proportion of events that are a gain") +
    ggtitle(title, subtitle = subtitle) +
    theme_light() +
    theme(
      text            = element_text(size = 10),
      plot.title      = element_text(size = 16, face = "bold"),
      plot.subtitle   = element_text(size = 9, face = "italic", color = "gray30"),
      legend.position = "right"
    )

  if (length(unique(df$Comparison)) > 1) p <- p + facet_wrap(~Comparison)

  return(p)
}
