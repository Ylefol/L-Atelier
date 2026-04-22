# ==============================================================================
# AETHER - ATAC-seq QC visualisation
# ==============================================================================


#' Plot TSS Enrichment Profile and Per-sample Scores
#'
#' @description
#' Produces two complementary plots from a \code{hades_tss_enrichment} object
#' returned by \code{\link{HADES_tss_enrichment}}:
#' \enumerate{
#'   \item \strong{Profile plot}: normalised mean signal profile across all TSSs
#'     per sample (lines), with a vertical dashed line at the TSS and a
#'     horizontal dotted line at y = 1 (background level). A well-enriched
#'     library shows a sharp peak at the centre.
#'   \item \strong{Score plot}: per-sample TSS enrichment score bar chart with
#'     an optional horizontal threshold line.
#' }
#'
#' @param tss_result A \code{hades_tss_enrichment} object from
#'   \code{\link{HADES_tss_enrichment}}.
#' @param color_by Character. Column used to colour both plots. Can be
#'   \code{"sample"} (each sample a unique colour), \code{"pass"} (pass vs fail
#'   relative to the stored threshold), or the name of any metadata column
#'   carried through from the sample sheet (e.g. \code{"condition"},
#'   \code{"batch"}). Falls back to \code{"sample"} with a message if the
#'   requested column is not present in \code{tss_result$scores}.
#' @param profile_alpha Numeric in (0, 1]. Line transparency in the profile
#'   plot. Default \code{0.85}.
#' @param show_threshold Logical. Draw a horizontal reference line at the
#'   stored threshold on the score plot. Default \code{TRUE}.
#' @param title Character or \code{NULL}. Title for the profile plot. If
#'   \code{NULL} a default title is used.
#'
#' @return A named list with two \code{ggplot} objects:
#'   \describe{
#'     \item{\code{$profile}}{Normalised signal profiles across the TSS window.}
#'     \item{\code{$score}}{Per-sample TSS enrichment score bar chart.}
#'   }
#'   Combine with \code{patchwork::wrap_plots(result, ncol = 1)} if desired.
#'
#' @seealso \code{\link{HADES_tss_enrichment}}
#' @export
AETHER_plot_tss_enrichment <- function(
  tss_result,
  color_by       = "condition",
  profile_alpha  = 0.85,
  show_threshold = TRUE,
  title          = NULL
) {

  if (!inherits(tss_result, "hades_tss_enrichment"))
    stop("'tss_result' must be a hades_tss_enrichment object from ",
         "HADES_tss_enrichment().", call. = FALSE)

  params   <- tss_result$params
  scores   <- tss_result$scores
  profiles <- tss_result$profiles

  # Resolve colour column — accept any column present in scores
  reserved <- c("sample", "pass")
  if (color_by == "sample") {
    color_col <- "sample_id"
  } else if (color_by == "pass") {
    color_col <- "pass"
  } else if (color_by %in% colnames(scores)) {
    color_col <- color_by
  } else {
    cat("[AETHER] color_by column '", color_by,
        "' not found in scores — colouring by sample.\n", sep = "")
    color_col <- "sample_id"
  }

  # X-axis positions in bp relative to TSS
  x_pos <- seq(-params$window, params$window, length.out = params$n_bins)

  # --- Build long data.frame for profile plot ----------------------------------
  # Repeat every metadata column from scores across the bins for each sample
  meta_cols <- setdiff(colnames(scores), c("sample_id", "tss_enrichment_score"))

  profile_df <- do.call(rbind, lapply(scores$sample_id, function(sid) {
    prof <- profiles[[sid]]
    row  <- scores[scores$sample_id == sid, , drop = FALSE]
    df   <- data.frame(
      sample_id    = sid,
      bin_position = x_pos,
      norm_signal  = prof,
      stringsAsFactors = FALSE
    )
    for (col in meta_cols)
      df[[col]] <- row[[col]]
    df
  }))

  # Convert pass to readable labels for colouring
  if (color_col == "pass")
    profile_df$pass <- ifelse(profile_df$pass, "Pass", "Fail")

  # --- Profile plot ------------------------------------------------------------
  half_kb <- params$window / 1000

  p_profile <- ggplot(
    profile_df,
    aes(x      = .data$bin_position,
        y      = .data$norm_signal,
        colour = .data[[color_col]],
        group  = .data$sample_id)
  ) +
    geom_line(alpha = profile_alpha, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dotted",
               colour = "grey70", linewidth = 0.4) +
    scale_x_continuous(
      breaks = c(-params$window, -params$window / 2L, 0L,
                  params$window / 2L,  params$window),
      labels = c(paste0("-", half_kb, "kb"),
                 paste0("-", half_kb / 2, "kb"),
                 "TSS",
                 paste0("+", half_kb / 2, "kb"),
                 paste0("+", half_kb, "kb"))
    ) +
    labs(
      x      = "Position relative to TSS (bp)",
      y      = "Normalised signal",
      colour = color_col,
      title  = if (!is.null(title)) title else "TSS Enrichment Profile"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(face = "bold", size = 11)
    )

  # --- Score bar plot ----------------------------------------------------------
  scores_plt <- scores
  # Carry through pass labels if colouring by pass
  if (color_col == "pass")
    scores_plt$pass <- ifelse(scores_plt$pass, "Pass", "Fail")

  # Fix sample order to match profile plot (input order)
  scores_plt$sample_id <- factor(scores_plt$sample_id,
                                  levels = scores_plt$sample_id)

  p_score <- ggplot(
    scores_plt,
    aes(x    = .data$sample_id,
        y    = .data$tss_enrichment_score,
        fill = .data[[color_col]])
  ) +
    geom_col(width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x     = NULL,
      y     = "TSS enrichment score",
      fill  = color_col,
      title = "Per-sample TSS enrichment score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = 11)
    )

  if (isTRUE(show_threshold) && !is.na(params$threshold)) {
    p_score <- p_score +
      geom_hline(yintercept = params$threshold, linetype = "dashed",
                 colour = "firebrick", linewidth = 0.5) +
      annotate("text",
               x      = Inf,
               y      = params$threshold,
               label  = paste0("threshold = ", params$threshold),
               hjust  = 1.05,
               vjust  = -0.4,
               size   = 3,
               colour = "firebrick")
  }

  list(profile = p_profile, score = p_score)
}
