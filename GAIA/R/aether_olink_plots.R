#' Olink QC Plots
#'
#' @description Generates QC outputs for an \code{olink_data} object:
#' (1) per-plate NPX boxplot, (2) SampleQC pass/fail summary (table), and
#' optionally (3) a histogram of AssayQC warn fractions (plot).
#'
#' @param olink_data An \code{olink_data} object from \code{ELEUTHIA_load_olink()}.
#' @param plate_col Character. Column in \code{$sample_meta} (or \code{$data})
#'   to group the distribution plot by (typically PlateID). Falls back to
#'   PlateID → SampleID if not found. Default: \code{"PlateID"}.
#' @param plot_warn_proteins Logical. Include a histogram of AssayQC warn
#'   fractions across all proteins, with a vertical line at
#'   \code{warn_threshold}. Default: FALSE.
#' @param warn_threshold Numeric (0–1). Reference line drawn on the warn
#'   fraction histogram. Default: 0.1.
#' @param font_size Numeric. Base font size for all plots. Default: 9.
#' @param title Character or NULL. Optional prefix prepended to each plot title.
#'   Default: NULL.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{npx_distributions}{ggplot. Boxplot of NPX values per plate (or group).}
#'   \item{sample_qc_table}{data.frame. PASS/FAIL counts per SampleType (or
#'     overall if SampleType is absent). More legible than a bar chart when most
#'     samples pass.}
#'   \item{warn_proteins}{(Only if \code{plot_warn_proteins = TRUE}.) ggplot.
#'     Histogram of warn_fraction across all proteins, with a dashed reference
#'     line at \code{warn_threshold}.}
#' }
#'
#' @details
#' SampleQC data is drawn from \code{$sample_meta}, which retains all samples
#' (biological + controls) regardless of the \code{keep_controls} setting used
#' at load time. This gives a complete picture of plate-level QC.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
#'
#' # Basic QC plots
#' qc <- AETHER_plot_olink_qc(ol)
#' print(qc$npx_distributions)
#' print(qc$sample_qc_table)
#'
#' # Include warn protein barplot
#' qc <- AETHER_plot_olink_qc(ol, plot_warn_proteins = TRUE, warn_threshold = 0.05)
#' print(qc$warn_proteins)
#' }
AETHER_plot_olink_qc <- function(olink_data,
                                  plate_col          = "PlateID",
                                  plot_warn_proteins = FALSE,
                                  warn_threshold     = 0.1,
                                  font_size          = 9,
                                  title              = NULL) {

  if (!inherits(olink_data, "olink_data")) {
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  }

  long_df <- olink_data$data
  smeta   <- olink_data$sample_meta
  ameta   <- olink_data$assay_meta
  sid_col <- olink_data$params$sample_col

  title_pfx <- if (!is.null(title)) paste0(title, " \u2014 ") else ""

  plots <- list()

  # ---------------------------------------------------------------------------
  # 1. Per-plate NPX distributions (boxplot)
  # ---------------------------------------------------------------------------

  # Determine grouping column for x-axis.
  # Checks $data first (legacy), then $sample_meta; joins into long_df if needed.
  if (plate_col %in% colnames(long_df)) {
    x_col <- plate_col
    x_lab <- plate_col
  } else if (plate_col %in% colnames(smeta)) {
    x_col   <- plate_col
    x_lab   <- plate_col
    pjoin   <- smeta[, c(sid_col, plate_col), drop = FALSE]
    long_df <- merge(long_df, pjoin, by = sid_col, all.x = TRUE, sort = FALSE)
  } else if ("PlateID" %in% colnames(smeta)) {
    x_col <- "PlateID"
    x_lab <- "PlateID"
    if (plate_col != "PlateID") {
      warning("plate_col '", plate_col,
              "' not found; falling back to 'PlateID'.")
    }
    pjoin   <- smeta[, c(sid_col, "PlateID"), drop = FALSE]
    long_df <- merge(long_df, pjoin, by = sid_col, all.x = TRUE, sort = FALSE)
  } else {
    x_col <- sid_col
    x_lab <- "Sample"
    warning("No plate column found; showing per-sample distributions.")
  }

  p_dist <- ggplot2::ggplot(long_df,
                             ggplot2::aes(x = .data[[x_col]],
                                          y = .data[["NPX"]])) +
    ggplot2::geom_boxplot(outlier.size = 0.4, fill = "grey90") +
    ggplot2::labs(
      title = paste0(title_pfx, "NPX distributions by ", x_lab),
      x     = x_lab,
      y     = "NPX (log2)"
    ) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  plots$npx_distributions <- p_dist

  # ---------------------------------------------------------------------------
  # 2. SampleQC summary table
  # A data.frame is more legible than a stacked bar when most samples pass.
  # Uses $sample_meta (all samples including controls).
  # ---------------------------------------------------------------------------

  if ("SampleQC" %in% colnames(smeta)) {
    if ("SampleType" %in% colnames(smeta)) {
      qc_tbl <- as.data.frame(
        table(SampleType = smeta$SampleType, SampleQC = smeta$SampleQC),
        stringsAsFactors = FALSE
      )
      colnames(qc_tbl)[colnames(qc_tbl) == "Freq"] <- "Count"
    } else {
      qc_tbl <- as.data.frame(
        table(SampleQC = smeta$SampleQC),
        stringsAsFactors = FALSE
      )
      colnames(qc_tbl)[colnames(qc_tbl) == "Freq"] <- "Count"
    }
    plots$sample_qc_table <- qc_tbl
  }

  # ---------------------------------------------------------------------------
  # 3. AssayQC warn fraction histogram (optional)
  # Shows the distribution of warn_fraction across all proteins, with a
  # reference line at warn_threshold. More diagnostic than a ranked barplot.
  # ---------------------------------------------------------------------------

  if (plot_warn_proteins &&
      !is.null(ameta) &&
      "warn_fraction" %in% colnames(ameta)) {

    warn_df <- ameta[!is.na(ameta$warn_fraction), ]

    n_above <- sum(warn_df$warn_fraction >= warn_threshold)

    p_warn <- ggplot2::ggplot(warn_df,
                               ggplot2::aes(x = .data[["warn_fraction"]])) +
      ggplot2::geom_histogram(binwidth = 0.05,
                              boundary = 0,
                              fill     = "#e08214",
                              colour   = "white") +
      ggplot2::geom_vline(xintercept = warn_threshold,
                          linetype   = "dashed",
                          colour     = "grey30") +
      ggplot2::annotate("text",
                        x     = warn_threshold + 0.01,
                        y     = Inf,
                        label = paste0(n_above, " proteins \u2265 threshold"),
                        hjust = 0,
                        vjust = 1.5,
                        size  = font_size / 3) +
      ggplot2::scale_x_continuous(
        labels = function(x) paste0(round(x * 100), "%"),
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        title = paste0(title_pfx, "AssayQC warn fraction distribution"),
        x     = "Fraction of samples with WARN",
        y     = "Number of proteins"
      ) +
      ggplot2::theme_bw(base_size = font_size)

    plots$warn_proteins <- p_warn
  }

  return(plots)
}
