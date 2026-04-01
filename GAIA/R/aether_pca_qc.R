#' PC-Metadata Association Heatmap
#'
#' @description Heatmap of -log10(p-values) from a PC-metadata association
#' test. Rows are principal components, columns are metadata variables. The
#' colour intensity indicates strength of association; a dashed line marks
#' the significance threshold.
#'
#' @param assoc_result An \code{artemis_pc_assoc} object from
#'   \code{ARTEMIS_pc_metadata_association()}.
#' @param p_threshold Numeric. Significance threshold for the cell border
#'   annotation. Default: 0.05.
#' @param top_n Integer or NULL. If set, retains only the \code{top_n} most
#'   significant variables (lowest p-value) per PC, then shows the union across
#'   all PCs. Useful when many variables are significant. NULL = show all
#'   tested variables. Default: NULL.
#' @param order_vars Logical. Order metadata variables (columns) by mean
#'   -log10(p) descending so the most associated variables appear first.
#'   Default: TRUE.
#' @param max_log10p Numeric or NULL. Cap the colour scale at this value.
#'   NULL = use the data maximum. Default: NULL.
#' @param color_low Character. Colour for low -log10(p) (no association).
#'   Default: "white".
#' @param color_high Character. Colour for high -log10(p) (strong association).
#'   Default: "#b2182b" (red).
#' @param show_values Logical. Overlay -log10(p) values as text in each cell.
#'   Default: FALSE.
#' @param font_size Numeric. Base font size. Default: 9.
#' @param title Character or NULL. Plot title. Default: "PC-Metadata Association".
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' assoc <- ARTEMIS_pc_metadata_association(ol$wide, ol$sample_meta,
#'                                           sample_col = "SampleID")
#' p <- AETHER_plot_pc_association(assoc)
#' print(p)
#' }
AETHER_plot_pc_association <- function(assoc_result,
                                        p_threshold = 0.05,
                                        top_n       = NULL,
                                        order_vars  = TRUE,
                                        max_log10p  = NULL,
                                        color_low   = "white",
                                        color_high  = "#b2182b",
                                        show_values = FALSE,
                                        font_size   = 9,
                                        title       = "PC-Metadata Association") {

  if (!inherits(assoc_result, "artemis_pc_assoc")) {
    stop("assoc_result must be an artemis_pc_assoc object from ",
         "ARTEMIS_pc_metadata_association().")
  }

  pmat <- assoc_result$pvalues

  # Remove skipped (all-NA) columns
  keep_cols <- colSums(!is.na(pmat)) > 0
  pmat      <- pmat[, keep_cols, drop = FALSE]

  if (ncol(pmat) == 0)
    stop("No testable variables in assoc_result.")

  # Optionally restrict to the top N variables per PC (union across all PCs)
  if (!is.null(top_n)) {
    top_n    <- as.integer(top_n)
    selected <- unique(unlist(lapply(seq_len(nrow(pmat)), function(i) {
      p_row <- pmat[i, ]
      valid <- which(!is.na(p_row))
      if (length(valid) == 0L) return(character(0L))
      top_idx <- valid[order(p_row[valid])][seq_len(min(top_n, length(valid)))]
      colnames(pmat)[top_idx]
    })))
    pmat <- pmat[, selected, drop = FALSE]
    cat("top_n = ", top_n, ": showing ", ncol(pmat),
            " unique variables across ", nrow(pmat), " PCs")
  }

  # Append top_n to title if set
  if (!is.null(top_n))
    title <- paste0(title, " \u2014 top ", top_n, " per PC")

  # -log10 transform (NA stays NA)
  log10p_mat        <- -log10(pmat)
  log10p_mat[is.na(log10p_mat)] <- 0

  # Optional: order columns by mean association strength
  if (order_vars) {
    col_order <- order(colMeans(log10p_mat, na.rm = TRUE), decreasing = TRUE)
    log10p_mat <- log10p_mat[, col_order, drop = FALSE]
  }

  # Significance threshold line value
  thresh_line <- -log10(p_threshold)

  # Cap colour scale
  cap_val <- if (!is.null(max_log10p)) max_log10p else max(log10p_mat, na.rm = TRUE)
  cap_val <- max(cap_val, thresh_line + 0.1)  # always show threshold

  # Build long data.frame for ggplot
  pc_names  <- rownames(log10p_mat)
  var_names <- colnames(log10p_mat)

  # Factor levels: PCs in order (PC1 at top), variables in association order
  pc_factor  <- factor(rep(pc_names,  times = ncol(log10p_mat)),
                        levels = rev(pc_names))   # rev so PC1 is at top
  var_factor <- factor(rep(var_names, each  = nrow(log10p_mat)),
                        levels = var_names)

  # Variance explained labels for y-axis
  var_exp_pct <- round(100 * assoc_result$var_explained[
    match(pc_names, paste0("PC", seq_along(assoc_result$var_explained)))
  ], 1)
  pc_labels <- paste0(pc_names, " (", var_exp_pct, "%)")
  names(pc_labels) <- pc_names

  plot_df <- data.frame(
    PC       = pc_factor,
    Variable = var_factor,
    log10p   = pmin(as.numeric(log10p_mat), cap_val),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data[["Variable"]],
                                     y = .data[["PC"]],
                                     fill = .data[["log10p"]])) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(
      low      = color_low,
      high     = color_high,
      limits   = c(0, cap_val),
      name     = expression(-log[10](p))
    ) +
    ggplot2::geom_hline(yintercept = seq_len(length(pc_names)) - 0.5,
                         colour = "grey85", linewidth = 0.2) +
    ggplot2::scale_y_discrete(labels = pc_labels) +
    ggplot2::labs(
      title = title,
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "grey60")
    )

  # Annotate cells that cross the significance threshold
  sig_mask <- as.numeric(log10p_mat) >= thresh_line
  if (any(sig_mask, na.rm = TRUE)) {
    sig_df <- plot_df[sig_mask & !is.na(sig_mask), ]
    p <- p + ggplot2::geom_tile(
      data    = sig_df,
      mapping = ggplot2::aes(x = .data[["Variable"]],
                             y = .data[["PC"]]),
      fill    = NA,
      colour  = "black",
      linewidth = 0.6
    )
  }

  # Optional: overlay numeric values
  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(.data[["log10p"]] > 0,
                                  round(.data[["log10p"]], 1), "")),
      size = font_size * 0.25
    )
  }

  return(p)
}
