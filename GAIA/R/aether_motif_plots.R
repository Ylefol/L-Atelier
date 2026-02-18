# GAIA/Aether/motif_plots.R
# Visualization functions for HOMER motif enrichment results
#
# Works with homer_motif and homer_motif_batch objects from
# APOLLO_homer_motif_enrichment() or APOLLO_load_homer_results().



#' Dotplot of HOMER Motif Enrichment Results
#'
#' @description Creates a dotplot showing enriched transcription factor motifs.
#' For single peak sets, shows one column of dots. For batch results (multiple
#' peak sets), shows motifs on the y-axis and peak sets on the x-axis.
#'
#' @param homer_result A homer_motif or homer_motif_batch object.
#' @param top_n Integer. Number of top motifs to show per peak set (default = 15).
#' @param color_by Character. Variable for color scale: "q_value" (default,
#'   shown as -log10) or "p_value".
#' @param size_by Character. Variable for dot size: "pct_target" (default)
#'   or "n_target".
#' @param sets Character vector. For batch results, specific peak sets to include.
#'   If NULL, all sets are shown.
#' @param q_thresh Numeric. Only show motifs significant in at least one set
#'   (default = 0.05). Set to 1 to show all.
#' @param title Character. Plot title. If NULL, auto-generated.
#' @param font_size Numeric. Base font size (default = 9).
#' @param max_name_length Integer. Maximum characters for motif names (default = 40).
#' @param color_low Character. Color for low significance (default = "#2166ac" blue).
#' @param color_high Character. Color for high significance (default = "#b2182b" red).
#' @param dot_range Numeric vector of length 2. Min and max dot sizes (default = c(2, 8)).
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_motif_enrichment <- function(homer_result,
                                          top_n = 15,
                                          color_by = "q_value",
                                          size_by = "pct_target",
                                          sets = NULL,
                                          q_thresh = 0.05,
                                          title = NULL,
                                          font_size = 9,
                                          max_name_length = 40,
                                          color_low = "#2166ac",
                                          color_high = "#b2182b",
                                          dot_range = c(2, 8)) {

  # --- Extract data ---
  if (inherits(homer_result, "homer_motif_batch")) {
    df <- homer_result$combined
    is_batch <- TRUE
  } else if (inherits(homer_result, "homer_motif")) {
    df <- homer_result$known
    if (!is.null(df)) df$peak_set <- "peaks"
    is_batch <- FALSE
  } else {
    stop("homer_result must be a homer_motif or homer_motif_batch object.",
         call. = FALSE)
  }

  if (is.null(df) || nrow(df) == 0) {
    warning("No motif results to plot")
    return(ggplot() + theme_void() +
             labs(title = title %||% "Motif Enrichment",
                  subtitle = "No motif results available"))
  }

  # --- Filter by sets ---
  if (!is.null(sets) && is_batch) {
    df <- df[df$peak_set %in% sets, ]
  }

  # --- Filter by significance ---
  # Keep motifs significant in at least one set
  sig_motifs <- unique(df$motif_name[!is.na(df$q_value) & df$q_value <= q_thresh])
  if (length(sig_motifs) == 0) {
    warning("No significant motifs at q < ", q_thresh)
    return(ggplot() + theme_void() +
             labs(title = title %||% "Motif Enrichment",
                  subtitle = paste("No significant motifs at q <", q_thresh)))
  }
  df <- df[df$motif_name %in% sig_motifs, ]

  # --- Select top N per set ---
  df <- do.call(rbind, lapply(split(df, df$peak_set), function(x) {
    x <- x[order(x$q_value), ]
    head(x, top_n)
  }))
  rownames(df) <- NULL

  # --- Truncate motif family names for display ---
  df$motif_label <- ifelse(
    nchar(df$motif_family) > max_name_length,
    paste0(substr(df$motif_family, 1, max_name_length - 3), "..."),
    df$motif_family
  )

  # --- Order motifs by best q-value ---
  motif_stats <- aggregate(q_value ~ motif_label, data = df,
                            FUN = function(x) min(x, na.rm = TRUE))
  motif_order <- motif_stats$motif_label[order(motif_stats$q_value, decreasing = TRUE)]
  df$motif_label <- factor(df$motif_label, levels = motif_order)

  # --- Set up color variable ---
  if (color_by == "p_value") {
    df$color_var <- -log10(df$p_value)
    color_label <- "-log10(p-value)"
  } else {
    df$color_var <- -log10(df$q_value)
    color_label <- "-log10(q-value)"
  }

  # --- Set up size variable ---
  if (size_by == "n_target") {
    df$size_var <- df$n_target
    size_label <- "Target Count"
  } else {
    df$size_var <- df$pct_target
    size_label <- "% Targets"
  }

  # --- Title ---
  if (is.null(title)) {
    if (is_batch) {
      title <- "Motif Enrichment Comparison"
    } else {
      title <- "HOMER Known Motif Enrichment"
    }
  }

  # --- Build plot ---
  if (is_batch) {
    # Multi-set: motifs on y-axis, peak sets on x-axis
    p <- ggplot(df, aes(x = peak_set, y = motif_label)) +
      geom_point(aes(size = size_var, color = color_var)) +
      labs(x = NULL)
  } else {
    # Single set: motifs on y-axis, significance on x-axis
    p <- ggplot(df, aes(x = color_var, y = motif_label)) +
      geom_point(aes(size = size_var, color = color_var)) +
      labs(x = color_label)
  }

  p <- p +
    scale_color_gradient(low = color_low, high = color_high, name = color_label) +
    scale_size_continuous(range = dot_range, name = size_label) +
    labs(
      title = title,
      subtitle = paste("Top", top_n, "motifs per set (q <", q_thresh, ")"),
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = font_size + 4),
      plot.subtitle = element_text(size = font_size + 1, color = "gray40"),
      axis.text.y = element_text(size = font_size),
      axis.text.x = element_text(size = font_size + 1, angle = 45, hjust = 1,
                                  face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}


#' Heatmap of Motif Enrichment Across Peak Sets
#'
#' @description Creates a heatmap showing -log10(q-value) for motifs across
#' multiple peak sets. Useful for comparing which TF motifs are enriched in
#' different conditions or clusters.
#'
#' @param homer_batch A homer_motif_batch object.
#' @param top_n Integer. Number of top motifs to include (default = 20).
#'   Selected by best q-value across all sets.
#' @param motifs Character vector. Specific motif family names to include.
#'   Overrides top_n if provided.
#' @param cluster_sets Logical. Hierarchically cluster columns (peak sets).
#'   Default FALSE (preserves input order).
#' @param cluster_motifs Logical. Hierarchically cluster rows (motifs).
#'   Default TRUE.
#' @param q_thresh Numeric. Significance threshold for motif selection
#'   (default = 0.05). Motifs must be significant in at least one set.
#' @param colors Character vector of length 3. Low, mid, high colors for
#'   the heatmap. Default: white to red.
#' @param title Character. Plot title (default = "Motif Enrichment Comparison").
#' @param show_values Logical. Show -log10(q) values in cells. Default FALSE.
#' @param ...
#'   Additional arguments passed to pheatmap.
#'
#' @return A pheatmap object (invisibly).
#'
#' @export
AETHER_plot_motif_comparison <- function(homer_batch,
                                          top_n = 20,
                                          motifs = NULL,
                                          cluster_sets = FALSE,
                                          cluster_motifs = TRUE,
                                          q_thresh = 0.05,
                                          colors = NULL,
                                          title = "Motif Enrichment Comparison",
                                          show_values = FALSE,
                                          ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')",
         call. = FALSE)
  }

  if (!inherits(homer_batch, "homer_motif_batch")) {
    stop("homer_batch must be a homer_motif_batch object.", call. = FALSE)
  }

  df <- homer_batch$combined
  if (is.null(df) || nrow(df) == 0) {
    stop("No combined motif results available.", call. = FALSE)
  }

  # --- Select motifs ---
  if (!is.null(motifs)) {
    # Use specified motifs
    df <- df[df$motif_family %in% motifs, ]
  } else {
    # Select top N by best q-value across sets
    # First filter to significant in at least one set
    sig_motifs <- unique(df$motif_family[!is.na(df$q_value) & df$q_value <= q_thresh])
    df <- df[df$motif_family %in% sig_motifs, ]

    if (nrow(df) > 0) {
      best_q <- aggregate(q_value ~ motif_family, data = df,
                           FUN = function(x) min(x, na.rm = TRUE))
      best_q <- best_q[order(best_q$q_value), ]
      selected <- head(best_q$motif_family, top_n)
      df <- df[df$motif_family %in% selected, ]
    }
  }

  if (nrow(df) == 0) {
    stop("No motifs remaining after filtering.", call. = FALSE)
  }

  # --- Build matrix: motifs (rows) x peak_sets (columns) ---
  # Value: -log10(q-value), capped for visualization
  set_names <- unique(df$peak_set)
  motif_names <- unique(df$motif_family)

  mat <- matrix(0, nrow = length(motif_names), ncol = length(set_names),
                dimnames = list(motif_names, set_names))

  for (i in seq_len(nrow(df))) {
    row_name <- df$motif_family[i]
    col_name <- df$peak_set[i]
    q_val <- df$q_value[i]
    if (!is.na(q_val) && q_val > 0) {
      mat[row_name, col_name] <- -log10(q_val)
    }
  }

  # Cap extreme values for visualization
  max_val <- max(mat, na.rm = TRUE)
  cap <- min(max_val, 50)  # Cap at -log10(1e-50)
  mat[mat > cap] <- cap

  # --- Color palette ---
  if (is.null(colors)) {
    colors <- c("white", "#fee0d2", "#b2182b")
  }

  color_fun <- colorRampPalette(colors)(100)

  # --- Display values ---
  display_numbers <- show_values
  number_format <- "%.1f"

  # --- Generate heatmap ---
  p <- pheatmap::pheatmap(
    mat,
    main             = title,
    color            = color_fun,
    cluster_rows     = cluster_motifs,
    cluster_cols     = cluster_sets,
    display_numbers  = display_numbers,
    number_format    = number_format,
    fontsize         = 10,
    fontsize_row     = 9,
    fontsize_col     = 10,
    angle_col        = 45,
    border_color     = "gray90",
    legend_labels    = "-log10(q-value)",
    ...
  )

  invisible(p)
}
