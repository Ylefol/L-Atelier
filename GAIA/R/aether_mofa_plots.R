# GAIA/Aether/mofa_plots.R
# Visualization functions for MOFA2 multi-omics factor analysis results
#
# Works with hephaestus_mofa objects from HEPHAESTUS_run_mofa() /
# HEPHAESTUS_load_mofa_model().


# ==============================================================================
# VARIANCE EXPLAINED HEATMAP
# ==============================================================================

#' Plot variance explained per view and factor
#'
#' Creates a heatmap showing how much variance each MOFA factor explains in
#' each omics view.
#'
#' @param mofa_result A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
#'   HEPHAESTUS_load_mofa_model().
#' @param show_values Logical. Show variance explained (%) as cell text.
#'   Default: TRUE.
#' @param colors Two-color low/high gradient. Default: white to dark blue.
#' @param title Plot title. Default: "Variance Explained per View and Factor".
#'
#' @return A ggplot object. Faceted by group if the model has more than one.
#'
#' @examples
#' \dontrun{
#' result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
#' AETHER_plot_mofa_variance_explained(result)
#'
#' }
#' @export
AETHER_plot_mofa_variance_explained <- function(mofa_result,
                                                  show_values = TRUE,
                                                  colors = NULL,
                                                  title = "Variance Explained per View and Factor") {

  if (!inherits(mofa_result, "hephaestus_mofa")) {
    stop("mofa_result must be a hephaestus_mofa object")
  }

  ve_df <- mofa_result$variance_explained

  # Order factors numerically (Factor1, Factor2, ..., Factor10) rather than
  # alphabetically (Factor1, Factor10, Factor2, ...)
  factor_levels <- unique(ve_df$factor)
  factor_levels <- factor_levels[order(as.integer(gsub("[^0-9]", "", factor_levels)))]
  ve_df$factor <- factor(ve_df$factor, levels = factor_levels)

  if (is.null(colors)) {
    colors <- c("white", "#08519C")
  }

  p <- ggplot(ve_df, aes(x = .data[["factor"]], y = .data[["view"]],
                         fill = .data[["variance_explained"]])) +
    geom_tile(color = "grey85") +
    scale_fill_gradientn(colors = colors, name = "Variance\nexplained (%)") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )

  if (show_values) {
    p <- p + geom_text(aes(label = sprintf("%.1f", .data[["variance_explained"]])), size = 3)
  }

  if (length(unique(ve_df$group)) > 1) {
    p <- p + facet_wrap(~group)
  }

  p
}


# ==============================================================================
# FACTOR SCATTER
# ==============================================================================

#' Plot MOFA factor scores for two factors
#'
#' Creates a scatter plot of sample scores on two MOFA factors, optionally
#' colored by an external metadata variable (e.g. clinical group, timepoint).
#'
#' @param mofa_result A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
#'   HEPHAESTUS_load_mofa_model().
#' @param factors Length-2 numeric vector, which factor indices to plot on
#'   x/y. Default: c(1, 2).
#' @param color_by Optional named vector/factor, sample ID -> value, used to
#'   color points (e.g. a clinical metadata column keyed by sample ID). Not
#'   part of the hephaestus_mofa object itself -- supply it separately.
#' @param color_label Legend title when color_by is supplied. Default: "Group".
#' @param point_size Point size. Default: 2.
#' @param title Plot title. If NULL, auto-generated from the factor names.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
#' AETHER_plot_mofa_factors(result, factors = c(1, 2))
#'
#' # Colored by clinical group (named vector keyed by sample ID)
#' AETHER_plot_mofa_factors(result, factors = c(1, 2),
#'                           color_by = clinical_group, color_label = "Group")
#'
#' }
#' @export
AETHER_plot_mofa_factors <- function(mofa_result,
                                      factors = c(1, 2),
                                      color_by = NULL,
                                      color_label = "Group",
                                      point_size = 2,
                                      title = NULL) {

  if (!inherits(mofa_result, "hephaestus_mofa")) {
    stop("mofa_result must be a hephaestus_mofa object")
  }
  if (length(factors) != 2) {
    stop("factors must be a length-2 vector of factor indices to plot")
  }

  factor_names <- paste0("Factor", factors)
  fac_df <- mofa_result$factors
  missing_factors <- setdiff(factor_names, unique(fac_df$factor))
  if (length(missing_factors) > 0) {
    stop("Factor(s) not found in mofa_result: ", paste(missing_factors, collapse = ", "))
  }

  wide_df <- stats::reshape(
    fac_df[fac_df$factor %in% factor_names, c("sample", "factor", "value")],
    idvar = "sample", timevar = "factor", direction = "wide"
  )
  colnames(wide_df) <- gsub("^value\\.", "", colnames(wide_df))

  if (!is.null(color_by)) {
    wide_df$color_value <- color_by[wide_df$sample]
    n_missing <- sum(is.na(wide_df$color_value))
    if (n_missing > 0) {
      warning(n_missing, " sample(s) have no color_by value and will be shown as NA.")
    }
  }

  if (is.null(title)) {
    title <- paste0("MOFA Factors: ", factor_names[1], " vs ", factor_names[2])
  }

  p <- ggplot(wide_df, aes(x = .data[[factor_names[1]]], y = .data[[factor_names[2]]]))

  if (!is.null(color_by)) {
    p <- p + geom_point(aes(color = .data[["color_value"]]), size = point_size, alpha = 0.8) +
      labs(color = color_label)
  } else {
    p <- p + geom_point(size = point_size, alpha = 0.8, color = "#0072B2")
  }

  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    labs(x = factor_names[1], y = factor_names[2], title = title) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  p
}


# ==============================================================================
# TOP FEATURE WEIGHTS
# ==============================================================================

#' Plot top-loading features for one view/factor
#'
#' Creates a bar plot of the highest-magnitude feature weights for one view
#' on one MOFA factor.
#'
#' @param mofa_result A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
#'   HEPHAESTUS_load_mofa_model().
#' @param view Character. View name (must match a name in the X list passed
#'   to HEPHAESTUS_run_mofa()).
#' @param factor Numeric factor index (e.g. 1) or character factor name
#'   (e.g. "Factor1").
#' @param top_n Number of top-loading features to show, ranked by absolute
#'   weight. Default: 20.
#' @param title Plot title. If NULL, auto-generated.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
#' AETHER_plot_mofa_weights(result, view = "rna", factor = 1)
#'
#' }
#' @export
AETHER_plot_mofa_weights <- function(mofa_result,
                                      view,
                                      factor,
                                      top_n = 20,
                                      title = NULL) {

  if (!inherits(mofa_result, "hephaestus_mofa")) {
    stop("mofa_result must be a hephaestus_mofa object")
  }

  factor_name <- if (is.numeric(factor)) paste0("Factor", factor) else factor

  w_df <- mofa_result$weights
  if (!view %in% unique(w_df$view)) {
    stop("view '", view, "' not found in mofa_result. Available views: ",
         paste(unique(w_df$view), collapse = ", "))
  }
  if (!factor_name %in% unique(w_df$factor)) {
    stop("factor '", factor_name, "' not found in mofa_result. Available factors: ",
         paste(unique(w_df$factor), collapse = ", "))
  }

  sub_df <- w_df[w_df$view == view & w_df$factor == factor_name, ]
  sub_df <- sub_df[order(-abs(sub_df$value)), ]
  sub_df <- head(sub_df, top_n)
  # base::factor() explicitly -- the `factor` parameter above shadows it
  sub_df$feature <- base::factor(sub_df$feature, levels = rev(sub_df$feature))
  sub_df$direction <- ifelse(sub_df$value >= 0, "Positive", "Negative")

  if (is.null(title)) {
    title <- paste0(view, " -- top ", nrow(sub_df), " weights on ", factor_name)
  }

  p <- ggplot(sub_df, aes(x = .data[["feature"]], y = .data[["value"]], fill = .data[["direction"]])) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Positive" = "#D55E00", "Negative" = "#0072B2"), name = NULL) +
    labs(x = NULL, y = "Weight", title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )

  p
}
