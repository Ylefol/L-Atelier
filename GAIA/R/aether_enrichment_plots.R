###############################################################################
########### Enrichment Result Visualizations ###########
###############################################################################


#' Multi-Module Dotplot for gprofiler2 Enrichment Results
#'
#' @description Creates a dotplot showing top enriched terms across multiple
#' modules/gene lists from gprofiler2 results. Modules are shown on the X-axis,
#' terms on the Y-axis, with dot size representing odds ratio (default) or
#' gene count, and color representing significance.
#'
#' @param gost_result A gost_enrichment object from APOLLO_enrich_gost(), or
#'   a data.frame with the required columns (typically the $combined element).
#' @param source Character. Which source to plot (e.g., "GO:BP", "KEGG", "REAC").
#'   Required - only one source per plot for clarity.
#' @param top_n Integer. Number of top terms to show per module (default = 10).
#' @param modules Character vector. Specific modules to include. If NULL (default),
#'   includes all modules.
#' @param color_by Character. Variable for color scale: "p_value" (default),
#'   "odds_ratio", "precision", or "recall".
#' @param size_by Character. Variable for dot size: "odds_ratio" (default),
#'   "intersection_size", "precision", or "recall".
#' @param order_by Character. How to order terms on Y-axis: "p_value" (default),
#'   "odds_ratio", "intersection_size", or "term_name".
#' @param title Character. Plot title. If NULL, auto-generated from source.
#' @param font_size Numeric. Base font size (default = 9).
#' @param max_term_length Integer. Maximum characters for term names (default = 50).
#' @param color_low Character. Color for low values (default = "#2166ac" blue).
#' @param color_high Character. Color for high values (default = "#b2182b" red).
#' @param dot_range Numeric vector of length 2. Min and max dot sizes (default = c(2, 8)).
#'
#' @return A ggplot object.
#'
#' @details
#' The plot shows:
#' - X-axis: Modules (gene lists)
#' - Y-axis: Enriched terms (top N per module, union across modules)
#' - Dot size: Odds ratio (default) - strength of association
#' - Dot color: -log10(p-value) - significance
#'
#' Odds ratio is calculated as: (k/n) / (M/N) simplified, or more precisely
#' using the 2x2 contingency table where OR > 1 indicates enrichment.
#'
#' Terms are selected as top N per module, then displayed as a union. A term
#' may appear significant in one module but not another - dots only appear
#' where the term was significant in that module.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot GO:BP results for all modules
#' p <- AETHER_plot_gost_dotplot(gost_results, source = "GO:BP")
#'
#' # KEGG pathways, top 5 per module, color by odds ratio
#' p <- AETHER_plot_gost_dotplot(gost_results, source = "KEGG", top_n = 5,
#'                                color_by = "odds_ratio")
#'
#' # Specific modules only
#' p <- AETHER_plot_gost_dotplot(gost_results, source = "GO:BP",
#'                                modules = c("blue", "turquoise", "brown"))
#'
#' }
AETHER_plot_gost_dotplot <- function(gost_result,
                                      source,
                                      top_n = 10,
                                      modules = NULL,
                                      color_by = "p_value",
                                      size_by = "odds_ratio",
                                      order_by = "p_value",
                                      title = NULL,
                                      font_size = 9,
                                      max_term_length = 50,
                                      color_low = "#2166ac",
                                      color_high = "#b2182b",
                                      dot_range = c(2, 8)) {

  # ---------------------------------------------------------------------------
  # Extract data
  # ---------------------------------------------------------------------------
  if (inherits(gost_result, "gost_enrichment")) {
    df <- gost_result$combined
  } else if (is.data.frame(gost_result)) {
    df <- gost_result
  } else {
    stop("gost_result must be a gost_enrichment object or data.frame")
  }

  if (nrow(df) == 0) {
    warning("No enrichment results to plot")
    return(ggplot() + theme_void() +
             labs(title = title %||% source, subtitle = "No significant terms found"))
  }

  # ---------------------------------------------------------------------------
  # Filter by source
  # ---------------------------------------------------------------------------
  if (missing(source) || is.null(source)) {
    stop("source parameter is required. Choose one of: ",
         paste(unique(df$source), collapse = ", "))
  }

  df <- df[df$source == source, ]

  if (nrow(df) == 0) {
    warning("No results for source: ", source)
    return(ggplot() + theme_void() +
             labs(title = title %||% source,
                  subtitle = paste("No significant terms for", source)))
  }

  # ---------------------------------------------------------------------------
  # Filter by modules
  # ---------------------------------------------------------------------------
  if (!is.null(modules)) {
    df <- df[df$module %in% modules, ]
    if (nrow(df) == 0) {
      stop("No results for specified modules")
    }
  }

  # ---------------------------------------------------------------------------
  # Check required columns
  # ---------------------------------------------------------------------------
  required_cols <- c("term_name", "p_value", "intersection_size", "module")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Calculate odds ratio
  # ---------------------------------------------------------------------------
  # OR = (k/n) / ((M-k)/(N-n)) where:
  # k = intersection_size, n = query_size, M = term_size, N = effective_domain_size
  if (all(c("query_size", "term_size", "effective_domain_size") %in% colnames(df))) {
    k <- df$intersection_size
    n <- df$query_size
    M <- df$term_size
    N <- df$effective_domain_size

    # 2x2 table: a=k, b=n-k, c=M-k, d=N-M-n+k
    a <- k
    b <- n - k
    c <- M - k
    d <- N - M - n + k

    # Add small constant to avoid division by zero
    df$odds_ratio <- (a * d) / (b * c + 1e-10)
  } else {
    # Fallback: use precision/recall ratio as proxy
    if (all(c("precision", "recall") %in% colnames(df))) {
      df$odds_ratio <- df$precision / (df$recall + 1e-10)
    } else {
      df$odds_ratio <- df$intersection_size
      warning("Could not calculate odds ratio - using intersection_size instead")
    }
  }

  # Calculate precision and recall if not present
  if (!"precision" %in% colnames(df) && "query_size" %in% colnames(df)) {
    df$precision <- df$intersection_size / df$query_size
  }
  if (!"recall" %in% colnames(df) && "term_size" %in% colnames(df)) {
    df$recall <- df$intersection_size / df$term_size
  }

  # ---------------------------------------------------------------------------
  # Select top N per module
  # ---------------------------------------------------------------------------
  df <- df[order(df$p_value), ]
  df <- do.call(rbind, lapply(split(df, df$module), function(x) {
    head(x, top_n)
  }))
  rownames(df) <- NULL

  if (nrow(df) == 0) {
    warning("No terms remaining after filtering")
    return(ggplot() + theme_void() + labs(title = title %||% source))
  }

  # ---------------------------------------------------------------------------
  # Wrap long term names across multiple lines
  # ---------------------------------------------------------------------------
  df$term_short <- sapply(df$term_name, function(name) {
    paste(strwrap(name, width = max_term_length), collapse = "\n")
  }, USE.NAMES = FALSE)

  # ---------------------------------------------------------------------------
  # Create module labels with annotation coverage percentage
  # ---------------------------------------------------------------------------
  # query_size = genes recognized in this specific database/source
  # module_sizes (from metadata) = total genes in the original gene list
  # Percentage shows how much of the module is annotated in the source
  total_sizes <- NULL
  if (inherits(gost_result, "gost_enrichment") &&
      !is.null(gost_result$metadata$module_sizes)) {
    total_sizes <- gost_result$metadata$module_sizes
  }

  if ("query_size" %in% colnames(df) && !is.null(total_sizes)) {
    query_sizes <- aggregate(query_size ~ module, data = df, FUN = function(x) x[1])
    module_labels <- sapply(query_sizes$module, function(mod) {
      total <- total_sizes[[mod]]
      qs <- query_sizes$query_size[query_sizes$module == mod]
      if (!is.null(total) && total > 0) {
        pct <- round(qs / total * 100, 1)
        paste0(mod, " (", pct, "%)")
      } else {
        mod
      }
    })
    module_labels <- setNames(module_labels, query_sizes$module)
    df$module_label <- module_labels[as.character(df$module)]
    df$module_label <- factor(df$module_label, levels = module_labels[unique(df$module)])
  } else {
    df$module_label <- df$module
  }

  # ---------------------------------------------------------------------------
  # Order terms for Y-axis
  # ---------------------------------------------------------------------------
  # Get unique terms and their best (minimum) p-value across modules
  term_stats <- aggregate(
    cbind(p_value, odds_ratio, intersection_size) ~ term_short,
    data = df,
    FUN = function(x) if (is.numeric(x)) min(x) else x[which.min(x)]
  )

  if (order_by == "p_value") {
    term_order <- term_stats$term_short[order(term_stats$p_value, decreasing = TRUE)]
  } else if (order_by == "odds_ratio") {
    term_order <- term_stats$term_short[order(term_stats$odds_ratio)]
  } else if (order_by == "intersection_size") {
    term_order <- term_stats$term_short[order(term_stats$intersection_size)]
  } else if (order_by == "term_name") {
    term_order <- sort(unique(df$term_short), decreasing = TRUE)
  } else {
    term_order <- term_stats$term_short[order(term_stats$p_value, decreasing = TRUE)]
  }

  df$term_short <- factor(df$term_short, levels = term_order)

  # ---------------------------------------------------------------------------
  # Set up size variable
  # ---------------------------------------------------------------------------
  if (size_by == "intersection_size") {
    df$size_var <- df$intersection_size
    size_label <- "Gene Count"
  } else if (size_by == "precision" && "precision" %in% colnames(df)) {
    df$size_var <- df$precision
    size_label <- "Precision"
  } else if (size_by == "recall" && "recall" %in% colnames(df)) {
    df$size_var <- df$recall
    size_label <- "Recall"
  } else {
    # Default: odds_ratio
    df$size_var <- df$odds_ratio
    size_label <- "Odds Ratio"
  }

  # ---------------------------------------------------------------------------
  # Set up color variable
  # ---------------------------------------------------------------------------
  if (color_by == "odds_ratio") {
    df$color_var <- log2(df$odds_ratio)
    color_label <- "log2(Odds Ratio)"
    color_scale <- scale_color_gradient(low = color_low, high = color_high,
                                         name = color_label)
  } else if (color_by == "precision" && "precision" %in% colnames(df)) {
    df$color_var <- df$precision
    color_label <- "Precision"
    color_scale <- scale_color_gradient(low = color_low, high = color_high,
                                         name = color_label)
  } else if (color_by == "recall" && "recall" %in% colnames(df)) {
    df$color_var <- df$recall
    color_label <- "Recall"
    color_scale <- scale_color_gradient(low = color_low, high = color_high,
                                         name = color_label)
  } else {
    # Default: p_value with -log10 transform
    df$color_var <- -log10(df$p_value)
    color_label <- "-log10(p-value)"
    color_scale <- scale_color_gradient(low = color_low, high = color_high,
                                         name = color_label)
  }

  # ---------------------------------------------------------------------------
  # Build title
  # ---------------------------------------------------------------------------
  if (is.null(title)) {
    source_names <- c(
      "GO:BP" = "GO Biological Process",
      "GO:MF" = "GO Molecular Function",
      "GO:CC" = "GO Cellular Component",
      "KEGG" = "KEGG Pathways",
      "REAC" = "Reactome Pathways",
      "WP" = "WikiPathways"
    )
    title <- source_names[source]
    if (is.na(title)) title <- source
  }

  # ---------------------------------------------------------------------------
  # Build plot - modules on X-axis, terms on Y-axis
  # ---------------------------------------------------------------------------
  p <- ggplot(df, aes(x = module_label, y = term_short)) +
    geom_point(aes(size = size_var, color = color_var)) +
    color_scale +
    scale_size_continuous(range = dot_range, name = size_label) +
    labs(
      title = title,
      subtitle = paste("Top", top_n, "terms per module"),
      x = NULL,
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


#' Plot Full (Unfiltered) g:GOSt Manhattan Plot
#'
#' @description Wraps \code{gprofiler2::gostplot()} to render the traditional,
#' interactive g:GOSt Manhattan-style plot for a single gene list/module,
#' showing every evaluated term — not just those passing the significance
#' threshold. Complements the significance-filtered results from
#' \code{APOLLO_enrich_gost()} so a cluster's biology can still be inspected
#' visually when no term clears multiple-testing correction.
#'
#' @param gost_result A single raw gost result object — one element of
#'   \code{enrich_result$results} (e.g. \code{enrich_result$results[["C1"]]})
#'   as returned by \code{APOLLO_enrich_gost()}. This already holds the full,
#'   unfiltered term set.
#' @param capped Logical. Cap the -log10(p-value) y-axis at 16 for readability
#'   (gprofiler2 default behavior). Default: TRUE.
#' @param save_html Character or NULL. Path to save as a standalone
#'   interactive HTML file. Default: NULL.
#'
#' @return A plotly htmlwidget object.
#'
#' @details
#' Interactive only — gprofiler2's Manhattan plot is plotly-native and has no
#' meaningful static (ggplot2) equivalent; hovering each point reveals the
#' term name, source, and p-value, which is what makes the non-significant
#' terms interpretable.
#'
#' @examples
#' \dontrun{
#' enrich <- APOLLO_enrich_gost(cluster_genes)
#' AETHER_plot_gost_full(enrich$results[["C1"]], save_html = "C1_gostplot_full.html")
#'
#' }
#' @export
AETHER_plot_gost_full <- function(gost_result, capped = TRUE, save_html = NULL) {

  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required. Install with: install.packages('gprofiler2')")
  }

  if (is.null(gost_result) || is.null(gost_result$result) || nrow(gost_result$result) == 0) {
    stop("'gost_result' has no terms to plot.")
  }

  p <- gprofiler2::gostplot(gost_result, capped = capped, interactive = TRUE)

  if (!is.null(save_html)) {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
      warning("htmlwidgets package needed for HTML export.")
    } else {
      htmlwidgets::saveWidget(p, file = save_html, selfcontained = TRUE)
      cat("[AETHER] Saved full gostplot to:", save_html, "\n")
    }
  }

  return(p)
}


#' Calculate Odds Ratio from Enrichment Results
#'
#' @description Helper function to calculate odds ratio from GeneRatio and
#' BgRatio columns in clusterProfiler enrichment results.
#'
#' @param gene_ratio Character vector. GeneRatio values (e.g., "10/50").
#' @param bg_ratio Character vector. BgRatio values (e.g., "100/10000").
#'
#' @return Numeric vector of odds ratios.
#'
#' @details
#' Odds Ratio = (k/(n-k)) / ((M-k)/(N-M-n+k))
#' Where:
#' - k = genes in query that hit the pathway (from GeneRatio numerator)
#' - n = total genes in query (from GeneRatio denominator)
#' - M = total genes in pathway (from BgRatio numerator)
#' - N = total background genes (from BgRatio denominator)
#'
#' @keywords internal
#'
.calculate_odds_ratio <- function(gene_ratio, bg_ratio) {
  # Parse GeneRatio: "k/n" -> k hits out of n query genes
  gene_parts <- strsplit(gene_ratio, "/")
  k <- as.numeric(sapply(gene_parts, `[`, 1))  # hits

n <- as.numeric(sapply(gene_parts, `[`, 2))  # query size

  # Parse BgRatio: "M/N" -> M pathway genes out of N background
  bg_parts <- strsplit(bg_ratio, "/")
  M <- as.numeric(sapply(bg_parts, `[`, 1))  # pathway size
  N <- as.numeric(sapply(bg_parts, `[`, 2))  # background size

  # 2x2 contingency table:
  #                 In Pathway    Not in Pathway
  # In Query            a              b
  # Not in Query        c              d
  #
  # a = k, b = n - k, c = M - k, d = N - M - n + k

  a <- k
  b <- n - k
  c <- M - k
  d <- N - M - n + k

  # Odds ratio = (a*d) / (b*c)
  # Add small constant to avoid division by zero
  odds_ratio <- (a * d) / (b * c + 1e-10)

  return(odds_ratio)
}


#' Dotplot for Enrichment Results
#'
#' @description Creates a dotplot visualization for GO/KEGG enrichment results.
#' Supports coloring by p-value, adjusted p-value, or odds ratio.
#'
#' @param enrich_result An enrichResult object from clusterProfiler (output of
#'   APOLLO_enrich_go or APOLLO_enrich_kegg), or a data.frame with the required
#'   columns.
#' @param show_category Integer. Number of top categories to display (default = 20).
#' @param color_by Character. Variable for color scale: "p.adjust" (default),
#'   "pvalue", or "odds_ratio".
#' @param size_by Character. Variable for dot size: "Count" (default) or
#'   "GeneRatio".
#' @param order_by Character. Variable to order terms: "p.adjust" (default),
#'   "pvalue", "odds_ratio", or "Count".
#' @param title Character. Plot title (default = "Enrichment Analysis").
#' @param x_label Character. X-axis label. If NULL, auto-generated based on
#'   size_by parameter.
#' @param font_size Numeric. Base font size for term labels (default = 10).
#' @param color_low Character. Color for low values (default = "#2166ac" blue).
#' @param color_high Character. Color for high values (default = "#b2182b" red).
#' @param color_midpoint Numeric or NULL. Midpoint for diverging color scale.
#'   If NULL, uses sequential scale.
#'
#' @return A ggplot object.
#'
#' @details
#' The dotplot shows:
#' - Y-axis: enriched terms (GO terms or KEGG pathways)
#' - X-axis: gene ratio or count
#' - Dot size: number of genes (or gene ratio)
#' - Dot color: significance (p-value/q-value) or effect size (odds ratio)
#'
#' **Why use Odds Ratio for color?**
#' P-values conflate effect size with sample size. A pathway with OR=10 and
#' p=0.01 may be more biologically relevant than OR=1.5 with p=1e-10.
#' Odds ratio directly measures the strength of association:
#' - OR = 1: no association
#' - OR > 1: enrichment (higher = stronger)
#' - OR < 1: depletion
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard dotplot (color by adjusted p-value)
#' p <- AETHER_plot_enrichment_dotplot(go_results)
#'
#' # Color by odds ratio for effect size
#' p <- AETHER_plot_enrichment_dotplot(go_results, color_by = "odds_ratio")
#'
#' # Show more categories, order by odds ratio
#' p <- AETHER_plot_enrichment_dotplot(
#'   go_results,
#'   show_category = 30,
#'   color_by = "odds_ratio",
#'   order_by = "odds_ratio"
#' )
#'
#' }
AETHER_plot_enrichment_dotplot <- function(enrich_result,
                                            show_category = 20,
                                            color_by = "p.adjust",
                                            size_by = "Count",
                                            order_by = "p.adjust",
                                            title = "Enrichment Analysis",
                                            x_label = NULL,
                                            font_size = 10,
                                            color_low = "#2166ac",
                                            color_high = "#b2182b",
                                            color_midpoint = NULL) {

  # ---------------------------------------------------------------------------
  # Extract data from enrichResult object or use data.frame directly
  # ---------------------------------------------------------------------------
  if (inherits(enrich_result, "enrichResult")) {
    df <- as.data.frame(enrich_result)
  } else if (is.data.frame(enrich_result)) {
    df <- enrich_result
  } else {
    stop("enrich_result must be an enrichResult object or data.frame")
  }

  if (nrow(df) == 0) {
    warning("No enrichment results to plot")
    return(ggplot() + theme_void() +
             labs(title = title, subtitle = "No significant terms found"))
  }

  # ---------------------------------------------------------------------------
  # Check required columns
  # ---------------------------------------------------------------------------
  required_cols <- c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Calculate odds ratio
  # ---------------------------------------------------------------------------
  df$odds_ratio <- .calculate_odds_ratio(df$GeneRatio, df$BgRatio)

  # Parse GeneRatio for numeric use
  gene_parts <- strsplit(df$GeneRatio, "/")
  df$GeneRatio_numeric <- as.numeric(sapply(gene_parts, `[`, 1)) /
                          as.numeric(sapply(gene_parts, `[`, 2))

  # ---------------------------------------------------------------------------
  # Validate color_by, size_by, order_by parameters
  # ---------------------------------------------------------------------------
  valid_color <- c("p.adjust", "pvalue", "odds_ratio")
  if (!color_by %in% valid_color) {
    stop("color_by must be one of: ", paste(valid_color, collapse = ", "))
  }

  valid_size <- c("Count", "GeneRatio")
  if (!size_by %in% valid_size) {
    stop("size_by must be one of: ", paste(valid_size, collapse = ", "))
  }

  valid_order <- c("p.adjust", "pvalue", "odds_ratio", "Count")
  if (!order_by %in% valid_order) {
    stop("order_by must be one of: ", paste(valid_order, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Order and subset data
  # ---------------------------------------------------------------------------
  # Determine sort direction (ascending for p-values, descending for OR/Count)
  if (order_by %in% c("p.adjust", "pvalue")) {
    df <- df[order(df[[order_by]]), ]
  } else {
    df <- df[order(df[[order_by]], decreasing = TRUE), ]
  }

  # Take top categories
  if (nrow(df) > show_category) {
    df <- df[1:show_category, ]
  }

  # Truncate long descriptions for readability
  df$Description_short <- ifelse(
    nchar(df$Description) > 50,
    paste0(substr(df$Description, 1, 47), "..."),
    df$Description
  )

  # Order factor levels for plotting (reverse so top is at top of plot)
  df$Description_short <- factor(df$Description_short,
                                  levels = rev(df$Description_short))

  # ---------------------------------------------------------------------------
  # Set up size variable
  # ---------------------------------------------------------------------------
  if (size_by == "GeneRatio") {
    df$size_var <- df$GeneRatio_numeric
    size_label <- "Gene Ratio"
  } else {
    df$size_var <- df$Count
    size_label <- "Gene Count"
  }

  # ---------------------------------------------------------------------------
  # Set up color variable and scale
  # ---------------------------------------------------------------------------
  if (color_by == "odds_ratio") {
    df$color_var <- df$odds_ratio
    color_label <- "Odds Ratio"

    # For odds ratio, use log2 transform for better visualization
    # and set midpoint at 1 (no enrichment) if using diverging scale
    df$color_var_plot <- log2(df$color_var)

    if (is.null(color_midpoint)) {
      # Sequential scale from low to high OR
      color_scale <- scale_color_gradient(
        low = color_low,
        high = color_high,
        name = paste0(color_label, "\n(log2)")
      )
    } else {
      color_scale <- scale_color_gradient2(
        low = color_low,
        mid = "white",
        high = color_high,
        midpoint = log2(color_midpoint),
        name = paste0(color_label, "\n(log2)")
      )
    }
  } else {
    # For p-values, use -log10 transform
    df$color_var <- df[[color_by]]
    df$color_var_plot <- -log10(df$color_var)
    color_label <- ifelse(color_by == "p.adjust", "-log10(q-value)", "-log10(p-value)")

    color_scale <- scale_color_gradient(
      low = color_low,
      high = color_high,
      name = color_label
    )
  }

  # ---------------------------------------------------------------------------
  # Set up x-axis label
  # ---------------------------------------------------------------------------
  if (is.null(x_label)) {
    x_label <- ifelse(size_by == "GeneRatio", "Gene Ratio", "Gene Count")
  }

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot(df, aes(x = size_var, y = Description_short)) +
    geom_point(aes(size = size_var, color = color_var_plot)) +
    color_scale +
    scale_size_continuous(range = c(3, 10), name = size_label) +
    labs(
      title = title,
      x = x_label,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = font_size + 4),
      axis.text.y = element_text(size = font_size),
      axis.text.x = element_text(size = font_size - 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  # Add subtitle with summary stats
  n_total <- nrow(as.data.frame(enrich_result))
  subtitle <- sprintf("Showing top %d of %d terms (ordered by %s)",
                      nrow(df), n_total, order_by)
  p <- p + labs(subtitle = subtitle)

  return(p)
}


#' Bar Plot for Enrichment Results
#'
#' @description Creates a horizontal bar plot for GO/KEGG enrichment results.
#' Simpler alternative to dotplot, showing either gene count or odds ratio.
#'
#' @param enrich_result An enrichResult object or data.frame.
#' @param show_category Integer. Number of top categories to display (default = 15).
#' @param x_var Character. Variable for bar length: "Count", "odds_ratio", or
#'   "GeneRatio" (default = "Count").
#' @param fill_by Character. Variable for bar fill: "p.adjust" (default),
#'   "pvalue", or a fixed color string.
#' @param order_by Character. Variable to order terms (default = same as x_var).
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#'
#' @export
#'
AETHER_plot_enrichment_bar <- function(enrich_result,
                                        show_category = 15,
                                        x_var = "Count",
                                        fill_by = "p.adjust",
                                        order_by = NULL,
                                        title = "Enrichment Analysis") {

  # Extract data
  if (inherits(enrich_result, "enrichResult")) {
    df <- as.data.frame(enrich_result)
  } else if (is.data.frame(enrich_result)) {
    df <- enrich_result
  } else {
    stop("enrich_result must be an enrichResult object or data.frame")
  }

  if (nrow(df) == 0) {
    warning("No enrichment results to plot")
    return(ggplot() + theme_void() +
             labs(title = title, subtitle = "No significant terms found"))
  }

  # Calculate odds ratio
  df$odds_ratio <- .calculate_odds_ratio(df$GeneRatio, df$BgRatio)

  # Parse GeneRatio
  gene_parts <- strsplit(df$GeneRatio, "/")
  df$GeneRatio_numeric <- as.numeric(sapply(gene_parts, `[`, 1)) /
                          as.numeric(sapply(gene_parts, `[`, 2))

  # Set order_by default
  if (is.null(order_by)) {
    order_by <- x_var
  }

  # Validate x_var
  if (!x_var %in% c("Count", "odds_ratio", "GeneRatio")) {
    stop("x_var must be 'Count', 'odds_ratio', or 'GeneRatio'")
  }

  # Set up x variable
  if (x_var == "GeneRatio") {
    df$x_val <- df$GeneRatio_numeric
    x_label <- "Gene Ratio"
  } else if (x_var == "odds_ratio") {
    df$x_val <- df$odds_ratio
    x_label <- "Odds Ratio"
  } else {
    df$x_val <- df$Count
    x_label <- "Gene Count"
  }

  # Order data
  if (order_by %in% c("p.adjust", "pvalue")) {
    df <- df[order(df[[order_by]]), ]
  } else if (order_by == "GeneRatio") {
    df <- df[order(df$GeneRatio_numeric, decreasing = TRUE), ]
  } else if (order_by == "odds_ratio") {
    df <- df[order(df$odds_ratio, decreasing = TRUE), ]
  } else {
    df <- df[order(df[[order_by]], decreasing = TRUE), ]
  }

  # Subset
  if (nrow(df) > show_category) {
    df <- df[1:show_category, ]
  }

  # Truncate descriptions
  df$Description_short <- ifelse(
    nchar(df$Description) > 50,
    paste0(substr(df$Description, 1, 47), "..."),
    df$Description
  )

  df$Description_short <- factor(df$Description_short,
                                  levels = rev(df$Description_short))

  # Build plot
  if (fill_by %in% c("p.adjust", "pvalue")) {
    df$fill_var <- -log10(df[[fill_by]])
    fill_label <- ifelse(fill_by == "p.adjust", "-log10(q-value)", "-log10(p-value)")

    p <- ggplot(df, aes(x = x_val, y = Description_short, fill = fill_var)) +
      geom_col() +
      scale_fill_gradient(low = "#deebf7", high = "#2166ac", name = fill_label)
  } else {
    # Fixed color
    p <- ggplot(df, aes(x = x_val, y = Description_short)) +
      geom_col(fill = fill_by)
  }

  p <- p +
    labs(
      title = title,
      x = x_label,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  return(p)
}


#' Genomic Annotation Distribution Bar Plot
#'
#' @description Creates a stacked horizontal bar plot showing the distribution
#' of genomic features (Promoter, Exon, Intron, Intergenic, etc.) for one or
#' more sets of annotated peaks. Useful for comparing annotation profiles
#' between different datasets (e.g., ATAC vs ChIP peaks).
#'
#' @param annotated_list A named list of data.frames from APOLLO_annotate_peaks(),
#'   or a single data.frame. Each data.frame must have an 'annotation_simple'
#'   or 'annotation' column.
#' @param use_simple Logical. If TRUE (default), uses 'annotation_simple' column.
#'   If FALSE, uses full 'annotation' column.
#' @param show_percentage Logical. If TRUE (default), shows percentages.
#'   If FALSE, shows raw counts.
#' @param colors Named character vector of colors for each annotation category.
#'   If NULL, uses a default color scheme.
#' @param title Character. Plot title (default = "Genomic Feature Distribution").
#' @param bar_height Numeric. Height of each bar (default = 0.7).
#' @param show_labels Logical. If TRUE (default), shows percentage labels on bars.
#' @param min_label_pct Numeric. Minimum percentage to show label (default = 5).
#'   Prevents cluttering small segments with labels.
#'
#' @return A ggplot object.
#'
#' @details
#' The plot shows one horizontal stacked bar per dataset in annotated_list.
#' Each segment's width is proportional to the count/percentage of peaks
#' in that genomic category.
#'
#' Default annotation categories (from annotation_simple):
#' - Promoter: regions within TSS window
#' - 5' UTR, 3' UTR: untranslated regions
#' - Exon: exonic regions
#' - Intron: intronic regions
#' - Downstream: downstream of genes
#' - Intergenic: between genes
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single dataset
#' p <- AETHER_plot_annotation_bar(atac_annotated)
#'
#' # Compare ATAC vs ChIP annotation profiles
#' annotated_list <- list(
#'   "ATAC peaks" = atac_annotated,
#'   "ChIP peaks" = chip_annotated
#' )
#' p <- AETHER_plot_annotation_bar(annotated_list)
#'
#' # Custom colors
#' my_colors <- c(Promoter = "#e41a1c", Exon = "#377eb8", Intron = "#4daf4a")
#' p <- AETHER_plot_annotation_bar(annotated_list, colors = my_colors)
#'
#' }
AETHER_plot_annotation_bar <- function(annotated_list,
                                        use_simple = TRUE,
                                        show_percentage = TRUE,
                                        colors = NULL,
                                        title = "Genomic Feature Distribution",
                                        bar_height = 0.7,
                                        show_labels = TRUE,
                                        min_label_pct = 5) {

  # ---------------------------------------------------------------------------
  # Handle single data.frame input
  # ---------------------------------------------------------------------------
  if (is.data.frame(annotated_list)) {
    annotated_list <- list("Peaks" = annotated_list)
  }

  if (!is.list(annotated_list) || length(annotated_list) == 0) {
    stop("annotated_list must be a named list of data.frames or a single data.frame")
  }

  # Ensure list is named
  if (is.null(names(annotated_list))) {
    names(annotated_list) <- paste0("Dataset_", seq_along(annotated_list))
  }

  # ---------------------------------------------------------------------------
  # Determine annotation column to use
  # ---------------------------------------------------------------------------
  anno_col <- if (use_simple) "annotation_simple" else "annotation"

  # ---------------------------------------------------------------------------
  # Build combined data frame
  # ---------------------------------------------------------------------------
  combined_df <- do.call(rbind, lapply(names(annotated_list), function(name) {
    df <- annotated_list[[name]]

    # Check for annotation column
    if (!anno_col %in% colnames(df)) {
      if (anno_col == "annotation_simple" && "annotation" %in% colnames(df)) {
        # Try to create simplified annotation
        df$annotation_simple <- sapply(df$annotation, function(x) {
          if (grepl("Promoter", x)) return("Promoter")
          if (grepl("5' UTR", x)) return("5' UTR")
          if (grepl("3' UTR", x)) return("3' UTR")
          if (grepl("Exon", x)) return("Exon")
          if (grepl("Intron", x)) return("Intron")
          if (grepl("Downstream", x)) return("Downstream")
          if (grepl("Intergenic", x)) return("Intergenic")
          return("Other")
        })
      } else {
        stop("Data.frame '", name, "' missing '", anno_col, "' column. ",
             "Run APOLLO_annotate_peaks() first.")
      }
    }

    # Count annotations
    counts <- as.data.frame(table(df[[anno_col]]), stringsAsFactors = FALSE)
    colnames(counts) <- c("Feature", "Count")
    counts$Dataset <- name
    counts$Total <- nrow(df)
    counts$Percentage <- 100 * counts$Count / counts$Total

    return(counts)
  }))

  # ---------------------------------------------------------------------------
  # Set up feature order (consistent across datasets)
  # ---------------------------------------------------------------------------
  # Define preferred order for genomic features
  preferred_order <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR",
                       "Downstream", "Intergenic", "Other")

  all_features <- unique(combined_df$Feature)
  # Order: preferred features first (in order), then any others alphabetically
  feature_order <- c(
    intersect(preferred_order, all_features),
    setdiff(all_features, preferred_order)
  )

  combined_df$Feature <- factor(combined_df$Feature, levels = feature_order)

  # Set dataset order (preserve input order) and add counts to labels
  totals <- aggregate(Total ~ Dataset, combined_df, function(x) x[1])
  dataset_labels <- setNames(
    paste0(totals$Dataset, " (", format(totals$Total, big.mark = ","), ")"),
    totals$Dataset
  )
  combined_df$Dataset_label <- dataset_labels[as.character(combined_df$Dataset)]
  combined_df$Dataset_label <- factor(
    combined_df$Dataset_label,
    levels = rev(dataset_labels[names(annotated_list)])
  )

  # ---------------------------------------------------------------------------
  # Set up colors
  # ---------------------------------------------------------------------------
  if (is.null(colors)) {
    # Default ChIPseeker-inspired color scheme
    colors <- c(
      "Promoter" = "#F8766D",      # Red/coral
      "5' UTR" = "#CD9600",        # Gold
      "Exon" = "#7CAE00",          # Green
      "Intron" = "#00BE67",        # Teal-green
      "3' UTR" = "#00BFC4",        # Cyan
      "Downstream" = "#00A9FF",    # Light blue
      "Intergenic" = "#C77CFF",    # Purple
      "Other" = "#878787"          # Gray
    )
  }

  # Ensure all features have colors
  missing_colors <- setdiff(all_features, names(colors))
  if (length(missing_colors) > 0) {
    # Generate colors for missing features
    extra_colors <- scales::hue_pal()(length(missing_colors))
    names(extra_colors) <- missing_colors
    colors <- c(colors, extra_colors)
  }

  # ---------------------------------------------------------------------------
  # Create label text (positioning handled by position_stack)
  # ---------------------------------------------------------------------------
  if (show_labels) {
    combined_df$label_text <- ifelse(
      combined_df$Percentage >= min_label_pct,
      sprintf("%.1f%%", combined_df$Percentage),
      ""
    )
  }

  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  n_features <- length(feature_order)
  legend_nrow <- ceiling(n_features / 4L)   # max 4 items per row

  y_var <- if (show_percentage) "Percentage" else "Count"

  p <- ggplot(combined_df, aes(x = .data[[y_var]], y = Dataset_label, fill = Feature)) +
    geom_bar(stat = "identity", position = "stack", width = bar_height,
             color = "white", linewidth = 0.3) +
    scale_fill_manual(values = colors, name = "Genomic\nFeature") +
    labs(
      title = title,
      x = if (show_percentage) "Percentage (%)" else "Count",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(nrow = legend_nrow))

  # Add percentage labels using position_stack for automatic alignment
  # IMPORTANT: Use full combined_df (not filtered) so position_stack sees all segments
  # label_text is already "" for small segments, so they render as invisible
  if (show_labels && show_percentage) {
    p <- p + geom_text(
      data = combined_df,
      aes(label = label_text),
      position = position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3
    )
  }

  # Cap x-axis at 100 for percentage view
  if (show_percentage) {
    p <- p + coord_cartesian(xlim = c(0, 100))
  }

  return(p)
}


###############################################################################
########### GO Treemap via Semantic Similarity Reduction ###########
###############################################################################

#' Internal: map gprofiler2 organism to OrgDb package name
#' @noRd
.aether_gost_orgdb <- function(organism) {
  map <- c(
    hsapiens      = "org.Hs.eg.db",
    mmusculus     = "org.Mm.eg.db",
    rnorvegicus   = "org.Rn.eg.db",
    drerio        = "org.Dr.eg.db",
    dmelanogaster = "org.Dm.eg.db",
    celegans      = "org.Ce.eg.db",
    scerevisiae   = "org.Sc.sgd.db",
    sscrofa       = "org.Ss.eg.db",
    btaurus       = "org.Bt.eg.db",
    cfamiliaris   = "org.Cf.eg.db",
    ggallus       = "org.Gg.eg.db"
  )
  pkg <- map[tolower(trimws(organism))]
  if (is.na(pkg))
    stop("Organism '", organism, "' not recognised. ",
         "Supported gprofiler2 organisms: ", paste(names(map), collapse = ", "),
         ".\nIf your organism is not listed, open an issue or pass the OrgDb package ",
         "name directly to rrvgo::calculateSimMatrix().")
  unname(pkg)
}


#' Internal: generate n categorical colors for module labeling
#' @noRd
.aether_n_colors <- function(n) {
  base_cols <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
  )
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  grDevices::colorRampPalette(base_cols)(n)
}

# Internal: override "up"/"down" module colours to red/blue respectively.
# Applied after auto-assignment so directional comparisons use intuitive colours.
# Case-insensitive: matches "up", "Up", "UP", "down", "Down", "DOWN".
.aether_updown_override <- function(mod_colors) {
  nms_lower <- tolower(names(mod_colors))
  if ("up"   %in% nms_lower) mod_colors[nms_lower == "up"]   <- "#E41A1C"  # red
  if ("down" %in% nms_lower) mod_colors[nms_lower == "down"] <- "#377EB8"  # blue
  mod_colors
}


#' GO Term Treemap via Semantic Similarity Reduction
#'
#' @description Reduces a set of enriched GO terms to representative clusters
#' using semantic similarity (via \pkg{rrvgo}), then renders a treemap where
#' each tile represents one GO term coloured by the module/gene-list that found
#' it most significantly.  Tiles are grouped into larger semantic clusters
#' (e.g., "immune response") based on GO hierarchy.
#'
#' @param enrichment_result A \code{gost_enrichment} object from
#'   \code{\link{APOLLO_enrich_gost}} or \code{\link{DEMETER_load_enrichment}},
#'   or a plain data.frame with at least the columns \code{term_id},
#'   \code{source}, \code{p_value}, and \code{module}.
#' @param orgdb An OrgDb object (e.g. \code{org.Hs.eg.db}) or a package name
#'   string (e.g. \code{"org.Hs.eg.db"}) passed to
#'   \code{rrvgo::calculateSimMatrix()} and \code{rrvgo::reduceSimMatrix()}.
#'   Required — load the package first with
#'   \code{library(org.Hs.eg.db)} (human) or
#'   \code{library(org.Mm.eg.db)} (mouse), then pass the object.
#' @param ont Character vector. GO ontology/ontologies to plot.
#'   One or more of \code{"BP"}, \code{"MF"}, \code{"CC"}.
#'   Multiple values produce one panel per ontology in a single figure.
#'   Default: \code{"BP"}.
#' @param top_n Integer. Top N GO terms per module (by p-value) pooled before
#'   semantic reduction.  Higher values include more terms but may slow down
#'   \code{rrvgo::calculateSimMatrix()}.  Default: \code{20L}.
#' @param threshold Numeric (0–1). Similarity threshold for
#'   \code{rrvgo::reduceSimMatrix()}.  Higher values → fewer, broader clusters;
#'   lower values → more, finer clusters.  Default: \code{0.7}.
#' @param method Character. Semantic similarity measure passed to
#'   \code{rrvgo::calculateSimMatrix()}.  One of \code{"Rel"} (default),
#'   \code{"Wang"}, \code{"Lin"}, \code{"Resnik"}, \code{"Jiang"}.
#' @param show_shared Logical. If \code{FALSE} (default), each GO term is
#'   assigned to the single module with the best (lowest) p-value for that
#'   term ("winner takes all").  If \code{TRUE}, terms found in multiple modules
#'   are assigned to a combined group labelled \code{"C1-C2"} (module names
#'   joined by \code{"-"}, sorted alphabetically) and given a distinct color
#'   automatically; useful for spotting convergent biology across modules.
#' @param show_group_labels Logical. If \code{TRUE} (default), each large
#'   parent-cluster tile shows the module group attribution beneath the semantic
#'   cluster name, e.g. \code{"immune response\n(C7)"}.  The group shown is the
#'   one assigned to the cluster's representative (highest-scoring) term.
#'   Set to \code{FALSE} to display the semantic cluster name only.
#' @param min_label_scale Numeric (0–1). Minimum text scaling fraction passed
#'   to \code{treemap}'s \code{lowerbound.cex.labels}.  When a tile is too
#'   small to fit text at the full \code{fontsize.labels}, treemap scales the
#'   text down; if the required scale falls below this value the label is
#'   suppressed entirely.  Default \code{0.1} is more permissive than
#'   \code{treemap}'s own default of \code{0.4}, so more labels appear in
#'   small tiles.  Set to \code{0} to force all labels regardless of tile size
#'   (may produce very small text).  Opening a larger graphics device before
#'   calling this function also helps, as larger tiles require less scaling.
#' @param palette Named character vector. Module name → hex color mapping.
#'   If \code{NULL} (default), colors are auto-assigned from an internal
#'   20-color categorical palette.  Unknown modules fall back to auto-generated
#'   colors.  The combined/shared group always uses \code{"#AAAAAA"} regardless
#'   of the palette.
#' @param title Character. Main plot title.  If \code{NULL} (default), the
#'   ontology label is used (e.g., "GO Biological Process").  When multiple
#'   ontologies are requested, the ontology label is appended automatically.
#' @param verbose Logical. Print progress messages.  Default: \code{TRUE}.
#'
#' @return Invisibly returns a named list of \code{reducedTerms} data.frames
#'   (one per ontology, named by ontology code: "BP", "MF", "CC").  The
#'   primary output is the treemap rendered to the current graphics device.
#'   Save with \code{png()} / \code{pdf()} wrappers before calling.
#'
#' @details
#' **Workflow:**
#' 1. Filters the combined enrichment table to the requested GO source(s).
#' 2. Takes the top \code{top_n} terms per module by p-value and pools them.
#' 3. For each unique GO term, uses the best (lowest) p-value across modules
#'    as the rrvgo score (\code{-log10(p)}).
#' 4. Computes a semantic similarity matrix via
#'    \code{rrvgo::calculateSimMatrix()} (requires the relevant OrgDb package).
#' 5. Reduces to representative clusters via \code{rrvgo::reduceSimMatrix()}.
#' 6. Assigns each term a module group (winner or shared, per \code{show_shared}).
#' 7. Renders via \code{treemap::treemap()} with parent clusters as outer
#'    groupings and individual terms as inner tiles coloured by module group.
#'
#' **OrgDb auto-detection:** The organism is read from
#' \code{enrichment_result$metadata$organism} (set automatically by
#' \code{APOLLO_enrich_gost()}).  Supported organisms and their OrgDb packages:
#' \code{hsapiens} → \code{org.Hs.eg.db};
#' \code{mmusculus} → \code{org.Mm.eg.db};
#' \code{rnorvegicus} → \code{org.Rn.eg.db};
#' \code{drerio} → \code{org.Dr.eg.db};
#' \code{dmelanogaster} → \code{org.Dm.eg.db};
#' \code{celegans} → \code{org.Ce.eg.db};
#' \code{scerevisiae} → \code{org.Sc.sgd.db};
#' \code{sscrofa} → \code{org.Ss.eg.db};
#' \code{btaurus} → \code{org.Bt.eg.db};
#' \code{cfamiliaris} → \code{org.Cf.eg.db};
#' \code{ggallus} → \code{org.Gg.eg.db}.
#' The corresponding OrgDb package must be installed.
#'
#' **Required packages:** \pkg{rrvgo} (Bioconductor) and \pkg{treemap} (CRAN),
#' plus the appropriate OrgDb package for the organism.
#'
#' @seealso \code{\link{APOLLO_enrich_gost}}, \code{\link{DEMETER_load_enrichment}}
#'
#' @examples
#' \dontrun{
#' library(org.Hs.eg.db)  # human; use org.Mm.eg.db for mouse
#'
#' # Basic GO:BP treemap
#' AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db)
#'
#' # Multiple ontologies side by side
#' AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, ont = c("BP", "MF"))
#'
#' # Show shared terms across modules in grey
#' AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, show_shared = TRUE)
#'
#' # Custom module palette + save to PNG
#' my_pal <- c(C1 = "#E41A1C", C2 = "#377EB8", C3 = "#4DAF4A")
#' png("treemap_BP.png", width = 10, height = 8, units = "in", res = 300)
#' AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, ont = "BP", palette = my_pal)
#' dev.off()
#' }
#' @export
AETHER_plot_go_treemap <- function(enrichment_result,
                                    orgdb,
                                    ont               = "BP",
                                    top_n             = 20L,
                                    threshold         = 0.7,
                                    method            = "Rel",
                                    show_shared       = FALSE,
                                    show_group_labels = TRUE,
                                    min_label_scale   = 0.1,
                                    palette           = NULL,
                                    title             = NULL,
                                    verbose           = TRUE) {

  # ---------------------------------------------------------------------------
  # Package checks
  # ---------------------------------------------------------------------------
  if (!requireNamespace("rrvgo", quietly = TRUE))
    stop("Package 'rrvgo' is required. Install with: BiocManager::install('rrvgo')")
  if (!requireNamespace("treemap", quietly = TRUE))
    stop("Package 'treemap' is required. Install with: install.packages('treemap')")

  # ---------------------------------------------------------------------------
  # Extract data
  # ---------------------------------------------------------------------------
  if (inherits(enrichment_result, "gost_enrichment")) {
    df <- enrichment_result$combined
  } else if (is.data.frame(enrichment_result)) {
    df <- enrichment_result
  } else {
    stop("enrichment_result must be a gost_enrichment object or data.frame")
  }

  if (nrow(df) == 0L) stop("No enrichment results to plot")

  # ---------------------------------------------------------------------------
  # Validate arguments
  # ---------------------------------------------------------------------------
  ont    <- match.arg(ont, c("BP", "MF", "CC"), several.ok = TRUE)
  method <- match.arg(method, c("Rel", "Wang", "Lin", "Resnik", "Jiang"))
  top_n  <- as.integer(top_n)

  # ---------------------------------------------------------------------------
  # Build base color palette for single modules.
  # Shared groups (if show_shared = TRUE) are assigned additional colors inside
  # the loop, extending this pool so all groups remain visually distinct.
  # ---------------------------------------------------------------------------
  all_modules <- sort(unique(as.character(df$module)))
  if (is.null(palette)) {
    mod_colors <- setNames(.aether_n_colors(length(all_modules)), all_modules)
    mod_colors <- .aether_updown_override(mod_colors)
  } else {
    mod_colors <- palette
    extra_mods <- setdiff(all_modules, names(mod_colors))
    if (length(extra_mods) > 0L) {
      extra_cols <- .aether_n_colors(length(extra_mods))
      mod_colors <- c(mod_colors, setNames(extra_cols, extra_mods))
    }
  }

  # ---------------------------------------------------------------------------
  # Ontology metadata
  # ---------------------------------------------------------------------------
  ont_source <- c(BP = "GO:BP", MF = "GO:MF", CC = "GO:CC")
  ont_label  <- c(
    BP = "GO Biological Process",
    MF = "GO Molecular Function",
    CC = "GO Cellular Component"
  )

  # ---------------------------------------------------------------------------
  # Filter to GO terms only
  # ---------------------------------------------------------------------------
  df_go_all <- df[grepl("^GO:", df$source), ]
  if (nrow(df_go_all) == 0L)
    stop("No GO terms found in enrichment results. ",
         "Confirm GO:BP, GO:MF, or GO:CC were included as sources.")

  # ---------------------------------------------------------------------------
  # Multi-panel layout when multiple ontologies requested
  # ---------------------------------------------------------------------------
  n_ont <- length(ont)
  if (n_ont > 1L) {
    old_par <- graphics::par(mfrow = c(1L, n_ont))
    on.exit(graphics::par(old_par), add = TRUE)
  }

  reduced_list <- list()

  for (o in ont) {

    df_ont <- df_go_all[df_go_all$source == ont_source[o], ]
    if (nrow(df_ont) == 0L) {
      if (verbose) cat("[AETHER] No ", ont_source[o], " terms found — skipping.")
      next
    }

    # -------------------------------------------------------------------------
    # Pool top_n terms per module
    # -------------------------------------------------------------------------
    df_ont  <- df_ont[order(df_ont$p_value), ]
    df_pool <- do.call(rbind, lapply(split(df_ont, df_ont$module),
                                     function(x) utils::head(x, top_n)))
    rownames(df_pool) <- NULL

    # -------------------------------------------------------------------------
    # Best p-value per unique GO term (rrvgo score)
    # -------------------------------------------------------------------------
    go_ids <- unique(df_pool$term_id)
    best_p <- tapply(df_pool$p_value, df_pool$term_id, min)
    scores <- setNames(-log10(as.numeric(best_p[go_ids])), go_ids)

    # -------------------------------------------------------------------------
    # Semantic similarity matrix + reduction
    # -------------------------------------------------------------------------
    if (verbose)
      cat(sprintf("[AETHER] GO %s: computing semantic similarity for %d terms...\n",
                  o, length(go_ids)))

    sim_mat <- tryCatch(
      rrvgo::calculateSimMatrix(go_ids, orgdb = orgdb,
                                ont = o, method = method),
      error = function(e) {
        warning("calculateSimMatrix failed for GO:", o,
                " — ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(sim_mat)) next

    # Align to terms present in OrgDb
    common_ids <- intersect(rownames(sim_mat), names(scores))
    if (length(common_ids) < 2L) {
      if (verbose)
        cat("Too few GO:", o, " terms present in OrgDb after filtering — skipping.")
      next
    }
    sim_mat  <- sim_mat[common_ids, common_ids]
    scores   <- scores[common_ids]
    df_pool  <- df_pool[df_pool$term_id %in% common_ids, ]

    red <- rrvgo::reduceSimMatrix(sim_mat, scores = scores,
                                  threshold = threshold,
                                  orgdb = orgdb)

    # -------------------------------------------------------------------------
    # Assign module group to each reduced term
    # -------------------------------------------------------------------------
    red$module_group <- vapply(red$go, function(go_id) {
      rows <- df_pool[df_pool$term_id == go_id, ]
      mods <- unique(rows$module)
      if (length(mods) == 0L) return("unknown")
      if (length(mods) == 1L || !show_shared) {
        # Winner takes all: module with best p-value for this term
        as.character(rows$module[which.min(rows$p_value)])
      } else {
        paste(sort(mods), collapse = "-")
      }
    }, FUN.VALUE = character(1L))

    # -------------------------------------------------------------------------
    # Assign distinct colors to any shared groups not yet in mod_colors.
    # Shared groups get colors continuing from index (n_single + 1) in the
    # same .aether_n_colors() pool so all groups remain visually distinct.
    # treemap::treemap() assigns palette colors POSITIONALLY by the
    # alphabetically-sorted factor levels of vColor — sort all_groups to match.
    # -------------------------------------------------------------------------
    shared_groups <- sort(setdiff(unique(red$module_group), names(mod_colors)))
    if (length(shared_groups) > 0L) {
      n_existing   <- length(mod_colors)
      all_extended <- .aether_n_colors(n_existing + length(shared_groups))
      new_cols     <- all_extended[(n_existing + 1L):(n_existing + length(shared_groups))]
      mod_colors   <- c(mod_colors, setNames(new_cols, shared_groups))
    }

    all_groups   <- sort(unique(red$module_group))
    group_colors <- mod_colors[all_groups]

    # -------------------------------------------------------------------------
    # Optionally append the cluster representative's module group to the
    # parentTerm label so each big tile shows both the semantic cluster name
    # and which module/group it was found in.
    # -------------------------------------------------------------------------
    if (show_group_labels) {
      rep_group  <- setNames(red$module_group, red$go)
      red$parentTerm <- paste0(red$parentTerm, "\n(",
                               rep_group[red$parent], ")")
    }

    # -------------------------------------------------------------------------
    # Wrap term names for tile readability
    # -------------------------------------------------------------------------
    red$term_label <- vapply(red$term, function(x) {
      paste(strwrap(x, width = 22L), collapse = "\n")
    }, FUN.VALUE = character(1L))

    # -------------------------------------------------------------------------
    # Render treemap
    # -------------------------------------------------------------------------
    plot_title <- if (!is.null(title)) {
      if (n_ont > 1L) paste(title, "\u2014", ont_label[o]) else title
    } else {
      ont_label[o]
    }

    if (verbose) cat("[AETHER] Rendering treemap:", plot_title, "\n")

    treemap::treemap(
      red,
      index                  = c("parentTerm", "term_label"),
      vSize                  = "score",
      vColor                 = "module_group",
      type                   = "categorical",
      palette                = group_colors,
      title                  = plot_title,
      fontsize.title         = 14L,
      fontsize.labels        = c(11L, 8L),
      fontface.labels        = c("bold", "plain"),
      fontcolor.labels       = c("white", "white"),
      align.labels           = list(c("center", "top"), c("center", "center")),
      overlap.labels         = 0.5,
      border.col             = c("white", "white"),
      border.lwds            = c(2, 0.5),
      bg.labels              = 0L,
      lowerbound.cex.labels  = min_label_scale,
      aspRatio               = 1,
      title.legend           = "Module",
      position.legend        = "bottom"
    )

    reduced_list[[o]] <- red
  }

  invisible(reduced_list)
}


###############################################################################
########### GO DAG Plot ###########
###############################################################################

#' Internal: return the GO.db parents environment for the given ontology
#' @noRd
.go_parents_env <- function(ont) {
  switch(ont,
    BP = GO.db::GOBPPARENTS,
    MF = GO.db::GOMFPARENTS,
    CC = GO.db::GOCCPARENTS
  )
}

#' Internal: look up a GO term's display name; falls back to GO ID on failure
#' @noRd
.go_term_label <- function(go_id) {
  tryCatch({
    gt <- GO.db::GOTERM[[go_id]]
    if (is.null(gt)) return(go_id)
    GO.db::Term(gt)
  }, error = function(e) go_id)
}


#' GO DAG Plot coloured by module
#'
#' @description Builds the GO Directed Acyclic Graph (DAG) induced by a set of
#' enriched GO terms, traverses upward toward the root to include ancestor
#' context nodes, and renders the result as a hierarchical network plot using
#' \pkg{ggraph}.  Enriched nodes are coloured by module of origin; ancestor
#' (non-enriched) context nodes are shown in grey.
#'
#' @param enrichment_result A \code{gost_enrichment} object from
#'   \code{\link{APOLLO_enrich_gost}} / \code{\link{DEMETER_load_enrichment}},
#'   or a plain data.frame with columns \code{term_id}, \code{source},
#'   \code{p_value}, \code{module}.
#' @param ont Character. GO ontology: \code{"BP"} (default), \code{"MF"}, or
#'   \code{"CC"}.  Only one ontology per call.
#' @param top_n Integer. Top N enriched GO terms per module (by p-value) to
#'   include before building the DAG.  Default: \code{20L}.
#' @param min_ancestor_freq Integer. Minimum number of enriched terms an
#'   ancestor node must be an ancestor of in order to be retained.  Enriched
#'   terms themselves are always kept.  Lowering this value shows more
#'   context; raising it focuses on shared hubs only.  Default: \code{2L}.
#' @param max_depth Integer or \code{NULL}.  Maximum number of levels to
#'   traverse upward from enriched terms.  \code{NULL} traverses all the way
#'   to the ontology root.  Default: \code{4L}.
#' @param show_shared Logical.  If \code{FALSE} (default), each enriched term
#'   is coloured by its single best-p module.  If \code{TRUE}, terms found
#'   in multiple modules get a combined label (e.g., \code{"C1-C2"}) and a
#'   distinct auto-assigned colour.
#' @param edge_types Character vector.  GO relationship types to include as
#'   edges.  Any combination of \code{"is_a"}, \code{"part_of"},
#'   \code{"regulates"}, \code{"positively_regulates"},
#'   \code{"negatively_regulates"}.  Default: all five.
#' @param node_size_range Numeric vector of length 2.  Min and max point sizes
#'   for enriched nodes (scaled by \code{-log10(p)}).  Default: \code{c(3, 10)}.
#' @param palette Named character vector.  Module name → hex color mapping.
#'   \code{NULL} (default) auto-assigns from internal palette.
#' @param title Character.  Plot title.  \code{NULL} (default) auto-generates.
#' @param verbose Logical.  Print progress messages.  Default: \code{TRUE}.
#'
#' @return A \code{ggplot} / \code{ggraph} object.  Save with
#'   \code{ggplot2::ggsave()}.
#'
#' @details
#' **Algorithm:**
#' 1. Filters \code{enrichment_result$combined} to the requested ontology and
#'    takes the top \code{top_n} terms per module by p-value.
#' 2. Traverses the GO DAG upward (child → parent) from enriched terms using
#'    \code{GO.db::GOBPPARENTS} (or MF/CC), up to \code{max_depth} levels.
#' 3. For each ancestor node, counts how many enriched terms it is an ancestor
#'    of.  Ancestors below \code{min_ancestor_freq} are removed.
#' 4. Renders the filtered DAG with \pkg{ggraph} using the Sugiyama
#'    (layered/hierarchical) layout.  Enriched nodes are coloured by module
#'    group; ancestor nodes are grey.  Edge line type encodes GO relationship.
#'
#' **Required packages:** \pkg{GO.db} (Bioconductor), \pkg{ggraph},
#' \pkg{igraph} (both already in Suggests).
#'
#' @seealso \code{\link{APOLLO_enrich_gost}}, \code{\link{AETHER_plot_go_treemap}}
#'
#' @examples
#' \dontrun{
#' # Basic GO:BP DAG
#' p <- AETHER_plot_go_dag(gost_res, ont = "BP")
#' ggplot2::ggsave("go_dag_bp.png", p, width = 14, height = 10)
#'
#' # Tighter focus: only high-frequency hubs, shallow traversal
#' p <- AETHER_plot_go_dag(gost_res, ont = "BP",
#'                          min_ancestor_freq = 3, max_depth = 3)
#'
#' # is_a edges only (cleaner)
#' p <- AETHER_plot_go_dag(gost_res, ont = "BP", edge_types = "is_a")
#' }
#' @export
AETHER_plot_go_dag <- function(enrichment_result,
                                ont               = "BP",
                                top_n             = 20L,
                                min_ancestor_freq = 2L,
                                max_depth         = 4L,
                                show_shared       = FALSE,
                                edge_types        = c("is_a", "part_of", "regulates",
                                                      "positively_regulates",
                                                      "negatively_regulates"),
                                node_size_range   = c(3, 10),
                                palette           = NULL,
                                title             = NULL,
                                verbose           = TRUE) {

  # ---------------------------------------------------------------------------
  # Package checks
  # ---------------------------------------------------------------------------
  for (pkg in c("GO.db", "ggraph", "igraph")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with: BiocManager::install('", pkg, "')")
  }

  # ---------------------------------------------------------------------------
  # Extract and validate
  # ---------------------------------------------------------------------------
  if (inherits(enrichment_result, "gost_enrichment")) {
    df <- enrichment_result$combined
  } else if (is.data.frame(enrichment_result)) {
    df <- enrichment_result
  } else {
    stop("enrichment_result must be a gost_enrichment object or data.frame")
  }
  if (nrow(df) == 0L) stop("No enrichment results to plot")

  ont               <- match.arg(ont, c("BP", "MF", "CC"))
  top_n             <- as.integer(top_n)
  min_ancestor_freq <- as.integer(min_ancestor_freq)
  max_depth         <- if (!is.null(max_depth)) as.integer(max_depth) else NULL

  ont_source <- c(BP = "GO:BP", MF = "GO:MF", CC = "GO:CC")
  df_go <- df[df$source == ont_source[ont], ]
  if (nrow(df_go) == 0L)
    stop("No ", ont_source[ont], " terms found in enrichment results.")

  # ---------------------------------------------------------------------------
  # Pool top_n per module
  # ---------------------------------------------------------------------------
  df_go   <- df_go[order(df_go$p_value), ]
  df_pool <- do.call(rbind, lapply(split(df_go, df_go$module),
                                   function(x) utils::head(x, top_n)))
  rownames(df_pool) <- NULL
  enriched_ids <- unique(df_pool$term_id)
  if (verbose) cat(sprintf("[AETHER] GO %s DAG: %d enriched terms\n", ont, length(enriched_ids)))

  # ---------------------------------------------------------------------------
  # Module group + color assignment (same logic as treemap)
  # ---------------------------------------------------------------------------
  all_modules <- sort(unique(as.character(df$module)))
  if (is.null(palette)) {
    mod_colors <- setNames(.aether_n_colors(length(all_modules)), all_modules)
    mod_colors <- .aether_updown_override(mod_colors)
  } else {
    mod_colors <- palette
    extra_mods <- setdiff(all_modules, names(mod_colors))
    if (length(extra_mods) > 0L)
      mod_colors <- c(mod_colors,
                      setNames(.aether_n_colors(length(extra_mods)), extra_mods))
  }

  term_group <- vapply(enriched_ids, function(go_id) {
    rows <- df_pool[df_pool$term_id == go_id, ]
    mods <- unique(rows$module)
    if (length(mods) == 0L) return("unknown")
    if (length(mods) == 1L || !show_shared)
      as.character(rows$module[which.min(rows$p_value)])
    else
      paste(sort(mods), collapse = "-")
  }, FUN.VALUE = character(1L))

  shared_groups <- sort(setdiff(unique(term_group), names(mod_colors)))
  if (length(shared_groups) > 0L) {
    n_ex     <- length(mod_colors)
    new_cols <- .aether_n_colors(n_ex + length(shared_groups))
    mod_colors <- c(mod_colors,
                    setNames(new_cols[(n_ex + 1L):(n_ex + length(shared_groups))],
                             shared_groups))
  }

  best_p <- tapply(df_pool$p_value, df_pool$term_id, min)

  # ---------------------------------------------------------------------------
  # BFS upward through the GO DAG from enriched terms
  # Tracks, for each ancestor, the set of enriched terms that "found" it
  # (used later for min_ancestor_freq filtering).
  # ---------------------------------------------------------------------------
  # Pre-convert the bimap to a plain R list once — avoids S4 dispatch issues
  # with [[]] on Go3AnnDbBimap objects and is faster for repeated lookups.
  # Format: names(entry) = relationship types; values = parent GO IDs.
  if (verbose) cat("[AETHER] Loading GO", ont, "parent mappings...\n")
  parents_list  <- as.list(.go_parents_env(ont))

  all_node_ids  <- enriched_ids
  depth_map     <- setNames(rep(0L, length(enriched_ids)), enriched_ids)
  # ancestor_seeds: node_id -> set of enriched IDs that are descendants
  anc_seeds     <- setNames(as.list(enriched_ids), enriched_ids)
  edge_rows     <- list()
  frontier      <- enriched_ids

  while (length(frontier) > 0L) {
    new_frontier <- character(0)
    for (child_id in frontier) {
      child_depth <- depth_map[[child_id]]
      if (!is.null(max_depth) && child_depth >= max_depth) next

      parents_raw <- parents_list[[child_id]]
      if (is.null(parents_raw) || length(parents_raw) == 0L) next

      # names = relationship types; values = parent GO IDs
      rel_types  <- names(parents_raw)
      parent_ids <- as.character(parents_raw)

      for (i in seq_along(parent_ids)) {
        pid <- parent_ids[[i]]
        rel <- rel_types[[i]]
        if (pid == "all") next  # GO root sentinel
        # GO.db stores "isa" (no underscore); normalise to "is_a"
        rel <- sub("^isa$", "is_a", rel)
        if (!rel %in% edge_types) next  # filtered relationship type

        edge_rows[[length(edge_rows) + 1L]] <- c(pid, child_id, rel)

        child_seeds <- anc_seeds[[child_id]]
        if (!pid %in% all_node_ids) {
          all_node_ids    <- c(all_node_ids, pid)
          depth_map[[pid]] <- child_depth + 1L
          anc_seeds[[pid]] <- child_seeds
          new_frontier     <- c(new_frontier, pid)
        } else {
          # Merge enriched seeds for this ancestor (union)
          anc_seeds[[pid]] <- union(anc_seeds[[pid]], child_seeds)
        }
      }
    }
    frontier <- unique(new_frontier)
  }

  if (length(edge_rows) == 0L) {
    warning("No edges found in GO DAG. Check GO.db is installed and GO IDs are valid.")
    return(invisible(NULL))
  }

  edges_df           <- as.data.frame(do.call(rbind, edge_rows), stringsAsFactors = FALSE)
  colnames(edges_df) <- c("from", "to", "rel_type")
  edges_df           <- unique(edges_df)

  # ---------------------------------------------------------------------------
  # Apply min_ancestor_freq: drop non-enriched ancestors that connect to fewer
  # than min_ancestor_freq enriched terms. Enriched nodes are always kept.
  # ---------------------------------------------------------------------------
  anc_freq <- vapply(all_node_ids, function(nid) {
    if (nid %in% enriched_ids) return(.Machine$integer.max)
    length(anc_seeds[[nid]])
  }, FUN.VALUE = integer(1L))

  keep_nodes <- names(anc_freq)[anc_freq >= min_ancestor_freq]
  if (length(keep_nodes) < 2L) {
    warning("min_ancestor_freq = ", min_ancestor_freq,
            " left fewer than 2 nodes. Lowering to keep all nodes.")
    keep_nodes <- all_node_ids
  }

  edges_df     <- edges_df[edges_df$from %in% keep_nodes &
                             edges_df$to   %in% keep_nodes, ]
  all_node_ids <- unique(c(enriched_ids,
                           intersect(all_node_ids, keep_nodes)))

  if (nrow(edges_df) == 0L) {
    warning("No edges remain after filtering. Try lowering min_ancestor_freq.")
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Node attribute table
  # ---------------------------------------------------------------------------
  enr_name_map <- setNames(df_pool$term_name, df_pool$term_id)
  enr_name_map <- enr_name_map[!duplicated(names(enr_name_map))]

  raw_labels <- vapply(all_node_ids, function(nid) {
    if (nid %in% names(enr_name_map))
      enr_name_map[[nid]]
    else
      .go_term_label(nid)
  }, FUN.VALUE = character(1L))

  node_labels <- vapply(raw_labels, function(x) {
    paste(strwrap(x, width = 20L), collapse = "\n")
  }, FUN.VALUE = character(1L))

  node_group <- vapply(all_node_ids, function(nid) {
    if (nid %in% names(term_group)) term_group[[nid]] else "ancestor"
  }, FUN.VALUE = character(1L))

  # Node size: scaled -log10(p) for enriched; fixed small for ancestors
  node_size_raw <- vapply(all_node_ids, function(nid) {
    if (nid %in% names(best_p)) -log10(as.numeric(best_p[[nid]])) else NA_real_
  }, FUN.VALUE = numeric(1L))

  enr_idx  <- !is.na(node_size_raw)
  enr_vals <- node_size_raw[enr_idx]
  rng      <- range(enr_vals, na.rm = TRUE)
  if (diff(rng) > 0) {
    node_size_raw[enr_idx] <- node_size_range[1] +
      (enr_vals - rng[1]) / diff(rng) * diff(node_size_range)
  } else {
    node_size_raw[enr_idx] <- mean(node_size_range)
  }
  node_size_raw[!enr_idx] <- node_size_range[1] * 0.6  # ancestors: smaller

  node_df <- data.frame(
    name       = all_node_ids,
    label      = node_labels,
    node_group = node_group,
    node_size  = node_size_raw,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Build igraph (edges_df cols: from, to, rel_type)
  # ---------------------------------------------------------------------------
  g <- igraph::graph_from_data_frame(
    d        = edges_df,
    directed = TRUE,
    vertices = node_df
  )

  # ---------------------------------------------------------------------------
  # Color and linetype scales
  # ---------------------------------------------------------------------------
  present_groups <- sort(unique(node_group))
  enr_groups     <- present_groups[present_groups != "ancestor"]
  fill_colors    <- c(
    setNames(vapply(enr_groups, function(g_) mod_colors[[g_]],
                    FUN.VALUE = character(1L)),
             enr_groups),
    ancestor = "#CCCCCC"
  )

  # Edge colour palette: distinct colours per relationship type
  all_edge_rel  <- sort(unique(edges_df$rel_type))
  edge_col_pool <- c("is_a"                   = "#888888",
                     part_of                  = "#4477AA",
                     regulates                = "#EE6677",
                     positively_regulates     = "#228833",
                     negatively_regulates     = "#CC3311")
  # Fall back for any unexpected types
  fallback_cols <- grDevices::hcl.colors(length(all_edge_rel), "Dark 2")
  edge_colors   <- stats::setNames(
    vapply(seq_along(all_edge_rel), function(i) {
      rel <- all_edge_rel[[i]]
      if (rel %in% names(edge_col_pool)) edge_col_pool[[rel]] else fallback_cols[[i]]
    }, character(1L)),
    all_edge_rel
  )

  # ---------------------------------------------------------------------------
  # Render
  # ---------------------------------------------------------------------------
  plot_title <- title %||% paste0(
    "GO ", c(BP = "Biological Process",
             MF = "Molecular Function",
             CC = "Cellular Component")[ont], " DAG"
  )

  if (verbose) cat("[AETHER] Rendering GO DAG (", igraph::vcount(g), "nodes,",
                   igraph::ecount(g), "edges)\n")

  p <- ggraph::ggraph(g, layout = "sugiyama") +
    ggraph::geom_edge_link(
      ggplot2::aes(color = rel_type),
      arrow     = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
      end_cap   = ggraph::circle(3, "mm"),
      linewidth = 0.4
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = node_group, size = node_size),
      shape  = 21,
      color  = "white",
      stroke = 0.5
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      size   = 2.5,
      repel  = TRUE,
      family = "sans"
    ) +
    ggraph::scale_edge_color_manual(values = edge_colors, name = "Relationship") +
    ggplot2::scale_fill_manual(
      values = fill_colors,
      name   = "Module",
      breaks = enr_groups
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  return(p)
}


# =============================================================================
# AETHER_plot_enrichment_map
# =============================================================================

#' Enrichment Map (Term-Term Network)
#'
#' @description
#' Visualises enriched pathway/GO terms as a network where nodes are terms and
#' edges connect terms that share a substantial fraction of their gene sets
#' (Jaccard similarity \eqn{\ge} \code{jaccard_threshold}).  Node size encodes
#' the best \eqn{-\log_{10}(p)} across all modules; node fill encodes the
#' primary module (lowest p-value).  Edge width is proportional to Jaccard
#' similarity.  A force-directed layout naturally clusters related terms into
#' "super-groups", making it easy to spot shared biology across modules.
#'
#' Terms with no edges above the threshold still appear as isolated nodes —
#' they represent module-specific biology not shared with other enriched terms.
#'
#' @details
#' Jaccard similarity is computed from the \code{intersection} gene lists in
#' the \code{gost_enrichment} object (requires \code{evcodes = TRUE} in
#' \code{APOLLO_enrich_gost}).  The union of intersection gene sets across all
#' modules is used per term, so similarity reflects shared biology rather than
#' query-specific overlap.
#'
#' Requires \pkg{ggraph} and \pkg{igraph} (both already in Suggests).
#'
#' @param enrichment_result A \code{gost_enrichment} object from
#'   \code{APOLLO_enrich_gost()} or a plain data.frame with columns
#'   \code{term_id}, \code{term_name}, \code{source}, \code{module},
#'   \code{p_value}, and optionally \code{intersection}.
#' @param source Character(1).  Database to display.  Default \code{"REAC"}.
#' @param top_n Integer.  Top N terms per module (by p-value) pooled before
#'   deduplication.  Default \code{20L}.  Reduce to 10–15 for cleaner graphs.
#' @param jaccard_threshold Numeric in \code{[0, 1]}.  Minimum Jaccard
#'   similarity to draw an edge.  Default \code{0.2}.  Increase to reduce
#'   edge density; decrease if too few edges appear.
#' @param layout Character.  igraph/ggraph layout algorithm.  \code{"fr"}
#'   (Fruchterman-Reingold, default) and \code{"kk"} (Kamada-Kawai) both work
#'   well; \code{"fr"} tends to produce rounder, more separated clusters.
#' @param node_size_range Numeric(2).  Min and max point size for nodes.
#'   Default \code{c(3, 12)}.
#' @param edge_width_range Numeric(2).  Min and max edge linewidth.
#'   Default \code{c(0.3, 2)}.
#' @param max_label_width Integer.  Maximum label characters before truncation.
#'   Default \code{40L}.
#' @param label_size Numeric.  \code{ggrepel} text size.  Default \code{2.5}.
#' @param palette Named character vector of colours, one per module.
#'   \code{NULL} (default) auto-assigns from the internal 20-colour palette.
#' @param title Character.  Plot title.  \code{NULL} auto-generates.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return A \code{ggplot} object (ggsave-compatible).
#'
#' @examples
#' \dontrun{
#' p <- AETHER_plot_enrichment_map(enrich, source = "REAC", top_n = 20)
#' ggplot2::ggsave("emap.png", p, width = 12, height = 10)
#'
#' # Fewer terms, tighter edges
#' p <- AETHER_plot_enrichment_map(enrich, top_n = 10, jaccard_threshold = 0.3)
#' }
#' @export
AETHER_plot_enrichment_map <- function(enrichment_result,
                                        source            = "REAC",
                                        top_n             = 20L,
                                        jaccard_threshold = 0.2,
                                        layout            = "fr",
                                        node_size_range   = c(3, 12),
                                        edge_width_range  = c(0.3, 2),
                                        max_label_width   = 40L,
                                        label_size        = 2.5,
                                        palette           = NULL,
                                        title             = NULL,
                                        verbose           = TRUE) {

  # ---------------------------------------------------------------------------
  # Package checks
  # ---------------------------------------------------------------------------
  for (pkg in c("ggraph", "igraph")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')")
  }

  # ---------------------------------------------------------------------------
  # Input normalisation
  # ---------------------------------------------------------------------------
  if (inherits(enrichment_result, "gost_enrichment")) {
    combined <- enrichment_result$combined
  } else if (is.data.frame(enrichment_result)) {
    combined <- enrichment_result
  } else {
    stop("enrichment_result must be a gost_enrichment object or data.frame")
  }

  req_cols <- c("term_id", "term_name", "source", "module", "p_value")
  missing  <- setdiff(req_cols, colnames(combined))
  if (length(missing) > 0L)
    stop("Missing required columns: ", paste(missing, collapse = ", "))

  if (!source %in% combined$source) {
    avail <- paste(sort(unique(combined$source)), collapse = ", ")
    stop("Source '", source, "' not found. Available: ", avail)
  }

  df      <- combined[combined$source == source, , drop = FALSE]
  modules <- sort(unique(df$module))
  if (length(modules) == 0L)
    stop("No modules found for source '", source, "'")

  # ---------------------------------------------------------------------------
  # Module colours
  # ---------------------------------------------------------------------------
  if (!is.null(palette)) {
    if (is.null(names(palette)))
      stop("palette must be a named vector")
    mod_colors <- palette[modules]
  } else {
    mod_colors <- stats::setNames(.aether_n_colors(length(modules)), modules)
    mod_colors <- .aether_updown_override(mod_colors)
  }

  # ---------------------------------------------------------------------------
  # Select top N per module, deduplicate
  # ---------------------------------------------------------------------------
  top_df <- do.call(rbind, lapply(modules, function(m) {
    sub <- df[df$module == m, , drop = FALSE]
    sub <- sub[order(sub$p_value), , drop = FALSE]
    utils::head(sub, as.integer(top_n))
  }))

  all_term_ids <- unique(top_df$term_id)
  n_terms      <- length(all_term_ids)
  if (verbose) cat(source, "emap:", n_terms, "unique terms across",
                   length(modules), "modules\n")

  # ---------------------------------------------------------------------------
  # Gene sets per term (union of intersections across modules)
  # ---------------------------------------------------------------------------
  has_evcodes <- "intersection" %in% colnames(top_df) &&
                 any(nchar(as.character(top_df$intersection)) > 0, na.rm = TRUE)

  if (!has_evcodes)
    stop("No 'intersection' column found. ",
         "AETHER_plot_enrichment_map requires evcodes=TRUE in APOLLO_enrich_gost.")

  gene_sets <- lapply(all_term_ids, function(tid) {
    rows  <- top_df[top_df$term_id == tid, , drop = FALSE]
    genes <- unique(unlist(strsplit(as.character(rows$intersection), ",")))
    trimws(genes[nchar(trimws(genes)) > 0])
  })
  names(gene_sets) <- all_term_ids

  # ---------------------------------------------------------------------------
  # Pairwise Jaccard → edge list
  # ---------------------------------------------------------------------------
  if (verbose) cat("[AETHER] Computing pairwise Jaccard similarity...\n")

  edge_rows <- list()
  for (i in seq_len(n_terms - 1L)) {
    a <- gene_sets[[i]]
    for (j in seq(i + 1L, n_terms)) {
      b   <- gene_sets[[j]]
      uni <- length(union(a, b))
      if (uni == 0L) next
      jacc <- length(intersect(a, b)) / uni
      if (jacc >= jaccard_threshold)
        edge_rows[[length(edge_rows) + 1L]] <-
          c(all_term_ids[[i]], all_term_ids[[j]], jacc)
    }
  }

  if (length(edge_rows) == 0L) {
    warning("No edges found at jaccard_threshold = ", jaccard_threshold,
            ". Try lowering the threshold.")
    # Still proceed — isolated-node graph is informative
    edge_df <- data.frame(from = character(0), to = character(0),
                           weight = numeric(0), stringsAsFactors = FALSE)
  } else {
    edge_mat <- do.call(rbind, edge_rows)
    edge_df  <- data.frame(from   = edge_mat[, 1],
                            to     = edge_mat[, 2],
                            weight = as.numeric(edge_mat[, 3]),
                            stringsAsFactors = FALSE)
  }

  if (verbose) cat(nrow(edge_df), "edges above Jaccard threshold",
                   jaccard_threshold, "\n")

  # ---------------------------------------------------------------------------
  # Node attributes
  # ---------------------------------------------------------------------------
  # Best -log10(p) per term (across all modules)
  best_nlp <- vapply(all_term_ids, function(tid) {
    sub <- top_df[top_df$term_id == tid, , drop = FALSE]
    -log10(min(sub$p_value))
  }, numeric(1L))

  # Primary module (lowest p)
  primary_module <- vapply(all_term_ids, function(tid) {
    sub <- top_df[top_df$term_id == tid, , drop = FALSE]
    sub$module[which.min(sub$p_value)]
  }, character(1L))

  # Number of modules that found this term
  n_mods <- vapply(all_term_ids, function(tid) {
    length(unique(top_df$module[top_df$term_id == tid]))
  }, integer(1L))

  # Truncated labels
  term_labels <- vapply(all_term_ids, function(tid) {
    nm <- top_df$term_name[top_df$term_id == tid][1L]
    if (is.na(nm) || nm == "") return(tid)
    w <- as.integer(max_label_width)
    if (nchar(nm) > w) paste0(substr(nm, 1L, w - 3L), "...") else nm
  }, character(1L))

  node_df <- data.frame(
    name           = all_term_ids,
    label          = term_labels,
    primary_module = primary_module,
    best_nlp       = best_nlp,
    n_modules      = n_mods,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Build igraph
  # ---------------------------------------------------------------------------
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)

  # ---------------------------------------------------------------------------
  # Render
  # ---------------------------------------------------------------------------
  plot_title <- title %||% paste0(source, " Enrichment Map")

  # Node border: thicker stroke for multi-module terms
  stroke_vals <- ifelse(igraph::V(g)$n_modules > 1L, 1.2, 0.4)

  if (verbose) cat("[AETHER] Rendering enrichment map (", igraph::vcount(g), "nodes,",
                   igraph::ecount(g), "edges)...\n")

  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      color = "grey30",
      alpha = 0.6
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = primary_module, size = best_nlp),
      shape  = 21,
      color  = "white",
      stroke = stroke_vals
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      size   = label_size,
      repel  = TRUE,
      family = "sans",
      color  = "grey20"
    ) +
    ggraph::scale_edge_width_continuous(
      range = edge_width_range,
      name  = "Jaccard"
    ) +
    ggplot2::scale_fill_manual(
      values = unname(mod_colors[modules]),
      labels = modules,
      name   = "Module"
    ) +
    ggplot2::scale_size_continuous(range = node_size_range, guide = "none") +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right",
      plot.margin     = ggplot2::margin(10, 10, 10, 10)
    )

  return(p)
}


###############################################################################
########### GSEA Result Visualizations (fgsea) ###########
###############################################################################

#' Word-wrap a pathway label at a given width, treating "_" as a break point
#'
#' MSigDB-style pathway names (e.g. "GOBP_CHROMOSOME_SEGREGATION") have no
#' spaces, so \code{strwrap()} on its own has nothing to wrap at. Underscores
#' are substituted for spaces, wrapped with \code{strwrap()}, then any
#' surviving space (i.e. one \code{strwrap()} didn't consume as a line break)
#' is restored to "_" -- so the only underscore actually lost is the one at
#' each chosen break point, same as how a normal word-wrap drops the space it
#' breaks on.
#' @noRd
.wrap_gsea_label <- function(label, width) {
  if (nchar(label) <= width) return(label)
  spaced  <- gsub("_", " ", label, fixed = TRUE)
  wrapped <- paste(strwrap(spaced, width = width), collapse = "\n")
  gsub(" ", "_", wrapped, fixed = TRUE)
}


#' NES Dotplot for GSEA Results
#'
#' @description Creates a dotplot showing the top enriched and depleted gene
#' sets from a GSEA result (\code{APOLLO_gsea()}), ranked by normalized
#' enrichment score (NES).
#'
#' @param gsea_result An \code{apollo_gsea} object from \code{APOLLO_gsea()},
#'   or a data.frame with columns \code{pathway}, \code{NES}, \code{padj},
#'   \code{size} (typically \code{$results} or \code{$significant}).
#' @param top_n Integer. Total number of gene sets to show, split between top
#'   enriched (NES > 0) and top depleted (NES < 0) by padj. Default: 20.
#' @param title Character. Plot title. Default: "GSEA".
#' @param font_size Numeric. Base font size for pathway labels. Default: 8.
#' @param max_label_length Integer. Pathway names longer than this are
#'   word-wrapped onto multiple lines (never truncated/"..."-ed). Default: 55.
#'
#' @return A ggplot object, or \code{NULL} (with a warning) if there are no
#'   results to plot.
#'
#' @details
#' Dot size = number of genes in the leading edge (\code{size} column); dot
#' color = -log10(padj). A dashed vertical line marks NES = 0.
#'
#' Long pathway names are wrapped (via \code{strwrap()}), not truncated --
#' MSigDB-style names (e.g. \code{"GOBP_CHROMOSOME_SEGREGATION"}) have no
#' spaces to wrap on, so underscores are treated as break points the same way
#' spaces would be, with the underscore at the chosen break consumed by the
#' line break (matching how a normal word-wrap drops the space it breaks on).
#'
#' @examples
#' \dontrun{
#' gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))
#' p <- AETHER_plot_gsea_dotplot(gsea_result, top_n = 30)
#' }
#' @export
AETHER_plot_gsea_dotplot <- function(gsea_result,
                                      top_n = 20,
                                      title = "GSEA",
                                      font_size = 8,
                                      max_label_length = 55) {

  if (inherits(gsea_result, "apollo_gsea")) {
    df <- gsea_result$results
  } else if (is.data.frame(gsea_result)) {
    df <- gsea_result
  } else {
    stop("gsea_result must be an apollo_gsea object or data.frame")
  }

  required_cols <- c("pathway", "NES", "padj", "size")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df[!is.na(df$padj) & !is.na(df$NES), , drop = FALSE]
  if (nrow(df) == 0) {
    warning("No GSEA results to plot")
    return(NULL)
  }

  n_each  <- ceiling(top_n / 2)
  pos_df  <- df[df$NES > 0, , drop = FALSE]
  neg_df  <- df[df$NES < 0, , drop = FALSE]
  pos_top <- head(pos_df[order(pos_df$padj), ], n_each)
  neg_top <- head(neg_df[order(neg_df$padj), ], n_each)
  plot_df <- rbind(pos_top, neg_top)

  if (nrow(plot_df) == 0) {
    warning("No gene sets remaining after enriched/depleted split")
    return(NULL)
  }

  plot_df$label <- vapply(plot_df$pathway, .wrap_gsea_label, character(1), width = max_label_length)
  plot_df       <- plot_df[order(plot_df$NES), , drop = FALSE]
  plot_df$label <- factor(plot_df$label, levels = unique(plot_df$label))

  p <- ggplot(plot_df, aes(x = NES, y = label, size = size, color = -log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_color_gradient(low = "#deebf7", high = "#08519c", name = "-log10(padj)") +
    scale_size_continuous(name = "Gene set\noverlap", range = c(2, 8)) +
    labs(
      title = title,
      subtitle = paste0("Top ", nrow(pos_top), " enriched / ",
                        nrow(neg_top), " depleted (by padj)"),
      x = "NES", y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = font_size + 4),
      plot.subtitle = element_text(size = font_size + 1, color = "gray40"),
      axis.text.y = element_text(size = font_size),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}


#' Running-Score (Mountain) Plot for a Single GSEA Pathway
#'
#' @description Wraps \code{fgsea::plotEnrichment()} to draw the running
#' enrichment score curve for one gene set from a GSEA result
#' (\code{APOLLO_gsea()}).
#'
#' @param gsea_result An \code{apollo_gsea} object from \code{APOLLO_gsea()}.
#'   Must retain \code{$gene_sets} and \code{$ranked_genes} (both included by
#'   default).
#' @param pathway Character. Name of the gene set to plot (must be present in
#'   \code{gsea_result$gene_sets}).
#' @param title Character or \code{NULL}. Plot title. \code{NULL} (default)
#'   uses the pathway name.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' gsea_result <- APOLLO_gsea(ranked, collection = "H")
#' p <- AETHER_plot_gsea_enrichment(gsea_result, "HALLMARK_INFLAMMATORY_RESPONSE")
#' }
#' @export
AETHER_plot_gsea_enrichment <- function(gsea_result, pathway, title = NULL) {

  if (!inherits(gsea_result, "apollo_gsea")) {
    stop("gsea_result must be an apollo_gsea object from APOLLO_gsea()")
  }

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
  }

  if (!pathway %in% names(gsea_result$gene_sets)) {
    stop("pathway '", pathway, "' not found in gsea_result$gene_sets. ",
         "Use one of the names in gsea_result$results$pathway that passed ",
         "min_size/max_size filtering.")
  }

  p <- fgsea::plotEnrichment(
    pathway = gsea_result$gene_sets[[pathway]],
    stats   = gsea_result$ranked_genes
  )

  p <- p + ggtitle(title %||% pathway) + theme_minimal()

  return(p)
}

