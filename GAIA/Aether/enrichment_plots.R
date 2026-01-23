library(ggplot2)

###############################################################################
########### Enrichment Result Visualizations ###########
###############################################################################


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

  # Set dataset order (preserve input order)
  combined_df$Dataset <- factor(combined_df$Dataset, levels = rev(names(annotated_list)))

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
  y_var <- if (show_percentage) "Percentage" else "Count"

  p <- ggplot(combined_df, aes(x = .data[[y_var]], y = Dataset, fill = Feature)) +
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
    guides(fill = guide_legend(nrow = 1))

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

  # Add count annotations on right side
  totals <- aggregate(Count ~ Dataset, combined_df, sum)
  totals$label <- paste0("n=", format(totals$Count, big.mark = ","))

  if (show_percentage) {
    p <- p +
      geom_text(
        data = totals,
        aes(x = 102, y = Dataset, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3.5,
        color = "gray30"
      ) +
      coord_cartesian(xlim = c(0, 115), clip = "off")
  }

  return(p)
}
