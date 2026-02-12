# GAIA/Aether/wgcna_plots.R
# Visualization functions for WGCNA analysis results
#
# Works with objects from ARTEMIS WGCNA functions.

library(ggplot2)
library(WGCNA)

# ==============================================================================
# POWER SELECTION PLOT
# ==============================================================================

#' Plot soft threshold power selection diagnostics
#'
#' Creates diagnostic plots for selecting the soft threshold power in WGCNA.
#' Shows scale-free topology fit and mean connectivity across powers.
#'
#' @param power_result A wgcna_power object from ARTEMIS_wgcna_pick_power(),
#'   or a data.frame of fit indices with columns: Power, SFT.R.sq, slope, mean.k.
#' @param r2_cutoff R-squared threshold line to display. Default: 0.80.
#' @param highlight_power Power value to highlight. If NULL and power_result is
#'   a wgcna_power object, uses the selected power.
#' @param return_plots Logical. If TRUE, returns list of ggplot objects instead
#'   of plotting. Default: FALSE.
#'
#' @return If return_plots = TRUE, returns list with scale_free and connectivity plots.
#'   Otherwise, displays plots and returns invisible NULL.
#'
#' @examples
#' power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
#' AETHER_plot_wgcna_power(power_result)
#'
#' @export
AETHER_plot_wgcna_power <- function(power_result,
                                     r2_cutoff = 0.80,
                                     highlight_power = NULL,
                                     return_plots = FALSE) {

  # Extract fit indices

  if (inherits(power_result, "wgcna_power")) {
    fit_indices <- power_result$fit_indices
    if (is.null(highlight_power)) {
      highlight_power <- power_result$power
    }
  } else {
    fit_indices <- as.data.frame(power_result)
  }

  # Calculate signed R² if not present
  if (!"signed_R2" %in% colnames(fit_indices)) {
    fit_indices$signed_R2 <- -sign(fit_indices[, 3]) * fit_indices[, 2]
  }

  # Prepare highlight data
  highlight_df <- NULL
  if (!is.null(highlight_power) && !is.na(highlight_power)) {
    idx <- which(fit_indices$Power == highlight_power)
    if (length(idx) > 0) {
      highlight_df <- fit_indices[idx, , drop = FALSE]
    }
  }

  # Plot 1: Scale-free topology fit
  p1 <- ggplot(fit_indices, aes(x = Power, y = signed_R2)) +
    geom_text(aes(label = Power), color = "red", size = 4) +
    geom_hline(yintercept = r2_cutoff, linetype = "dashed", color = "red") +
    labs(
      x = "Soft Threshold (power)",
      y = expression("Scale Free Topology Model Fit (signed R"^2*")"),
      title = "Scale Independence"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  if (!is.null(highlight_df)) {
    p1 <- p1 + geom_point(data = highlight_df, aes(x = Power, y = signed_R2),
                          size = 5, color = "blue", shape = 19)
  }

  # Plot 2: Mean connectivity
  p2 <- ggplot(fit_indices, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), color = "red", size = 4) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  if (!is.null(highlight_df)) {
    p2 <- p2 + geom_point(data = highlight_df, aes(x = Power, y = mean.k.),
                          size = 5, color = "blue", shape = 19)
  }

  if (return_plots) {
    return(list(scale_free = p1, connectivity = p2))
  }

  # Display side by side
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  invisible(NULL)
}


# ==============================================================================
# MODULE DENDROGRAM
# ==============================================================================

#' Plot module dendrogram with colors
#'
#' Creates a dendrogram of genes colored by module assignment.
#'
#' @param modules A wgcna_modules object from ARTEMIS_wgcna_detect_modules().
#' @param block Which block's dendrogram to plot (for large datasets split into
#'   blocks). Default: 1.
#' @param main Plot title. Default: "Gene dendrogram and module colors".
#' @param ... Additional arguments passed to plotDendroAndColors().
#'
#' @return Invisible NULL. Plots to current device.
#'
#' @examples
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
#' AETHER_plot_wgcna_dendrogram(modules)
#'
#' @export
AETHER_plot_wgcna_dendrogram <- function(modules,
                                          block = 1,
                                          main = "Gene dendrogram and module colors",
                                          ...) {

  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules object")
  }

  if (is.null(modules$dendrograms) || length(modules$dendrograms) < block) {
    stop("Dendrogram not available for block ", block)
  }

  # Get display colors for genes in this block
  block_genes <- modules$net$blockGenes[[block]]
  block_module_names <- modules$module_names[block_genes]
  block_colors <- modules$module_colors[block_module_names]

  plotDendroAndColors(
    modules$dendrograms[[block]],
    unname(block_colors),
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = main,
    ...
  )

  invisible(NULL)
}


# ==============================================================================
# MODULE-TRAIT HEATMAP
# ==============================================================================

#' Plot module-trait correlation heatmap
#'
#' Creates a heatmap showing correlations between module eigengenes and traits,
#' with correlation values and p-values displayed in cells.
#'
#' @param trait_cor A wgcna_trait_cor object from ARTEMIS_wgcna_module_traits().
#' @param use_padj Logical. Use adjusted p-values for display. Default: FALSE
#'   (uses raw p-values like original WGCNA).
#' @param show_values Logical. Show correlation and p-value in cells. Default: TRUE.
#' @param colors Color palette for heatmap. Default: blueWhiteRed(50).
#' @param main Plot title. Default: "Module-trait relationships".
#' @param cex.text Text size for cell values. Default: 0.7.
#' @param ... Additional arguments passed to labeledHeatmap().
#'
#' @return Invisible NULL. Plots to current device.
#'
#' @examples
#' trait_cor <- ARTEMIS_wgcna_module_traits(modules, wgcna_data)
#' AETHER_plot_module_trait_heatmap(trait_cor)
#'
#' @export
AETHER_plot_module_trait_heatmap <- function(trait_cor,
                                              use_padj = FALSE,
                                              show_values = TRUE,
                                              colors = NULL,
                                              main = "Module-trait relationships",
                                              cex.text = 0.7,
                                              ...) {

  if (!inherits(trait_cor, "wgcna_trait_cor")) {
    stop("trait_cor must be a wgcna_trait_cor object")
  }

  cor_matrix <- trait_cor$cor_matrix
  pvalue_matrix <- if (use_padj) trait_cor$padj_matrix else trait_cor$pvalue_matrix

  # Default colors
  if (is.null(colors)) {
    colors <- blueWhiteRed(50)
  }

  # Create text matrix
  if (show_values) {
    text_matrix <- paste(
      signif(cor_matrix, 2), "\n(",
      signif(pvalue_matrix, 1), ")",
      sep = ""
    )
    dim(text_matrix) <- dim(cor_matrix)
  } else {
    text_matrix <- NULL
  }

  # Set margins
  par(mar = c(6, 10, 3, 3))

  labeledHeatmap(
    Matrix = cor_matrix,
    xLabels = colnames(cor_matrix),
    yLabels = rownames(cor_matrix),
    ySymbols = rownames(cor_matrix),
    colorLabels = FALSE,
    colors = colors,
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text = cex.text,
    zlim = c(-1, 1),
    main = main,
    ...
  )

  invisible(NULL)
}


#' Plot module-trait heatmap using ggplot2
#'
#' Alternative ggplot2-based heatmap for module-trait correlations.
#' More customizable than the base R version.
#'
#' @param trait_cor A wgcna_trait_cor object from ARTEMIS_wgcna_module_traits().
#' @param use_padj Logical. Use adjusted p-values. Default: FALSE.
#' @param sig_threshold P-value threshold for significance stars. Default: 0.05.
#' @param show_values What to show in cells: "cor", "pval", "both", or "stars".
#'   Default: "both".
#' @param pval_format Format for p-values: "scientific" (e.g., 1.2e-03) or
#'   "decimal" (e.g., 0.001). Default: "decimal".
#' @param pval_digits Integer. Decimal places for p-values (decimal format only).
#'   Default: 3.
#' @param colors Color palette. Default: blue-white-red gradient.
#' @param title Plot title.
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_module_trait_heatmap_gg <- function(trait_cor,
                                                 use_padj = FALSE,
                                                 sig_threshold = 0.05,
                                                 show_values = "both",
                                                 pval_format = "decimal",
                                                 pval_digits = 3,
                                                 colors = NULL,
                                                 title = "Module-trait relationships") {

  if (!inherits(trait_cor, "wgcna_trait_cor")) {
    stop("trait_cor must be a wgcna_trait_cor object")
  }

  cor_matrix <- trait_cor$cor_matrix
  pvalue_matrix <- if (use_padj) trait_cor$padj_matrix else trait_cor$pvalue_matrix

  # Convert to long format
  cor_df <- expand.grid(
    Module = rownames(cor_matrix),
    Trait = colnames(cor_matrix),
    stringsAsFactors = FALSE
  )
  cor_df$Correlation <- as.vector(cor_matrix)
  cor_df$Pvalue <- as.vector(pvalue_matrix)
  cor_df$Significant <- cor_df$Pvalue < sig_threshold

  # Format p-values based on user preference
  if (pval_format == "decimal") {
    pval_fmt <- paste0("%.", pval_digits, "f")
    pval_str <- sprintf(pval_fmt, cor_df$Pvalue)
  } else {
    pval_str <- sprintf("%.1e", cor_df$Pvalue)
  }

  # Create label
  cor_df$Label <- switch(
    show_values,
    "cor" = sprintf("%.2f", cor_df$Correlation),
    "pval" = pval_str,
    "both" = paste0(sprintf("%.2f", cor_df$Correlation), "\n(", pval_str, ")"),
    "stars" = ifelse(cor_df$Pvalue < 0.001, "***",
                     ifelse(cor_df$Pvalue < 0.01, "**",
                            ifelse(cor_df$Pvalue < 0.05, "*", "")))
  )

  # Default colors
  if (is.null(colors)) {
    colors <- c("#2166AC", "white", "#B2182B")  # Blue-white-red
  }

  p <- ggplot(cor_df, aes(x = Trait, y = Module, fill = Correlation)) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = Label), size = 3) +
    scale_fill_gradient2(
      low = colors[1], mid = colors[2], high = colors[3],
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )

  return(p)
}


# ==============================================================================
# SAMPLE DENDROGRAM WITH TRAITS
# ==============================================================================

#' Plot sample dendrogram with trait colors
#'
#' Creates a sample clustering dendrogram with trait values shown as colored
#' bars underneath.
#'
#' @param cluster_result A wgcna_cluster object from ARTEMIS_wgcna_cluster_samples(),
#'   or an hclust object.
#' @param traits Data.frame of traits (samples as rows) or wgcna_data object.
#' @param main Plot title.
#' @param ... Additional arguments passed to plotDendroAndColors().
#'
#' @return Invisible NULL. Plots to current device.
#'
#' @examples
#' clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data)
#' AETHER_plot_sample_dendrogram(clust, wgcna_data)
#'
#' @export
AETHER_plot_sample_dendrogram <- function(cluster_result,
                                           traits = NULL,
                                           main = "Sample dendrogram and trait heatmap",
                                           ...) {

  # Extract dendrogram

  if (inherits(cluster_result, "wgcna_cluster")) {
    sample_tree <- cluster_result$sample_tree
    if (is.null(traits) && !is.null(cluster_result$datTraits)) {
      traits <- cluster_result$datTraits
    }
  } else if (inherits(cluster_result, "hclust")) {
    sample_tree <- cluster_result
  } else {
    stop("cluster_result must be a wgcna_cluster or hclust object")
  }

  # Extract traits
  if (inherits(traits, "wgcna_data")) {
    traits <- traits$datTraits
  }

  if (is.null(traits)) {
    # Plot without trait colors
    plot(sample_tree, main = main, sub = "", xlab = "")
  } else {
    # Ensure traits are numeric
    traits_numeric <- as.data.frame(lapply(traits, function(x) {
      if (is.numeric(x)) x else as.numeric(as.factor(x))
    }))
    rownames(traits_numeric) <- rownames(traits)

    # Convert to colors
    trait_colors <- numbers2colors(traits_numeric, signed = FALSE)

    plotDendroAndColors(
      sample_tree,
      trait_colors,
      groupLabels = colnames(traits),
      main = main,
      ...
    )
  }

  invisible(NULL)
}


# ==============================================================================
# GENE SIGNIFICANCE PLOTS
# ==============================================================================

#' Plot gene significance vs module membership
#'
#' Creates a scatter plot of gene significance (GS) vs module membership (MM)
#' for genes in a specific module. Useful for identifying hub genes.
#'
#' @param gene_sig A wgcna_gene_sig object from ARTEMIS_wgcna_gene_significance().
#' @param module Module color to plot.
#' @param trait_name Trait name for gene significance. If NULL, uses first trait.
#' @param mm_threshold Module membership threshold for highlighting. Default: 0.8.
#' @param gs_threshold Gene significance threshold for highlighting. Default: 0.2.
#' @param n_label Number of top genes to label. Default: 10.
#' @param module_color Hex color for hub gene points. If NULL, uses default green.
#' @param title Plot title. If NULL, auto-generated.
#'
#' @return A ggplot object.
#'
#' @examples
#' gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)
#' AETHER_plot_gs_vs_mm(gene_sig, module = "blue", trait_name = "severity")
#'
#' @export
AETHER_plot_gs_vs_mm <- function(gene_sig,
                                  module,
                                  trait_name = NULL,
                                  mm_threshold = 0.8,
                                  gs_threshold = 0.2,
                                  n_label = 10,
                                  module_color = NULL,
                                  title = NULL) {

  if (!inherits(gene_sig, "wgcna_gene_sig")) {
    stop("gene_sig must be a wgcna_gene_sig object")
  }

  # Get trait name
  if (is.null(trait_name)) {
    trait_name <- gene_sig$trait_names[1]
  }

  # Get data for this module
  gene_info <- gene_sig$gene_info
  module_genes <- gene_info[gene_info$module == module, ]

  if (nrow(module_genes) == 0) {
    stop("No genes found in module: ", module)
  }

  # Column names
  gs_col <- paste0("GS.", trait_name)
  mm_col <- paste0("MM.", module)

  if (!gs_col %in% colnames(module_genes)) {
    stop("Trait '", trait_name, "' not found")
  }
  if (!mm_col %in% colnames(module_genes)) {
    stop("MM column for module '", module, "' not found")
  }

  # Prepare data
  plot_df <- data.frame(
    gene = module_genes$gene,
    GS = abs(module_genes[[gs_col]]),
    MM = abs(module_genes[[mm_col]]),
    stringsAsFactors = FALSE
  )

  # Mark hub genes
  plot_df$is_hub <- plot_df$MM >= mm_threshold & plot_df$GS >= gs_threshold

  # Get genes to label (top by GS * MM)
  plot_df$score <- plot_df$GS * plot_df$MM
  plot_df <- plot_df[order(-plot_df$score), ]
  genes_to_label <- head(plot_df$gene, n_label)
  plot_df$label <- ifelse(plot_df$gene %in% genes_to_label, plot_df$gene, "")

  # Calculate correlation
  cor_val <- cor(plot_df$GS, plot_df$MM, use = "complete.obs")

  # Title
  if (is.null(title)) {
    title <- sprintf("Module: %s | Trait: %s | cor = %.2f", module, trait_name, cor_val)
  }

  # Resolve display color for hub points
  if (is.null(module_color)) {
    module_color <- "#4DAF4A"  # default green
  }

  # Plot
  p <- ggplot(plot_df, aes(x = MM, y = GS)) +
    geom_point(aes(color = is_hub), alpha = 0.6, size = 2) +
    geom_vline(xintercept = mm_threshold, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = gs_threshold, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = module_color),
                       guide = "none") +
    labs(
      x = paste0("Module Membership (|MM.", module, "|)"),
      y = paste0("Gene Significance (|GS.", trait_name, "|)"),
      title = title
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  # Add labels
  if (n_label > 0) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20
    )
  }

  return(p)
}


# ==============================================================================
# HUB GENE VISUALIZATION
# ==============================================================================

#' Plot hub genes across modules
#'
#' Bar plot showing the number of hub genes per module, with top hub gene labeled.
#'
#' @param hubs A wgcna_hubs object from ARTEMIS_wgcna_hub_genes().
#' @param top_n Number of modules to show. Default: NULL (all).
#' @param show_top_gene Logical. Label bars with top hub gene name. Default: TRUE.
#' @param module_colors Named character vector of module colors (module_name -> hex).
#'   If NULL, generates colors automatically.
#' @param title Plot title.
#'
#' @return A ggplot object.
#'
#' @export
AETHER_plot_hub_summary <- function(hubs,
                                     top_n = NULL,
                                     show_top_gene = TRUE,
                                     module_colors = NULL,
                                     title = "Hub genes per module") {

  if (!inherits(hubs, "wgcna_hubs")) {
    stop("hubs must be a wgcna_hubs object")
  }

  hub_summary <- hubs$hub_summary
  hub_summary <- hub_summary[hub_summary$n_hubs > 0, ]

  if (!is.null(top_n)) {
    hub_summary <- head(hub_summary, top_n)
  }

  # Order by n_hubs
  hub_summary$module <- factor(hub_summary$module,
                               levels = hub_summary$module[order(hub_summary$n_hubs)])

  # Resolve module colors
  if (is.null(module_colors)) {
    # Generate colors for each module present
    mod_levels <- levels(hub_summary$module)
    module_colors <- setNames(
      grDevices::hcl.colors(length(mod_levels), palette = "Dark 3"),
      mod_levels
    )
  }

  p <- ggplot(hub_summary, aes(x = module, y = n_hubs, fill = module)) +
    geom_col() +
    scale_fill_manual(values = module_colors) +
    coord_flip() +
    labs(x = NULL, y = "Number of hub genes", title = title) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )

  if (show_top_gene) {
    p <- p + geom_text(aes(label = top_hub), hjust = -0.1, size = 3)
  }

  return(p)
}
