# GAIA/Aether/timeseries_plots.R
# Visualization functions for time series analysis and PART clustering results
#
# Works with objects from ARTEMIS_part() and ARTEMIS_timeseries_*() functions.



# ==============================================================================
# PART HEATMAP (ComplexHeatmap - TiSA style)
# ==============================================================================

#' Plot PART Clustering Heatmap
#'
#' Creates a heatmap of PART clustered data using ComplexHeatmap. Matches the
#' TiSA package visualization style with:
#' \itemize{
#'   \item Row annotation showing cluster color blocks
#'   \item Top annotation with group blocks and (optionally) timepoint colors
#'   \item Column splitting by group
#'   \item Row splitting by cluster
#'   \item Custom legends for Z-score, clusters, groups, and timepoints
#' }
#'
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part()}.
#' @param sample_info Data.frame with sample metadata. Rownames must match
#'   colnames of the PART data matrix. Required column: group.
#'   Optional: timepoint (for time series data).
#' @param group_col Character. Column name in sample_info for group. Default: "group".
#' @param time_col Character or NULL. Column name in sample_info for timepoint.
#'   Default: NULL. If NULL or column not found, timepoint annotation is skipped
#'   (useful for non-time-series data).
#' @param group_colors Named vector of colors for groups. Names should match
#'   group levels. Default: NULL (auto-generated).
#' @param time_colors Named vector of colors for timepoints. Default: NULL
#'   (uses yellow-orange-red gradient).
#' @param show_row_names Logical. Show gene names. Default: FALSE.
#' @param row_names_side Character. Side for row names: "left" or "right".
#'   Default: "right".
#' @param show_column_labels Logical. Whether to show column (sample) labels
#'   on the heatmap. Default: TRUE.
#' @param save_path Character. If provided, saves the heatmap to this path
#'   (supports .png, .pdf, .svg). Default: NULL (returns plot object).
#' @param width Numeric. Width in inches for saved plot. Default: 10.
#' @param height Numeric. Height in inches for saved plot. Default: 12.
#'
#' @return If save_path is NULL, returns the ComplexHeatmap draw object.
#'   Otherwise saves to file and returns invisibly.
#'
#' @details
#' This function is ported from TiSA's PART_heat_map() function and uses
#' ComplexHeatmap for high-quality visualization. The heatmap displays:
#' \itemize{
#'   \item Genes as rows, ordered and split by PART cluster
#'   \item Samples as columns, split by group and ordered by timepoint (if available)
#'   \item Z-score color scale (blue-white-red)
#'   \item Cluster color blocks on the left
#'   \item Group and timepoint annotations on top (timepoint optional)
#' }
#'
#' For non-time-series data (standard DEA), set \code{time_col = NULL} or leave
#' it as default. The function will create a heatmap with only group annotations.
#'
#' @examples
#' \dontrun{
#' # Basic usage with time series data
#' AETHER_plot_part_heatmap(part_result, sample_info, time_col = "timepoint")
#'
#' # Non-time-series data (no timepoint annotation)
#' AETHER_plot_part_heatmap(part_result, sample_info)
#'
#' # With custom colors
#' AETHER_plot_part_heatmap(
#'   part_result, sample_info,
#'   group_colors = c("Control" = "blue", "Treatment" = "red")
#' )
#'
#' # Save to file
#' AETHER_plot_part_heatmap(part_result, sample_info, save_path = "heatmap.png")
#'
#' }
#' @export
AETHER_plot_part_heatmap <- function(part_result,
                                      sample_info,
                                      group_col = "group",
                                      time_col = NULL,
                                      group_colors = NULL,
                                      time_colors = NULL,
                                      show_row_names = FALSE,
                                      show_column_labels = TRUE,
                                      row_names_side = "right",
                                      save_path = NULL,
                                      width = 10,
                                      height = 12) {

  # --- Validate inputs ---
  if (!inherits(part_result, "artemis_part")) {
    stop("'part_result' must be an artemis_part object from ARTEMIS_part()")
  }

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install with:\n",
         "  install.packages('circlize')")
  }

  requireNamespace("grid", quietly = TRUE)

  if (!group_col %in% colnames(sample_info)) {
    stop("'", group_col, "' column not found in sample_info")
  }

  # time_col is optional (NULL for non-time-series data)
  has_timepoints <- !is.null(time_col) && time_col %in% colnames(sample_info)
  if (!is.null(time_col) && !has_timepoints) {
    warning("'", time_col, "' column not found in sample_info. ",
            "Skipping timepoint annotation.")
  }

  # --- Prepare matrix and align samples ---
  mat <- part_result$data
  common_samples <- intersect(colnames(mat), rownames(sample_info))
  if (length(common_samples) == 0) {
    stop("No matching samples between part_result and sample_info rownames")
  }

  # Order samples by group (then timepoint if available)
  sample_info <- sample_info[common_samples, , drop = FALSE]
  if (has_timepoints) {
    sample_order <- order(sample_info[[group_col]], sample_info[[time_col]])
  } else {
    sample_order <- order(sample_info[[group_col]])
  }
  sample_info <- sample_info[sample_order, , drop = FALSE]
  mat <- mat[, rownames(sample_info), drop = FALSE]

  # --- Extract cluster info ---
  clusters <- part_result$cluster_map$cluster[match(rownames(mat), part_result$cluster_map$gene)]
  cluster_levels <- unique(clusters)  # Already ordered by PART
  cluster_colors <- part_result$cluster_colors

  # --- Generate group colors if not provided ---
  groups <- sample_info[[group_col]]
  if (!is.null(group_colors)) {
    # Use order from group_colors as authoritative group ordering
    group_levels <- names(group_colors)
    # Keep only groups actually present in the data
    group_levels <- group_levels[group_levels %in% unique(groups)]
  } else {
    group_levels <- unique(groups)
    group_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                   "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    group_colors <- setNames(group_pal[seq_along(group_levels)], group_levels)
  }

  # --- Generate timepoint colors if not provided (only if timepoints exist) ---
  if (has_timepoints) {
    timepoints <- sample_info[[time_col]]
    time_levels <- unique(timepoints)
    if (is.null(time_colors)) {
      # Yellow-orange-red gradient (TiSA style)
      time_pal <- c('#ffeda0', '#fed976', '#feb24c', '#fd8d3c',
                    '#fc4e2a', '#e31a1c', '#bd0026', '#800026')
      if (length(time_levels) <= length(time_pal)) {
        time_colors <- setNames(time_pal[seq_along(time_levels)], time_levels)
      } else {
        time_colors <- setNames(
          grDevices::colorRampPalette(time_pal)(length(time_levels)),
          time_levels
        )
      }
    }
  }

  # --- Create column split factor ---
  # Split by group, with samples ordered by timepoint within each group
  col_split <- factor(groups, levels = group_levels)

  # --- Create row annotation (cluster blocks) ---
  row_annot <- ComplexHeatmap::rowAnnotation(
    cluster = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = cluster_colors[cluster_levels]),
      labels = NULL
    ),
    show_annotation_name = FALSE
  )

  # --- Create top annotation ---
  # With timepoints: group blocks + timepoint colors
  # Without timepoints: group blocks only
  # Build anno_block arguments (conditionally include labels)
  block_args <- list(
    gp = grid::gpar(fill = group_colors[group_levels])
  )
  if (show_column_labels) {
    block_args$labels <- group_levels
    block_args$labels_gp <- grid::gpar(fontsize = 10, fontface = "bold")
  }
  group_block <- do.call(ComplexHeatmap::anno_block, block_args)

  if (has_timepoints) {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      group = group_block,
      timepoint = timepoints,
      col = list(timepoint = time_colors),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        timepoint = list(title = "Timepoint")
      )
    )
  } else {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      group = group_block,
      show_annotation_name = FALSE
    )
  }

  # --- Calculate Z-score range for color scale ---
  target_z <- ceiling(max(abs(range(mat, na.rm = TRUE))))

  # --- Create legends ---
  cluster_lgd <- ComplexHeatmap::Legend(
    labels = cluster_levels,
    title = "Cluster",
    legend_gp = grid::gpar(fill = cluster_colors[cluster_levels])
  )

  group_lgd <- ComplexHeatmap::Legend(
    labels = names(group_colors),
    title = "Group",
    legend_gp = grid::gpar(fill = unname(group_colors))
  )

  # --- Build heatmap ---
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Z-score",
    col = circlize::colorRamp2(
      c(-target_z, 0, target_z),
      c("blue", "white", "red")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(clusters, levels = cluster_levels),
    column_split = col_split,
    row_gap = grid::unit(0.5, "mm"),
    column_gap = grid::unit(2, "mm"),
    left_annotation = row_annot,
    top_annotation = top_annot,
    show_row_names = show_row_names,
    show_column_names = FALSE,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = 6),
    row_title = NULL,
    column_title = NULL,
    show_heatmap_legend = TRUE,
    border = FALSE,
    use_raster = nrow(mat) > 500  # Use raster for large matrices
  )

  # --- Draw or save ---
  if (is.null(save_path)) {
    # Return drawable object
    # Need to wrap in pdf(NULL) to capture the draw object
    grDevices::pdf(NULL)
    ht_drawn <- ComplexHeatmap::draw(
      ht,
      annotation_legend_list = list(group_lgd, cluster_lgd),
      merge_legend = TRUE
    )
    grDevices::dev.off()
    return(ht_drawn)
  } else {
    # Save to file
    ext <- tolower(tools::file_ext(save_path))
    if (ext == "png") {
      grDevices::png(save_path, width = width, height = height, units = "in", res = 150)
    } else if (ext == "pdf") {
      grDevices::pdf(save_path, width = width, height = height)
    } else if (ext == "svg") {
      grDevices::svg(save_path, width = width, height = height)
    } else {
      stop("Unsupported file format: ", ext, ". Use .png, .pdf, or .svg")
    }

    ComplexHeatmap::draw(
      ht,
      annotation_legend_list = list(group_lgd, cluster_lgd),
      merge_legend = TRUE
    )
    grDevices::dev.off()

    if (ext != "pdf") {
      message("Heatmap saved to: ", save_path)
    }
    return(invisible(NULL))
  }
}


# ==============================================================================
# CLUSTER TRAJECTORY PLOT (TiSA style)
# ==============================================================================

#' Plot Cluster Expression Trajectories
#'
#' Creates faceted line plots showing expression trajectories per PART cluster.
#' Matches TiSA's visualization style where each panel shows one cluster for one
#' group, with individual gene lines (thin, colored) and a mean trajectory line
#' (thick, gray).
#'
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part()}.
#' @param sample_info Data.frame with sample metadata. Must have rownames
#'   matching colnames of the PART data matrix, and columns for timepoint and
#'   group.
#' @param time_col Character. Column in sample_info for timepoint. Default:
#'   "timepoint".
#' @param group_col Character. Column in sample_info for group. Default:
#'   "group".
#' @param norm_counts Optional matrix of normalized counts (genes x samples).
#'   If NULL, uses the z-scored data from part_result. Default: NULL.
#' @param clusters Character vector of clusters to plot (e.g., c("C1", "C3")).
#'   NULL = all non-outlier clusters. Default: NULL.
#' @param scale_features Logical. Apply scale feature sum transformation
#'   (value / row_sum) as in TiSA. Recommended when using norm_counts.
#'   Default: FALSE.
#' @param ncol Integer. Number of columns in facet grid. Default: 4.
#' @param colors Named character vector of colors for groups. NULL = auto.
#'   Default: NULL.
#' @param alpha Numeric. Transparency for individual gene lines. Default: 0.4.
#' @param title_size Numeric. Font size for facet titles. Default: 10.
#'
#' @return A ggplot object.
#'
#' @details
#' This function is ported from TiSA's \code{plot_cluster_traj()} function.
#' Key features:
#' \itemize{
#'   \item Faceted by cluster AND group (e.g., "C1 - 50 genes | IgM")
#'   \item Individual gene trajectories as thin colored lines
#'   \item Mean cluster trajectory as thick gray line
#'   \item X-axis shows timepoints, Y-axis shows expression
#' }
#'
#' @examples
#' \dontrun{
#' AETHER_plot_cluster_trajectories(part_result, sample_info)
#'
#' }
#' @export
AETHER_plot_cluster_trajectories <- function(part_result,
                                              sample_info,
                                              time_col = "timepoint",
                                              group_col = "group",
                                              norm_counts = NULL,
                                              clusters = NULL,
                                              scale_features = FALSE,
                                              ncol = 4,
                                              colors = NULL,
                                              alpha = 0.4,
                                              title_size = 10) {

  if (!inherits(part_result, "artemis_part")) {
    stop("'part_result' must be an artemis_part object")
  }

  if (!time_col %in% colnames(sample_info)) {
    stop("'", time_col, "' not found in sample_info columns")
  }
  if (!group_col %in% colnames(sample_info)) {
    stop("'", group_col, "' not found in sample_info columns")
  }

  # Use z-scored data from PART result or provided norm_counts
  if (is.null(norm_counts)) {
    mat <- part_result$data
  } else {
    # Subset norm_counts to genes in cluster_map
    mat <- norm_counts[rownames(norm_counts) %in% part_result$cluster_map$gene, ,
                       drop = FALSE]
  }

  cmap <- part_result$cluster_map

  # Align samples
  common_samples <- intersect(colnames(mat), rownames(sample_info))
  if (length(common_samples) == 0) {
    stop("No matching sample names between PART data and sample_info")
  }
  mat <- mat[, common_samples, drop = FALSE]
  sample_info <- sample_info[common_samples, , drop = FALSE]

  # Filter clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(cmap$cluster))
    clusters <- clusters[clusters != "C0"]  # Exclude outliers
  }

  # Get genes in selected clusters
  gene_cluster <- setNames(cmap$cluster, cmap$gene)
  genes_in_clusters <- names(gene_cluster)[gene_cluster %in% clusters]
  mat <- mat[genes_in_clusters, , drop = FALSE]

  # --- Calculate trajectory data (TiSA style) ---
  # For each group-timepoint, calculate mean expression across replicates
  groups <- unique(sample_info[[group_col]])
  timepoints <- sort(unique(sample_info[[time_col]]))

  # Build group-timepoint means matrix
  gt_means <- list()
  for (grp in groups) {
    for (tp in timepoints) {
      samples_gt <- rownames(sample_info)[
        sample_info[[group_col]] == grp & sample_info[[time_col]] == tp
      ]
      if (length(samples_gt) > 0) {
        if (length(samples_gt) > 1) {
          gt_means[[paste0(grp, "_", tp)]] <- rowMeans(mat[, samples_gt, drop = FALSE])
        } else {
          gt_means[[paste0(grp, "_", tp)]] <- mat[, samples_gt]
        }
      }
    }
  }
  gt_mat <- do.call(cbind, gt_means)

  # Apply scale feature sum if requested (TiSA style normalization)
  if (scale_features) {
    row_sums <- rowSums(gt_mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    gt_mat <- sweep(gt_mat, 1, row_sums, "/")
  }

  # --- Build long-format trajectory data ---
  ts_data <- data.frame(
    gene_id = rep(rownames(gt_mat), ncol(gt_mat)),
    variable = rep(colnames(gt_mat), each = nrow(gt_mat)),
    trans_mean = as.vector(gt_mat),
    stringsAsFactors = FALSE
  )

  # Parse group and timepoint from variable name (group_timepoint)
  # Handle cases where group name might contain underscore
  # Split from the end to get timepoint
  ts_data$timepoint <- sapply(strsplit(ts_data$variable, "_"), function(x) x[length(x)])
  ts_data$group <- sapply(strsplit(ts_data$variable, "_"), function(x) {
    paste(x[-length(x)], collapse = "_")
  })

  ts_data$timepoint <- as.numeric(ts_data$timepoint)

  # Add cluster info
  ts_data$cluster <- gene_cluster[ts_data$gene_id]

  # Create cluster counts
  cluster_counts <- table(gene_cluster[gene_cluster %in% clusters])

  # Create facet labels: "C1 - 50 genes | IgM"
  ts_data$labels <- paste0(
    ts_data$cluster, " - ",
    cluster_counts[as.character(ts_data$cluster)], " genes | ",
    ts_data$group
  )

  # Order labels by cluster then group
  label_order <- unique(ts_data$labels[order(ts_data$cluster, ts_data$group)])
  ts_data$labels <- factor(ts_data$labels, levels = label_order)

  # --- Calculate mean trajectory per cluster-group-timepoint ---
  mean_data <- stats::aggregate(
    trans_mean ~ cluster + group + timepoint + labels,
    data = ts_data,
    FUN = mean
  )

  # --- Set up colors ---
  if (is.null(colors)) {
    color_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                   "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    colors <- setNames(color_pal[seq_along(groups)], groups)
  }

  # --- Build plot (TiSA style) ---
  p <- ggplot(ts_data, aes(x = timepoint, y = trans_mean, color = group)) +
    # Individual gene lines (thin, colored, semi-transparent)
    geom_line(aes(group = gene_id), alpha = alpha, linewidth = 0.3) +
    geom_point(size = 0.5, alpha = alpha) +
    # Mean trajectory (thick gray line)
    geom_line(
      data = mean_data,
      aes(x = timepoint, y = trans_mean, group = group),
      linewidth = 1.5,
      color = "grey50",
      inherit.aes = FALSE
    ) +
    # Styling
    scale_color_manual(values = colors) +
    scale_x_continuous(expand = c(0.02, 0)) +
    facet_wrap(~ labels, scales = "free_x", ncol = ncol) +
    labs(
      x = "Timepoint",
      y = if (scale_features) "Scaled expression" else "Expression (z-score)",
      color = "Group"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = title_size, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  return(p)
}


# ==============================================================================
# CLUSTER MEANS OVERLAY
# ==============================================================================

#' Plot Cluster Mean Trajectories (Overlay)
#'
#' Creates a single panel showing mean expression trajectories for all clusters,
#' colored by cluster. Provides a compact overview of cluster behavior across
#' timepoints. When groups are present, lines are distinguished by linetype.
#'
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part()}.
#' @param sample_info Data.frame with sample metadata. Must have rownames
#'   matching colnames of the PART data matrix.
#' @param time_col Character. Column in sample_info for timepoint. Default:
#'   "timepoint".
#' @param group_col Character or NULL. Column in sample_info for group. If NULL,
#'   averages across all samples per timepoint. Default: "group".
#' @param clusters Character vector of clusters to plot. NULL = all non-outlier.
#'   Default: NULL.
#' @param colors Named character vector of colors for clusters. NULL = uses
#'   PART cluster_colors. Default: NULL.
#' @param title Character or NULL. Plot title. Default: NULL (auto-generated).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' AETHER_plot_cluster_means(part_result, sample_info)
#'
#' }
#' @export
AETHER_plot_cluster_means <- function(part_result,
                                       sample_info,
                                       time_col = "timepoint",
                                       group_col = "group",
                                       clusters = NULL,
                                       colors = NULL,
                                       title = NULL) {

  if (!inherits(part_result, "artemis_part")) {
    stop("'part_result' must be an artemis_part object")
  }

  if (!time_col %in% colnames(sample_info)) {
    stop("'", time_col, "' not found in sample_info columns")
  }

  mat <- part_result$data
  cmap <- part_result$cluster_map

  # Align samples
  common_samples <- intersect(colnames(mat), rownames(sample_info))
  if (length(common_samples) == 0) {
    stop("No matching sample names between PART data and sample_info")
  }
  mat <- mat[, common_samples, drop = FALSE]
  sample_info <- sample_info[common_samples, , drop = FALSE]

  # Filter clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(cmap$cluster))
    clusters <- clusters[clusters != "C0"]
  }

  gene_cluster <- setNames(cmap$cluster, cmap$gene)
  has_group <- !is.null(group_col) && group_col %in% colnames(sample_info)

  # Compute cluster mean per timepoint (and optionally per group)
  agg_rows <- list()

  for (cl in clusters) {
    cl_genes <- names(gene_cluster)[gene_cluster == cl]
    cl_mat <- mat[cl_genes, , drop = FALSE]

    for (samp in colnames(cl_mat)) {
      agg_rows[[length(agg_rows) + 1]] <- data.frame(
        cluster = cl,
        sample = samp,
        timepoint = as.character(sample_info[samp, time_col]),
        group = if (has_group) as.character(sample_info[samp, group_col]) else "All",
        mean_expr = mean(cl_mat[, samp], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, agg_rows)

  # Average across samples within timepoint x group x cluster
  agg_final <- stats::aggregate(
    mean_expr ~ cluster + timepoint + group,
    data = plot_df,
    FUN = mean
  )

  # Ensure timepoint ordering
  tp_vals <- unique(agg_final$timepoint)
  tp_order <- tryCatch(sort(as.numeric(as.character(tp_vals))),
                        warning = function(w) sort(tp_vals))
  if (is.numeric(tp_order)) {
    tp_labels <- tp_vals[order(as.numeric(as.character(tp_vals)))]
  } else {
    tp_labels <- sort(tp_vals)
  }
  agg_final$timepoint <- factor(agg_final$timepoint, levels = tp_labels)
  agg_final$tp_num <- as.numeric(agg_final$timepoint)

  # Default colors from PART result
  if (is.null(colors)) {
    colors <- part_result$cluster_colors[clusters]
  }

  # Build plot
  if (has_group && length(unique(agg_final$group)) > 1) {
    p <- ggplot(agg_final, aes(x = tp_num, y = mean_expr,
                                color = cluster, linetype = group,
                                group = interaction(cluster, group))) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      labs(linetype = "Group")
  } else {
    p <- ggplot(agg_final, aes(x = tp_num, y = mean_expr,
                                color = cluster, group = cluster)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5)
  }

  p <- p +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seq_along(tp_labels), labels = tp_labels) +
    labs(
      x = "Timepoint",
      y = "Mean Expression (z-score)",
      color = "Cluster",
      title = if (is.null(title)) "Cluster Mean Trajectories" else title
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  return(p)
}


# ==============================================================================
# TIME SERIES DEA SUMMARY PLOT
# ==============================================================================

#' Plot Time Series DEA Summary
#'
#' Creates a bar chart showing the number of up- and down-regulated genes per
#' comparison from a time series DEA result. Up-regulated genes are shown above
#' the axis, down-regulated below.
#'
#' @param ts_de_result An \code{artemis_ts_de} object from
#'   \code{ARTEMIS_timeseries_conditional()} or
#'   \code{ARTEMIS_timeseries_temporal()}.
#' @param colors Named character vector of length 2: c(up = ..., down = ...).
#'   Default: red for up, blue for down.
#' @param title Character or NULL. Plot title. Default: NULL (auto-generated).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' cond_de <- ARTEMIS_timeseries_conditional(counts, targets, "Ctrl", "Treat")
#' AETHER_plot_timeseries_summary(cond_de)
#'
#' }
#' @export
AETHER_plot_timeseries_summary <- function(ts_de_result,
                                            colors = NULL,
                                            title = NULL) {

  if (!inherits(ts_de_result, "artemis_ts_de")) {
    stop("'ts_de_result' must be an artemis_ts_de object")
  }

  summary_df <- ts_de_result$summary

  # Default colors
  if (is.null(colors)) {
    colors <- c(up = "#B2182B", down = "#2166AC")
  }

  # Build label column for x-axis
  if (ts_de_result$type == "conditional") {
    summary_df$label <- paste0("TP", summary_df$timepoint)
  } else {
    # Temporal: show TP_exp vs TP_ref
    if ("experiment_name" %in% colnames(summary_df)) {
      summary_df$label <- summary_df$experiment_name
    } else {
      summary_df$label <- paste0("TP", summary_df$tp_experiment,
                                  " vs TP", summary_df$tp_reference)
    }
  }

  # Reshape for stacked bars
  plot_df <- data.frame(
    label = rep(summary_df$label, 2),
    direction = rep(c("Up", "Down"), each = nrow(summary_df)),
    count = c(summary_df$n_sig_up, -summary_df$n_sig_down),
    abs_count = c(summary_df$n_sig_up, summary_df$n_sig_down),
    stringsAsFactors = FALSE
  )
  plot_df$label <- factor(plot_df$label, levels = summary_df$label)
  plot_df$direction <- factor(plot_df$direction, levels = c("Up", "Down"))

  # Auto title
  if (is.null(title)) {
    if (ts_de_result$type == "conditional") {
      title <- paste("DE Summary:", ts_de_result$comparison)
    } else {
      title <- "Temporal DE Summary"
    }
  }

  p <- ggplot(plot_df, aes(x = label, y = count, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_text(aes(label = abs_count,
                   vjust = ifelse(direction == "Up", -0.5, 1.5)),
              size = 3.5) +
    scale_fill_manual(values = c("Up" = colors[["up"]],
                                  "Down" = colors[["down"]])) +
    labs(
      x = "Comparison",
      y = "Number of DE Genes",
      fill = "Direction",
      title = title
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  return(p)
}
