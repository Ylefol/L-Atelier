# GAIA/Aether/pca_plots.R
# PCA visualization for normalized count data


#' PCA plot of normalized counts
#'
#' Performs PCA on normalized count data and creates a scatter plot colored
#' by group with optional shape mapping for a second variable (e.g., timepoint).
#'
#' @param counts Normalized count matrix (genes x samples), an artemis_norm
#'   object, or an artemis_ts_norm object. Genes as rows, samples as columns.
#' @param sample_info Data.frame with sample metadata. Rownames or a column
#'   must match the column names of the count matrix. If counts is an
#'   artemis_norm/artemis_ts_norm object, this is extracted automatically
#'   (but can be overridden).
#' @param group_col Character. Column in sample_info to use for point color.
#'   Default: "group".
#' @param shape_col Character or NULL. Column in sample_info to use for point
#'   shape. Default: NULL (all points same shape).
#' @param dims Integer vector of length 2. Which PCs to plot. Default: c(1, 2).
#' @param ntop Integer. Number of most variable genes to use for PCA. Set to
#'   NULL to use all genes. Default: 500.
#' @param log_transform Logical. Log2-transform counts before PCA (recommended
#'   for raw normalized counts). Adds a pseudocount of 1. Default: TRUE.
#' @param colors Named character vector of colors for groups, or NULL for
#'   default ggplot2 colors. Names should match levels in group_col.
#' @param point_size Numeric. Size of points. Default: 3.
#' @param label_samples Logical or "repel". TRUE for direct labels, "repel" for
#'   non-overlapping labels (requires ggrepel). Default: FALSE.
#' @param title Character or NULL. Plot title. Default: "PCA".
#' @param verbose Logical. Print PCA summary. Default: TRUE.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' # From artemis_norm object
#' p <- AETHER_plot_pca(norm_data, group_col = "group")
#'
#' # With timepoint shapes
#' p <- AETHER_plot_pca(norm_data, group_col = "group", shape_col = "timepoint")
#'
#' # From matrix + sample sheet
#' p <- AETHER_plot_pca(norm_counts, sample_info = targets,
#'                       group_col = "condition", shape_col = "batch")
#'
#' @export
AETHER_plot_pca <- function(counts,
                             sample_info = NULL,
                             group_col = "group",
                             shape_col = NULL,
                             dims = c(1, 2),
                             ntop = 500,
                             log_transform = TRUE,
                             colors = NULL,
                             point_size = 3,
                             label_samples = FALSE,
                             title = "PCA",
                             verbose = TRUE) {

  # --- Extract data from S3 objects ---
  if (inherits(counts, "artemis_norm")) {
    if (is.null(sample_info)) sample_info <- counts$targets
    counts <- counts$norm_counts
  } else if (inherits(counts, "artemis_ts_norm")) {
    if (is.null(sample_info)) sample_info <- counts$targets
    counts <- counts$norm_counts
  }

  # --- Validate inputs ---
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  if (is.null(sample_info)) {
    stop("'sample_info' is required when 'counts' is a plain matrix")
  }

  if (!group_col %in% colnames(sample_info)) {
    stop("group_col '", group_col, "' not found in sample_info")
  }

  if (!is.null(shape_col) && !shape_col %in% colnames(sample_info)) {
    stop("shape_col '", shape_col, "' not found in sample_info")
  }

  # Match sample order between counts and sample_info
  shared <- intersect(colnames(counts), rownames(sample_info))
  if (length(shared) == 0) {
    stop("No matching sample names between count columns and sample_info rownames")
  }
  if (length(shared) < ncol(counts) && verbose) {
    cat("Note:", ncol(counts) - length(shared), "samples in counts not found in sample_info, dropped\n")
  }
  counts <- counts[, shared, drop = FALSE]
  sample_info <- sample_info[shared, , drop = FALSE]

  # --- Prepare counts ---
  # Remove zero-variance genes
  gene_vars <- apply(counts, 1, var)
  counts <- counts[gene_vars > 0, , drop = FALSE]

  if (verbose) cat("Genes with non-zero variance:", nrow(counts), "\n")

  if (log_transform) {
    counts <- log2(counts + 1)
  }

  # Select top variable genes
  if (!is.null(ntop) && ntop < nrow(counts)) {
    gene_vars <- apply(counts, 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:ntop]
    counts <- counts[top_genes, , drop = FALSE]
    if (verbose) cat("Using top", ntop, "most variable genes\n")
  }

  # --- Run PCA ---
  pca <- prcomp(t(counts), center = TRUE, scale. = FALSE)

  # Variance explained
  var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  if (verbose) {
    cat("Variance explained:\n")
    top_pcs <- min(5, length(var_pct))
    for (i in seq_len(top_pcs)) {
      cat("  PC", i, ":", var_pct[i], "%\n")
    }
  }

  # --- Build plot data ---
  pc1 <- dims[1]
  pc2 <- dims[2]

  plot_df <- data.frame(
    PC_x = pca$x[, pc1],
    PC_y = pca$x[, pc2],
    group = factor(sample_info[[group_col]]),
    sample = rownames(sample_info),
    stringsAsFactors = FALSE
  )

  if (!is.null(shape_col)) {
    plot_df$shape_var <- factor(sample_info[[shape_col]])
  }

  # --- Plot ---
  if (!is.null(shape_col)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = PC_x, y = PC_y, color = group, shape = shape_var
    ))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = PC_x, y = PC_y, color = group
    ))
  }

  p <- p +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(
      x = paste0("PC", pc1, " (", var_pct[pc1], "%)"),
      y = paste0("PC", pc2, " (", var_pct[pc2], "%)"),
      color = group_col,
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(shape_col)) {
    p <- p + ggplot2::labs(shape = shape_col)
  }

  # Custom colors
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  # Labels
  if (!identical(label_samples, FALSE)) {
    if (identical(label_samples, "repel") && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),
        size = 2.5, max.overlaps = 20, show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sample),
        size = 2.5, vjust = -0.8, show.legend = FALSE
      )
    }
  }

  return(p)
}
