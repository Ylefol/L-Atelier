# ==============================================================================
# ASPIS — Marker gene visualisation
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Dot plot and heatmap for exploring cluster marker genes returned by
# KERAUNOS_find_markers() or user-supplied gene lists.
#
# Public functions:
#   ASPIS_plot_marker_dotplot()   — size = % expressing, fill = scaled expr
#   ASPIS_plot_marker_heatmap()   — cluster-averaged expression, faceted by
#                                   source cluster
#   ASPIS_plot_db_heatmap()       — marker database heatmap; rows = database
#                                   genes faceted by cell type, cols = clusters;
#                                   missing genes shown as grey NA tiles
#
# Internal helpers:
#   .aspis_extract_marker_genes()
#   .aspis_compute_dot_stats()
#   .aspis_compute_cluster_avg()
# ==============================================================================


# ==============================================================================
# Internal helpers
# ==============================================================================

# Resolve genes + source cluster from a keraunos_markers object or a character
# vector.  Returns data.frame(gene, source_cluster).
.aspis_extract_marker_genes <- function(markers, top_n) {
  if (inherits(markers, "keraunos_markers")) {
    top <- markers$top
    if (is.null(top) || nrow(top) == 0L)
      stop("keraunos_markers$top is empty. Re-run KERAUNOS_find_markers() ",
           "with a lower fdr_threshold.", call. = FALSE)
    genes_df <- do.call(rbind, lapply(
      split(top, top$cluster),
      function(df) {
        df <- df[order(df$Top, df$FDR), ]
        head(df[, c("gene", "cluster"), drop = FALSE], top_n)
      }
    ))
    rownames(genes_df) <- NULL
    names(genes_df)[names(genes_df) == "cluster"] <- "source_cluster"
    genes_df
  } else if (is.character(markers)) {
    data.frame(gene = markers, source_cluster = NA_character_,
               stringsAsFactors = FALSE)
  } else {
    stop("'markers' must be a keraunos_markers object or a character vector.",
         call. = FALSE)
  }
}

# Per-cluster dot statistics: mean expression and % expressing.
# Returns data.frame(gene, cluster, mean_expr, pct_expr).
.aspis_compute_dot_stats <- function(sce, genes, cluster_col, assay_name,
                                      min_expr) {
  clusters  <- as.character(colData(sce)[[cluster_col]])
  cl_levels <- sort(unique(clusters))
  mat       <- as.matrix(assay(sce, assay_name)[genes, , drop = FALSE])

  do.call(rbind, lapply(cl_levels, function(cl) {
    idx <- clusters == cl
    sub <- mat[, idx, drop = FALSE]
    data.frame(
      gene      = genes,
      cluster   = cl,
      mean_expr = rowMeans(sub),
      pct_expr  = rowMeans(sub > min_expr) * 100L,
      stringsAsFactors = FALSE
    )
  }))
}

# Cluster-averaged expression matrix (genes × clusters).
.aspis_compute_cluster_avg <- function(sce, genes, cluster_col, assay_name) {
  clusters  <- as.character(colData(sce)[[cluster_col]])
  cl_levels <- sort(unique(clusters))
  mat       <- as.matrix(assay(sce, assay_name)[genes, , drop = FALSE])

  avg <- vapply(cl_levels, function(cl) {
    rowMeans(mat[, clusters == cl, drop = FALSE])
  }, numeric(length(genes)))
  rownames(avg) <- genes
  colnames(avg) <- cl_levels
  avg
}

# Inline squish: clamp x to [range[1], range[2]] — avoids importing scales.
.aspis_squish <- function(x, range) pmin(pmax(x, range[1]), range[2])


# ==============================================================================
# ASPIS_plot_marker_dotplot
# ==============================================================================

#' Marker gene dot plot
#'
#' Visualises cluster marker genes as a dot plot where dot \strong{size}
#' encodes the fraction of cells expressing each gene and dot \strong{colour}
#' encodes scaled (or raw) mean expression.  Genes are ordered by source
#' cluster and dashed horizontal lines separate marker groups.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or
#'   other) assay.
#' @param markers A \code{keraunos_markers} object returned by
#'   \code{\link{KERAUNOS_find_markers}}, or a character vector of gene names.
#' @param cluster_col Character. \code{colData} column containing cluster
#'   labels.  Default \code{"cluster"}.
#' @param assay_name Character. Assay to pull expression from.
#'   Default \code{"logcounts"}.
#' @param top_n Integer. Maximum marker genes per cluster shown when
#'   \code{markers} is a \code{keraunos_markers} object.  Default \code{5L}.
#' @param scale Logical. Z-score mean expression per gene across clusters
#'   before plotting.  Default \code{TRUE}.
#' @param min_expr Numeric. Expression threshold used to count a cell as
#'   "expressing" the gene.  Default \code{0}.
#' @param dot_max_size Numeric. Maximum dot diameter (ggplot size units).
#'   Default \code{6}.
#' @param palette Character vector of colours for the fill gradient.
#'   \code{NULL} (default) uses a blue–white–red diverging palette.
#' @param title Character. Plot title.  Default \code{"Marker genes"}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_marker_dotplot <- function(sce,
                                       markers,
                                       cluster_col  = "cluster",
                                       assay_name   = "logcounts",
                                       top_n        = 5L,
                                       scale        = TRUE,
                                       min_expr     = 0,
                                       dot_max_size = 6,
                                       palette      = NULL,
                                       title        = "Marker genes") {

  if (!cluster_col %in% names(colData(sce)))
    stop("'", cluster_col, "' not found in colData(sce).", call. = FALSE)
  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.", call. = FALSE)

  genes_df <- .aspis_extract_marker_genes(markers, as.integer(top_n))
  genes    <- unique(genes_df$gene)

  missing <- setdiff(genes, rownames(sce))
  if (length(missing) > 0L)
    warning("Genes not found in SCE and skipped: ",
            paste(missing, collapse = ", "), call. = FALSE)
  genes <- intersect(genes, rownames(sce))
  if (length(genes) == 0L)
    stop("No valid genes to plot.", call. = FALSE)

  dot_df <- .aspis_compute_dot_stats(sce, genes, cluster_col, assay_name,
                                      min_expr)

  if (isTRUE(scale)) {
    dot_df <- do.call(rbind, lapply(split(dot_df, dot_df$gene), function(gdf) {
      m <- mean(gdf$mean_expr, na.rm = TRUE)
      s <- sd(gdf$mean_expr,   na.rm = TRUE)
      gdf$expr_plot <- if (is.na(s) || s == 0) rep(0, nrow(gdf)) else
                         (gdf$mean_expr - m) / s
      gdf
    }))
    fill_name <- "Scaled\nexpression"
    fill_lims <- c(-2.5, 2.5)
  } else {
    dot_df$expr_plot <- dot_df$mean_expr
    fill_name <- "Mean\nexpression"
    fill_lims <- range(dot_df$mean_expr, na.rm = TRUE)
  }

  # Ordering: preserve source-cluster order from genes_df
  gene_order      <- unique(genes_df$gene[genes_df$gene %in% genes])
  dot_df$gene     <- factor(dot_df$gene,    levels = rev(gene_order))
  dot_df$cluster  <- factor(dot_df$cluster, levels = sort(unique(dot_df$cluster)))

  # Dashed separator lines between source-cluster gene groups
  has_source <- !all(is.na(genes_df$source_cluster))
  if (has_source) {
    source_map <- setNames(genes_df$source_cluster, genes_df$gene)
    sources    <- source_map[gene_order]
    cl_sizes   <- rle(sources)$lengths
    n_genes    <- length(gene_order)
    boundaries <- cumsum(cl_sizes)[-length(cl_sizes)]
    hline_y    <- n_genes - boundaries + 0.5
  } else {
    hline_y <- numeric(0L)
  }

  fill_scale <- if (is.null(palette)) {
    scale_fill_gradient2(
      low      = "#3361A5",
      mid      = "white",
      high     = "#A31D1D",
      midpoint = 0,
      limits   = fill_lims,
      oob      = .aspis_squish,
      name     = fill_name
    )
  } else {
    scale_fill_gradientn(colours = palette, name = fill_name)
  }

  p <- ggplot(dot_df, aes(x = cluster, y = gene,
                            size = pct_expr, fill = expr_plot)) +
    geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
    fill_scale +
    scale_size_continuous(
      range  = c(0.5, dot_max_size),
      name   = "% expressing",
      limits = c(0, 100)
    ) +
    labs(title = title, x = cluster_col, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.ticks       = element_blank()
    )

  if (length(hline_y) > 0L)
    p <- p + geom_hline(yintercept = hline_y, colour = "grey40",
                         linewidth = 0.6, linetype = "dashed")

  p
}


# ==============================================================================
# ASPIS_plot_marker_heatmap
# ==============================================================================

#' Marker gene heatmap
#'
#' Visualises cluster-averaged expression of marker genes as a tiled heatmap.
#' When \code{markers} is a \code{keraunos_markers} object, genes are grouped
#' into facets by their source cluster so marker identity is immediately clear.
#' Columns are clusters; rows are genes; colour encodes z-scored (or raw) mean
#' expression.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or
#'   other) assay.
#' @param markers A \code{keraunos_markers} object returned by
#'   \code{\link{KERAUNOS_find_markers}}, or a character vector of gene names.
#' @param cluster_col Character. \code{colData} column containing cluster
#'   labels.  Default \code{"cluster"}.
#' @param assay_name Character. Assay to pull expression from.
#'   Default \code{"logcounts"}.
#' @param top_n Integer. Maximum marker genes per cluster shown when
#'   \code{markers} is a \code{keraunos_markers} object.  Default \code{5L}.
#' @param scale Logical. Z-score mean expression per gene across clusters.
#'   Default \code{TRUE}.
#' @param palette Character vector of colours for the fill gradient.
#'   \code{NULL} (default) uses a blue–white–red diverging palette.
#' @param show_gene_labels Logical. Show gene names on the y axis.
#'   Default \code{TRUE}.
#' @param title Character. Plot title.  Default \code{"Marker gene heatmap"}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_marker_heatmap <- function(sce,
                                       markers,
                                       cluster_col      = "cluster",
                                       assay_name       = "logcounts",
                                       top_n            = 5L,
                                       scale            = TRUE,
                                       palette          = NULL,
                                       show_gene_labels = TRUE,
                                       title            = "Marker gene heatmap") {

  if (!cluster_col %in% names(colData(sce)))
    stop("'", cluster_col, "' not found in colData(sce).", call. = FALSE)
  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.", call. = FALSE)

  genes_df <- .aspis_extract_marker_genes(markers, as.integer(top_n))
  genes    <- unique(genes_df$gene)

  missing <- setdiff(genes, rownames(sce))
  if (length(missing) > 0L)
    warning("Genes not found in SCE and skipped: ",
            paste(missing, collapse = ", "), call. = FALSE)
  genes <- intersect(genes, rownames(sce))
  if (length(genes) == 0L)
    stop("No valid genes to plot.", call. = FALSE)

  avg_mat <- .aspis_compute_cluster_avg(sce, genes, cluster_col, assay_name)

  if (isTRUE(scale)) {
    avg_mat           <- t(scale(t(avg_mat)))
    avg_mat[is.nan(avg_mat)] <- 0
    fill_name <- "Scaled\nexpression"
    fill_lims <- c(-2.5, 2.5)
  } else {
    fill_name <- "Mean\nexpression"
    fill_lims <- range(avg_mat, na.rm = TRUE)
  }

  # Gene and cluster ordering
  gene_order <- unique(genes_df$gene[genes_df$gene %in% genes])
  cl_order   <- sort(unique(as.character(colData(sce)[[cluster_col]])))

  # Long format
  long_df <- do.call(rbind, lapply(gene_order, function(g) {
    data.frame(
      gene    = g,
      cluster = cl_order,
      value   = avg_mat[g, cl_order],
      stringsAsFactors = FALSE
    )
  }))
  long_df$gene    <- factor(long_df$gene,    levels = rev(gene_order))
  long_df$cluster <- factor(long_df$cluster, levels = cl_order)

  # Source cluster facets (keraunos_markers input only)
  has_source <- !all(is.na(genes_df$source_cluster))
  if (has_source) {
    source_map  <- setNames(genes_df$source_cluster, genes_df$gene)
    sc_order    <- unique(source_map[gene_order])  # order of first appearance
    long_df$source_cluster <- factor(
      source_map[as.character(long_df$gene)],
      levels = rev(sc_order)   # rev → first cluster appears at top
    )
  }

  fill_scale <- if (is.null(palette)) {
    scale_fill_gradient2(
      low      = "#3361A5",
      mid      = "white",
      high     = "#A31D1D",
      midpoint = 0,
      limits   = fill_lims,
      oob      = .aspis_squish,
      name     = fill_name
    )
  } else {
    scale_fill_gradientn(colours = palette, name = fill_name)
  }

  p <- ggplot(long_df, aes(x = cluster, y = gene, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    fill_scale +
    labs(title = title, x = cluster_col, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks       = element_blank(),
      panel.spacing    = unit(2, "pt"),
      strip.text.y     = element_text(angle = 0, hjust = 0, size = 9),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )

  if (has_source)
    p <- p + facet_grid(source_cluster ~ ., scales = "free_y",
                         space = "free_y", switch = "y")

  if (!isTRUE(show_gene_labels))
    p <- p + theme(axis.text.y = element_blank())

  p
}


# ==============================================================================
# ASPIS_plot_db_heatmap
# ==============================================================================

#' Marker database heatmap
#'
#' Visualises cluster-averaged expression of genes drawn from a user-supplied
#' marker database.  Rows are genes grouped into facets by cell type; columns
#' are clusters.  Genes absent from the SCE are retained as grey \code{NA}
#' tiles so coverage gaps are immediately visible.  When a \code{neg_markers}
#' column is present, positive and negative markers form nested sub-facets
#' within each cell type.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} (or other)
#'   assay.
#' @param markers A \code{data.frame} with columns:
#'   \describe{
#'     \item{cell_type}{Cell type name used for row faceting.}
#'     \item{markers}{Positive marker genes, comma-separated string.}
#'     \item{neg_markers}{(Optional) Negative marker genes, comma-separated.
#'       When present, positive and negative genes are shown in nested
#'       sub-facets.}
#'   }
#' @param cluster_col Character. \code{colData} column containing cluster
#'   labels used for column grouping.  Default \code{"cluster"}.
#' @param assay_name Character. Assay to pull expression from.
#'   Default \code{"logcounts"}.
#' @param scale Logical. Z-score mean expression per gene across clusters.
#'   Default \code{TRUE}.
#' @param na_color Character. Fill colour for genes absent from the SCE.
#'   Default \code{"grey80"}.
#' @param palette Character vector of colours for the fill gradient.
#'   \code{NULL} (default) uses a blue-white-red diverging palette.
#' @param show_gene_labels Logical. Show gene names on the y axis.
#'   Default \code{TRUE}.
#' @param title Character. Plot title.  Default \code{"Marker database heatmap"}.
#'
#' @return A \code{ggplot} object.
#' @export
ASPIS_plot_db_heatmap <- function(sce,
                                   markers,
                                   cluster_col      = "cluster",
                                   assay_name       = "logcounts",
                                   scale            = TRUE,
                                   na_color         = "grey80",
                                   palette          = NULL,
                                   show_gene_labels = TRUE,
                                   title            = "Marker database heatmap") {

  # ── Input validation ─────────────────────────────────────────────────────────
  if (!cluster_col %in% names(colData(sce)))
    stop("'", cluster_col, "' not found in colData(sce).", call. = FALSE)
  if (!assay_name %in% assayNames(sce))
    stop("Assay '", assay_name, "' not found.", call. = FALSE)
  if (!is.data.frame(markers))
    stop("'markers' must be a data.frame.", call. = FALSE)
  if (!all(c("cell_type", "markers") %in% names(markers)))
    stop("'markers' must have columns 'cell_type' and 'markers'.", call. = FALSE)

  has_neg <- "neg_markers" %in% names(markers)

  # ── Parse comma-separated gene strings ───────────────────────────────────────
  .parse_genes <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
    trimws(unlist(strsplit(as.character(x), ",")))
  }

  # Build ordered gene metadata: cell_type order preserved from data.frame;
  # positive genes listed before negative within each cell type
  gene_meta <- do.call(rbind, lapply(seq_len(nrow(markers)), function(i) {
    ct  <- markers$cell_type[i]
    pos <- .parse_genes(markers$markers[i])
    neg <- if (has_neg) .parse_genes(markers$neg_markers[i]) else character(0)
    rows <- NULL
    if (length(pos) > 0L)
      rows <- rbind(rows, data.frame(gene = pos, cell_type = ct,
                                     marker_type = "positive",
                                     stringsAsFactors = FALSE))
    if (length(neg) > 0L)
      rows <- rbind(rows, data.frame(gene = neg, cell_type = ct,
                                     marker_type = "negative",
                                     stringsAsFactors = FALSE))
    rows
  }))

  if (is.null(gene_meta) || nrow(gene_meta) == 0L)
    stop("No genes parsed from 'markers'.", call. = FALSE)

  # Deduplicate within each cell type only: same gene can appear in multiple
  # cell type facets but not twice within the same facet
  gene_meta <- gene_meta[!duplicated(paste(gene_meta$gene, gene_meta$cell_type)), , drop = FALSE]

  all_genes_uniq <- unique(gene_meta$gene)   # for expression lookup
  all_genes      <- gene_meta$gene           # display order, may repeat across cell types

  # ── Report presence / absence ─────────────────────────────────────────────────
  present <- intersect(all_genes_uniq, rownames(sce))
  missing <- setdiff(all_genes_uniq, rownames(sce))

  cat("[ASPIS] Marker database heatmap\n")
  cat("    Total marker genes : ", length(all_genes_uniq), "\n", sep = "")
  cat("    Found in SCE       : ", length(present), "\n", sep = "")
  cat("    Missing            : ", length(missing), "\n", sep = "")
  if (length(missing) > 0L)
    cat("    Missing genes      : ", paste(missing, collapse = ", "), "\n", sep = "")

  if (length(present) == 0L)
    stop("None of the marker genes were found in the SCE.", call. = FALSE)

  # ── Cluster-averaged expression ───────────────────────────────────────────────
  cl_levels   <- sort(unique(as.character(colData(sce)[[cluster_col]])))
  avg_present <- .aspis_compute_cluster_avg(sce, present, cluster_col, assay_name)

  if (isTRUE(scale)) {
    avg_present              <- t(scale(t(avg_present)))
    avg_present[is.nan(avg_present)] <- 0
    fill_name <- "Scaled\nexpression"
    fill_lims <- c(-2.5, 2.5)
  } else {
    fill_name <- "Mean\nexpression"
    fill_lims <- range(avg_present, na.rm = TRUE)
  }

  # Build avg_mat keyed by unique gene names (NA rows for missing genes)
  if (length(missing) > 0L) {
    na_mat  <- matrix(NA_real_, nrow = length(missing), ncol = length(cl_levels),
                      dimnames = list(missing, cl_levels))
    avg_mat <- rbind(avg_present, na_mat)
  } else {
    avg_mat <- avg_present
  }

  # ── Build long data.frame ─────────────────────────────────────────────────────
  # Iterate gene_meta rows so each (gene, cell_type) pair gets its own cluster rows
  long_df <- do.call(rbind, lapply(seq_len(nrow(gene_meta)), function(i) {
    g <- gene_meta$gene[i]
    data.frame(
      gene        = g,
      cell_type   = gene_meta$cell_type[i],
      marker_type = gene_meta$marker_type[i],
      cluster     = cl_levels,
      value       = avg_mat[g, cl_levels],
      stringsAsFactors = FALSE
    )
  }))

  # Gene order: rev() so first gene in all_genes appears at top of plot
  long_df$gene    <- factor(long_df$gene,    levels = rev(unique(all_genes)))
  long_df$cluster <- factor(long_df$cluster, levels = cl_levels)

  # Build facet label: combined "cell_type - marker_type" when neg present,
  # plain cell_type otherwise — avoids double strip and the gap between them
  if (has_neg) {
    long_df$facet_label <- paste(long_df$cell_type, long_df$marker_type, sep = " - ")
    facet_order         <- unique(paste(gene_meta$cell_type, gene_meta$marker_type,
                                        sep = " - "))
    long_df$facet_label <- factor(long_df$facet_label, levels = rev(facet_order))
  } else {
    ct_order            <- unique(gene_meta$cell_type)
    long_df$facet_label <- factor(long_df$cell_type, levels = rev(ct_order))
  }

  # ── Fill scale ────────────────────────────────────────────────────────────────
  fill_scale <- if (is.null(palette)) {
    scale_fill_gradient2(
      low      = "#3361A5",
      mid      = "white",
      high     = "#A31D1D",
      midpoint = 0,
      limits   = fill_lims,
      oob      = .aspis_squish,
      na.value = na_color,
      name     = fill_name
    )
  } else {
    scale_fill_gradientn(colours = palette, na.value = na_color, name = fill_name)
  }

  # ── Build plot ────────────────────────────────────────────────────────────────
  p <- ggplot(long_df, aes(x = cluster, y = gene, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    fill_scale +
    labs(title = title, x = cluster_col, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks       = element_blank(),
      panel.spacing    = unit(2, "pt"),
      strip.text.y     = element_text(angle = 0, hjust = 0, size = 9),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )

  p <- p + facet_grid(facet_label ~ ., scales = "free_y", space = "free_y")

  if (!isTRUE(show_gene_labels))
    p <- p + theme(axis.text.y = element_blank())

  p
}
