#' @keywords internal
"_PACKAGE"

# ============================================================================
# Base R imports
# ============================================================================

#' @importFrom stats median mad IQR quantile sd var p.adjust wilcox.test model.matrix na.omit setNames cor aggregate
#' @importFrom utils head tail read.csv read.table write.csv capture.output sessionInfo
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf png rgb
#' @importFrom methods is as new

# ============================================================================
# ggplot2 imports (Depends — always attached)
# ============================================================================

#' @importFrom ggplot2 ggplot aes geom_bar geom_boxplot geom_col geom_errorbar geom_hline geom_jitter geom_line geom_point geom_text geom_tile
#' @importFrom ggplot2 geom_violin geom_vline annotate coord_flip element_blank element_line element_rect element_text facet_grid facet_wrap
#' @importFrom ggplot2 guide_legend guides labs margin scale_color_gradient scale_color_gradient2 scale_color_manual
#' @importFrom ggplot2 scale_fill_gradient scale_fill_gradient2 scale_fill_gradientn scale_fill_manual scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme theme_bw theme_minimal theme_void unit ggsave

# ============================================================================
# rlang (tidy eval)
# ============================================================================

#' @importFrom rlang .data

# ============================================================================
# SingleCellExperiment / Bioconductor core
# ============================================================================

#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim reducedDim<- reducedDims reducedDimNames colLabels colLabels<-
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assay assay<- assayNames
#' @importFrom BiocGenerics counts counts<- normalize
#' @importFrom S4Vectors DataFrame metadata metadata<-

# ============================================================================
# scran / scater / scDblFinder / SingleR
# ============================================================================

#' @importFrom BiocNeighbors findKNN
#' @importFrom cluster silhouette
#' @importFrom scran computeSumFactors modelGeneVar getTopHVGs findMarkers clusterCells buildSNNGraph
#' @importFrom scater logNormCounts runPCA runUMAP runTSNE plotReducedDim
#' @importFrom scuttle addPerCellQCMetrics perCellQCMetrics perCellQCFilters
#' @importFrom scDblFinder scDblFinder
#' @importFrom SingleR SingleR
#' @importFrom basilisk BasiliskEnvironment basiliskRun

# ============================================================================
# Parallel computing
# ============================================================================

#' @importFrom BiocParallel MulticoreParam SerialParam bpparam

# ============================================================================
# Matrix
# ============================================================================

#' @importFrom Matrix sparseMatrix colSums rowSums colMeans rowMeans

# ============================================================================
# Global variables (NSE / ggplot aes)
# ============================================================================

utils::globalVariables(c(
  "cell", "cluster", "label", "group", "sample", "condition",
  "n_counts", "n_genes", "pct_mt", "pct_ribo", "doublet_score",
  "x", "y", "value", "variable", "feature", "expression",
  "cell_type", "score", "p_value", "log2FC", "direction",
  "k", "metric_label", "n_clusters", "modularity", "mean_sil",
  "dim1", "dim2", "colour_val",
  "n_neighbors", "min_dist", "perplexity", "raw_val", "metric",
  "panel_label", "method", "resolution", "xintercept", "min_cl_size"
))
