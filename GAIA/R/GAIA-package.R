#' @keywords internal
"_PACKAGE"

# ============================================================================
# Base R imports (stats, utils, grDevices, graphics)
# ============================================================================

#' @importFrom stats IQR aggregate as.formula binomial chisq.test coef
#'   complete.cases contr.poly contrasts<- cor.test cov cutree density
#'   dist end gaussian glm hclust kruskal.test mad median model.matrix
#'   na.omit p.adjust pnorm prcomp qnorm quantile reshape rnorm runif sd
#'   setNames start t.test var wilcox.test
#' @importFrom utils capture.output head read.csv read.delim read.table
#'   sessionInfo tail write.csv
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf png rgb
#' @importFrom graphics abline par points text

# ============================================================================
# ggplot2 imports (listed in Imports)
# ============================================================================

#' @importFrom ggplot2 ggplot aes geom_area geom_bar geom_boxplot geom_col
#'   geom_errorbar geom_hline geom_jitter geom_line geom_path geom_point
#'   geom_polygon geom_raster geom_rect geom_segment geom_smooth geom_text
#'   geom_tile geom_violin geom_vline annotate arrow coord_cartesian
#'   coord_fixed coord_flip
#'   element_blank element_line element_rect element_text expansion facet_grid
#'   facet_wrap guide_legend guides labs margin position_stack
#'   scale_color_discrete scale_color_gradient scale_color_gradient2
#'   scale_color_manual scale_fill_gradient scale_fill_gradient2
#'   scale_fill_manual scale_linetype_manual scale_size_continuous
#'   scale_x_continuous scale_y_continuous sec_axis stat_ellipse theme
#'   theme_bw theme_minimal theme_void unit ggsave

# ============================================================================
# ggplot2 .data pronoun for tidy eval
# ============================================================================

#' @importFrom rlang .data

# ============================================================================
# DESeq2 imports (listed in Imports)
# ============================================================================

#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix estimateDispersions
#'   estimateSizeFactors nbinomWaldTest results sizeFactors counts
#'   design "design<-"
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom GenomicRanges seqinfo seqnames GRanges start end coverage
#'   tileGenome binnedAverage distanceToNearest findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import.bw export.bw BigWigFile

# ============================================================================
# WGCNA imports (listed in Imports)
# ============================================================================

#' @importFrom WGCNA blockwiseModules blueWhiteRed corPvalueStudent
#'   cutreeStatic goodSamplesGenes labeledHeatmap labels2colors
#'   moduleEigengenes numbers2colors orderMEs pickSoftThreshold
#'   plotDendroAndColors

# ============================================================================
# Global variables used in NSE / data.table / ggplot aes()
# These are declared to avoid R CMD check NOTEs about "no visible binding"
# ============================================================================

utils::globalVariables(c(
  # ggplot aesthetics / column references
  "activity", "method", "group", "proportion", "cell_type", "cluster",
  "direction", "label", "variable", "contribution", "dimension",
  "estimate", "ci_lower", "ci_upper", "color_group", "label_x",
  "label_hjust", "feature_index", "weight", "is_nonzero", "grid_type",
  "best_score", "cv_score", "param_value", "size_var", "color_var",
  "color_var_plot", "fill_var", "x_val", "Description_short",
  "label_text", "module_label", "term_short", "motif_label", "peak_set",
  "label_var", "sil_width", "timepoint", "trans_mean", "gene_id",
  "tp_num", "mean_expr", "log2FC", "log_avg_signal", "signal", "region",
  "abs_count", "count", "variance_percent", "cumulative_percent",
  "Dim1", "Dim2", "x", "y", "contrib", "type", "xend", "yend",
  "xmin", "xmax", "level", "PC_x", "PC_y", "shape_var", "PC1", "PC2",
  "batch", "k", "silhouette_avg", "lambda", "cve", "cvse",
  "mean_activity", "min_activity", "max_activity",
  "Dataset_label", "Feature", "Correlation", "Trait", "Module", "Label",
  "Power", "signed_R2", "mean.k.", "GS", "MM", "is_hub", "module",
  "n_hubs", "top_hub", "x_plot",

  # data.table symbols
  "chr", "start", "end", ":=", ".N", "bin_id", "overlap_bp",
  "i.end", "i.start", "region_idx",

  # Other NSE / internal
  "Xtilde_1", "CIBERSORT", "cor",

  # GenomeInfoDb (Suggested package)
  "seqlevels"
))
