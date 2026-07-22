#' Test for differences in cell-type proportions across samples (propeller)
#'
#' Wraps \code{speckle::propeller()} to test whether cell-type (cluster)
#' proportions differ between two or more experimental groups, using
#' sample-level replication. Unlike differential expression, this only reads
#' \code{colData(sce)} -- no assay matrix is touched, so the function runs in
#' a fraction of a second regardless of dataset size.
#'
#' propeller compares proportions at the level of biological samples (one
#' proportion vector per \code{sample_col} level), not individual cells --
#' treating cells as the unit of replication would pseudoreplicate. It
#' variance-stabilises sample-level proportions (arcsine-square-root or
#' logit transform), then fits an empirical-Bayes-moderated t-test (two
#' groups) or F-test/ANOVA (more than two groups) via \code{limma}, borrowing
#' variance information across clusters the same way limma borrows
#' information across genes.
#'
#' As with any test on few biological replicates, statistical power is
#' limited when there are only 2-3 samples per group -- propeller's variance
#' moderation makes this more reliable than an unmoderated per-cluster test,
#' but does not manufacture power the replicate count doesn't support; treat
#' marginal p-values from very small designs cautiously.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param cluster_col Character. \code{colData(sce)} column naming the
#'   cluster/cell-type assignment. Default \code{"cluster"}.
#' @param sample_col Character. \code{colData(sce)} column naming the
#'   biological replicate/sample each cell belongs to. Required -- this is
#'   the unit proportions are computed over.
#' @param group_col Character. \code{colData(sce)} column naming the
#'   experimental group/condition to compare (e.g. genotype). Must be
#'   constant within each \code{sample_col} level.
#' @param cluster_colors Named character vector, cluster label ->
#'   hex colour, one entry for every level of \code{cluster_col} present in
#'   \code{sce}. Optional -- if supplied, carried through in \code{$params}
#'   so \code{\link{ASPIS_plot_composition}} uses these colours instead of
#'   assigning its own, keeping cluster colours consistent with however
#'   they're coloured elsewhere (e.g. on a UMAP). Default \code{NULL} (the
#'   plot functions fall back to their own default palette).
#' @param group_colors Named character vector, group label -> hex colour,
#'   one entry for every level of \code{group_col} present in \code{sce}.
#'   Same purpose as \code{cluster_colors}, consumed by
#'   \code{\link{ASPIS_plot_propeller}} and \code{\link{ASPIS_plot_composition}}'s
#'   facet strips. Default \code{NULL}.
#' @param transform Character. Variance-stabilising transform applied to
#'   proportions before testing: \code{"logit"} (default, matches
#'   \code{speckle::propeller()}'s own default) or \code{"asin"}
#'   (arcsine-square-root).
#' @param trend Logical. Fit a mean-variance trend on the transformed
#'   proportions (passed to \code{limma::eBayes()}). Default \code{FALSE},
#'   matching the \code{propeller} default.
#' @param robust Logical. Use robust empirical Bayes shrinkage of variances.
#'   Default \code{TRUE}, matching the \code{propeller} default. Note
#'   \code{propeller} itself forces this to \code{FALSE} internally when
#'   fewer than 3 clusters are tested.
#' @param fdr_threshold Numeric. FDR cutoff used to populate
#'   \code{$significant}. Default \code{0.05}.
#' @param verbose Logical. Print \code{[KERAUNOS]} progress output, including
#'   a low-replicate warning when any group has fewer than 3 samples. Default
#'   \code{TRUE}.
#'
#' @return A classed list \code{keraunos_propeller}:
#'   \describe{
#'     \item{results}{Tidy data.frame, one row per cluster: \code{cluster},
#'       \code{baseline_prop}, per-group mean proportion
#'       (\code{PropMean.<group>}), \code{PropRatio}/\code{Tstatistic} (two
#'       groups) or \code{Fstatistic} (more than two groups), \code{P.Value},
#'       \code{FDR}. Ordered by \code{P.Value}. \code{baseline_prop} is
#'       normalised here to a single consistent name -- speckle itself calls
#'       it \code{BaselineProp.Freq} in the two-group case and
#'       \code{BaselineProp} in the >2-group case.}
#'     \item{significant}{Subset of \code{results} with \code{FDR <
#'       fdr_threshold}.}
#'     \item{proportions}{Tidy data.frame, one row per sample x cluster:
#'       \code{sample}, \code{cluster}, \code{proportion}, \code{group}. Feeds
#'       \code{\link{ASPIS_plot_propeller}} so the plot never needs \code{sce}
#'       directly.}
#'     \item{summary}{One-row data.frame: \code{n_clusters}, \code{n_sig}.}
#'     \item{params}{List of parameters used, including per-group sample
#'       counts (\code{n_per_group}) and any supplied \code{cluster_colors}/
#'       \code{group_colors}.}
#'   }
#' @export
KERAUNOS_propeller_proportions <- function(sce, cluster_col = "cluster",
                                            sample_col, group_col,
                                            cluster_colors = NULL,
                                            group_colors = NULL,
                                            transform = c("logit", "asin"),
                                            trend = FALSE, robust = TRUE,
                                            fdr_threshold = 0.05,
                                            verbose = TRUE) {

  if (!inherits(sce, "SingleCellExperiment"))
    stop("'sce' must be a SingleCellExperiment.", call. = FALSE)
  if (missing(sample_col) || is.null(sample_col))
    stop("'sample_col' is required -- name the colData(sce) column identifying biological replicates.",
         call. = FALSE)
  if (missing(group_col) || is.null(group_col))
    stop("'group_col' is required -- name the colData(sce) column identifying the groups to compare.",
         call. = FALSE)
  transform <- match.arg(transform)

  needed <- c(cluster_col, sample_col, group_col)
  missing_cols <- setdiff(needed, names(colData(sce)))
  if (length(missing_cols) > 0L)
    stop("colData(sce) is missing: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)

  if (!requireNamespace("speckle", quietly = TRUE))
    stop("Package 'speckle' is required for KERAUNOS_propeller_proportions(). ",
         "Install it with BiocManager::install(\"speckle\").", call. = FALSE)

  clusters <- as.character(colData(sce)[[cluster_col]])
  samples  <- as.character(colData(sce)[[sample_col]])
  groups   <- as.character(colData(sce)[[group_col]])

  cluster_colors <- .keraunos_validate_colors(cluster_colors, unique(clusters), "cluster_colors")
  group_colors   <- .keraunos_validate_colors(group_colors, unique(groups), "group_colors")

  sample_group <- unique(data.frame(sample = samples, group = groups))
  dup_samples <- sample_group$sample[duplicated(sample_group$sample)]
  if (length(dup_samples) > 0L)
    stop("'", sample_col, "' values map to more than one '", group_col,
         "' level, so it is not a per-sample covariate: ",
         paste(unique(dup_samples), collapse = ", "), call. = FALSE)

  n_per_group <- table(sample_group$group)

  if (isTRUE(verbose)) {
    cat("[KERAUNOS] Cell-type Proportion Testing (speckle::propeller)\n")
    cat("    Clusters:", length(unique(clusters)), "\n")
    cat("    Samples:", length(unique(samples)), "across", length(n_per_group), "groups\n")
    for (g in names(n_per_group))
      cat("      ", g, ": ", n_per_group[[g]], " sample(s)\n", sep = "")
    if (any(n_per_group < 3L))
      cat("    Note: at least one group has < 3 samples -- variance moderation is\n",
          "    less reliable here and power will be low; treat marginal results cautiously.\n", sep = "")
    cat("    Transform:", transform, "| trend:", trend, "| robust:", robust, "\n\n")
  }

  raw <- speckle::propeller(clusters = clusters, sample = samples, group = groups,
                             transform = transform, trend = trend, robust = robust)

  # speckle::propeller() shapes its baseline-proportion column differently
  # depending on which branch ran internally: the 2-group t-test path emits
  # "BaselineProp.clusters" + "BaselineProp.Freq", but the >2-group ANOVA path
  # emits only a plain "BaselineProp" column and puts cluster identity solely
  # in rownames(raw) instead. rownames(raw) holds the cluster labels in BOTH
  # branches, so use that as the single source of truth rather than relying
  # on a column name that only exists in one of the two branches.
  results <- as.data.frame(raw, stringsAsFactors = FALSE)
  results$cluster <- rownames(raw)
  results$BaselineProp.clusters <- NULL
  names(results)[names(results) == "BaselineProp.Freq"] <- "baseline_prop"
  names(results)[names(results) == "BaselineProp"] <- "baseline_prop"
  results <- results[, c("cluster", setdiff(names(results), "cluster")), drop = FALSE]
  results <- results[order(results$P.Value), , drop = FALSE]
  rownames(results) <- results$cluster

  significant <- results[!is.na(results$FDR) & results$FDR < fdr_threshold, , drop = FALSE]

  summary_df <- data.frame(n_clusters = nrow(results), n_sig = nrow(significant))

  # Per-sample, per-cluster proportions -- not produced by speckle::propeller()
  # itself (which only returns group-level test statistics), but needed by
  # ASPIS_plot_propeller() to show the actual per-sample points behind the
  # test. Computed here (not in ASPIS) so the plot function never has to
  # touch sce/colData directly -- ASPIS renders, KERAUNOS computes.
  prop_tab <- prop.table(table(sample = samples, cluster = clusters), margin = 1)
  proportions <- as.data.frame(prop_tab, stringsAsFactors = FALSE)
  names(proportions)[names(proportions) == "Freq"] <- "proportion"
  proportions$sample  <- as.character(proportions$sample)
  proportions$cluster <- as.character(proportions$cluster)
  proportions$group <- sample_group$group[match(proportions$sample, sample_group$sample)]

  out <- list(
    results     = results,
    significant = significant,
    proportions = proportions,
    summary     = summary_df,
    params      = list(cluster_col = cluster_col, sample_col = sample_col,
                        group_col = group_col, transform = transform,
                        trend = trend, robust = robust,
                        fdr_threshold = fdr_threshold,
                        n_per_group = as.list(n_per_group),
                        cluster_colors = cluster_colors,
                        group_colors = group_colors)
  )
  class(out) <- "keraunos_propeller"

  if (isTRUE(verbose)) {
    cat("[KERAUNOS] --- Summary ---\n")
    cat("    Tested:", summary_df$n_clusters, "clusters\n")
    cat("    Significant (FDR <", fdr_threshold, "):", summary_df$n_sig, "\n\n")
  }

  out
}

# Errors if `colors` doesn't name every value in `levels_present` -- a
# partial map would silently reintroduce the exact per-plot colour
# inconsistency this parameter exists to prevent. Returns NULL unchanged
# (the "use default palette" signal to the ASPIS plot functions).
.keraunos_validate_colors <- function(colors, levels_present, arg_name) {
  if (is.null(colors)) return(NULL)
  if (is.null(names(colors)) || any(!nzchar(names(colors))))
    stop("'", arg_name, "' must be a named character vector (name = category label, value = hex colour).",
         call. = FALSE)
  missing_levels <- setdiff(levels_present, names(colors))
  if (length(missing_levels) > 0L)
    stop("'", arg_name, "' is missing colours for: ", paste(missing_levels, collapse = ", "),
         call. = FALSE)
  colors
}

#' @export
print.keraunos_propeller <- function(x, ...) {
  cat("[KERAUNOS] Cell-type proportion test (propeller)\n")
  cat("    Comparing:", x$params$group_col, "| Samples per group:\n")
  for (g in names(x$params$n_per_group))
    cat("      ", g, ": ", x$params$n_per_group[[g]], "\n", sep = "")
  cat("    Clusters tested:", x$summary$n_clusters, "\n")
  cat("    Significant (FDR <", x$params$fdr_threshold, "):", x$summary$n_sig, "\n")
  if (nrow(x$significant) > 0L) {
    cat("\n")
    print(x$significant)
  }
  invisible(x)
}
