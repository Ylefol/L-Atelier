#' Artemis - Peak Annotation Distribution Comparison
#'
#' @description Statistical comparison of genomic feature distributions
#' between two groups of samples, using per-sample proportions.


###############################################################################
########### Peak Annotation Distribution Comparison ###########
###############################################################################

#' Compare Peak Annotation Distributions Between Two Groups
#'
#' @description
#' For each sample, computes the percentage of peaks falling in each genomic
#' feature category (Promoter, Intron, etc.), then tests whether those
#' percentages differ between two groups using a Wilcoxon rank-sum or
#' two-sample t-test. Multiple-testing correction is applied across features.
#'
#' This per-sample approach correctly treats each biological replicate as
#' independent, unlike pooling all peaks per group (which ignores replication).
#'
#' @param annotated_list Named list of data.frames from \code{APOLLO_annotate_peaks()}.
#'   Each element corresponds to one sample; the name is used as the sample
#'   identifier and must match a value in \code{metadata$sample_id} (or
#'   \code{rownames(metadata)}).
#' @param metadata Data.frame with at least columns \code{sample_id} and
#'   \code{group_col}. If \code{sample_id} is absent, \code{rownames(metadata)}
#'   is used.
#' @param group_col Character. Column in \code{metadata} containing group labels.
#'   Exactly two unique values are required. Default: \code{"group"}.
#' @param feature_col Character. Column in each annotated data.frame to use as
#'   the feature category. Default: \code{"annotation_simple"}. If absent and
#'   \code{"annotation"} exists, a simplified version is derived automatically.
#' @param method Character. Statistical test to compare group proportions.
#'   One of \code{"wilcoxon"} (default; more robust with small n) or
#'   \code{"t.test"}.
#' @param p_adj_method Character. Multiple-testing correction method passed to
#'   \code{p.adjust()}. Default: \code{"BH"} (Benjamini-Hochberg).
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return A named list with class \code{"artemis_annotation_test"} containing:
#'   \describe{
#'     \item{results}{Data.frame of per-feature test results: feature,
#'       group mean counts, mean difference, test statistic, p.value, p.adj.}
#'     \item{counts}{Long-format data.frame of per-sample peak counts
#'       with group labels attached (for plotting).}
#'     \item{groups}{Character vector of the two group labels.}
#'     \item{group_col}{The group column name used.}
#'     \item{feature_col}{The feature column name used.}
#'     \item{method}{The test method used.}
#'   }
#'
#' @details
#' \strong{Why per-sample counts?}
#' Each sample's per-feature peak count is treated as the unit of inference.
#' This preserves absolute peak numbers, which may be biologically meaningful
#' when total peak counts differ between groups.
#'
#' \strong{Confounding note:} If samples differ in total peak count, raw counts
#' can be confounded. This is not a concern when peak libraries have been
#' balanced (e.g., via spike-in normalization) within each batch.
#'
#' \strong{Power note:} With few replicates (n = 2-3 per group), statistical
#' power is limited regardless of method. Results should be interpreted
#' alongside effect sizes (mean difference) rather than p-values alone.
#'
#' @seealso \code{\link{AETHER_plot_annotation_comparison}} for visualization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Annotate peaks per sample
#' peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' annotated_list <- lapply(peak_list, APOLLO_annotate_peaks, txdb = txdb)
#'
#' # Build minimal metadata
#' meta <- sample_sheet[, c("sample_id", "group")]
#'
#' # Run comparison
#' test_result <- ARTEMIS_compare_annotation_distribution(
#'   annotated_list = annotated_list,
#'   metadata       = meta,
#'   group_col      = "group"
#' )
#'
#' # Visualise
#' p <- AETHER_plot_annotation_comparison(test_result)
#' }
ARTEMIS_compare_annotation_distribution <- function(annotated_list,
                                                     metadata,
                                                     group_col    = "group",
                                                     feature_col  = "annotation_simple",
                                                     method       = "wilcoxon",
                                                     p_adj_method = "BH",
                                                     verbose      = TRUE) {

  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  if (is.data.frame(annotated_list))
    stop("annotated_list must be a named list of data.frames (one per sample), ",
         "not a single data.frame.", call. = FALSE)
  if (!is.list(annotated_list) || is.null(names(annotated_list)))
    stop("annotated_list must be a named list.", call. = FALSE)

  method <- match.arg(method, c("wilcoxon", "t.test"))

  if (!"sample_id" %in% colnames(metadata)) {
    rn <- rownames(metadata)
    if (!is.null(rn) && !all(rn == seq_len(nrow(metadata))))
      metadata$sample_id <- rn
    else
      stop("metadata must have a 'sample_id' column or informative rownames.", call. = FALSE)
  }
  if (!group_col %in% colnames(metadata))
    stop("group_col '", group_col, "' not found in metadata.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Compute per-sample feature proportions
  # ---------------------------------------------------------------------------
  if (verbose) cat("[ARTEMIS] Computing per-sample peak annotation counts...\n")

  count_list <- lapply(names(annotated_list), function(sid) {
    df <- annotated_list[[sid]]

    if (!feature_col %in% colnames(df)) {
      if (feature_col == "annotation_simple" && "annotation" %in% colnames(df)) {
        df[[feature_col]] <- sapply(df$annotation, function(x) {
          if (grepl("Promoter",   x)) return("Promoter")
          if (grepl("5' UTR",     x)) return("5' UTR")
          if (grepl("3' UTR",     x)) return("3' UTR")
          if (grepl("Exon",       x)) return("Exon")
          if (grepl("Intron",     x)) return("Intron")
          if (grepl("Downstream", x)) return("Downstream")
          if (grepl("Intergenic", x)) return("Intergenic")
          return("Other")
        })
      } else {
        stop("Column '", feature_col, "' not found for sample '", sid, "'.", call. = FALSE)
      }
    }

    tab <- table(df[[feature_col]])
    data.frame(
      sample_id = sid,
      feature   = names(tab),
      count     = as.numeric(tab),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  prop_df <- do.call(rbind, count_list)

  # Fill zeros for features absent in some samples
  all_features <- unique(prop_df$feature)
  all_samples  <- names(annotated_list)
  full_grid    <- expand.grid(sample_id  = all_samples,
                              feature    = all_features,
                              stringsAsFactors = FALSE)
  prop_df <- merge(full_grid, prop_df, by = c("sample_id", "feature"), all.x = TRUE)
  prop_df$count[is.na(prop_df$count)] <- 0L

  # Attach group labels
  prop_df <- merge(prop_df,
                   metadata[, c("sample_id", group_col), drop = FALSE],
                   by = "sample_id", all.x = TRUE)

  groups <- sort(unique(prop_df[[group_col]]))
  if (length(groups) != 2L)
    stop("Exactly 2 groups are required. Found: ",
         paste(groups, collapse = ", "), call. = FALSE)

  if (verbose) {
    cat("    Samples            : ", length(all_samples), "\n", sep = "")
    cat("    Groups             : ", paste(groups, collapse = " vs "), "\n", sep = "")
    cat("    Feature categories : ", length(all_features), "\n", sep = "")
    cat("    Test               : ", method, " + ", p_adj_method, " correction\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # Per-feature statistical test
  # ---------------------------------------------------------------------------
  if (verbose) cat("[ARTEMIS] Testing annotation proportions...\n")

  results <- do.call(rbind, lapply(all_features, function(feat) {
    sub_df <- prop_df[prop_df$feature == feat, ]
    g1 <- sub_df$count[sub_df[[group_col]] == groups[1]]
    g2 <- sub_df$count[sub_df[[group_col]] == groups[2]]

    res <- tryCatch({
      if (method == "wilcoxon") {
        wilcox.test(g1, g2, exact = FALSE)
      } else {
        t.test(g1, g2)
      }
    }, error = function(e) list(statistic = NA_real_, p.value = NA_real_))

    data.frame(
      feature   = feat,
      mean_g1   = mean(g1, na.rm = TRUE),
      mean_g2   = mean(g2, na.rm = TRUE),
      diff      = mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE),
      statistic = as.numeric(res$statistic),
      p.value   = res$p.value,
      stringsAsFactors = FALSE
    )
  }))

  colnames(results)[colnames(results) == "mean_g1"] <- paste0("mean_", groups[1])
  colnames(results)[colnames(results) == "mean_g2"] <- paste0("mean_", groups[2])
  colnames(results)[colnames(results) == "diff"]    <-
    paste0("diff_", groups[2], "_minus_", groups[1])

  results$p.adj     <- p.adjust(results$p.value, method = p_adj_method)
  results           <- results[order(results$p.value), ]
  rownames(results) <- NULL

  if (verbose) {
    cat("[ARTEMIS] Results summary:\n")
    sig <- results$feature[!is.na(results$p.adj) & results$p.adj < 0.05]
    if (length(sig) > 0L)
      cat("    Significant (p.adj < 0.05) : ", paste(sig, collapse = ", "), "\n", sep = "")
    else
      cat("    No features significant at p.adj < 0.05\n")
    cat("    Note: test is on raw counts; valid when peak libraries are balanced\n")
    cat("          (e.g., spike-in normalized) within the batch.\n")
    cat("\n")
    print(results[, c("feature",
                      paste0("mean_", groups[1]),
                      paste0("mean_", groups[2]),
                      "p.value", "p.adj")],
          digits = 3, row.names = FALSE)
    cat("\n")
  }

  list(
    results     = results,
    counts      = prop_df,
    groups      = groups,
    group_col   = group_col,
    feature_col = feature_col,
    method      = method
  )
}
