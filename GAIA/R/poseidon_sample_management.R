###############################################################################
########### Sample/Feature Management for quant_result Objects ##############
###############################################################################
# Operations on the $counts + $targets ("quant_result") structure used by
# the ATAC-seq/ChIP-seq/RNA-seq quantification pipelines.

#' Filter Low-Count Features (Peaks/Genes)
#'
#' @description Filters features (peaks, genes, etc.) with insufficient counts
#' across samples. Commonly used to remove unreliable features before
#' differential analysis or integration.
#'
#' @param quant_result A list containing at minimum:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item annotation: Data.frame with feature annotations (optional)
#'   }
#' @param min_count Integer. Minimum count threshold (default = 10).
#' @param min_samples Integer. Minimum number of samples that must meet
#'   min_count threshold (default = 2).
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return The input list with filtered counts and annotation.
#'
#' @details
#' A feature is retained if at least \code{min_samples} samples have
#' counts >= \code{min_count}. This filtering removes features that are:
#' \itemize{
#'   \item Not detected in most samples
#'   \item Too low-count for reliable statistical analysis
#'   \item Likely to introduce noise in downstream analyses
#' }
#'
#' Common thresholds:
#' \itemize{
#'   \item ATAC-seq peaks: min_count = 10, min_samples = 3
#'   \item RNA-seq genes: min_count = 10, min_samples = 2-3
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter peaks with at least 10 counts in at least 3 samples
#' counts_filtered <- POSEIDON_filter_low_counts(counts, min_count = 10, min_samples = 3)
#'
#' }
POSEIDON_filter_low_counts <- function(quant_result,
                                        min_count = 10,
                                        min_samples = 2,
                                        verbose = TRUE) {

  counts <- quant_result$counts

  # Count how many samples meet threshold for each feature
  n_above_threshold <- rowSums(counts >= min_count)

  # Filter

  keep <- n_above_threshold >= min_samples

  n_before <- nrow(counts)
  n_after <- sum(keep)

  if (verbose) {
    cat("[POSEIDON] Filtering low-count features:\n")
    cat("  Before:", n_before, "features\n")
    cat("  After:", n_after, "features\n")
    cat("  Removed:", n_before - n_after, "features\n")
    cat("  (min_count =", min_count, ", min_samples =", min_samples, ")\n")
  }

  # Update result
  quant_result$counts <- counts[keep, , drop = FALSE]

  # Update annotation if present
  if (!is.null(quant_result$annotation)) {
    quant_result$annotation <- quant_result$annotation[keep, , drop = FALSE]
  }

  return(quant_result)
}


#' Filter Samples by Group
#'
#' @description Filters a count dataset to include only samples belonging
#' to specified group(s). Works with the standard quant_result structure.
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata
#'   }
#' @param group_col Character string. Column name in targets containing
#'   group information (default = "group").
#' @param groups Character vector. Group value(s) to keep.
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered quant_result with only samples from specified group(s).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter for WT samples only
#' wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")
#'
#' # Filter for multiple groups
#' subset_data <- POSEIDON_filter_by_group(data, groups = c("WT", "Control"))
#'
#' }
POSEIDON_filter_by_group <- function(quant_result,
                                      group_col = "group",
                                      groups,
                                      verbose = TRUE) {

  targets <- quant_result$targets
  counts <- quant_result$counts

  # Check group column exists
 if (!group_col %in% colnames(targets)) {
    stop("Group column '", group_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  # Find samples in specified groups
  keep_samples <- targets[[group_col]] %in% groups
  n_before <- nrow(targets)
  n_after <- sum(keep_samples)

  if (n_after == 0) {
    stop("No samples found in group(s): ", paste(groups, collapse = ", "))
  }

  if (verbose) {
    cat("[POSEIDON] Filtering by group:\n")
    cat("  Keeping groups:", paste(groups, collapse = ", "), "\n")
    cat("  Samples before:", n_before, "\n")
    cat("  Samples after:", n_after, "\n")
  }

  # Filter counts and targets
  quant_result$counts <- counts[, keep_samples, drop = FALSE]
  quant_result$targets <- targets[keep_samples, , drop = FALSE]

  # Update annotation if present (no change needed, features stay the same)

  return(quant_result)
}


#' Filter Samples by Metadata Value
#'
#' @description General function to filter samples by any metadata column.
#' More flexible than filter_by_group.
#'
#' @param quant_result A list containing counts and targets.
#' @param column Character string. Column name in targets to filter by.
#' @param values Vector. Value(s) to keep.
#' @param verbose Logical. Print filtering summary (default = TRUE).
#'
#' @return Filtered quant_result.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter by batch
#' batch1_data <- POSEIDON_filter_by_metadata(data, column = "batch", values = 1)
#'
#' # Filter by biological replicate
#' bio1_data <- POSEIDON_filter_by_metadata(data, column = "bio_rep", values = c(1, 2))
#'
#' }
POSEIDON_filter_by_metadata <- function(quant_result,
                                         column,
                                         values,
                                         verbose = TRUE) {

  targets <- quant_result$targets
  counts <- quant_result$counts

  if (!column %in% colnames(targets)) {
    stop("Column '", column, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  keep_samples <- targets[[column]] %in% values
  n_before <- nrow(targets)
  n_after <- sum(keep_samples)

  if (n_after == 0) {
    stop("No samples found with ", column, " in: ", paste(values, collapse = ", "))
  }

  if (verbose) {
    cat("[POSEIDON] Filtering by", column, ":\n")
    cat("  Keeping values:", paste(values, collapse = ", "), "\n")
    cat("  Samples before:", n_before, "\n")
    cat("  Samples after:", n_after, "\n")
  }

  quant_result$counts <- counts[, keep_samples, drop = FALSE]
  quant_result$targets <- targets[keep_samples, , drop = FALSE]

  return(quant_result)
}


#' Match Samples Across Datasets
#'
#' @description Identifies and filters datasets to include only samples
#' that are present in all datasets, based on a matching column (e.g., bio_rep).
#'
#' @param data_list Named list of quant_result objects (each with counts and targets).
#' @param match_col Character string. Column in targets to match on
#'   (default = "bio_rep").
#' @param verbose Logical. Print matching summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{datasets}{Named list of filtered quant_results with matched samples}
#'   \item{matched_values}{Vector of matched values (e.g., bio_rep IDs)}
#'   \item{sample_map}{Data.frame showing sample correspondence across datasets}
#' }
#'
#' @details
#' This function finds the intersection of values in match_col across all
#' datasets and filters each dataset to include only samples with those values.
#'
#' Note: This does NOT reorder samples. Use POSEIDON_align_samples() after
#' matching to ensure consistent sample ordering across datasets.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_list <- list(
#'   ATAC = atac_counts,
#'   ChIP = chip_counts,
#'   RNA = rna_data
#' )
#' matched <- POSEIDON_match_samples(data_list, match_col = "bio_rep")
#'
#' }
POSEIDON_match_samples <- function(data_list,
                                    match_col = "bio_rep",
                                    verbose = TRUE) {

  if (!is.list(data_list) || length(data_list) < 2) {
    stop("data_list must be a list with at least 2 datasets")
  }

  dataset_names <- names(data_list)
  if (is.null(dataset_names)) {
    dataset_names <- paste0("Dataset_", seq_along(data_list))
    names(data_list) <- dataset_names
  }

  # Get values of match_col for each dataset
  values_per_dataset <- lapply(data_list, function(d) {
    if (!match_col %in% colnames(d$targets)) {
      stop("Column '", match_col, "' not found in targets")
    }
    unique(d$targets[[match_col]])
  })

  if (verbose) {
    cat("[POSEIDON] Sample matching by:", match_col, "\n\n")
    cat("[POSEIDON] Values per dataset:\n")
    for (nm in dataset_names) {
      cat("  ", nm, ":", paste(sort(values_per_dataset[[nm]]), collapse = ", "), "\n")
    }
  }

  # Find intersection
  matched_values <- Reduce(intersect, values_per_dataset)

  if (length(matched_values) == 0) {
    stop("No matching values found across all datasets")
  }

  if (verbose) {
    cat("\n[POSEIDON] Matched values:", paste(sort(matched_values), collapse = ", "), "\n")
    cat("[POSEIDON] Number of matched", match_col, ":", length(matched_values), "\n\n")
  }

  # Filter each dataset to matched values
  filtered_list <- list()

  for (nm in dataset_names) {
    filtered_list[[nm]] <- POSEIDON_filter_by_metadata(
      data_list[[nm]],
      column = match_col,
      values = matched_values,
      verbose = verbose
    )
    if (verbose) cat("\n")
  }

  # Create sample map showing correspondence
  sample_map <- data.frame(match_value = matched_values)
  colnames(sample_map) <- match_col

  for (nm in dataset_names) {
    targets <- filtered_list[[nm]]$targets
    # Get sample counts per match value
    samples_per_value <- sapply(matched_values, function(v) {
      sum(targets[[match_col]] == v)
    })
    sample_map[[paste0(nm, "_n")]] <- samples_per_value
  }

  if (verbose) {
    cat("[POSEIDON] Sample map:\n")
    print(sample_map)
  }

  return(list(
    datasets = filtered_list,
    matched_values = matched_values,
    sample_map = sample_map
  ))
}


#' Average Technical Replicates
#'
#' @description Aggregates technical replicates by averaging counts within
#' each biological replicate. Produces one sample per unique combination of
#' grouping variables (excluding tech_rep).
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata
#'   }
#' @param bio_rep_col Character string. Column identifying biological replicates
#'   (default = "bio_rep").
#' @param tech_rep_col Character string. Column identifying technical replicates
#'   to average over (default = "tech_rep").
#' @param group_by Character vector. Additional columns to group by when
#'   averaging. Samples are averaged within each unique combination of
#'   bio_rep_col + group_by columns. Default includes "group" and "batch"
#'   if present in targets.
#' @param method Character string. Aggregation method: "mean" (default) or
#'   "median".
#' @param verbose Logical. Print summary (default = TRUE).
#'
#' @return A quant_result with:
#'   \itemize{
#'     \item counts: Averaged count matrix (features x biological samples)
#'     \item targets: Updated metadata with one row per biological sample
#'     \item n_tech_reps: Number of technical replicates averaged per sample
#'   }
#'
#' @details
#' This function averages technical replicates to produce one observation per
#' biological sample. This is appropriate when:
#' \itemize{
#'   \item Preparing data for integration methods requiring matched samples
#'   \item Reducing technical noise while preserving biological variation
#'   \item Creating a cleaner dataset for correlation-based analyses
#' }
#'
#' The function automatically detects which columns to use for grouping:
#' \itemize{
#'   \item Always groups by bio_rep_col
#'   \item Includes "group" if present (preserves experimental groups)
#'   \item Includes "batch" if present (keeps batches separate)
#'   \item Additional columns can be specified via group_by parameter
#' }
#'
#' Can be used before or after filtering by group - the function respects
#' whatever samples are present in the input.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic averaging
#' averaged <- POSEIDON_average_tech_reps(atac_counts)
#'
#' # Average after filtering for WT
#' wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")
#' wt_averaged <- POSEIDON_average_tech_reps(wt_data)
#'
#' # Average before filtering (groups will be preserved)
#' averaged <- POSEIDON_average_tech_reps(atac_counts)
#' wt_averaged <- POSEIDON_filter_by_group(averaged, groups = "WT")
#'
#' # Custom grouping
#' averaged <- POSEIDON_average_tech_reps(data, group_by = c("group", "batch", "treatment"))
#'
#' }
POSEIDON_average_tech_reps <- function(quant_result,
                                        bio_rep_col = "bio_rep",
                                        tech_rep_col = "tech_rep",
                                        group_by = NULL,
                                        method = "mean",
                                        verbose = TRUE) {

  counts <- quant_result$counts
  targets <- quant_result$targets

  # Validate columns exist
  if (!bio_rep_col %in% colnames(targets)) {
    stop("bio_rep_col '", bio_rep_col, "' not found in targets")
  }

  if (!tech_rep_col %in% colnames(targets)) {
    stop("tech_rep_col '", tech_rep_col, "' not found in targets")
  }

  # Determine grouping columns
  # Always include bio_rep_col, plus any standard columns that exist
  standard_group_cols <- c("group", "batch")
  auto_group_by <- standard_group_cols[standard_group_cols %in% colnames(targets)]

  if (is.null(group_by)) {
    group_cols <- c(bio_rep_col, auto_group_by)
  } else {
    group_cols <- unique(c(bio_rep_col, group_by))
  }

  # Remove tech_rep_col from grouping if accidentally included
  group_cols <- setdiff(group_cols, tech_rep_col)

  if (verbose) {
    cat("[POSEIDON] Averaging technical replicates:\n")
    cat("  Grouping by:", paste(group_cols, collapse = ", "), "\n")
    cat("  Averaging over:", tech_rep_col, "\n")
    cat("  Method:", method, "\n")
  }

  # Create grouping key for each sample
  group_key <- apply(targets[, group_cols, drop = FALSE], 1, paste, collapse = "_")

  # Get unique groups
  unique_groups <- unique(group_key)
  n_groups <- length(unique_groups)

  if (verbose) {
    cat("  Samples before:", ncol(counts), "\n")
    cat("  Samples after:", n_groups, "\n")
  }

  # Aggregate counts
  agg_func <- if (method == "mean") rowMeans else function(x) apply(x, 1, median)

  agg_counts <- matrix(
    nrow = nrow(counts),
    ncol = n_groups,
    dimnames = list(rownames(counts), NULL)
  )

  agg_targets <- data.frame(matrix(
    nrow = n_groups,
    ncol = length(group_cols),
    dimnames = list(NULL, group_cols)
  ), stringsAsFactors = FALSE)

  n_tech_reps <- numeric(n_groups)

  for (i in seq_along(unique_groups)) {
    grp <- unique_groups[i]
    sample_idx <- which(group_key == grp)

    # Average counts
    if (length(sample_idx) == 1) {
      agg_counts[, i] <- counts[, sample_idx]
    } else {
      agg_counts[, i] <- agg_func(counts[, sample_idx, drop = FALSE])
    }

    # Take first row's metadata for grouping columns
    agg_targets[i, ] <- targets[sample_idx[1], group_cols]

    # Record number of tech reps
    n_tech_reps[i] <- length(sample_idx)
  }

  # Create sample IDs for averaged data
  new_sample_ids <- apply(agg_targets, 1, paste, collapse = "_")
  colnames(agg_counts) <- new_sample_ids
  rownames(agg_targets) <- new_sample_ids
  agg_targets$sample_id <- new_sample_ids
  agg_targets$n_tech_reps <- n_tech_reps

  if (verbose) {
    cat("\n[POSEIDON] Technical replicates per biological sample:\n")
    tech_rep_summary <- table(n_tech_reps)
    for (n in names(tech_rep_summary)) {
      cat("  ", tech_rep_summary[n], "sample(s) with", n, "tech rep(s)\n")
    }
  }

  # Build result
  result <- list(
    counts = agg_counts,
    targets = agg_targets
  )

  # Preserve annotation if present
  if (!is.null(quant_result$annotation)) {
    result$annotation <- quant_result$annotation
  }

  # Add metadata
  result$averaged <- TRUE
  result$averaging_method <- method
  result$group_cols <- group_cols

  return(result)
}
