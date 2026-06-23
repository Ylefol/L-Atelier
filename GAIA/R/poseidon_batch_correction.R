###############################################################################
########################### Batch Correction #################################
###############################################################################
# ComBat/limma-based batch correction for two distinct data structures:
#   - POSEIDON_correct_batch_counts()    - $counts + $targets (quant_result)
#   - POSEIDON_correct_batch_massspec()  - $wide + $sample_meta (massspec_data)
# Plus the accompanying before/after PCA visualization for quant_result.

#' Batch Correction for Count Data
#'
#' @description Removes batch effects from count data using ComBat or limma.
#' Preserves biological variation while removing technical batch effects.
#'
#' For mass spec data (\code{massspec_data} objects, already log2-scale with
#' NAs), use \code{POSEIDON_correct_batch_massspec()} instead -- it does not
#' share this function's count-data assumptions (raw counts, log2(x+1)
#' transform, no missingness).
#'
#' @param quant_result A list containing:
#'   \itemize{
#'     \item counts: Matrix of counts (features x samples)
#'     \item targets: Data.frame with sample metadata including batch information
#'   }
#' @param batch_col Character string. Column name in targets containing batch
#'   information (default = "batch").
#' @param group_col Character string. Optional column name for biological groups
#'   to preserve during correction (default = "group").
#' @param method Character string. Batch correction method:
#'   "combat" (default) or "limma".
#' @param log_transform Logical. Should data be log-transformed before ComBat?
#'   ComBat assumes approximately normal data. If TRUE, applies log2(counts + 1)
#'   before correction and back-transforms after (default = TRUE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return The input list with batch-corrected counts matrix.
#'
#' @details
#' Two methods are available:
#'
#' \strong{ComBat} (sva package):
#' \itemize{
#'   \item Empirical Bayes framework for batch correction
#'   \item Recommended for most applications
#'   \item Can preserve biological groups during correction
#'   \item Assumes approximately normal data (use log_transform = TRUE for counts)
#' }
#'
#' \strong{limma::removeBatchEffect}:
#' \itemize{
#'   \item Linear model-based batch removal
#'   \item Faster than ComBat
#'   \item Also preserves biological groups
#'   \item Works well with log-transformed data
#' }
#'
#' For count data (ATAC-seq, RNA-seq), log transformation is recommended
#' before batch correction since these methods assume approximately normal
#' distributions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic batch correction
#' counts_corrected <- POSEIDON_correct_batch_counts(atac_counts)
#'
#' # Preserve group differences while correcting batch
#' counts_corrected <- POSEIDON_correct_batch_counts(atac_counts,
#'                                             batch_col = "batch",
#'                                             group_col = "group")
#'
#' # Use limma instead of ComBat
#' counts_corrected <- POSEIDON_correct_batch_counts(atac_counts, method = "limma")
#'
#' }
POSEIDON_correct_batch_counts <- function(quant_result,
                                    batch_col = "batch",
                                    group_col = "group",
                                    method = "combat",
                                    log_transform = TRUE,
                                    verbose = TRUE) {

  # Validate method

  method <- tolower(method)
  if (!method %in% c("combat", "limma")) {
    stop("method must be either 'combat' or 'limma'")
  }

  # Check required packages
  if (method == "combat" && !requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBat batch correction. ",
         "Install with: BiocManager::install('sva')")
  }

  if (method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for limma batch correction. ",
         "Install with: BiocManager::install('limma')")
  }

  # Extract data
  counts <- quant_result$counts
  targets <- quant_result$targets

  # Check batch column exists
  if (!batch_col %in% colnames(targets)) {
    stop("Batch column '", batch_col, "' not found in targets. ",
         "Available columns: ", paste(colnames(targets), collapse = ", "))
  }

  batch <- targets[[batch_col]]

  # Check we have multiple batches
  n_batches <- length(unique(batch))
  if (n_batches < 2) {
    warning("Only one batch found. No batch correction needed.")
    return(quant_result)
  }

  if (verbose) {
    cat("[POSEIDON] Batch correction using", method, "\n")
    cat("  Batches:", paste(unique(batch), collapse = ", "), "\n")
    cat("  Samples per batch:\n")
    batch_counts <- table(batch)
    for (b in names(batch_counts)) {
      cat("    Batch", b, ":", batch_counts[b], "samples\n")
    }
  }

  # Check for group to preserve
  mod <- NULL
  if (!is.null(group_col) && group_col %in% colnames(targets)) {
    group <- factor(targets[[group_col]])
    mod <- model.matrix(~ group)
    if (verbose) {
      cat("  Preserving group:", group_col, "\n")
      cat("    Groups:", paste(levels(group), collapse = ", "), "\n")
    }
  }

  # Log transform if requested
  if (log_transform) {
    if (verbose) cat("  Log-transforming data (log2(x + 1))...\n")
    counts_input <- log2(counts + 1)
  } else {
    counts_input <- counts
  }

  # Apply batch correction
  if (verbose) cat("  Running batch correction...\n")

  if (method == "combat") {
    # ComBat batch correction
    counts_corrected <- sva::ComBat(
      dat = counts_input,
      batch = batch,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
  } else {
    # limma batch correction
    if (!is.null(mod)) {
      # Remove batch effect while preserving group
      design <- mod
      counts_corrected <- limma::removeBatchEffect(
        counts_input,
        batch = batch,
        design = design
      )
    } else {
      counts_corrected <- limma::removeBatchEffect(
        counts_input,
        batch = batch
      )
    }
  }

  # Back-transform if we log-transformed
  if (log_transform) {
    if (verbose) cat("  Back-transforming (2^x - 1)...\n")
    counts_corrected <- 2^counts_corrected - 1
    # Ensure no negative values (can happen with correction)
    counts_corrected[counts_corrected < 0] <- 0
  }

  if (verbose) {
    cat("[POSEIDON] Batch correction complete.\n")
    cat("  Output dimensions:", nrow(counts_corrected), "x",
        ncol(counts_corrected), "\n")
  }

  # Update result
  quant_result$counts <- counts_corrected
  quant_result$batch_corrected <- TRUE
  quant_result$batch_method <- method

  return(quant_result)
}


#' Batch Correction for Mass Spectrometry Data
#'
#' @description Removes batch effects from a \code{massspec_data} object's
#' log2 intensity matrix (\code{$wide}) using ComBat or
#' \code{limma::removeBatchEffect}. Unlike \code{POSEIDON_correct_batch_counts()},
#' this operates on data that is already log2-transformed and contains NAs --
#' no log/back-transform is applied, and missingness is handled via a
#' disposable imputed working copy used only to make the underlying algorithm
#' runnable.
#'
#' @param massspec_data A \code{massspec_data} object (from
#'   \code{ELEUTHIA_load_massspec()}), typically after
#'   \code{HADES_filter_massspec()}/\code{HADES_filter_massspec_proteins()}
#'   and \code{POSEIDON_normalize_massspec()}.
#' @param batch_col Character string. Column in \code{$sample_meta} containing
#'   batch information (default = \code{"PlateID"}).
#' @param group_col Character string or NULL. Optional column in
#'   \code{$sample_meta} for a biological covariate to preserve during
#'   correction (default = NULL, no covariate preserved). Samples with NA in
#'   this column (e.g. POOL samples, which have no biological group) are
#'   assigned their own \code{"no_group"} level rather than being dropped, so
#'   every sample is still corrected.
#' @param method Character string. \code{"combat"} (default) or
#'   \code{"limma"}.
#' @param impute_method Character string. How the disposable working copy
#'   used to make ComBat/limma runnable is filled in:
#'   \itemize{
#'     \item \code{"batchmean"} (default) -- each missing cell is filled with
#'       its protein's mean WITHIN its own batch (observed cells only).
#'       Falls back to the protein's global mean for any protein/batch
#'       combination with zero observed values in that batch.
#'     \item \code{"global"} -- each missing cell is filled with its
#'       protein's mean across ALL samples, regardless of batch (the
#'       original behavior of this function).
#'   }
#'   Comparing the two on the SEPSOMICS dataset (see
#'   \code{projects/SEPSOMICS/massSpec_imputation_comparison.R}) showed
#'   \code{"global"} under-corrects proteins with high missingness in a
#'   given batch -- the placeholder pulls that batch's estimate toward the
#'   grand mean instead of its own, diluting the estimated batch effect.
#'   \code{"batchmean"} does not have this failure mode, but that was
#'   established on one dataset/filtering setup; \code{"global"} is kept
#'   available in case a different experiment's missingness pattern behaves
#'   differently.
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return The input \code{massspec_data} object with \code{$wide}
#'   batch-corrected. Values that were NA before correction remain NA in the
#'   output -- imputation is used only internally to make ComBat/limma
#'   runnable; the imputed values themselves are never written back.
#'
#' @details
#' DIA mass spec retains substantial missingness even after filtering, and
#' both ComBat and \code{limma::removeBatchEffect} require a complete matrix.
#' A temporary copy of \code{$wide} is imputed (per \code{impute_method})
#' purely so the algorithm can run; the correction is computed on that copy,
#' but only the originally non-NA cells of the result are written back.
#' Cells that were NA before correction stay NA.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_filter_massspec(ms_raw, keep_types = c("SAMPLE", "POOL"))
#' ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide)
#' ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID")
#'
#' # Preserve a biological covariate during correction
#' ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID", group_col = "Group")
#'
#' # Use the original global-mean imputation instead
#' ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID", impute_method = "global")
#' }
POSEIDON_correct_batch_massspec <- function(massspec_data,
                                             batch_col = "PlateID",
                                             group_col = NULL,
                                             method    = "combat",
                                             impute_method = "batchmean",
                                             verbose   = TRUE) {

  if (!inherits(massspec_data, "massspec_data"))
    stop("massspec_data must be a massspec_data object from ELEUTHIA_load_massspec().")

  method <- tolower(method)
  if (!method %in% c("combat", "limma"))
    stop("method must be either 'combat' or 'limma'")

  impute_method <- tolower(impute_method)
  if (!impute_method %in% c("batchmean", "global"))
    stop("impute_method must be either 'batchmean' or 'global'")

  if (method == "combat" && !requireNamespace("sva", quietly = TRUE))
    stop("Package 'sva' is required for ComBat batch correction. ",
         "Install with: BiocManager::install('sva')")
  if (method == "limma" && !requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' is required for limma batch correction. ",
         "Install with: BiocManager::install('limma')")

  wide    <- massspec_data$wide
  smeta   <- massspec_data$sample_meta
  sid_col <- massspec_data$params$sample_col
  if (is.null(sid_col)) sid_col <- "SampleID"

  if (!batch_col %in% colnames(smeta))
    stop("Batch column '", batch_col, "' not found in sample_meta. ",
         "Available columns: ", paste(colnames(smeta), collapse = ", "))

  # Align sample_meta rows to the matrix's column order
  smeta <- smeta[match(colnames(wide), smeta[[sid_col]]), ]

  batch     <- smeta[[batch_col]]
  n_batches <- length(unique(batch))
  if (n_batches < 2) {
    warning("Only one batch found. No batch correction performed.")
    return(massspec_data)
  }

  if (verbose) {
    cat("[POSEIDON] Batch correction (mass spec) using", method, "\n")
    cat("    Batches:", paste(unique(batch), collapse = ", "), "\n")
    cat("    Samples per batch:\n")
    batch_counts <- table(batch)
    for (b in names(batch_counts))
      cat("        ", b, ": ", batch_counts[b], " samples\n", sep = "")
  }

  # Build covariate model matrix, if requested. NA (e.g. POOL samples with no
  # biological group) get their own level rather than being dropped, so the
  # design matrix still covers every sample in `wide`.
  mod <- NULL
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(smeta))
      stop("Group column '", group_col, "' not found in sample_meta.")
    group_raw <- as.character(smeta[[group_col]])
    group_raw[is.na(group_raw)] <- "no_group"
    group <- factor(group_raw)
    mod <- model.matrix(~ group)
    if (verbose) {
      cat("    Preserving group:", group_col, "\n")
      cat("        Groups:", paste(levels(group), collapse = ", "), "\n")
    }
  }

  # ComBat/limma require a complete matrix -- impute a disposable working
  # copy purely to make the algorithm runnable. Corrected values are written
  # back only at originally non-NA positions; real NAs stay NA.
  na_idx <- which(is.na(wide), arr.ind = TRUE)

  if (impute_method == "global") {
    wide_imputed <- wide
    if (nrow(na_idx) > 0) {
      protein_means <- rowMeans(wide, na.rm = TRUE)
      wide_imputed[na_idx] <- protein_means[na_idx[, 1]]
    }
  } else {
    # "batchmean": fill each missing cell with its protein's mean WITHIN its
    # own batch (observed cells only). Falls back to the protein's global
    # mean for any protein/batch combination with zero observed values.
    wide_imputed  <- wide
    protein_means <- rowMeans(wide, na.rm = TRUE)
    n_fallback    <- 0
    for (b in unique(batch)) {
      cols_b <- which(batch == b)
      sub    <- wide[, cols_b, drop = FALSE]
      batch_means  <- rowMeans(sub, na.rm = TRUE)
      fallback_idx <- is.na(batch_means)
      n_fallback   <- n_fallback + sum(fallback_idx)
      batch_means[fallback_idx] <- protein_means[fallback_idx]
      na_in_block <- is.na(sub)
      if (any(na_in_block)) {
        sub[na_in_block] <- batch_means[row(sub)[na_in_block]]
        wide_imputed[, cols_b] <- sub
      }
    }
    if (verbose && n_fallback > 0)
      cat("    ", n_fallback, " protein/batch combinations had zero observed ",
          "values in that batch (fell back to global mean)\n", sep = "")
  }

  if (verbose)
    cat("    Imputed", nrow(na_idx), "missing values for correction using '",
        impute_method, "' (", round(100 * nrow(na_idx) / length(wide), 1),
        "% of matrix); NAs restored after correction\n", sep = "")

  if (verbose) cat("    Running batch correction...\n")

  if (method == "combat") {
    corrected_full <- sva::ComBat(
      dat = wide_imputed,
      batch = batch,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
  } else {
    if (!is.null(mod)) {
      corrected_full <- limma::removeBatchEffect(wide_imputed, batch = batch, design = mod)
    } else {
      corrected_full <- limma::removeBatchEffect(wide_imputed, batch = batch)
    }
  }

  # Restore real NAs; keep corrected values only where data was observed.
  corrected <- wide
  observed  <- !is.na(wide)
  corrected[observed] <- corrected_full[observed]

  if (verbose) {
    cat("[POSEIDON] Batch correction complete.\n")
    cat("    Output dimensions:", nrow(corrected), "x", ncol(corrected), "\n")
  }

  massspec_data$wide <- corrected
  massspec_data$batch_corrected <- TRUE
  massspec_data$batch_method    <- method

  return(massspec_data)
}


#' Visualize Batch Effects (Before/After Correction)
#'
#' @description Creates PCA plots to visualize batch effects before and
#' optionally after batch correction.
#'
#' @param quant_result A list containing counts and targets.
#' @param batch_col Character string. Column name for batch (default = "batch").
#' @param group_col Character string. Column name for group (default = "group").
#' @param log_transform Logical. Log-transform counts for PCA (default = TRUE).
#' @param title Character string. Plot title (default = "PCA - Batch Effect").
#'
#' @return A ggplot object showing PCA colored by batch and shaped by group.
#'
#' @export
#'
POSEIDON_plot_batch_pca <- function(quant_result,
                                     batch_col = "batch",
                                     group_col = "group",
                                     log_transform = TRUE,
                                     title = "PCA - Batch Effect") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  counts <- quant_result$counts
  targets <- quant_result$targets

  # Log transform if requested
  if (log_transform) {
    counts <- log2(counts + 1)
  }

  # Run PCA on transposed data (samples as rows)
  pca_result <- prcomp(t(counts), scale. = TRUE, center = TRUE)

  # Create plot data
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample_id = rownames(pca_result$x)
  )

  # Add batch and group info
  pca_data$batch <- factor(targets[[batch_col]])

  if (!is.null(group_col) && group_col %in% colnames(targets)) {
    pca_data$group <- factor(targets[[group_col]])
  } else {
    pca_data$group <- factor("all")
  }

  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

  # Create plot
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2,
                                               color = batch,
                                               shape = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      color = "Batch",
      shape = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}
