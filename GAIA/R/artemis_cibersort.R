# GAIA/Artemis/cibersort.R
# CIBERSORT immune cell deconvolution wrapper
#
# Wraps the original CIBERSORT.R source code (Newman et al., 2015) for
# estimating immune cell type proportions from bulk gene expression data.
#
# The original CIBERSORT.R source code is NOT modified in any way.
# This wrapper sources it and provides a GAIA-compatible interface.
#
# Requires: CIBERSORT.R source file and a signature matrix (e.g., LM22.txt)
# Dependencies (loaded by CIBERSORT.R): e1071, parallel, preprocessCore


#' Run CIBERSORT immune cell deconvolution
#'
#' Estimates cell type proportions from bulk gene expression data using
#' nu-support vector regression. Wraps the original CIBERSORT source code
#' with a GAIA-compatible interface that accepts R objects directly.
#'
#' @param mixture Expression matrix (genes as rows, samples as columns) or a
#'   file path to a tab-delimited mixture file. Gene names should be in rownames
#'   (or first column for files). Can also be an artemis_norm object, in which
#'   case norm_counts is extracted automatically.
#' @param cibersort_path Path to the original CIBERSORT.R source file.
#' @param sig_matrix Path to the signature matrix file (e.g., LM22.txt for
#'   22 immune cell types).
#' @param output_dir Directory where the mixture file and CIBERSORT results
#'   will be saved. A "cibersort" subdirectory is created within it.
#'   Default: "data/cibersort".
#' @param mixture_name Character. Base name for the saved mixture file
#'   (without extension). Default: "mixture". The file is saved as
#'   `<output_dir>/cibersort/<mixture_name>.txt`.
#' @param perm Number of permutations for p-value calculation. Default 100
#' @param QN Logical. Apply quantile normalization. Default TRUE for microarray
#'   data. Set to FALSE for RNA-seq data.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return S3 object of class "artemis_cibersort" with components:
#'   \describe{
#'     \item{proportions}{Matrix of cell type proportions (samples x cell types).
#'       Rows sum to 1.}
#'     \item{p_values}{Numeric vector of per-sample p-values (NULL if perm=0)}
#'     \item{correlations}{Numeric vector of Pearson correlations per sample}
#'     \item{rmse}{Numeric vector of RMSE per sample}
#'     \item{raw}{The full raw CIBERSORT output matrix}
#'     \item{metadata}{List with parameters and run info}
#'   }
#'
#' @details
#' CIBERSORT uses nu-SVR to deconvolve bulk expression into cell type
#' proportions using a reference signature matrix. The standard LM22 matrix
#' contains signatures for 22 human immune cell types.
#'
#' For RNA-seq data, set QN=FALSE. Quantile normalization is designed for
#' microarray data and can distort RNA-seq distributions.
#'
#' The mixture file is saved to `<output_dir>/cibersort/<mixture_name>.txt`
#' for reproducibility and re-use. CIBERSORT's side-effect output
#' (CIBERSORT-Results.txt) is also written to the same directory.
#'
#' @examples
#' \dontrun{
#' result <- ARTEMIS_cibersort(
#'   mixture        = norm_counts,
#'   cibersort_path = "data/CIBERSORT.R",
#'   sig_matrix     = "data/LM22.txt",
#'   output_dir     = "results",
#'   mixture_name   = "my_experiment",
#'   perm           = 100,
#'   QN             = FALSE  # RNA-seq
#' )
#'
#' # View proportions
#' head(result$proportions)
#'
#' # From an artemis_norm object
#' result <- ARTEMIS_cibersort(
#'   mixture        = my_norm,
#'   cibersort_path = "data/CIBERSORT.R",
#'   sig_matrix     = "data/LM22.txt",
#'   output_dir     = "results",
#'   QN             = FALSE
#' )
#'
#' }
#' @export
ARTEMIS_cibersort <- function(mixture,
                               cibersort_path,
                               sig_matrix,
                               output_dir = "data",
                               mixture_name = "mixture",
                               perm = 100,
                               QN = TRUE,
                               verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  if (!file.exists(cibersort_path)) {
    stop("CIBERSORT.R source file not found: ", cibersort_path, call. = FALSE)
  }

  if (!file.exists(sig_matrix)) {
    stop("Signature matrix file not found: ", sig_matrix, call. = FALSE)
  }

  # Create cibersort output subdirectory
  cibersort_dir <- file.path(output_dir, "cibersort")
  if (!dir.exists(cibersort_dir)) {
    dir.create(cibersort_dir, recursive = TRUE)
  }

  # ---------------------------------------------------------------------------
  # Source CIBERSORT.R (defines CIBERSORT, CoreAlg, doPerm functions)
  # ---------------------------------------------------------------------------
  if (verbose) cat("Sourcing CIBERSORT.R...\n")
  source(cibersort_path, local = TRUE)

  # Verify the function was loaded
  if (!exists("CIBERSORT", inherits = FALSE)) {
    stop("CIBERSORT function not found after sourcing: ", cibersort_path,
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Handle mixture input
  # ---------------------------------------------------------------------------
  if (is.character(mixture) && length(mixture) == 1 && file.exists(mixture)) {
    # File path provided directly
    mixture_file <- normalizePath(mixture, mustWork = TRUE)
    if (verbose) cat("Using mixture file:", mixture_file, "\n")

  } else {
    # R object — extract matrix if needed
    if (inherits(mixture, "artemis_norm")) {
      if (verbose) cat("Extracting norm_counts from artemis_norm object\n")
      mix_mat <- mixture$norm_counts
    } else if (inherits(mixture, "artemis_ts_norm")) {
      if (verbose) cat("Extracting norm_counts from artemis_ts_norm object\n")
      mix_mat <- mixture$norm_counts
    } else if (is.matrix(mixture) || is.data.frame(mixture)) {
      mix_mat <- as.matrix(mixture)
    } else {
      stop("mixture must be a matrix, data.frame, artemis_norm object, or file path.",
           call. = FALSE)
    }

    if (is.null(rownames(mix_mat))) {
      stop("mixture matrix must have gene names as rownames.", call. = FALSE)
    }

    if (verbose) {
      cat("Mixture:", nrow(mix_mat), "genes x", ncol(mix_mat), "samples\n")
    }

    # Write mixture to cibersort output directory
    mixture_file <- file.path(cibersort_dir, paste0(mixture_name, ".txt"))
    utils::write.table(
      data.frame(GeneSymbol = rownames(mix_mat), mix_mat, check.names = FALSE),
      file = mixture_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    mixture_file <- normalizePath(mixture_file, mustWork = TRUE)

    if (verbose) cat("Mixture file saved:", mixture_file, "\n")
  }

  # ---------------------------------------------------------------------------
  # Run CIBERSORT
  # (CIBERSORT writes "CIBERSORT-Results.txt" as a side effect — capture it
  #  in the cibersort output directory)
  # ---------------------------------------------------------------------------
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(cibersort_dir)

  if (verbose) {
    cat("Running CIBERSORT (perm =", perm, ", QN =", QN, ")...\n")
  }

  start_time <- proc.time()

  raw_result <- CIBERSORT(
    sig_matrix   = sig_matrix,
    mixture_file = mixture_file,
    perm         = perm,
    QN           = QN
  )

  elapsed <- (proc.time() - start_time)["elapsed"]

  if (verbose) cat("CIBERSORT completed in", round(elapsed, 1), "seconds.\n")

  # ---------------------------------------------------------------------------
  # Parse results
  # ---------------------------------------------------------------------------
  result_cols <- colnames(raw_result)

  # Identify cell type columns (everything except P-value, Correlation, RMSE)
  meta_cols <- c("P-value", "Correlation", "RMSE")
  cell_type_cols <- setdiff(result_cols, meta_cols)

  proportions <- raw_result[, cell_type_cols, drop = FALSE]

  p_values <- NULL
  if ("P-value" %in% result_cols) {
    p_values <- raw_result[, "P-value"]
    names(p_values) <- rownames(raw_result)
  }

  correlations <- NULL
  if ("Correlation" %in% result_cols) {
    correlations <- raw_result[, "Correlation"]
    names(correlations) <- rownames(raw_result)
  }

  rmse <- NULL
  if ("RMSE" %in% result_cols) {
    rmse <- raw_result[, "RMSE"]
    names(rmse) <- rownames(raw_result)
  }

  if (verbose) {
    cat("Results: ", nrow(proportions), " samples x ",
        ncol(proportions), " cell types\n", sep = "")
    if (!is.null(correlations)) {
      cat("Mean correlation:", round(mean(correlations, na.rm = TRUE), 3), "\n")
      cat("Mean RMSE:", round(mean(rmse, na.rm = TRUE), 3), "\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Build result object
  # ---------------------------------------------------------------------------
  result <- list(
    proportions  = proportions,
    p_values     = p_values,
    correlations = correlations,
    rmse         = rmse,
    raw          = raw_result,
    metadata     = list(
      sig_matrix       = normalizePath(sig_matrix, mustWork = FALSE),
      mixture_file     = mixture_file,
      output_dir       = normalizePath(cibersort_dir, mustWork = FALSE),
      perm             = perm,
      QN               = QN,
      n_samples        = nrow(proportions),
      n_cell_types     = ncol(proportions),
      cell_types       = cell_type_cols,
      elapsed_seconds  = as.numeric(elapsed)
    )
  )
  class(result) <- c("artemis_cibersort", "list")

  return(result)
}


#' @method print artemis_cibersort
#' @export
print.artemis_cibersort <- function(x, ...) {
  cat("CIBERSORT Deconvolution Result\n")
  cat("  Samples:", x$metadata$n_samples, "\n")
  cat("  Cell types:", x$metadata$n_cell_types, "\n")
  cat("  Permutations:", x$metadata$perm, "\n")
  cat("  QN:", x$metadata$QN, "\n")

  if (!is.null(x$correlations)) {
    cat("  Mean correlation:", round(mean(x$correlations, na.rm = TRUE), 3), "\n")
    cat("  Mean RMSE:", round(mean(x$rmse, na.rm = TRUE), 3), "\n")
  }

  # Show top cell types by mean proportion
  mean_props <- sort(colMeans(x$proportions), decreasing = TRUE)
  top <- head(mean_props[mean_props > 0.01], 5)
  if (length(top) > 0) {
    cat("\n  Top cell types (mean proportion):\n")
    for (i in seq_along(top)) {
      cat("    ", names(top)[i], ": ", round(top[i] * 100, 1), "%\n", sep = "")
    }
  }

  invisible(x)
}
