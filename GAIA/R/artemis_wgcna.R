# GAIA/Artemis/wgcna.R
# WGCNA (Weighted Gene Co-expression Network Analysis) functions
#
# Core workflow functions for network-based co-expression analysis.
# Integrates with ELEUTHIA data loading functions.


# ==============================================================================
# DATA PREPARATION
# ==============================================================================

#' Prepare data for WGCNA analysis
#'
#' Formats expression data and traits for WGCNA. Accepts either pre-loaded
#' data objects (from ELEUTHIA functions) or file paths.
#'
#' @param counts Expression matrix/data.frame OR path to count file/directory.
#'   If matrix/data.frame: genes as rows, samples as columns (will be transposed).
#'   If path: will be loaded based on counts_format parameter.
#' @param traits Data.frame OR path to traits file.
#'   Must have sample identifiers that match count data.
#' @param counts_format How to interpret counts if it's a path:
#'   - "file": single CSV/TSV file with genes as rows
#'   - "directory": folder of individual count files (one per sample)
#'   - "auto": attempt to detect (default when counts is a path)
#'   Ignored if counts is already a matrix/data.frame.
#' @param sample_col Column name in traits containing sample identifiers.
#'   If NULL, uses rownames of traits. Default: NULL.
#' @param gene_col Column name/index for gene identifiers in count data.
#'   Default: 1 (first column) or rownames.
#' @param transpose Logical. If TRUE, transpose counts so samples are rows.
#'   Default: TRUE (standard gene x sample input becomes sample x gene).
#' @param remove_zero_variance Logical. Remove genes with zero variance.
#'   Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_data" containing:
#'   \item{datExpr}{Expression matrix with samples as rows, genes as columns}
#'   \item{datTraits}{Traits data.frame with samples as rows, matching datExpr}
#'   \item{gene_names}{Character vector of gene names}
#'   \item{sample_names}{Character vector of sample names}
#'   \item{n_genes}{Number of genes}
#'   \item{n_samples}{Number of samples}
#'   \item{removed_genes}{Genes removed due to zero variance (if any)}
#'   \item{removed_samples}{Samples removed due to missing data (if any)}
#'
#' @examples
#' \dontrun{
#' # From pre-loaded ELEUTHIA data
#' counts <- ELEUTHIA_load_rnaseq_from_sheet(sample_sheet)
#' traits <- sample_sheet[, c("sample_id", "age", "condition")]
#' wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits, sample_col = "sample_id")
#'
#' # From file paths
#' wgcna_data <- ARTEMIS_wgcna_prepare(
#'   counts = "path/to/counts.csv",
#'   traits = "path/to/metadata.csv",
#'   sample_col = "sample_id"
#' )
#'
#' }
#' @export
ARTEMIS_wgcna_prepare <- function(counts,
                                   traits,
                                   counts_format = "auto",
                                   sample_col = NULL,
                                   gene_col = 1,
                                   transpose = TRUE,
                                   remove_zero_variance = TRUE,
                                   verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] ARTEMIS WGCNA Data Preparation \n")

  # --------------------------------------------------------------------------

  # Load counts if path provided
  # --------------------------------------------------------------------------
  removed_genes <- character(0)
  removed_samples <- character(0)

  if (is.character(counts) && length(counts) == 1) {
    if (verbose) cat("    Loading counts from path:", counts, "\n")
    counts <- .wgcna_load_counts(counts, counts_format, gene_col, verbose)
  }

  # Ensure counts is a matrix
  if (is.data.frame(counts)) {
    # Check if first column is gene names (character)
    if (is.character(counts[[1]]) || is.factor(counts[[1]])) {
      rownames(counts) <- counts[[1]]
      counts <- counts[, -1, drop = FALSE]
    }
    counts <- as.matrix(counts)
  }

  # Ensure numeric (double precision - required by WGCNA)
  if (!is.numeric(counts)) {
    counts <- apply(counts, 2, as.numeric)
  }
  storage.mode(counts) <- "double"

  # --------------------------------------------------------------------------
  # Load traits if path provided
  # --------------------------------------------------------------------------
  if (is.character(traits) && length(traits) == 1 && file.exists(traits)) {
    if (verbose) cat("    Loading traits from path:", traits, "\n")
    traits <- read.csv(traits, stringsAsFactors = FALSE)
  }

  # --------------------------------------------------------------------------
  # Align sample identifiers
  # --------------------------------------------------------------------------
  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(traits)) {
      stop("sample_col '", sample_col, "' not found in traits data")
    }
    rownames(traits) <- traits[[sample_col]]
    # Remove sample_col from traits (it's now rownames)
    traits <- traits[, colnames(traits) != sample_col, drop = FALSE]
  }

  # Transpose counts if needed (genes x samples -> samples x genes)
  if (transpose) {
    datExpr <- t(counts)
  } else {
    datExpr <- counts
  }

  # Get sample names from both datasets
  expr_samples <- rownames(datExpr)
  trait_samples <- rownames(traits)

  if (is.null(expr_samples)) {
    stop("Expression data must have sample names (rownames after transpose)")
  }
  if (is.null(trait_samples)) {
    stop("Traits data must have sample names (rownames or sample_col)")
  }

  # Find common samples
  common_samples <- intersect(expr_samples, trait_samples)

  if (length(common_samples) == 0) {
    stop("No matching samples found between expression and traits data.\n",
         "Expression samples: ", paste(head(expr_samples, 5), collapse = ", "),
         ifelse(length(expr_samples) > 5, "...", ""), "\n",
         "Trait samples: ", paste(head(trait_samples, 5), collapse = ", "),
         ifelse(length(trait_samples) > 5, "...", ""))
  }

  # Report missing samples
  missing_from_expr <- setdiff(trait_samples, expr_samples)
  missing_from_traits <- setdiff(expr_samples, trait_samples)

  if (length(missing_from_expr) > 0) {
    if (verbose) {
      cat("    Samples in traits but not in expression data (removed):\n")
      cat("    ", paste(missing_from_expr, collapse = ", "), "\n")
    }
    removed_samples <- c(removed_samples, missing_from_expr)
  }

  if (length(missing_from_traits) > 0) {
    if (verbose) {
      cat("    Samples in expression but not in traits (removed):\n")
      cat("        ", paste(missing_from_traits, collapse = ", "), "\n")
    }
    removed_samples <- c(removed_samples, missing_from_traits)
  }

  # Subset and align
  datExpr <- datExpr[common_samples, , drop = FALSE]
  datTraits <- traits[common_samples, , drop = FALSE]

  # --------------------------------------------------------------------------
  # Convert traits to numeric where possible
  # --------------------------------------------------------------------------
  datTraits <- .wgcna_prepare_traits(datTraits, verbose)

  # --------------------------------------------------------------------------
  # Remove zero-variance genes
  # --------------------------------------------------------------------------
  if (remove_zero_variance) {
    gene_vars <- apply(datExpr, 2, var, na.rm = TRUE)
    zero_var_genes <- names(gene_vars)[gene_vars == 0 | is.na(gene_vars)]

    if (length(zero_var_genes) > 0) {
      if (verbose) {
        cat("    Removing", length(zero_var_genes), "genes with zero variance\n")
      }
      datExpr <- datExpr[, !colnames(datExpr) %in% zero_var_genes, drop = FALSE]
      removed_genes <- zero_var_genes
    }
  }

  # --------------------------------------------------------------------------
  # Build result object
  # --------------------------------------------------------------------------
  result <- list(
    datExpr = datExpr,
    datTraits = datTraits,
    gene_names = colnames(datExpr),
    sample_names = rownames(datExpr),
    n_genes = ncol(datExpr),
    n_samples = nrow(datExpr),
    removed_genes = removed_genes,
    removed_samples = removed_samples
  )
  class(result) <- c("wgcna_data", "list")

  if (verbose) {
    cat("\n[ARTEMIS] Data Summary ---\n")
    cat("    Samples:", result$n_samples, "\n")
    cat("    Genes/Features:", result$n_genes, "\n")
    cat("    Traits:", ncol(datTraits), "\n")
    if (length(removed_genes) > 0) {
      cat("    Genes removed (zero variance):", length(removed_genes), "\n")
    }
    if (length(removed_samples) > 0) {
      cat("    Samples removed (unmatched):", length(removed_samples), "\n")
    }
  }

  return(result)
}


# ==============================================================================
# INTERNAL HELPER FUNCTIONS
# ==============================================================================

#' Load counts from file or directory
#' @keywords internal
.wgcna_load_counts <- function(path, format, gene_col, verbose) {

  # Auto-detect format

  if (format == "auto") {
    if (dir.exists(path)) {
      format <- "directory"
    } else if (file.exists(path)) {
      format <- "file"
    } else {
      stop("Path does not exist: ", path)
    }
  }

  if (format == "file") {
    # Single file - detect delimiter
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      counts <- read.csv(path, row.names = NULL, stringsAsFactors = FALSE)
    } else {
      counts <- read.delim(path, row.names = NULL, stringsAsFactors = FALSE)
    }
    return(counts)

  } else if (format == "directory") {
    # Directory of count files
    count_files <- list.files(path, pattern = "\\.(txt|tsv|counts)$",
                              full.names = TRUE)
    if (length(count_files) == 0) {
      stop("No count files found in directory: ", path)
    }

    if (verbose) cat("[ARTEMIS] Found", length(count_files), "count files\n")

    # Load and merge
    merged <- NULL
    for (f in count_files) {
      sample_name <- tools::file_path_sans_ext(basename(f))
      # Remove common suffixes
      sample_name <- sub("\\.counts$", "", sample_name)

      temp <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
      colnames(temp) <- c("gene", sample_name)

      if (is.null(merged)) {
        merged <- temp
      } else {
        merged <- merge(merged, temp, by = "gene", all = TRUE)
      }
    }

    rownames(merged) <- merged$gene
    merged <- merged[, -1, drop = FALSE]
    return(as.matrix(merged))
  }

  stop("Unknown counts_format: ", format)
}


#' Prepare traits data - convert to numeric where appropriate
#' @keywords internal
.wgcna_prepare_traits <- function(traits, verbose) {

  for (col in colnames(traits)) {
    x <- traits[[col]]

    if (is.numeric(x)) {
      # Already numeric, keep as is
      next
    }

    if (is.character(x) || is.factor(x)) {
      # Try to convert to numeric
      x_numeric <- suppressWarnings(as.numeric(as.character(x)))

      if (!all(is.na(x_numeric))) {
        # Successful conversion (at least some values)
        traits[[col]] <- x_numeric
        if (verbose && any(is.na(x_numeric) & !is.na(x))) {
          cat("[ARTEMIS] Note: Some non-numeric values in '", col, "' converted to NA\n")
        }
      } else {
        # Categorical - convert to numeric factor codes
        x_factor <- as.factor(x)
        traits[[col]] <- as.numeric(x_factor)
        if (verbose) {
          cat("[ARTEMIS] Converted '", col, "' to numeric codes: ",
              paste(levels(x_factor), "=", seq_along(levels(x_factor)),
                    collapse = ", "), "\n")
        }
      }
    }
  }

  return(traits)
}


#' Print method for wgcna_data objects
#' @param x A wgcna_data object
#' @param ... Additional arguments (unused)
#' @method print wgcna_data
#' @export
print.wgcna_data <- function(x, ...) {
  cat("WGCNA Data Object\n")
  cat("------------------------------\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Genes:", x$n_genes, "\n")
  cat("Traits:", ncol(x$datTraits), "-",
      paste(colnames(x$datTraits), collapse = ", "), "\n")

  if (length(x$removed_genes) > 0) {
    cat("Removed genes:", length(x$removed_genes), "\n")
  }
  if (length(x$removed_samples) > 0) {
    cat("Removed samples:", length(x$removed_samples), "\n")
  }
}


# ==============================================================================
# QUALITY CONTROL
# ==============================================================================

#' Quality control for WGCNA data
#'
#' Checks for genes and samples with too many missing values or zero variance.
#' Wraps WGCNA's goodSamplesGenes() with additional reporting.
#'
#' @param wgcna_data A wgcna_data object from ARTEMIS_wgcna_prepare(), or
#'   a matrix with samples as rows and genes as columns.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return If wgcna_data object provided, returns updated object with QC applied.
#'   If matrix provided, returns filtered matrix.
#'   In both cases, attributes record what was removed.
#'
#' @details
#' Uses WGCNA::goodSamplesGenes() which checks for:
#' - Genes with too many missing values
#' - Samples with too many missing values
#' - Genes with zero variance
#'
#' @examples
#' \dontrun{
#' wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
#' wgcna_data <- ARTEMIS_wgcna_qc(wgcna_data)
#'
#' }
#' @export
ARTEMIS_wgcna_qc <- function(wgcna_data, verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Quality Control \n")


  # Extract expression matrix
  if (inherits(wgcna_data, "wgcna_data")) {
    datExpr <- wgcna_data$datExpr
    is_wgcna_obj <- TRUE
  } else if (is.matrix(wgcna_data) || is.data.frame(wgcna_data)) {
    datExpr <- as.matrix(wgcna_data)
    is_wgcna_obj <- FALSE
  } else {
    stop("Input must be a wgcna_data object or expression matrix")
  }

  # Run WGCNA QC
  gsg <- goodSamplesGenes(datExpr, verbose = ifelse(verbose, 3, 0))

  removed_genes_qc <- character(0)
  removed_samples_qc <- character(0)

  if (!gsg$allOK) {
    # Report what will be removed
    if (sum(!gsg$goodGenes) > 0) {
      removed_genes_qc <- colnames(datExpr)[!gsg$goodGenes]
      if (verbose) {
        cat("    Removing", length(removed_genes_qc), "genes that failed QC\n")
        if (length(removed_genes_qc) <= 10) {
          cat("        ", paste(removed_genes_qc, collapse = ", "), "\n")
        } else {
          cat("        ", paste(head(removed_genes_qc, 10), collapse = ", "), "...\n")
        }
      }
    }

    if (sum(!gsg$goodSamples) > 0) {
      removed_samples_qc <- rownames(datExpr)[!gsg$goodSamples]
      if (verbose) {
        cat("    Removing", length(removed_samples_qc), "samples that failed QC:\n")
        cat("        ", paste(removed_samples_qc, collapse = ", "), "\n")
      }
    }

    # Apply filtering
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  } else {
    if (verbose) cat("[ARTEMIS] All samples and genes passed QC checks\n")
  }

  if (verbose) {
    cat("\n[ARTEMIS] Post-QC Summary ---\n")
    cat("    Samples:", nrow(datExpr), "\n")
    cat("    Genes:", ncol(datExpr), "\n")
  }

  # Return appropriate object type
  if (is_wgcna_obj) {
    wgcna_data$datExpr <- datExpr
    wgcna_data$datTraits <- wgcna_data$datTraits[rownames(datExpr), , drop = FALSE]
    wgcna_data$gene_names <- colnames(datExpr)
    wgcna_data$sample_names <- rownames(datExpr)
    wgcna_data$n_genes <- ncol(datExpr)
    wgcna_data$n_samples <- nrow(datExpr)
    wgcna_data$removed_genes <- c(wgcna_data$removed_genes, removed_genes_qc)
    wgcna_data$removed_samples <- c(wgcna_data$removed_samples, removed_samples_qc)
    wgcna_data$qc_passed <- TRUE
    return(wgcna_data)
  } else {
    attr(datExpr, "removed_genes") <- removed_genes_qc
    attr(datExpr, "removed_samples") <- removed_samples_qc
    return(datExpr)
  }
}


# ==============================================================================
# SAMPLE CLUSTERING
# ==============================================================================

#' Cluster samples and detect outliers
#'
#' Performs hierarchical clustering of samples to visualize relationships and
#' identify potential outliers. Can optionally remove outlier samples based on
#' a height cutoff.
#'
#' @param wgcna_data A wgcna_data object or expression matrix (samples as rows).
#' @param method Clustering method for hclust(). Default: "average".
#'   Options: "average", "complete", "single", "ward.D", "ward.D2", etc.
#' @param cut_height Optional height cutoff for outlier removal. Samples that
#'   cluster above this height may be removed. Default: NULL (no removal).
#' @param min_cluster_size Minimum cluster size to keep when cutting.
#'   Default: 10.
#' @param plot Logical. Generate dendrogram plot. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_cluster" containing:
#'   \item{sample_tree}{hclust object}
#'   \item{datExpr}{Expression matrix (possibly with outliers removed)}
#'   \item{datTraits}{Traits data (if wgcna_data provided)}
#'   \item{outliers_removed}{Names of removed outlier samples}
#'   \item{method}{Clustering method used}
#'   \item{cut_height}{Height cutoff used (if any)}
#'
#' @examples
#' \dontrun{
#' wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
#' clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data)
#'
#' # With outlier removal
#' clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data, cut_height = 100)
#'
#' }
#' @export
ARTEMIS_wgcna_cluster_samples <- function(wgcna_data,
                                           method = "average",
                                           cut_height = NULL,
                                           min_cluster_size = 10,
                                           plot = TRUE,
                                           verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Sample Clustering \n")

  # Extract data
  if (inherits(wgcna_data, "wgcna_data")) {
    datExpr <- wgcna_data$datExpr
    datTraits <- wgcna_data$datTraits
  } else {
    datExpr <- as.matrix(wgcna_data)
    datTraits <- NULL
  }

  # Perform hierarchical clustering
  if (verbose) cat("    Clustering samples using method:", method, "\n")
  sample_tree <- hclust(dist(datExpr), method = method)

  outliers_removed <- character(0)

  # Optional outlier removal based on height cutoff
  if (!is.null(cut_height)) {
    if (verbose) cat("    Applying cut height:", cut_height, "\n")

    clust <- cutreeStatic(sample_tree, cutHeight = cut_height,
                          minSize = min_cluster_size)

    # Cluster 0 typically contains outliers (samples not in any cluster)
    # Keep samples in the largest cluster or those below the cut
    keep_samples <- clust != 0

    if (sum(!keep_samples) > 0) {
      outliers_removed <- rownames(datExpr)[!keep_samples]
      if (verbose) {
        cat("    Removing", length(outliers_removed), "outlier samples:\n")
        cat("        ", paste(outliers_removed, collapse = ", "), "\n")
      }

      datExpr <- datExpr[keep_samples, , drop = FALSE]
      if (!is.null(datTraits)) {
        datTraits <- datTraits[keep_samples, , drop = FALSE]
      }

      # Re-cluster without outliers
      sample_tree <- hclust(dist(datExpr), method = method)
    } else {
      if (verbose) cat("    No outliers detected at this cut height\n")
    }
  }

  # Generate plot if requested
  if (plot) {
    plot(sample_tree,
         main = "Sample clustering to detect outliers",
         sub = paste("Method:", method),
         xlab = "",
         cex.lab = 1.2,
         cex.axis = 1.2,
         cex.main = 1.4)

    if (!is.null(cut_height)) {
      abline(h = cut_height, col = "red", lty = 2)
    }
  }

  if (verbose) {
    cat("\n[ARTEMIS] Clustering Summary \n")
    cat("    Samples:", nrow(datExpr), "\n")
    cat("    Method:", method, "\n")
    if (length(outliers_removed) > 0) {
      cat("    Outliers removed:", length(outliers_removed), "\n")
    }
  }

  # Build result
  result <- list(
    sample_tree = sample_tree,
    datExpr = datExpr,
    datTraits = datTraits,
    outliers_removed = outliers_removed,
    method = method,
    cut_height = cut_height,
    n_samples = nrow(datExpr)
  )
  class(result) <- c("wgcna_cluster", "list")

  return(result)
}


# ==============================================================================
# SOFT THRESHOLD POWER SELECTION
# ==============================================================================

#' Select soft threshold power for WGCNA
#'
#' Determines the optimal soft threshold power for network construction.
#' Can automatically select power or return diagnostics for manual selection.
#'
#' @param wgcna_data A wgcna_data object, wgcna_cluster object, or expression
#'   matrix (samples as rows).
#' @param powers Numeric vector of powers to test. Default: c(1:10, seq(12, 30, by = 2)).
#' @param power_auto Logical. If TRUE, automatically select power. If FALSE,
#'   return diagnostics for manual selection. Default: TRUE.
#' @param r2_cutoff R-squared threshold for scale-free topology. Power is selected
#'   where R^2 first exceeds this value. Default: 0.80.
#' @param mean_k_cutoff Optional. If R^2 criterion yields high mean connectivity,
#'   may select power with lower connectivity. Default: NULL (disabled).
#' @param network_type Network type: "unsigned", "signed", or "signed hybrid".
#'   Default: "unsigned".
#' @param plot Logical. Generate diagnostic plots. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_power" containing:
#'   \item{power}{Selected power (or NA if auto-selection disabled/failed)}
#'   \item{power_estimate}{WGCNA's built-in power estimate}
#'   \item{fit_indices}{Data.frame of fit indices for all tested powers}
#'   \item{selection_method}{How power was selected}
#'   \item{r2_at_power}{R^2 value at selected power}
#'   \item{mean_k_at_power}{Mean connectivity at selected power}
#'   \item{recommendation}{Text recommendation}
#'
#' @details
#' Power selection uses multiple criteria in order:
#' 1. WGCNA's built-in powerEstimate (if available and meets R^2 threshold)
#' 2. First power where scale-free topology R^2 exceeds r2_cutoff
#' 3. Power with best R^2 if none exceed threshold (with warning)
#'
#' The diagnostic plots show:
#' - Scale-free topology fit (R^2) vs power
#' - Mean connectivity vs power
#'
#' @examples
#' \dontrun{
#' wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
#' power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
#' # Use: power_result$power
#'
#' # Manual selection (just show plots, don't auto-select)
#' power_result <- ARTEMIS_wgcna_pick_power(wgcna_data, power_auto = FALSE)
#'
#' }
#' @export
ARTEMIS_wgcna_pick_power <- function(wgcna_data,
                                      powers = NULL,
                                      power_auto = TRUE,
                                      r2_cutoff = 0.80,
                                      mean_k_cutoff = NULL,
                                      network_type = "unsigned",
                                      plot = TRUE,
                                      verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Power Selection \n")

  # Default powers
  if (is.null(powers)) {
    powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))
  }

  # Extract expression matrix
  if (inherits(wgcna_data, "wgcna_data")) {
    datExpr <- wgcna_data$datExpr
  } else if (inherits(wgcna_data, "wgcna_cluster")) {
    datExpr <- wgcna_data$datExpr
  } else {
    datExpr <- as.matrix(wgcna_data)
  }

  # Calculate soft threshold
  if (verbose) cat("    Testing powers:", paste(powers, collapse = ", "), "\n")
  if (verbose) cat("    Network type:", network_type, "\n")

  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = network_type,
    verbose = ifelse(verbose, 2, 0)
  )

  fit_indices <- sft$fitIndices

  # Calculate signed R^2 (scale-free topology fit)
  signed_r2 <- -sign(fit_indices[, 3]) * fit_indices[, 2]
  fit_indices$signed_R2 <- signed_r2

  # ----------------------------------------------------------
  # Power selection logic
  # ----------------------------------------------------------
  selected_power <- NA
  selection_method <- "none"

  if (power_auto) {
    # Method 1: WGCNA's built-in estimate
    if (!is.na(sft$powerEstimate)) {
      idx <- which(fit_indices$Power == sft$powerEstimate)
      if (length(idx) > 0 && signed_r2[idx] >= r2_cutoff) {
        selected_power <- sft$powerEstimate
        selection_method <- "WGCNA_powerEstimate"
        if (verbose) cat("    Using WGCNA's powerEstimate:", selected_power, "\n")
      }
    }

    # Method 2: First power above R^2 threshold
    if (is.na(selected_power)) {
      above_threshold <- which(signed_r2 >= r2_cutoff)
      if (length(above_threshold) > 0) {
        selected_power <- fit_indices$Power[above_threshold[1]]
        selection_method <- "first_above_r2_threshold"
        if (verbose) {
          cat("    Selected first power with R^2 >=", r2_cutoff, ":",
              selected_power, "\n")
        }
      }
    }

    # Method 3: Best available if none meet threshold
    if (is.na(selected_power)) {
      best_idx <- which.max(signed_r2)
      selected_power <- fit_indices$Power[best_idx]
      selection_method <- "best_available"
      if (verbose) {
        warning("No power achieved R^2 >= ", r2_cutoff,
                ". Using best available: ", selected_power,
                " (R^2 = ", round(signed_r2[best_idx], 3), ")")
      }
    }

    # Optional: check mean connectivity isn't too high
    if (!is.null(mean_k_cutoff) && !is.na(selected_power)) {
      idx <- which(fit_indices$Power == selected_power)
      mean_k <- fit_indices$mean.k.[idx]
      if (mean_k > mean_k_cutoff) {
        # Find power with acceptable connectivity
        acceptable <- which(fit_indices$mean.k. <= mean_k_cutoff &
                              signed_r2 >= r2_cutoff * 0.9)  # Slightly relaxed
        if (length(acceptable) > 0) {
          selected_power <- fit_indices$Power[acceptable[1]]
          selection_method <- "connectivity_adjusted"
          if (verbose) {
            cat("    Adjusted for connectivity. New power:", selected_power, "\n")
          }
        }
      }
    }
  }

  # Get metrics at selected power
  if (!is.na(selected_power)) {
    idx <- which(fit_indices$Power == selected_power)
    r2_at_power <- signed_r2[idx]
    mean_k_at_power <- fit_indices$mean.k.[idx]
  } else {
    r2_at_power <- NA
    mean_k_at_power <- NA
  }

  # ----------------------------------------------------------
  # Generate diagnostic plots
  # ----------------------------------------------------------
  if (plot) {
    par(mfrow = c(1, 2))
    cex1 <- 0.9

    # Plot 1: Scale-free topology fit
    plot(fit_indices$Power, signed_r2,
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit (signed R^2)",
         type = "n",
         main = "Scale independence")
    text(fit_indices$Power, signed_r2,
         labels = fit_indices$Power, cex = cex1, col = "red")
    abline(h = r2_cutoff, col = "red", lty = 2)

    # Mark selected power
    if (!is.na(selected_power)) {
      idx <- which(fit_indices$Power == selected_power)
      points(selected_power, signed_r2[idx], pch = 19, cex = 2, col = "blue")
    }

    # Plot 2: Mean connectivity
    plot(fit_indices$Power, fit_indices$mean.k.,
         xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity",
         type = "n",
         main = "Mean connectivity")
    text(fit_indices$Power, fit_indices$mean.k.,
         labels = fit_indices$Power, cex = cex1, col = "red")

    if (!is.na(selected_power)) {
      idx <- which(fit_indices$Power == selected_power)
      points(selected_power, fit_indices$mean.k.[idx],
             pch = 19, cex = 2, col = "blue")
    }

    par(mfrow = c(1, 1))
  }

  # Build recommendation text
  if (!is.na(selected_power)) {
    recommendation <- sprintf(
      "Recommended power: %d (R^2 = %.3f, mean connectivity = %.1f)",
      selected_power, r2_at_power, mean_k_at_power
    )
  } else {
    recommendation <- "Could not auto-select power. Please review plots and select manually."
  }

  if (verbose) {
    cat("\n[ARTEMIS] Power Selection Summary \n")
    cat("    ", recommendation, "\n")
    cat("     Selection method:", selection_method, "\n")
  }

  # Build result
  result <- list(
    power = selected_power,
    power_estimate = sft$powerEstimate,
    fit_indices = fit_indices,
    selection_method = selection_method,
    r2_at_power = r2_at_power,
    mean_k_at_power = mean_k_at_power,
    r2_cutoff = r2_cutoff,
    network_type = network_type,
    recommendation = recommendation
  )
  class(result) <- c("wgcna_power", "list")

  return(result)
}


#' Print method for wgcna_power objects
#' @param x A wgcna_power object
#' @param ... Additional arguments (unused)
#' @method print wgcna_power
#' @export
print.wgcna_power <- function(x, ...) {
  cat("WGCNA Power Selection\n")
  cat("------------------------------\n")
  cat("Selected power:", x$power, "\n")
  cat("Selection method:", x$selection_method, "\n")
  cat("R^2 at power:", round(x$r2_at_power, 3), "\n")
  cat("Mean connectivity:", round(x$mean_k_at_power, 1), "\n")
  cat("R^2 cutoff used:", x$r2_cutoff, "\n")
  cat("Network type:", x$network_type, "\n")
}


# ==============================================================================
# MODULE DETECTION
# ==============================================================================

#' Detect co-expression modules
#'
#' Identifies modules of co-expressed genes using WGCNA's blockwiseModules.
#' This is the core network construction and module detection step.
#'
#' @param wgcna_data A wgcna_data object, wgcna_cluster object, or expression
#'   matrix (samples as rows, genes as columns).
#' @param power Soft threshold power. Can be:
#'   - Numeric value (e.g., 6)
#'   - A wgcna_power object from ARTEMIS_wgcna_pick_power()
#'   - NULL to auto-select (calls ARTEMIS_wgcna_pick_power internally)
#' @param network_type Network type: "unsigned", "signed", or "signed hybrid".
#'   Default: "unsigned".
#' @param tom_type TOM type: "unsigned", "signed", "signed Nowick", "none".
#'   Default: "unsigned".
#' @param min_module_size Minimum number of genes in a module. Default: 30.
#' @param merge_cut_height Dendrogram cut height for merging similar modules.
#'   Lower values = more merging. Default: 0.25.
#' @param reassign_threshold Threshold for reassigning genes to modules.
#'   Default: 0 (no reassignment).
#' @param max_block_size Maximum block size for network calculation. Increase
#'   for more genes but requires more memory. Default: 20000.
#' @param n_threads Number of threads for parallel computation. Default: 1.
#' @param save_tom Logical. Save TOM matrix (large, but needed for some analyses).
#'   Default: FALSE.
#' @param tom_file_base Base name for TOM file if save_tom = TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_modules" containing:
#'   \item{module_names}{Named vector: gene -> module name (module_0, module_1, ...)}
#'   \item{module_labels}{Named vector: gene -> module number (0, 1, 2, ...)}
#'   \item{module_colors}{Named vector: module_name -> hex display color}
#'   \item{color_names}{Named vector: module_name -> WGCNA color name (for reference)}
#'   \item{module_eigengenes}{Data.frame of module eigengenes (MEmodule_1, ...)}
#'   \item{module_summary}{Data.frame: module, n_genes, color_name}
#'   \item{gene_module_df}{Data.frame: gene, module_name, module_label, color_name}
#'   \item{dendrograms}{List of gene dendrograms}
#'   \item{net}{Raw blockwiseModules output (for advanced use)}
#'   \item{power}{Power used}
#'   \item{datExpr}{Expression matrix used}
#'   \item{n_modules}{Number of modules detected (excluding module_0)}
#'
#' @details
#' The function wraps WGCNA::blockwiseModules() with sensible defaults and
#' returns results in a more accessible format.
#'
#' Module colors are assigned by WGCNA. The "grey" module contains genes that
#' couldn't be assigned to any module.
#'
#' @examples
#' \dontrun{
#' wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
#' power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = power_result)
#'
#' # Or with auto power selection
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = NULL)
#'
#' }
#' @export
ARTEMIS_wgcna_detect_modules <- function(wgcna_data,
                                          power,
                                          network_type = "unsigned",
                                          tom_type = "unsigned",
                                          min_module_size = 30,
                                          merge_cut_height = 0.25,
                                          reassign_threshold = 0,
                                          max_block_size = 20000,
                                          n_threads = 1,
                                          save_tom = FALSE,
                                          tom_file_base = "TOM",
                                          verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Module Detection \n")

  # --------------------------------------------------------------------------
  # Check for WGCNA cor() masking issue
  # --------------------------------------------------------------------------
  # WGCNA has its own cor() with extra arguments (weights.x, weights.y, cosine).
  # blockwiseModules passes these args internally; if stats::cor is found instead
  # of WGCNA::cor the call fails. We check the runtime search path (not GAIA's
  # namespace) and assign WGCNA::cor to .GlobalEnv when needed.
  wgcna_cor_ok <- tryCatch({
    cor_on_path <- get("cor", envir = .GlobalEnv, inherits = TRUE)
    identical(environment(cor_on_path), asNamespace("WGCNA"))
  }, error = function(e) FALSE)

  if (!wgcna_cor_ok) {
    if (verbose) {
      cat("    Note: WGCNA::cor() not found on search path. Assigning WGCNA::cor to global environment for blockwiseModules.\n")
    }
    assign("cor", WGCNA::cor, envir = .GlobalEnv)
  }

  # --------------------------------------------------------------------------
  # Extract expression matrix
  # --------------------------------------------------------------------------
  if (inherits(wgcna_data, "wgcna_data")) {
    datExpr <- wgcna_data$datExpr
  } else if (inherits(wgcna_data, "wgcna_cluster")) {
    datExpr <- wgcna_data$datExpr
  } else {
    datExpr <- as.matrix(wgcna_data)
  }

  # Ensure numeric (double, not integer - WGCNA requirement)
  # WGCNA's blockwiseModules requires double precision, not integer
  if (!is.numeric(datExpr)) {
    datExpr <- apply(datExpr, 2, as.numeric)
  }
  storage.mode(datExpr) <- "double"

  # --------------------------------------------------------------------------
  # Handle power parameter
  # --------------------------------------------------------------------------
  if (is.null(power)) {
    if (verbose) cat("    Power not specified, auto-selecting...\n")
    power_result <- ARTEMIS_wgcna_pick_power(datExpr,
                                              network_type = network_type,
                                              verbose = verbose)
    power_value <- power_result$power
  } else if (inherits(power, "wgcna_power")) {
    power_value <- power$power
  } else {
    power_value <- as.numeric(power)
  }

  if (is.na(power_value)) {
    stop("Could not determine power value. Please specify manually.")
  }

  if (verbose) cat("    Using power:", power_value, "\n")

  # --------------------------------------------------------------------------
  # Run blockwiseModules
  # --------------------------------------------------------------------------
  if (verbose) {
    cat("[ARTEMIS] Detecting modules...\n")
    cat("    Network type:", network_type, "\n")
    cat("    TOM type:", tom_type, "\n")
    cat("    Min module size:", min_module_size, "\n")
    cat("    Merge cut height:", merge_cut_height, "\n")
    cat("    Max block size:", max_block_size, "\n")
    if (n_threads > 1) cat("        Threads:", n_threads, "\n")
  }

  net <- blockwiseModules(
    datExpr,
    power = power_value,
    networkType = network_type,
    TOMType = tom_type,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    reassignThreshold = reassign_threshold,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = save_tom,
    saveTOMFileBase = tom_file_base,
    maxBlockSize = max_block_size,
    nThreads = n_threads,
    verbose = ifelse(verbose, 3, 0)
  )

  # --------------------------------------------------------------------------
  # Process results
  # --------------------------------------------------------------------------
  module_labels <- net$colors
  names(module_labels) <- colnames(datExpr)

  # WGCNA color names (needed internally for moduleEigengenes, kept for reference)
  wgcna_color_names <- labels2colors(net$colors)
  names(wgcna_color_names) <- colnames(datExpr)

  # Generic module names: module_0 (unassigned), module_1 (largest), etc.
  module_names <- paste0("module_", module_labels)
  names(module_names) <- colnames(datExpr)

  # Build label -> color_name mapping for eigengene renaming
  unique_labels <- sort(unique(module_labels))
  label_to_module <- setNames(paste0("module_", unique_labels), unique_labels)
  label_to_color <- setNames(
    labels2colors(unique_labels),
    unique_labels
  )

  # Display color map: module_name -> hex color
  unique_module_names <- paste0("module_", unique_labels)
  n_colors <- length(unique_module_names)
  display_colors <- grDevices::hcl.colors(max(n_colors * 3, 50), palette = "Dark 3")
  # Filter out near-white colors
  rgb_vals <- grDevices::col2rgb(display_colors)
  brightness <- colMeans(rgb_vals)
  display_colors <- display_colors[brightness <= 200]
  display_colors <- display_colors[seq_len(n_colors)]
  # module_0 (unassigned) gets grey
  if ("module_0" %in% unique_module_names) {
    display_colors[which(unique_module_names == "module_0")] <- "#AAAAAA"
  }
  module_color_map <- setNames(display_colors, unique_module_names)

  # color_names: module_name -> WGCNA color name (for reference)
  color_name_map <- setNames(label_to_color[as.character(unique_labels)], unique_module_names)

  # Get module eigengenes (requires WGCNA color names internally)
  MEs0 <- moduleEigengenes(datExpr, wgcna_color_names)$eigengenes
  MEs <- orderMEs(MEs0)

  # Rename eigengene columns: MEturquoise -> MEmodule_1
  color_to_module <- setNames(
    paste0("module_", unique_labels),
    labels2colors(unique_labels)
  )
  colnames(MEs) <- paste0("ME", color_to_module[gsub("^ME", "", colnames(MEs))])

  # Create gene-module data.frame
  gene_module_df <- data.frame(
    gene = colnames(datExpr),
    module_name = module_names,
    module_label = module_labels,
    color_name = wgcna_color_names,
    stringsAsFactors = FALSE
  )

  # Create module summary
  module_counts <- table(module_names)
  module_summary <- data.frame(
    module = names(module_counts),
    n_genes = as.integer(module_counts),
    color_name = color_name_map[names(module_counts)],
    stringsAsFactors = FALSE
  )
  module_summary <- module_summary[order(-module_summary$n_genes), ]
  rownames(module_summary) <- NULL

  # Count modules (excluding module_0/unassigned)
  n_modules <- length(unique(module_names[module_names != "module_0"]))

  if (verbose) {
    cat("\n[ARTEMIS] Module Detection Summary \n")
    cat("    Total genes:", ncol(datExpr), "\n")
    cat("    Modules detected:", n_modules, "(plus module_0/unassigned)\n")
    cat("    \nModule sizes:\n")
    print(module_summary, row.names = FALSE)
  }

  # --------------------------------------------------------------------------
  # Build result object
  # --------------------------------------------------------------------------
  result <- list(
    module_names = module_names,
    module_labels = module_labels,
    module_colors = module_color_map,
    color_names = color_name_map,
    module_eigengenes = MEs,
    module_summary = module_summary,
    gene_module_df = gene_module_df,
    dendrograms = net$dendrograms,
    net = net,
    power = power_value,
    network_type = network_type,
    datExpr = datExpr,
    n_modules = n_modules,
    n_genes = ncol(datExpr),
    n_samples = nrow(datExpr)
  )
  class(result) <- c("wgcna_modules", "list")

  return(result)
}


#' Print method for wgcna_modules objects
#' @param x A wgcna_modules object
#' @param ... Additional arguments (unused)
#' @method print wgcna_modules
#' @export
print.wgcna_modules <- function(x, ...) {
  cat("WGCNA Module Detection Results\n")
  cat("------------------------------\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Genes:", x$n_genes, "\n")
  cat("Modules:", x$n_modules, "(excluding grey)\n")
  cat("Power:", x$power, "\n")
  cat("Network type:", x$network_type, "\n")
  cat("\nTop modules by size:\n")
  print(head(x$module_summary, 10), row.names = FALSE)
}


# ==============================================================================
# MODULE-TRAIT RELATIONSHIPS
# ==============================================================================

#' Correlate modules with traits
#'
#' Calculates correlations between module eigengenes and sample traits.
#' Identifies which modules are associated with which clinical/phenotypic variables.
#'
#' @param modules A wgcna_modules object from ARTEMIS_wgcna_detect_modules().
#' @param traits Data.frame of traits, or wgcna_data object (uses datTraits).
#'   Must have samples as rows matching the module data.
#' @param cor_method Correlation method: "pearson" or "spearman". Default: "pearson".
#' @param p_adjust Method for p-value adjustment: "BH", "bonferroni", "none", etc.
#'   Default: "BH" (Benjamini-Hochberg).
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_trait_cor" containing:
#'   \item{cor_matrix}{Matrix of correlations (modules x traits)}
#'   \item{pvalue_matrix}{Matrix of p-values}
#'   \item{padj_matrix}{Matrix of adjusted p-values}
#'   \item{significant_associations}{Data.frame of significant module-trait pairs}
#'   \item{n_significant}{Number of significant associations}
#'
#' @examples
#' \dontrun{
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
#' trait_cor <- ARTEMIS_wgcna_module_traits(modules, wgcna_data)
#'
#' }
#' @export
ARTEMIS_wgcna_module_traits <- function(modules,
                                         traits,
                                         cor_method = "pearson",
                                         p_adjust = "BH",
                                         verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Module-Trait Correlation \n")

  # --------------------------------------------------------------------------
  # Extract data
  # --------------------------------------------------------------------------
  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules object from ARTEMIS_wgcna_detect_modules()")
  }

  MEs <- modules$module_eigengenes

  # Extract traits
  if (inherits(traits, "wgcna_data")) {
    datTraits <- traits$datTraits
  } else if (inherits(traits, "wgcna_cluster")) {
    datTraits <- traits$datTraits
  } else {
    datTraits <- as.data.frame(traits)
  }

  # Ensure numeric
  for (col in colnames(datTraits)) {
    if (!is.numeric(datTraits[[col]])) {
      datTraits[[col]] <- as.numeric(as.factor(datTraits[[col]]))
    }
  }

  # Align samples
  common_samples <- intersect(rownames(MEs), rownames(datTraits))
  if (length(common_samples) == 0) {
    stop("No matching samples between modules and traits")
  }

  if (length(common_samples) < nrow(MEs)) {
    if (verbose) {
      cat("    Using", length(common_samples), "of", nrow(MEs), "samples\n")
    }
  }

  MEs <- MEs[common_samples, , drop = FALSE]
  datTraits <- datTraits[common_samples, , drop = FALSE]

  n_samples <- nrow(MEs)

  # --------------------------------------------------------------------------
  # Calculate correlations
  # --------------------------------------------------------------------------
  if (verbose) cat("    Calculating correlations using method:", cor_method, "\n")

  if (cor_method == "pearson") {
    cor_matrix <- cor(MEs, datTraits, use = "pairwise.complete.obs")
    # Calculate p-values using WGCNA function
    pvalue_matrix <- corPvalueStudent(cor_matrix, n_samples)
  } else {
    # Spearman - calculate manually
    cor_matrix <- matrix(NA, ncol(MEs), ncol(datTraits))
    pvalue_matrix <- matrix(NA, ncol(MEs), ncol(datTraits))
    rownames(cor_matrix) <- rownames(pvalue_matrix) <- colnames(MEs)
    colnames(cor_matrix) <- colnames(pvalue_matrix) <- colnames(datTraits)

    for (i in seq_len(ncol(MEs))) {
      for (j in seq_len(ncol(datTraits))) {
        test <- cor.test(MEs[, i], datTraits[, j],
                         method = "spearman", exact = FALSE)
        cor_matrix[i, j] <- test$estimate
        pvalue_matrix[i, j] <- test$p.value
      }
    }
  }

  # --------------------------------------------------------------------------
  # Adjust p-values
  # --------------------------------------------------------------------------
  if (p_adjust != "none") {
    padj_vector <- p.adjust(as.vector(pvalue_matrix), method = p_adjust)
    padj_matrix <- matrix(padj_vector, nrow = nrow(pvalue_matrix))
    rownames(padj_matrix) <- rownames(pvalue_matrix)
    colnames(padj_matrix) <- colnames(pvalue_matrix)
  } else {
    padj_matrix <- pvalue_matrix
  }

  # --------------------------------------------------------------------------
  # Identify significant associations
  # --------------------------------------------------------------------------
  sig_threshold <- 0.05

  sig_list <- list()
  for (i in seq_len(nrow(cor_matrix))) {
    for (j in seq_len(ncol(cor_matrix))) {
      if (!is.na(padj_matrix[i, j]) && padj_matrix[i, j] < sig_threshold) {
        sig_list[[length(sig_list) + 1]] <- data.frame(
          module = rownames(cor_matrix)[i],
          trait = colnames(cor_matrix)[j],
          correlation = cor_matrix[i, j],
          pvalue = pvalue_matrix[i, j],
          padj = padj_matrix[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(sig_list) > 0) {
    significant_associations <- do.call(rbind, sig_list)
    significant_associations <- significant_associations[
      order(significant_associations$padj), ]
    rownames(significant_associations) <- NULL
  } else {
    significant_associations <- data.frame(
      module = character(0),
      trait = character(0),
      correlation = numeric(0),
      pvalue = numeric(0),
      padj = numeric(0)
    )
  }

  n_significant <- nrow(significant_associations)

  if (verbose) {
    cat("[ARTEMIS] Module-Trait Correlation Summary \n")
    cat("    Modules:", nrow(cor_matrix), "\n")
    cat("    Traits:", ncol(cor_matrix), "\n")
    cat("    Significant associations (padj < 0.05):", n_significant, "\n")

    if (n_significant > 0) {
      cat("    Top significant associations:\n")
      print(head(significant_associations, 10), row.names = FALSE)
    }
  }

  # --------------------------------------------------------------------------
  # Build result
  # --------------------------------------------------------------------------
  result <- list(
    cor_matrix = cor_matrix,
    pvalue_matrix = pvalue_matrix,
    padj_matrix = padj_matrix,
    significant_associations = significant_associations,
    n_significant = n_significant,
    cor_method = cor_method,
    p_adjust = p_adjust,
    n_samples = n_samples,
    module_eigengenes = MEs,
    traits = datTraits
  )
  class(result) <- c("wgcna_trait_cor", "list")

  return(result)
}


#' Print method for wgcna_trait_cor objects
#' @param x A wgcna_trait_cor object
#' @param ... Additional arguments (unused)
#' @method print wgcna_trait_cor
#' @export
print.wgcna_trait_cor <- function(x, ...) {
  cat("WGCNA Module-Trait Correlations\n")
  cat("-------------------------------\n")
  cat("Modules:", nrow(x$cor_matrix), "\n")
  cat("Traits:", ncol(x$cor_matrix), "\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Correlation method:", x$cor_method, "\n")
  cat("P-value adjustment:", x$p_adjust, "\n")
  cat("Significant associations:", x$n_significant, "\n")

  if (x$n_significant > 0) {
    cat("\nTop associations:\n")
    print(head(x$significant_associations, 5), row.names = FALSE)
  }
}


# ==============================================================================
# GENE SIGNIFICANCE
# ==============================================================================

#' Calculate gene significance for traits
#'
#' Computes gene significance (GS) and module membership (MM) for each gene.
#' Gene significance measures the correlation of each gene with a trait.
#' Module membership measures how well each gene fits its assigned module.
#'
#' @param modules A wgcna_modules object from ARTEMIS_wgcna_detect_modules().
#' @param traits Data.frame of traits, or wgcna_data object.
#' @param trait_name Name of specific trait to analyze. If NULL, analyzes all traits.
#' @param cor_method Correlation method: "pearson" or "spearman". Default: "pearson".
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_gene_sig" containing:
#'   \item{gene_info}{Data.frame with gene significance and module membership}
#'   \item{by_module}{List of data.frames, one per module}
#'   \item{top_genes_per_module}{Top significant genes for each module}
#'   \item{trait_names}{Names of traits analyzed}
#'
#' @details
#' For each gene, calculates:
#' - Gene Significance (GS): correlation between gene expression and trait
#' - Module Membership (MM): correlation between gene expression and module eigengene
#'
#' High GS + High MM = genes that are both important for the trait AND
#' central to their module (good biomarker candidates).
#'
#' @examples
#' \dontrun{
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
#' gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)
#'
#' # For specific trait
#' gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data,
#'                                              trait_name = "age")
#'
#' }
#' @export
ARTEMIS_wgcna_gene_significance <- function(modules,
                                             traits,
                                             trait_name = NULL,
                                             cor_method = "pearson",
                                             verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Gene Significance \n")

  # --------------------------------------------------------------------------
  # Extract data
  # --------------------------------------------------------------------------
  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules object")
  }

  datExpr <- modules$datExpr
  MEs <- modules$module_eigengenes
  module_names_vec <- modules$module_names

  # Extract traits
  if (inherits(traits, "wgcna_data")) {
    datTraits <- traits$datTraits
  } else if (inherits(traits, "wgcna_cluster")) {
    datTraits <- traits$datTraits
  } else {
    datTraits <- as.data.frame(traits)
  }

  # Subset to specific trait if requested
  if (!is.null(trait_name)) {
    if (!trait_name %in% colnames(datTraits)) {
      stop("Trait '", trait_name, "' not found. Available: ",
           paste(colnames(datTraits), collapse = ", "))
    }
    datTraits <- datTraits[, trait_name, drop = FALSE]
  }

  # Align samples
  common_samples <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr <- datExpr[common_samples, , drop = FALSE]
  datTraits <- datTraits[common_samples, , drop = FALSE]
  MEs <- MEs[common_samples, , drop = FALSE]

  n_samples <- nrow(datExpr)

  if (verbose) {
    cat("    Samples:", n_samples, "\n")
    cat("    Genes:", ncol(datExpr), "\n")
    cat("    Traits:", ncol(datTraits), "-",
        paste(colnames(datTraits), collapse = ", "), "\n")
  }

  # --------------------------------------------------------------------------
  # Calculate Module Membership (MM)
  # --------------------------------------------------------------------------
  if (verbose) cat("    Calculating module membership...\n")

  # MM = correlation of each gene with each module eigengene
  MM <- cor(datExpr, MEs, use = "pairwise.complete.obs",
            method = cor_method)
  colnames(MM) <- paste0("MM.", gsub("^ME", "", colnames(MEs)))

  # P-values for MM
  MM_pvalue <- corPvalueStudent(MM, n_samples)
  colnames(MM_pvalue) <- paste0("p.MM.", gsub("^ME", "", colnames(MEs)))

  # --------------------------------------------------------------------------
  # Calculate Gene Significance (GS)
  # --------------------------------------------------------------------------
  if (verbose) cat("    Calculating gene significance...\n")

  GS <- cor(datExpr, datTraits, use = "pairwise.complete.obs",
            method = cor_method)
  colnames(GS) <- paste0("GS.", colnames(datTraits))

  GS_pvalue <- corPvalueStudent(GS, n_samples)
  colnames(GS_pvalue) <- paste0("p.GS.", colnames(datTraits))

  # --------------------------------------------------------------------------
  # Build gene info table
  # --------------------------------------------------------------------------
  gene_info <- data.frame(
    gene = colnames(datExpr),
    module = module_names_vec,
    stringsAsFactors = FALSE
  )

  # Add GS columns
  gene_info <- cbind(gene_info, GS, GS_pvalue)

  # Add MM for each gene's own module
  gene_info$MM_own_module <- NA
  gene_info$p.MM_own_module <- NA
  for (i in seq_len(nrow(gene_info))) {
    mod_col <- paste0("MM.", gene_info$module[i])
    p_col <- paste0("p.MM.", gene_info$module[i])
    if (mod_col %in% colnames(MM)) {
      gene_info$MM_own_module[i] <- MM[i, mod_col]
      gene_info$p.MM_own_module[i] <- MM_pvalue[i, p_col]
    }
  }

  # Add all MM columns
  gene_info <- cbind(gene_info, MM, MM_pvalue)

  # --------------------------------------------------------------------------
  # Organize by module
  # --------------------------------------------------------------------------
  modules_list <- split(gene_info, gene_info$module)

  # Sort each module by GS (first trait) descending
  first_gs_col <- paste0("GS.", colnames(datTraits)[1])
  for (mod in names(modules_list)) {
    modules_list[[mod]] <- modules_list[[mod]][
      order(-abs(modules_list[[mod]][[first_gs_col]])), ]
  }

  # Get top genes per module
  top_n <- 10
  top_genes <- lapply(modules_list, function(df) {
    head(df[, c("gene", "module", first_gs_col,
                paste0("p.", first_gs_col), "MM_own_module")], top_n)
  })

  if (verbose) {
    cat("[ARTEMIS] Gene Significance Summary \n")
    cat("    Total genes analyzed:", nrow(gene_info), "\n")

    # Count significant genes per module (for first trait)
    p_col <- paste0("p.GS.", colnames(datTraits)[1])
    sig_genes <- gene_info[gene_info[[p_col]] < 0.05, ]
    cat("    Genes with significant GS (p < 0.05) for",
        colnames(datTraits)[1], ":", nrow(sig_genes), "\n")

    if (nrow(sig_genes) > 0) {
      sig_by_mod <- table(sig_genes$module)
      cat("    Significant genes by module:\n")
      print(sort(sig_by_mod, decreasing = TRUE))
    }
  }

  # --------------------------------------------------------------------------
  # Build result
  # --------------------------------------------------------------------------
  result <- list(
    gene_info = gene_info,
    by_module = modules_list,
    top_genes_per_module = top_genes,
    trait_names = colnames(datTraits),
    module_membership = MM,
    gene_significance = GS,
    n_samples = n_samples,
    cor_method = cor_method
  )
  class(result) <- c("wgcna_gene_sig", "list")

  return(result)
}


#' Print method for wgcna_gene_sig objects
#' @param x A wgcna_gene_sig object
#' @param ... Additional arguments (unused)
#' @method print wgcna_gene_sig
#' @export
print.wgcna_gene_sig <- function(x, ...) {
  cat("WGCNA Gene Significance Analysis\n")
  cat("---------------------------------\n")
  cat("Genes:", nrow(x$gene_info), "\n")
  cat("Traits:", paste(x$trait_names, collapse = ", "), "\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Correlation method:", x$cor_method, "\n")
  cat("Modules:", length(x$by_module), "\n")
}


# ==============================================================================
# HUB GENES
# ==============================================================================

#' Identify hub genes in modules
#'
#' Identifies hub genes - genes that are highly connected within their module
#' and/or highly correlated with traits of interest.
#'
#' @param modules A wgcna_modules object from ARTEMIS_wgcna_detect_modules().
#' @param gene_sig Optional. A wgcna_gene_sig object. If provided, uses GS
#'   to identify trait-relevant hubs.
#' @param trait_name Trait name to use for trait-relevant hub identification.
#'   Required if gene_sig is provided and has multiple traits.
#' @param n_top Number of top hub genes to return per module. Default: 10.
#'   Set to NULL to return all genes above mm_threshold (useful for enrichment
#'   analysis where a kME threshold is more appropriate than a fixed count).
#' @param mm_threshold Module membership threshold for hub genes. Default: 0.8.
#' @param gs_threshold Gene significance threshold (if using trait). Default: 0.2.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with class "wgcna_hubs" containing:
#'   \item{hub_genes}{Data.frame of all hub genes across modules}
#'   \item{by_module}{List of hub genes per module}
#'   \item{hub_summary}{Summary of hub genes per module}
#'   \item{criteria}{Criteria used for hub identification}
#'
#' @details
#' Hub genes are identified based on:
#' - Module Membership (MM): How well the gene's expression pattern matches
#'   the module eigengene. High MM = gene is central to module.
#' - Gene Significance (GS, optional): Correlation with trait of interest.
#'   High GS + High MM = biologically relevant hub.
#'
#' @examples
#' \dontrun{
#' modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
#' hubs <- ARTEMIS_wgcna_hub_genes(modules)
#'
#' # With trait information
#' gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)
#' hubs <- ARTEMIS_wgcna_hub_genes(modules, gene_sig, trait_name = "age")
#'
#' }
#' @export
ARTEMIS_wgcna_hub_genes <- function(modules,
                                     gene_sig = NULL,
                                     trait_name = NULL,
                                     n_top = 10,
                                     mm_threshold = 0.8,
                                     gs_threshold = 0.2,
                                     verbose = TRUE) {

  if (verbose) cat("[ARTEMIS] WGCNA Hub Gene Identification \n")

  # --------------------------------------------------------------------------
  # Extract data
  # --------------------------------------------------------------------------
  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules object")
  }

  datExpr <- modules$datExpr
  MEs <- modules$module_eigengenes
  module_names_vec <- modules$module_names

  # Calculate MM if gene_sig not provided
  if (is.null(gene_sig)) {
    if (verbose) cat("    Calculating module membership...\n")
    MM <- cor(datExpr, MEs, use = "pairwise.complete.obs")
    colnames(MM) <- gsub("^ME", "", colnames(MM))
    GS <- NULL
    trait_name <- NULL
  } else {
    if (!inherits(gene_sig, "wgcna_gene_sig")) {
      stop("gene_sig must be a wgcna_gene_sig object")
    }

    # Extract MM (remove "MM." prefix for matching)
    mm_cols <- grep("^MM\\.", colnames(gene_sig$gene_info), value = TRUE)
    MM <- as.matrix(gene_sig$gene_info[, mm_cols])
    colnames(MM) <- gsub("^MM\\.", "", colnames(MM))
    rownames(MM) <- gene_sig$gene_info$gene

    # Extract GS
    if (!is.null(trait_name)) {
      gs_col <- paste0("GS.", trait_name)
      if (!gs_col %in% colnames(gene_sig$gene_info)) {
        stop("Trait '", trait_name, "' not found in gene_sig")
      }
      GS <- gene_sig$gene_info[[gs_col]]
      names(GS) <- gene_sig$gene_info$gene
    } else if (length(gene_sig$trait_names) == 1) {
      trait_name <- gene_sig$trait_names[1]
      gs_col <- paste0("GS.", trait_name)
      GS <- gene_sig$gene_info[[gs_col]]
      names(GS) <- gene_sig$gene_info$gene
    } else {
      GS <- NULL
    }
  }

  if (verbose) {
    cat("    MM threshold:", mm_threshold, "\n")
    if (!is.null(GS)) {
      cat("    GS threshold:", gs_threshold, "(trait:", trait_name, ")\n")
    }
  }

  # --------------------------------------------------------------------------
  # Identify hub genes per module
  # --------------------------------------------------------------------------
  unique_modules <- unique(module_names_vec)
  unique_modules <- unique_modules[unique_modules != "module_0"]  # Exclude unassigned

  hub_list <- list()

  for (mod in unique_modules) {
    # Get genes in this module
    mod_genes <- names(module_names_vec)[module_names_vec == mod]

    # Get MM for these genes (for their own module)
    if (mod %in% colnames(MM)) {
      mod_mm <- MM[mod_genes, mod]
    } else {
      next
    }

    # Build data.frame
    mod_df <- data.frame(
      gene = mod_genes,
      module = mod,
      MM = mod_mm,
      stringsAsFactors = FALSE
    )

    # Add GS if available
    if (!is.null(GS)) {
      mod_df$GS <- GS[mod_genes]
      mod_df$abs_GS <- abs(mod_df$GS)
    }

    # Filter by thresholds
    mod_df$is_hub_mm <- abs(mod_df$MM) >= mm_threshold

    if (!is.null(GS)) {
      mod_df$is_hub_gs <- abs(mod_df$GS) >= gs_threshold
      mod_df$is_hub <- mod_df$is_hub_mm & mod_df$is_hub_gs

      # Sort by combined criteria
      mod_df <- mod_df[order(-abs(mod_df$MM) * abs(mod_df$GS)), ]
    } else {
      mod_df$is_hub <- mod_df$is_hub_mm
      mod_df <- mod_df[order(-abs(mod_df$MM)), ]
    }

    # Get top hubs (NULL = all genes above mm_threshold)
    hub_list[[mod]] <- if (is.null(n_top)) mod_df else head(mod_df, n_top)
  }

  # Combine all hub genes
  all_hubs <- do.call(rbind, hub_list)
  rownames(all_hubs) <- NULL

  # Filter to actual hubs
  actual_hubs <- all_hubs[all_hubs$is_hub, ]

  # Summary
  hub_summary <- data.frame(
    module = names(hub_list),
    n_genes = sapply(unique_modules, function(m) sum(module_names_vec == m)),
    n_hubs = sapply(hub_list, function(df) sum(df$is_hub)),
    top_hub = sapply(hub_list, function(df) df$gene[1]),
    stringsAsFactors = FALSE
  )
  hub_summary <- hub_summary[order(-hub_summary$n_hubs), ]
  rownames(hub_summary) <- NULL

  if (verbose) {
    cat("[ARTEMIS] Hub Gene Summary \n")
    cat("    Total hub genes identified:", nrow(actual_hubs), "\n")
    cat("    Hubs per module:\n")
    print(hub_summary, row.names = FALSE)
  }

  # --------------------------------------------------------------------------
  # Build result
  # --------------------------------------------------------------------------
  criteria <- list(
    mm_threshold = mm_threshold,
    gs_threshold = if (!is.null(GS)) gs_threshold else NA,
    trait_name = trait_name,
    n_top = n_top
  )

  result <- list(
    hub_genes = actual_hubs,
    top_per_module = all_hubs,
    by_module = hub_list,
    hub_summary = hub_summary,
    criteria = criteria
  )
  class(result) <- c("wgcna_hubs", "list")

  return(result)
}


#' Print method for wgcna_hubs objects
#' @param x A wgcna_hubs object
#' @param ... Additional arguments (unused)
#' @method print wgcna_hubs
#' @export
print.wgcna_hubs <- function(x, ...) {
  cat("WGCNA Hub Genes\n")
  cat("---------------\n")
  cat("Total hub genes:", nrow(x$hub_genes), "\n")
  cat("MM threshold:", x$criteria$mm_threshold, "\n")
  if (!is.na(x$criteria$gs_threshold)) {
    cat("GS threshold:", x$criteria$gs_threshold,
        "(trait:", x$criteria$trait_name, ")\n")
  }
  cat("\nTop modules by hub count:\n")
  print(head(x$hub_summary, 10), row.names = FALSE)
}
