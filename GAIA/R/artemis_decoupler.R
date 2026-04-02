# GAIA/Artemis/decoupler.R
# Activity inference using decoupleR
#
# Functions for inferring transcription factor and pathway activities
# from gene expression data using various statistical methods.
#
# Core workflow:
# 1. Get prior knowledge network (Apollo functions)
# 2. Run activity inference (this file)
# 3. Visualize results (Aether functions)


# ==============================================================================
# AVAILABLE METHODS
# ==============================================================================

# Define available methods and their descriptions
.DECOUPLER_METHODS <- list(
  ulm = list(
    name = "Univariate Linear Model",
    description = "Fits a linear model for each source. Fast and robust.",
    recommended = TRUE
  ),
  mlm = list(
    name = "Multivariate Linear Model",
    description = "Fits a single linear model with all sources. Accounts for co-regulation.",
    recommended = TRUE
  ),
  viper = list(
    name = "VIPER",
    description = "Virtual Inference of Protein-activity by Enriched Regulon analysis.",
    recommended = FALSE
  ),
  wsum = list(
    name = "Weighted Sum",
    description = "Weighted sum of target expression. Simple but effective.",
    recommended = FALSE
  ),
  wmean = list(
    name = "Weighted Mean",
    description = "Weighted mean of target expression.",
    recommended = FALSE
  ),
  udt = list(
    name = "Univariate Decision Tree",
    description = "Decision tree for each source independently.",
    recommended = FALSE
  ),
  mdt = list(
    name = "Multivariate Decision Tree",
    description = "Single decision tree with all sources.",
    recommended = FALSE
  ),
  ora = list(
    name = "Over-Representation Analysis",
    description = "Fisher's exact test for enrichment. Requires binary input.",
    recommended = FALSE
  ),
  gsva = list(
    name = "GSVA",
    description = "Gene Set Variation Analysis. Sample-wise enrichment scores.",
    recommended = FALSE
  ),
  aucell = list(
    name = "AUCell",
    description = "Area Under Curve for cell-level enrichment. Good for single-cell.",
    recommended = FALSE
  )
)

# ==============================================================================
# INTERNAL NORMALIZATION HELPERS
# ==============================================================================

.is_raw_counts <- function(mat) {
  sample_vals <- mat[seq_len(min(100, nrow(mat))), seq_len(min(10, ncol(mat))), drop = FALSE]
  all(sample_vals == floor(sample_vals), na.rm = TRUE)
}

.normalize_for_decoupler <- function(mat, norm_method, verbose = TRUE) {
  if (norm_method == "none") return(mat)
  if (norm_method == "log2cpm") {
    if (verbose) cat("[ARTEMIS] Normalizing: log2(CPM + 1)\n")
    lib_sizes <- colSums(mat)
    cpm <- sweep(mat, 2, lib_sizes, "/") * 1e6
    return(log2(cpm + 1))
  }
  stop("Unknown norm_method: '", norm_method, "'. Options: 'log2cpm', 'none'.")
}


#' List available decoupleR methods
#'
#' Prints information about available statistical methods for activity inference.
#'
#' @param recommended_only Logical. If TRUE, only show recommended methods.
#'   Default: FALSE.
#'
#' @return Character vector of method names (invisibly).
#'
#' @examples
#' \dontrun{
#' ARTEMIS_decoupler_list_methods()
#' ARTEMIS_decoupler_list_methods(recommended_only = TRUE)
#'
#' }
#' @export
ARTEMIS_decoupler_list_methods <- function(recommended_only = FALSE) {

  cat("[ARTEMIS] Available decoupleR Methods \n\n")

  methods_to_show <- .DECOUPLER_METHODS
  if (recommended_only) {
    methods_to_show <- Filter(function(x) x$recommended, methods_to_show)
  }

  for (method_name in names(methods_to_show)) {
    info <- methods_to_show[[method_name]]
    rec_tag <- if (info$recommended) " [RECOMMENDED]" else ""
    cat(sprintf("    %s%s\n", method_name, rec_tag))
    cat(sprintf("        %s\n", info$name))
    cat(sprintf("        %s\n\n", info$description))
  }

  if (!recommended_only) {
    cat("[ARTEMIS] Note: 'ulm' and 'mlm' are recommended as default methods.\n")
    cat("[ARTEMIS] Use recommended_only = TRUE to see only recommended methods.\n")
  }

  invisible(names(methods_to_show))
}


# ==============================================================================
# CORE ACTIVITY INFERENCE
# ==============================================================================

#' Run decoupleR activity inference
#'
#' Core function for inferring biological activities (TF or pathway) from
#' gene expression data using a specified statistical method.
#'
#' @param mat Numeric matrix of gene expression data (genes x samples), an
#'   \code{artemis_norm} object from \code{ARTEMIS_normalize_counts()}, or raw
#'   integer count matrix. Rownames must be gene identifiers matching the network.
#' @param network Prior knowledge network (from APOLLO_get_* functions or custom).
#'   Must have columns: 'source' (TF/pathway), 'target' (gene), and 'mor' or 'weight'.
#' @param method Character. Statistical method to use. REQUIRED - no default.
#'   Options: "ulm", "mlm", "viper", "wsum", "wmean", "udt", "mdt", "ora", "gsva", "aucell".
#' @param minsize Integer. Minimum number of targets per source to include.
#'   Default: 5.
#' @param norm_method Character. Normalization to apply to raw counts before inference.
#'   \code{"log2cpm"} (default): log2(CPM + 1). \code{"none"}: pass matrix through as-is
#'   (use when already normalized). If \code{mat} is an \code{artemis_norm} object,
#'   DESeq2-normalized counts are extracted and log2-transformed automatically,
#'   regardless of this parameter.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#' @param ... Additional arguments passed to the specific method function.
#'
#' @return A list with class "decoupler_result" containing:
#'   \item{activities}{Matrix of activity scores (sources x samples)}
#'   \item{pvalues}{Matrix of p-values (if available for method)}
#'   \item{results_long}{Full results in long format (from decoupleR)}
#'   \item{method}{Method used}
#'   \item{network_info}{Information about the network used}
#'   \item{n_sources}{Number of sources (TFs/pathways) with results}
#'   \item{n_samples}{Number of samples}
#'
#' @details
#' This function requires an explicit method choice to ensure users understand
#' which statistical approach they are using. For method recommendations, use
#' \code{ARTEMIS_decoupler_list_methods()}.
#'
#' DecoupleR methods expect log-normalized continuous expression values, not raw
#' counts. Use \code{norm_method = "log2cpm"} for raw counts, or pass an
#' \code{artemis_norm} object to reuse DESeq2 size factors from a prior normalization.
#'
#' @examples
#' \dontrun{
#' # Get network
#' network <- APOLLO_get_collectri()
#'
#' # Run ULM (recommended)
#' result <- ARTEMIS_run_decoupler(expr_matrix, network, method = "ulm")
#'
#' # Run MLM
#' result <- ARTEMIS_run_decoupler(expr_matrix, network, method = "mlm")
#'
#' }
#' @export
ARTEMIS_run_decoupler <- function(mat,
                                   network,
                                   method,
                                   minsize = 5,
                                   norm_method = "log2cpm",
                                   verbose = TRUE,
                                   ...) {

  # Validate method is provided
  if (missing(method)) {
    stop("'method' is required. Use ARTEMIS_decoupler_list_methods() to see available options.")
  }

  method <- tolower(method)
  valid_methods <- names(.DECOUPLER_METHODS)
  if (!method %in% valid_methods) {
    stop("Invalid method '", method, "'. Valid options: ",
         paste(valid_methods, collapse = ", "))
  }

  # Handle artemis_norm input
  if (inherits(mat, "artemis_norm")) {
    if (verbose) cat("[ARTEMIS] artemis_norm object detected — extracting log2(DESeq2-normalized + 1) counts.\n")
    mat <- log2(mat$norm_counts + 1)
    norm_method <- "none"
  }

  # Validate matrix
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("'mat' must be a matrix, data.frame, or artemis_norm object")
  }
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (is.null(rownames(mat))) {
    stop("'mat' must have rownames (gene identifiers)")
  }

  # Apply normalization or warn if raw counts passed with norm_method = "none"
  if (norm_method != "none") {
    if (verbose && .is_raw_counts(mat)) {
      cat("[ARTEMIS] Raw integer counts detected — applying", norm_method, "normalization.\n")
    }
    mat <- .normalize_for_decoupler(mat, norm_method, verbose)
  } else if (.is_raw_counts(mat)) {
    warning("Input appears to be raw integer counts but norm_method = 'none'. ",
            "DecoupleR expects normalized expression values. ",
            "Consider norm_method = 'log2cpm' or pre-normalizing with ARTEMIS_normalize_counts().",
            call. = FALSE)
  }

  # Validate network
  required_cols <- c("source", "target")
  if (!all(required_cols %in% colnames(network))) {
    stop("'network' must have columns: ", paste(required_cols, collapse = ", "))
  }
  if (!"mor" %in% colnames(network) && !"weight" %in% colnames(network)) {
    stop("'network' must have a 'mor' or 'weight' column")
  }

  if (verbose) {
    cat("[ARTEMIS] Activity Inference \n")
    cat("    Method:", .DECOUPLER_METHODS[[method]]$name, "(", method, ")\n")
    cat("    Input matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

    # Network info
    db_name <- attr(network, "database")
    if (!is.null(db_name)) {
      cat("[ARTEMIS] Network:", db_name, "\n")
    }
    cat("    Network sources:", length(unique(network$source)), "\n")
    cat("    Min targets per source:", minsize, "\n")
  }

  # Check gene overlap
  n_overlap <- sum(rownames(mat) %in% network$target)
  pct_overlap <- round(100 * n_overlap / nrow(mat), 1)
  if (verbose) {
    cat("    Gene overlap:", n_overlap, "/", nrow(mat), "(", pct_overlap, "%)\n")
  }
  if (n_overlap < 100) {
    warning("Low gene overlap (", n_overlap, "). Check that gene ID types match.")
  }

  if (verbose) cat("[ARTEMIS] Running", method, "...\n")

  # Run the appropriate method
  run_fn <- switch(method,
    ulm = decoupleR::run_ulm,
    mlm = decoupleR::run_mlm,
    viper = decoupleR::run_viper,
    wsum = decoupleR::run_wsum,
    wmean = decoupleR::run_wmean,
    udt = decoupleR::run_udt,
    mdt = decoupleR::run_mdt,
    ora = decoupleR::run_ora,
    gsva = decoupleR::run_gsva,
    aucell = decoupleR::run_aucell
  )

  # Run method
  results_long <- run_fn(
    mat = mat,
    network = network,
    minsize = minsize,
    ...
  )

  # Convert to wide format (activities matrix)
  activities <- .pivot_activities(results_long, value_col = "score")
  pvalues <- .pivot_activities(results_long, value_col = "p_value")

  # Build result object
  result <- list(
    activities = activities,
    pvalues = pvalues,
    results_long = results_long,
    method = method,
    method_name = .DECOUPLER_METHODS[[method]]$name,
    network_info = list(
      database = attr(network, "database"),
      organism = attr(network, "organism"),
      network_type = attr(network, "network_type"),
      n_sources = length(unique(network$source)),
      n_targets = length(unique(network$target))
    ),
    n_sources = nrow(activities),
    n_samples = ncol(activities),
    gene_overlap = n_overlap,
    minsize = minsize
  )

  class(result) <- c("decoupler_result", "list")

  if (verbose) {
    cat("[ARTEMIS] Done. Found activities for", result$n_sources, "sources.\n")
  }

  return(result)
}


#' Compare multiple decoupleR methods
#'
#' Runs multiple statistical methods on the same data and compares results.
#' Useful for assessing robustness and identifying consensus activities.
#'
#' @param mat Numeric matrix of gene expression data (genes x samples), or an
#'   \code{artemis_norm} object from \code{ARTEMIS_normalize_counts()}. Rownames
#'   must be gene identifiers matching the network.
#' @param network Prior knowledge network.
#' @param methods Character vector. Methods to compare. REQUIRED - no default.
#' @param minsize Integer. Minimum targets per source. Default: 5.
#' @param norm_method Character. Normalization to apply before running methods.
#'   \code{"log2cpm"} (default): log2(CPM + 1). \code{"none"}: pass through as-is.
#'   Normalization is applied once to the full matrix before any method runs to
#'   ensure consistency. If \code{mat} is an \code{artemis_norm} object, DESeq2
#'   normalized counts are extracted and log2-transformed automatically.
#' @param consensus Logical. Compute consensus scores across methods. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#' @param ... Additional arguments passed to each method.
#'
#' @return A list with class "decoupler_comparison" containing:
#'   \item{results}{Named list of decoupler_result objects (one per method)}
#'   \item{summary}{Data.frame summarizing results across methods}
#'   \item{correlations}{Method-method correlation matrix}
#'   \item{method_stats}{Per-method statistics}
#'   \item{consensus}{Consensus activity scores (if consensus = TRUE)}
#'   \item{methods}{Methods used}
#'   \item{n_methods}{Number of methods}
#'
#' @details
#' The summary table includes for each source (TF/pathway):
#' - Mean activity across methods
#' - Standard deviation across methods
#' - Number of methods where significant (p < 0.05)
#' - Agreement score (1 - normalized SD)
#'
#' This helps identify which activities are robust across different
#' statistical approaches vs. method-dependent.
#'
#' @examples
#' \dontrun{
#' network <- APOLLO_get_collectri()
#' comparison <- ARTEMIS_decoupler_compare_methods(
#'   expr_matrix, network,
#'   methods = c("ulm", "mlm", "wsum")
#' )
#'
#' # View summary
#' head(comparison$summary)
#'
#' # View method correlations
#' comparison$correlations
#'
#' }
#' @export
ARTEMIS_decoupler_compare_methods <- function(mat,
                                               network,
                                               methods,
                                               minsize = 5,
                                               norm_method = "log2cpm",
                                               consensus = TRUE,
                                               verbose = TRUE,
                                               ...) {

  # Validate methods
  if (missing(methods)) {
    stop("'methods' is required. Use ARTEMIS_decoupler_list_methods() to see options.")
  }

  methods <- tolower(methods)
  valid_methods <- names(.DECOUPLER_METHODS)
  invalid <- methods[!methods %in% valid_methods]
  if (length(invalid) > 0) {
    stop("Invalid method(s): ", paste(invalid, collapse = ", "),
         "\nValid options: ", paste(valid_methods, collapse = ", "))
  }

  if (length(methods) < 2) {
    stop("At least 2 methods required for comparison.")
  }

  # Handle artemis_norm input
  if (inherits(mat, "artemis_norm")) {
    if (verbose) cat("[ARTEMIS] artemis_norm object detected — extracting log2(DESeq2-normalized + 1) counts.\n")
    mat <- log2(mat$norm_counts + 1)
    norm_method <- "none"
  } else if (norm_method != "none") {
    if (verbose) cat("[ARTEMIS] Normalizing matrix once before running", length(methods), "methods...\n")
    if (verbose && .is_raw_counts(mat)) cat("[ARTEMIS] Raw integer counts detected — applying", norm_method, "normalization.\n")
    mat <- .normalize_for_decoupler(mat, norm_method, verbose = FALSE)
    norm_method <- "none"
  } else if (.is_raw_counts(mat)) {
    warning("Input appears to be raw integer counts but norm_method = 'none'. ",
            "DecoupleR expects normalized expression values. ",
            "Consider norm_method = 'log2cpm' or pre-normalizing with ARTEMIS_normalize_counts().",
            call. = FALSE)
  }

  if (verbose) {
    cat("[ARTEMIS] Method Comparison \n")
    cat("    Methods:", paste(methods, collapse = ", "), "\n")
    cat("    Input:", nrow(mat), "genes x", ncol(mat), "samples\n\n")
  }

  # Run each method (matrix already normalized above — pass norm_method = "none")
  results <- list()
  for (m in methods) {
    if (verbose) cat("[ARTEMIS] Running", m, "...\n")
    results[[m]] <- ARTEMIS_run_decoupler(
      mat = mat,
      network = network,
      method = m,
      minsize = minsize,
      norm_method = "none",
      verbose = FALSE,
      ...
    )
  }

  if (verbose) cat("[ARTEMIS] \nCompiling comparison...\n")

  # Find common sources across all methods
  common_sources <- Reduce(intersect, lapply(results, function(x) rownames(x$activities)))
  common_samples <- Reduce(intersect, lapply(results, function(x) colnames(x$activities)))

  if (length(common_sources) == 0) {
    stop("No common sources found across methods.")
  }

  if (verbose) {
    cat("    Common sources:", length(common_sources), "\n")
    cat("    Common samples:", length(common_samples), "\n")
  }

  # Extract aligned activity matrices
  activity_matrices <- lapply(results, function(x) {
    x$activities[common_sources, common_samples, drop = FALSE]
  })

  # Method-method correlations
  correlations <- .compute_method_correlations(activity_matrices)

  # Per-method stats
  method_stats <- .compute_method_stats(results, activity_matrices)

  # Build summary table
  summary_df <- .build_comparison_summary(results, activity_matrices, common_sources, common_samples)

  # Consensus scores (mean of z-scored activities)
  consensus_scores <- NULL
  if (consensus) {
    consensus_scores <- .compute_consensus(activity_matrices)
  }

  # Build result object
  comparison <- list(
    results = results,
    summary = summary_df,
    correlations = correlations,
    method_stats = method_stats,
    consensus = consensus_scores,
    methods = methods,
    n_methods = length(methods),
    common_sources = common_sources,
    common_samples = common_samples,
    network_info = results[[1]]$network_info
  )

  class(comparison) <- c("decoupler_comparison", "list")

  if (verbose) {
    cat("[ARTEMIS] Comparison Summary \n")
    cat("    Method correlations (mean):", round(mean(correlations[lower.tri(correlations)]), 3), "\n")
    cat("    Sources with high agreement (>0.7):", sum(summary_df$agreement_score > 0.7), "/",
        length(common_sources), "\n")
  }

  return(comparison)
}


# ==============================================================================
# CONVENIENCE WRAPPERS
# ==============================================================================

#' Infer transcription factor activity
#'
#' Convenience wrapper for TF activity inference. Fetches the appropriate
#' network and runs activity inference.
#'
#' @param mat Numeric matrix of gene expression data. Genes as rows, samples as columns.
#' @param method Character. Statistical method to use. REQUIRED - no default.
#' @param database Character. TF-target database to use.
#'   Options: "collectri" (default), "dorothea".
#' @param organism Character. Organism: "human" or "mouse". Default: "human".
#' @param dorothea_levels Character vector. If database = "dorothea", which
#'   confidence levels to include. Default: c("A", "B", "C", "D", "E").
#' @param minsize Integer. Minimum targets per TF. Default: 5.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @param ... Additional arguments passed to ARTEMIS_run_decoupler().
#'
#' @return A decoupler_result object with TF activity scores.
#'
#' @examples
#' \dontrun{
#' # Using CollecTRI (default)
#' tf_activity <- ARTEMIS_infer_tf_activity(expr_matrix, method = "ulm")
#'
#' # Using DoRothEA with high-confidence regulons
#' tf_activity <- ARTEMIS_infer_tf_activity(
#'   expr_matrix,
#'   method = "mlm",
#'   database = "dorothea",
#'   dorothea_levels = c("A", "B")
#' )
#'
#' }
#' @export
ARTEMIS_infer_tf_activity <- function(mat,
                                       method,
                                       database = "collectri",
                                       organism = "human",
                                       dorothea_levels = c("A", "B", "C", "D", "E"),
                                       minsize = 5,
                                       verbose = TRUE,
                                       ...) {

  if (missing(method)) {
    stop("'method' is required. Use ARTEMIS_decoupler_list_methods() to see options.")
  }

  database <- tolower(database)
  if (!database %in% c("collectri", "dorothea")) {
    stop("database must be 'collectri' or 'dorothea'")
  }

  if (verbose) {
    cat("[ARTEMIS] TF Activity Inference \n")
    cat("    Database:", database, "\n")
  }

  # Get network
  if (database == "collectri") {
    network <- APOLLO_get_collectri(organism = organism, verbose = verbose)
  } else {
    network <- APOLLO_get_dorothea(
      organism = organism,
      levels = dorothea_levels,
      verbose = verbose
    )
  }

  # Run inference
  result <- ARTEMIS_run_decoupler(
    mat = mat,
    network = network,
    method = method,
    minsize = minsize,
    verbose = verbose,
    ...
  )

  # Add TF-specific info
  result$analysis_type <- "tf_activity"

  return(result)
}


#' Infer pathway activity
#'
#' Convenience wrapper for pathway activity inference using PROGENy signatures.
#'
#' @param mat Numeric matrix of gene expression data. Genes as rows, samples as columns.
#' @param method Character. Statistical method to use. REQUIRED - no default.
#' @param database Character. Pathway database. Currently only "progeny" supported.
#'   Default: "progeny".
#' @param organism Character. Organism: "human" or "mouse". Default: "human".
#' @param top Integer. Number of top responsive genes per pathway. Default: 500.
#' @param minsize Integer. Minimum genes per pathway. Default: 5.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @param ... Additional arguments passed to ARTEMIS_run_decoupler().
#'
#' @return A decoupler_result object with pathway activity scores.
#'
#' @details
#' Note: This provides pathway activity via a different approach than
#' enrichment analysis (gprofiler2). PROGENy uses pathway-responsive genes
#' derived from perturbation experiments, while enrichment uses pathway
#' membership. Both approaches are complementary.
#'
#' @examples
#' \dontrun{
#' pathway_activity <- ARTEMIS_infer_pathway_activity(expr_matrix, method = "mlm")
#'
#' }
#' @export
ARTEMIS_infer_pathway_activity <- function(mat,
                                            method,
                                            database = "progeny",
                                            organism = "human",
                                            top = 500,
                                            minsize = 5,
                                            verbose = TRUE,
                                            ...) {

  if (missing(method)) {
    stop("'method' is required. Use ARTEMIS_decoupler_list_methods() to see options.")
  }

  database <- tolower(database)
  if (database != "progeny") {
    warning("Currently only 'progeny' is supported for pathway activity. Using progeny.")
    database <- "progeny"
  }

  if (verbose) {
    cat("[ARTEMIS] Pathway Activity Inference \n")
    cat("    Database:", database, "\n")
  }

  # Get network
  network <- APOLLO_get_progeny(organism = organism, top = top, verbose = verbose)

  # Run inference
  result <- ARTEMIS_run_decoupler(
    mat = mat,
    network = network,
    method = method,
    minsize = minsize,
    verbose = verbose,
    ...
  )

  # Add pathway-specific info
  result$analysis_type <- "pathway_activity"

  return(result)
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Pivot long results to wide activity matrix
#' @noRd
.pivot_activities <- function(results_long, value_col = "score") {

  if (!value_col %in% colnames(results_long)) {
    return(NULL)
  }

  # Get unique sources and conditions (samples)
  sources <- unique(results_long$source)
  conditions <- unique(results_long$condition)

  # Create matrix
  mat <- matrix(NA, nrow = length(sources), ncol = length(conditions),
                dimnames = list(sources, conditions))

  # Fill matrix
  for (i in seq_len(nrow(results_long))) {
    row <- results_long[i, ]
    mat[row$source, row$condition] <- row[[value_col]]
  }

  return(mat)
}


#' Compute method-method correlations
#' @noRd
.compute_method_correlations <- function(activity_matrices) {

  methods <- names(activity_matrices)
  n_methods <- length(methods)

  cor_mat <- matrix(1, nrow = n_methods, ncol = n_methods,
                    dimnames = list(methods, methods))

  for (i in seq_len(n_methods - 1)) {
    for (j in (i + 1):n_methods) {
      # Flatten matrices and correlate
      vec_i <- as.vector(activity_matrices[[i]])
      vec_j <- as.vector(activity_matrices[[j]])
      cor_val <- cor(vec_i, vec_j, use = "pairwise.complete.obs")
      cor_mat[i, j] <- cor_val
      cor_mat[j, i] <- cor_val
    }
  }

  return(cor_mat)
}


#' Compute per-method statistics
#' @noRd
.compute_method_stats <- function(results, activity_matrices) {

  stats <- data.frame(
    method = names(results),
    n_sources = sapply(results, function(x) x$n_sources),
    n_samples = sapply(results, function(x) x$n_samples),
    mean_activity = sapply(activity_matrices, function(x) mean(x, na.rm = TRUE)),
    sd_activity = sapply(activity_matrices, function(x) sd(x, na.rm = TRUE)),
    n_significant = sapply(results, function(x) {
      if (!is.null(x$pvalues)) {
        sum(x$pvalues < 0.05, na.rm = TRUE)
      } else {
        NA
      }
    }),
    stringsAsFactors = FALSE
  )

  return(stats)
}


#' Build comparison summary table
#' @noRd
.build_comparison_summary <- function(results, activity_matrices, sources, samples) {

  n_methods <- length(activity_matrices)
  methods <- names(activity_matrices)

  # For each source, compute stats across methods
  summary_list <- lapply(sources, function(src) {

    # Get activities for this source across all methods and samples
    activities <- sapply(activity_matrices, function(mat) mat[src, ])

    # Mean activity per method (across samples)
    mean_per_method <- colMeans(activities, na.rm = TRUE)

    # Get p-values if available
    pvals <- sapply(results, function(res) {
      if (!is.null(res$pvalues) && src %in% rownames(res$pvalues)) {
        min(res$pvalues[src, ], na.rm = TRUE)  # Use min p-value across samples
      } else {
        NA
      }
    })

    n_significant <- sum(pvals < 0.05, na.rm = TRUE)

    data.frame(
      source = src,
      mean_activity = mean(mean_per_method, na.rm = TRUE),
      sd_activity = sd(mean_per_method, na.rm = TRUE),
      min_activity = min(mean_per_method, na.rm = TRUE),
      max_activity = max(mean_per_method, na.rm = TRUE),
      n_methods_significant = n_significant,
      agreement_score = 1 - (sd(mean_per_method, na.rm = TRUE) /
                              (abs(mean(mean_per_method, na.rm = TRUE)) + 0.01)),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  # Clamp agreement score to [0, 1]
  summary_df$agreement_score <- pmax(0, pmin(1, summary_df$agreement_score))

  # Sort by mean activity
  summary_df <- summary_df[order(abs(summary_df$mean_activity), decreasing = TRUE), ]

  return(summary_df)
}


#' Compute consensus scores
#' @noRd
.compute_consensus <- function(activity_matrices) {

  # Z-score each method's activities
  z_matrices <- lapply(activity_matrices, function(mat) {
    (mat - mean(mat, na.rm = TRUE)) / sd(mat, na.rm = TRUE)
  })

  # Stack into 3D array and take mean
  sources <- rownames(z_matrices[[1]])
  samples <- colnames(z_matrices[[1]])
  n_methods <- length(z_matrices)

  # Initialize consensus matrix
  consensus <- matrix(0, nrow = length(sources), ncol = length(samples),
                      dimnames = list(sources, samples))

  for (mat in z_matrices) {
    consensus <- consensus + mat
  }
  consensus <- consensus / n_methods

  return(consensus)
}


# ==============================================================================
# PRINT METHODS
# ==============================================================================

#' @method print decoupler_result
#' @export
print.decoupler_result <- function(x, ...) {
  cat("decoupleR Activity Result\n")
  cat("------------------------------\n")
  cat("Method:", x$method_name, "(", x$method, ")\n")
  cat("Sources:", x$n_sources, "\n")
  cat("Samples:", x$n_samples, "\n")

  if (!is.null(x$network_info$database)) {
    cat("Network:", x$network_info$database, "\n")
  }

  cat("Gene overlap:", x$gene_overlap, "\n")

  if (!is.null(x$analysis_type)) {
    cat("Analysis type:", x$analysis_type, "\n")
  }

  cat("\nUse $activities for activity matrix, $pvalues for p-values.\n")
}


#' @method print decoupler_comparison
#' @export
print.decoupler_comparison <- function(x, ...) {
  cat("decoupleR Method Comparison\n")
  cat("---------------------------\n")
  cat("Methods:", paste(x$methods, collapse = ", "), "\n")
  cat("Common sources:", length(x$common_sources), "\n")
  cat("Common samples:", length(x$common_samples), "\n")

  cat("\nMethod correlations:\n")
  print(round(x$correlations, 3))

  cat("\nTop 5 sources by mean activity:\n")
  print(head(x$summary[, c("source", "mean_activity", "agreement_score")], 5))

  cat("\nUse $results for per-method results, $summary for full comparison,\n")
  cat("$consensus for consensus scores.\n")
}
