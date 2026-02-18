# GAIA/Demeter/result_loaders.R
# Functions for loading saved analysis results
#
# Helper functions to reload exported results (e.g., from ELEUTHIA_export_*)
# without re-running the full analysis pipeline.


# ==============================================================================
# INTERNAL UTILITIES
# ==============================================================================

#' Strip Ensembl version suffix from gene IDs
#'
#' Removes the version number after the dot (e.g., ENSG00000103888.18 -> ENSG00000103888).
#' Works on character vectors, data.frame columns, matrix rownames, and named vectors.
#'
#' @param x Character vector of gene IDs.
#' @return Character vector with version suffixes removed.
#' @noRd
.strip_ensembl_version <- function(x) {
  sub("\\.\\d+$", "", x)
}


#' Load WGCNA results from export directory
#'
#' Loads previously exported WGCNA results. Prioritizes RDS files (complete objects)
#' and falls back to CSV reconstruction if RDS not available.
#'
#' @param results_dir Path to the WGCNA export directory (output of ELEUTHIA_export_wgcna_results).
#' @param prefix Filename prefix used during export. Default: "wgcna".
#' @param what Character vector specifying which results to load. Options:
#'   "modules", "trait_cor", "gene_sig", "hubs", "enrichment", "all".
#'   Default: "all".
#' @param strip_version Logical. If TRUE, strips Ensembl version suffixes
#'   (e.g., ".12") from gene IDs. Default: FALSE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A named list containing the requested result objects. Each element
#'   will be the appropriate class (e.g., wgcna_modules, wgcna_trait_cor) if
#'   loaded from RDS, or a reconstructed object if loaded from CSV.
#'
#' @details
#' The function first looks for RDS files (e.g., wgcna_modules.rds) which contain
#' the complete R objects with all attributes and class information. If RDS files
#' are not found, it attempts to reconstruct the objects from CSV files.
#'
#' For full compatibility with ARTEMIS functions (e.g., ARTEMIS_wgcna_module_traits),
#' ensure that `save_rds = TRUE` was used during export.
#'
#' @examples
#' \dontrun{
#' # Load all WGCNA results
#' results <- DEMETER_load_wgcna("results/wgcna_output")
#'
#' # Use modules with decoupleR integration
#' tf_trait_cor <- ARTEMIS_wgcna_module_traits(
#'   modules = results$modules,
#'   traits = t(tf_activities)
#' )
#'
#' # Load only modules
#' results <- DEMETER_load_wgcna("results/wgcna_output", what = "modules")
#'
#' }
#' @export
DEMETER_load_wgcna <- function(results_dir,
                                prefix = "wgcna",
                                what = "all",
                                strip_version = FALSE,
                                verbose = TRUE) {

  if (!dir.exists(results_dir)) {
    stop("Directory not found: ", results_dir)
  }

  if (verbose) cat("=== Loading WGCNA Results ===\n")
  if (verbose) cat("Directory:", results_dir, "\n")

  # Parse 'what' argument
  valid_types <- c("modules", "trait_cor", "gene_sig", "hubs", "enrichment")
  if ("all" %in% what) {
    what <- valid_types
  } else {
    invalid <- what[!what %in% valid_types]
    if (length(invalid) > 0) {
      stop("Invalid result type(s): ", paste(invalid, collapse = ", "),
           "\nValid options: ", paste(valid_types, collapse = ", "), ", all")
    }
  }

  results <- list()

  # --------------------------------------------------------------------------
  # Load modules
  # --------------------------------------------------------------------------
  if ("modules" %in% what) {
    modules <- .load_wgcna_modules(results_dir, prefix, verbose)
    if (!is.null(modules)) {
      results$modules <- modules
    }
  }

  # --------------------------------------------------------------------------
  # Load trait correlations
  # --------------------------------------------------------------------------
  if ("trait_cor" %in% what) {
    trait_cor <- .load_wgcna_trait_cor(results_dir, prefix, verbose)
    if (!is.null(trait_cor)) {
      results$trait_cor <- trait_cor
    }
  }

  # --------------------------------------------------------------------------
  # Load gene significance
  # --------------------------------------------------------------------------
  if ("gene_sig" %in% what) {
    gene_sig <- .load_wgcna_gene_sig(results_dir, prefix, verbose)
    if (!is.null(gene_sig)) {
      results$gene_sig <- gene_sig
    }
  }

  # --------------------------------------------------------------------------
  # Load hub genes
  # --------------------------------------------------------------------------
  if ("hubs" %in% what) {
    hubs <- .load_wgcna_hubs(results_dir, prefix, verbose)
    if (!is.null(hubs)) {
      results$hubs <- hubs
    }
  }

  # --------------------------------------------------------------------------
  # Load enrichment
  # --------------------------------------------------------------------------
  if ("enrichment" %in% what) {
    enrichment <- .load_wgcna_enrichment(results_dir, prefix, verbose)
    if (!is.null(enrichment)) {
      results$enrichment <- enrichment
    }
  }

  # --- Strip Ensembl version suffixes if requested ---
  if (strip_version) {
    if (verbose) cat("Stripping Ensembl version suffixes from gene IDs...\n")

    if (!is.null(results$modules)) {
      m <- results$modules
      # gene_module_df
      if (!is.null(m$gene_module_df)) {
        m$gene_module_df$gene <- .strip_ensembl_version(m$gene_module_df$gene)
        rownames(m$gene_module_df) <- m$gene_module_df$gene
      }
      # module_names (named vector: gene -> module_name)
      if (!is.null(m$module_names) && !is.null(names(m$module_names))) {
        names(m$module_names) <- .strip_ensembl_version(names(m$module_names))
      }
      # module_labels (named vector: gene -> numeric label)
      if (!is.null(m$module_labels) && !is.null(names(m$module_labels))) {
        names(m$module_labels) <- .strip_ensembl_version(names(m$module_labels))
      }
      results$modules <- m
    }

    if (!is.null(results$gene_sig)) {
      gs <- results$gene_sig
      if (!is.null(gs$gene_info)) {
        gene_col <- intersect(c("gene", "gene_id", "feature_id"), colnames(gs$gene_info))[1]
        if (!is.na(gene_col)) {
          gs$gene_info[[gene_col]] <- .strip_ensembl_version(gs$gene_info[[gene_col]])
        }
      }
      results$gene_sig <- gs
    }

    if (!is.null(results$hubs)) {
      h <- results$hubs
      if (!is.null(h$hub_genes)) {
        gene_col <- intersect(c("gene", "gene_id", "feature_id"), colnames(h$hub_genes))[1]
        if (!is.na(gene_col)) {
          h$hub_genes[[gene_col]] <- .strip_ensembl_version(h$hub_genes[[gene_col]])
        }
      }
      results$hubs <- h
    }
  }

  if (verbose) {
    cat("\n--- Loaded ---\n")
    cat("Objects:", paste(names(results), collapse = ", "), "\n")
  }

  return(results)
}


# ==============================================================================
# INTERNAL HELPER FUNCTIONS
# ==============================================================================

#' Load modules (RDS or CSV reconstruction)
#' @noRd
.load_wgcna_modules <- function(results_dir, prefix, verbose) {

  # Try RDS first
  rds_file <- file.path(results_dir, paste0(prefix, "_modules.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("Loading modules from RDS...\n")
    modules <- readRDS(rds_file)
    if (verbose) cat("  Loaded wgcna_modules object\n")
    return(modules)
  }

  # Fall back to CSV reconstruction
  if (verbose) cat("Loading modules from CSV (RDS not found)...\n")

  # Required files
  me_file <- file.path(results_dir, paste0(prefix, "_module_eigengenes.csv"))
  gm_file <- file.path(results_dir, paste0(prefix, "_gene_modules.csv"))
  summary_file <- file.path(results_dir, paste0(prefix, "_module_summary.csv"))

  if (!file.exists(me_file)) {
    if (verbose) cat("  Module eigengenes CSV not found, skipping modules\n")
    return(NULL)
  }

  # Load eigengenes
  me_df <- read.csv(me_file, stringsAsFactors = FALSE)
  rownames(me_df) <- me_df$sample
  me_df$sample <- NULL
  module_eigengenes <- as.data.frame(me_df)

  # Load gene-module assignments
  gene_module_df <- NULL
  module_names <- NULL
  if (file.exists(gm_file)) {
    gene_module_df <- read.csv(gm_file, stringsAsFactors = FALSE)
    module_names <- setNames(gene_module_df$module_name, gene_module_df$gene)
  }

  # Load summary
  module_summary <- NULL
  if (file.exists(summary_file)) {
    module_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  }

  # Load color map if available
  color_map_file <- file.path(results_dir, paste0(prefix, "_color_map.csv"))
  module_colors <- NULL
  color_names <- NULL
  if (file.exists(color_map_file)) {
    cm <- read.csv(color_map_file, stringsAsFactors = FALSE)
    module_colors <- setNames(cm$display_color, cm$module)
    if ("color_name" %in% colnames(cm)) {
      color_names <- setNames(cm$color_name, cm$module)
    }
  } else if (!is.null(module_names)) {
    # Generate colors if no color map file
    unique_mods <- sort(unique(module_names))
    module_colors <- DEMETER_generate_colors(length(unique_mods))
    if ("module_0" %in% unique_mods) {
      module_colors[which(unique_mods == "module_0")] <- "#AAAAAA"
    }
    module_colors <- setNames(module_colors, unique_mods)
  }

  # Reconstruct module_labels from module_names
  module_labels <- NULL
  if (!is.null(module_names)) {
    module_labels <- as.integer(gsub("^module_", "", module_names))
    names(module_labels) <- names(module_names)
  }

  # Reconstruct object
  modules <- list(
    module_names = module_names,
    module_labels = module_labels,
    module_colors = module_colors,
    color_names = color_names,
    module_eigengenes = module_eigengenes,
    gene_module_df = gene_module_df,
    module_summary = module_summary,
    n_modules = if (!is.null(module_summary)) nrow(module_summary) else ncol(module_eigengenes),
    n_genes = if (!is.null(module_names)) length(module_names) else NA,
    .reconstructed = TRUE  # Flag that this was reconstructed from CSV
  )

  class(modules) <- c("wgcna_modules", "list")

  if (verbose) cat("  Reconstructed wgcna_modules from CSV\n")
  return(modules)
}


#' Load trait correlations (RDS or CSV reconstruction)
#' @noRd
.load_wgcna_trait_cor <- function(results_dir, prefix, verbose) {

  # Try RDS first
  rds_file <- file.path(results_dir, paste0(prefix, "_trait_cor.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("Loading trait correlations from RDS...\n")
    trait_cor <- readRDS(rds_file)
    if (verbose) cat("  Loaded wgcna_trait_cor object\n")
    return(trait_cor)
  }

  # Fall back to CSV reconstruction
  if (verbose) cat("Loading trait correlations from CSV...\n")

  cor_file <- file.path(results_dir, paste0(prefix, "_trait_correlations.csv"))
  pval_file <- file.path(results_dir, paste0(prefix, "_trait_pvalues.csv"))
  sig_file <- file.path(results_dir, paste0(prefix, "_significant_associations.csv"))

  if (!file.exists(cor_file)) {
    if (verbose) cat("  Correlation CSV not found, skipping trait_cor\n")
    return(NULL)
  }

  # Load correlations
  cor_df <- read.csv(cor_file, stringsAsFactors = FALSE)
  rownames(cor_df) <- cor_df$module
  cor_df$module <- NULL
  cor_matrix <- as.matrix(cor_df)

  # Load p-values
  pvalue_matrix <- NULL
  if (file.exists(pval_file)) {
    pval_df <- read.csv(pval_file, stringsAsFactors = FALSE)
    rownames(pval_df) <- pval_df$module
    pval_df$module <- NULL
    pvalue_matrix <- as.matrix(pval_df)
  }

  # Load significant associations
  significant_associations <- data.frame()
  if (file.exists(sig_file)) {
    significant_associations <- read.csv(sig_file, stringsAsFactors = FALSE)
  }

  # Reconstruct object
  trait_cor <- list(
    cor_matrix = cor_matrix,
    pvalue_matrix = pvalue_matrix,
    padj_matrix = pvalue_matrix,  # Assume same if not stored separately
    significant_associations = significant_associations,
    n_significant = nrow(significant_associations),
    .reconstructed = TRUE
  )

  class(trait_cor) <- c("wgcna_trait_cor", "list")

  if (verbose) cat("  Reconstructed wgcna_trait_cor from CSV\n")
  return(trait_cor)
}


#' Load gene significance (RDS or CSV reconstruction)
#' @noRd
.load_wgcna_gene_sig <- function(results_dir, prefix, verbose) {

  # Try RDS first
  rds_file <- file.path(results_dir, paste0(prefix, "_gene_sig.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("Loading gene significance from RDS...\n")
    gene_sig <- readRDS(rds_file)
    if (verbose) cat("  Loaded wgcna_gene_sig object\n")
    return(gene_sig)
  }

  # Fall back to CSV
  if (verbose) cat("Loading gene significance from CSV...\n")

  gs_file <- file.path(results_dir, paste0(prefix, "_gene_significance.csv"))
  if (!file.exists(gs_file)) {
    if (verbose) cat("  Gene significance CSV not found, skipping\n")
    return(NULL)
  }

  gene_info <- read.csv(gs_file, stringsAsFactors = FALSE)

  gene_sig <- list(
    gene_info = gene_info,
    .reconstructed = TRUE
  )

  class(gene_sig) <- c("wgcna_gene_sig", "list")

  if (verbose) cat("  Reconstructed wgcna_gene_sig from CSV\n")
  return(gene_sig)
}


#' Load hub genes (RDS or CSV reconstruction)
#' @noRd
.load_wgcna_hubs <- function(results_dir, prefix, verbose) {

  # Try RDS first
  rds_file <- file.path(results_dir, paste0(prefix, "_hubs.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("Loading hub genes from RDS...\n")
    hubs <- readRDS(rds_file)
    if (verbose) cat("  Loaded wgcna_hubs object\n")
    return(hubs)
  }

  # Fall back to CSV
  if (verbose) cat("Loading hub genes from CSV...\n")

  hub_file <- file.path(results_dir, paste0(prefix, "_hub_genes.csv"))
  summary_file <- file.path(results_dir, paste0(prefix, "_hub_summary.csv"))

  if (!file.exists(hub_file)) {
    if (verbose) cat("  Hub genes CSV not found, skipping\n")
    return(NULL)
  }

  hub_genes <- read.csv(hub_file, stringsAsFactors = FALSE)
  hub_summary <- NULL
  if (file.exists(summary_file)) {
    hub_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  }

  hubs <- list(
    hub_genes = hub_genes,
    hub_summary = hub_summary,
    .reconstructed = TRUE
  )

  class(hubs) <- c("wgcna_hubs", "list")

  if (verbose) cat("  Reconstructed wgcna_hubs from CSV\n")
  return(hubs)
}


#' Load enrichment results (RDS only - too complex for CSV reconstruction)
#' @noRd
.load_wgcna_enrichment <- function(results_dir, prefix, verbose) {

  # Check for enrichment subdirectory
  enrich_dir <- file.path(results_dir, paste0(prefix, "_enrichment"))

  # Try RDS in enrichment subdirectory
  rds_file <- file.path(enrich_dir, "enrichment.rds")
  if (file.exists(rds_file)) {
    if (verbose) cat("Loading enrichment from RDS...\n")
    enrichment <- readRDS(rds_file)
    if (verbose) cat("  Loaded gost_enrichment object\n")
    return(enrichment)
  }

  # Try loading combined CSV (partial reconstruction)
  combined_file <- file.path(enrich_dir, "enrichment_combined.csv")
  if (file.exists(combined_file)) {
    if (verbose) cat("Loading enrichment from CSV (limited reconstruction)...\n")
    combined <- read.csv(combined_file, stringsAsFactors = FALSE)

    enrichment <- list(
      combined = combined,
      summary = NULL,
      .reconstructed = TRUE,
      .note = "Loaded from CSV - some features may be limited"
    )

    class(enrichment) <- c("gost_enrichment", "list")
    if (verbose) cat("  Reconstructed gost_enrichment from CSV (limited)\n")
    return(enrichment)
  }

  if (verbose) cat("  Enrichment results not found\n")
  return(NULL)
}


#' Select top N sources from activity results
#'
#' Filters activity results to keep only the top N sources (TFs/pathways)
#' based on variance, mean absolute activity, or other criteria.
#'
#' @param result A decoupler_result object, or a matrix of activities
#'   (sources as rows, samples as columns).
#' @param n Integer. Number of top sources to keep. Default: 20.
#' @param by Character. Criterion for selecting top sources:
#'   - "variance": highest variance across samples (most variable)
#'   - "activity": highest mean absolute activity (strongest signals)
#'   - "max": highest maximum absolute activity
#'   Default: "variance".
#' @param verbose Logical. Print selection info. Default: TRUE.
#'
#' @return If input was decoupler_result, returns filtered decoupler_result.
#'   If input was matrix, returns filtered matrix.
#'
#' @details
#' Use this to reduce the number of TFs/pathways before:
#' - Plotting heatmaps (too many rows are unreadable
#' - Using as traits in WGCNA module-trait correlation
#' - Any analysis where you want to focus on the most informative sources
#'
#' @examples
#' \dontrun{
#' # Filter to top 20 most variable TFs
#' tf_top <- DEMETER_top_activities(tf_result, n = 20, by = "variance")
#'
#' # Use with WGCNA
#' tf_trait_cor <- ARTEMIS_wgcna_module_traits(
#'   modules = modules,
#'   traits = t(tf_top$activities)
#' )
#'
#' # Filter to top 30 by activity strength
#' tf_top <- DEMETER_top_activities(tf_result, n = 30, by = "activity")
#'
#' }
#' @export
DEMETER_top_activities <- function(result,
                                    n = 20,
                                    by = "variance",
                                    verbose = TRUE) {

  # Extract matrix

  if (inherits(result, "decoupler_result")) {
    mat <- result$activities
    is_decoupler <- TRUE
  } else if (is.matrix(result)) {
    mat <- result
    is_decoupler <- FALSE
  } else {
    stop("'result' must be a decoupler_result object or matrix")
  }

  # Validate 'by'
  valid_by <- c("variance", "activity", "max")
  if (!by %in% valid_by) {
    stop("'by' must be one of: ", paste(valid_by, collapse = ", "))
  }

  # Calculate ranking metric
  if (by == "variance") {
    scores <- apply(mat, 1, var, na.rm = TRUE)
    metric_name <- "variance"
  } else if (by == "activity") {
    scores <- apply(mat, 1, function(x) mean(abs(x), na.rm = TRUE))
    metric_name <- "mean |activity|"
  } else if (by == "max") {
    scores <- apply(mat, 1, function(x) max(abs(x), na.rm = TRUE))
    metric_name <- "max |activity|"
  }

  # Handle n > nrow
  n <- min(n, nrow(mat))

  # Select top N
  top_sources <- names(sort(scores, decreasing = TRUE))[1:n]

  if (verbose) {
    cat("Selected top", n, "sources by", metric_name, "\n")
    cat("  Range:", round(min(scores[top_sources]), 3), "to",
        round(max(scores[top_sources]), 3), "\n")
  }

  # Filter and return
  if (is_decoupler) {
    # Return filtered decoupler_result
    result$activities <- mat[top_sources, , drop = FALSE]
    if (!is.null(result$pvalues)) {
      result$pvalues <- result$pvalues[top_sources, , drop = FALSE]
    }
    result$n_sources <- n
    result$top_selection <- list(
      n = n,
      by = by,
      sources = top_sources
    )
    return(result)
  } else {
    return(mat[top_sources, , drop = FALSE])
  }
}


#' Load decoupleR activity results
#'
#' Loads previously saved decoupleR results from RDS files.
#' Use this to reload activity inference results without re-running the analysis.
#'
#' @param file_path Path to the RDS file containing decoupler results.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return The loaded decoupler_result or decoupler_comparison object.
#'
#' @examples
#' \dontrun{
#' # Save results (in your analysis script)
#' saveRDS(tf_result, "results/tf_activity.rds")
#'
#' # Load later
#' tf_result <- DEMETER_load_activity("results/tf_activity.rds")
#'
#' }
#' @export
DEMETER_load_activity <- function(file_path, verbose = TRUE) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  if (verbose) cat("Loading activity results from:", file_path, "\n")

  result <- readRDS(file_path)

  if (inherits(result, "decoupler_result")) {
    if (verbose) {
      cat("  Type: decoupler_result\n")
      cat("  Method:", result$method, "\n")
      cat("  Sources:", result$n_sources, "\n")
      cat("  Samples:", result$n_samples, "\n")
    }
  } else if (inherits(result, "decoupler_comparison")) {
    if (verbose) {
      cat("  Type: decoupler_comparison\n")
      cat("  Methods:", paste(result$methods, collapse = ", "), "\n")
      cat("  Sources:", length(result$common_sources), "\n")
    }
  } else {
    warning("Loaded object is not a recognized decoupler result type")
  }

  return(result)
}


#' Load and combine DEA gene lists from multiple experiments
#'
#' Scans a parent directory for DEA result subdirectories and collects
#' significant gene IDs across experiments. Useful for building a combined
#' gene selection for downstream analyses (e.g., PART clustering).
#'
#' @param dea_dir Path to the parent directory containing DEA result
#'   subdirectories (one per experiment/comparison).
#' @param source Character. Which file to read genes from:
#'   - "significant" (default): reads the pre-filtered significant results CSV
#'   - "selected": reads the selected genes CSV (gene_id column)
#'   - "all": reads the full results CSV and applies threshold filtering
#' @param prefix Character. Filename prefix used during export. Default: "dea".
#' @param experiments Character vector. Specific subdirectory names to include.
#'   If NULL (default), all subdirectories containing the expected files are used.
#' @param l2fc_thresh Numeric. Absolute log2 fold-change threshold (only used
#'   when source = "all"). Default: 1.0.
#' @param p_thresh Numeric. Adjusted p-value threshold (only used when
#'   source = "all"). Default: 0.05.
#' @param gene_col Character. Column name containing gene IDs. If NULL (default),
#'   auto-detected from: feature_id, gene_id, gene, ensembl_id.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Character vector of unique gene IDs (union across all experiments).
#'   Attributes:
#'   - "per_experiment": named list of gene vectors per experiment
#'   - "n_per_experiment": named integer vector of gene counts
#'   - "source": which file source was used
#'   - "parameters": list of thresholds used (if source = "all")
#'
#' @examples
#' \dontrun{
#' # Default: collect from pre-filtered significant files
#' genes <- DEMETER_load_dea_genes("results/DEA/")
#'
#' # Only specific experiments
#' genes <- DEMETER_load_dea_genes("results/DEA/", experiments = c("C5RO", "CS1AN"))
#'
#' # Custom thresholds on all results
#' genes <- DEMETER_load_dea_genes("results/DEA/", source = "all",
#'                                  l2fc_thresh = 1.5, p_thresh = 0.01)
#'
#' # Use as gene selection for PART
#' part_mat <- ARTEMIS_prepare_part_matrix(norm_data, genes = genes)
#'
#' }
#' @export
DEMETER_load_dea_genes <- function(dea_dir,
                                    source = "significant",
                                    prefix = "dea",
                                    experiments = NULL,
                                    l2fc_thresh = 1.0,
                                    p_thresh = 0.05,
                                    gene_col = NULL,
                                    verbose = TRUE) {

  if (!dir.exists(dea_dir)) {
    stop("Directory not found: ", dea_dir)
  }

  valid_sources <- c("significant", "selected", "all")
  if (!source %in% valid_sources) {
    stop("'source' must be one of: ", paste(valid_sources, collapse = ", "))
  }

  if (verbose) cat("=== Loading DEA Genes ===\n")
  if (verbose) cat("Directory:", dea_dir, "\n")
  if (verbose) cat("Source:", source, "\n")

  # ---- Identify experiment subdirectories ----
  all_dirs <- list.dirs(dea_dir, full.names = FALSE, recursive = FALSE)

  if (!is.null(experiments)) {
    missing <- experiments[!experiments %in% all_dirs]
    if (length(missing) > 0) {
      warning("Experiment directories not found: ", paste(missing, collapse = ", "))
    }
    all_dirs <- all_dirs[all_dirs %in% experiments]
  }

  if (length(all_dirs) == 0) {
    stop("No experiment subdirectories found in: ", dea_dir)
  }

  # ---- Determine target file per source type ----
  target_file <- switch(source,
    "significant" = paste0(prefix, "_significant.csv"),
    "selected"    = paste0(prefix, "_selected_genes.csv"),
    "all"         = paste0(prefix, "_all_results.csv")
  )

  # ---- Collect genes from each experiment ----
  per_experiment <- list()

  for (exp_name in all_dirs) {
    filepath <- file.path(dea_dir, exp_name, target_file)

    if (!file.exists(filepath)) {
      if (verbose) cat("  [SKIP]", exp_name, "- file not found:", target_file, "\n")
      next
    }

    df <- read.csv(filepath, stringsAsFactors = FALSE)

    # Auto-detect gene column
    col <- gene_col
    if (is.null(col)) {
      candidates <- c("feature_id", "gene_id", "gene", "ensembl_id", "symbol")
      col <- candidates[candidates %in% colnames(df)][1]
      if (is.na(col)) {
        # Fall back to first column
        col <- colnames(df)[1]
        if (verbose) cat("  [WARN]", exp_name, "- no standard gene column found, using:", col, "\n")
      }
    }

    if (!col %in% colnames(df)) {
      if (verbose) cat("  [SKIP]", exp_name, "- column '", col, "' not found\n")
      next
    }

    # Apply threshold filtering for 'all' source
    if (source == "all") {
      # Detect l2fc column
      lfc_col <- intersect(c("log2FoldChange", "logFC", "l2fc"), colnames(df))[1]
      p_col <- intersect(c("padj", "p_adjust", "FDR", "adj.P.Val", "pvalue"), colnames(df))[1]

      if (is.na(lfc_col) || is.na(p_col)) {
        if (verbose) cat("  [SKIP]", exp_name, "- could not find l2fc/p-value columns\n")
        next
      }

      df <- df[!is.na(df[[p_col]]) &
               abs(df[[lfc_col]]) >= l2fc_thresh &
               df[[p_col]] <= p_thresh, , drop = FALSE]
    }

    genes <- unique(df[[col]])
    genes <- genes[!is.na(genes) & genes != ""]
    per_experiment[[exp_name]] <- genes

    if (verbose) cat("  ", exp_name, ":", length(genes), "genes\n")
  }

  if (length(per_experiment) == 0) {
    stop("No genes found in any experiment. Check directory structure and file prefix.")
  }

  # ---- Combine (union) ----
  all_genes <- unique(unlist(per_experiment, use.names = FALSE))

  if (verbose) {
    cat("\n--- Summary ---\n")
    cat("Experiments loaded:", length(per_experiment), "/", length(all_dirs), "\n")
    cat("Total unique genes:", length(all_genes), "\n")
  }

  # Attach metadata as attributes
  attr(all_genes, "per_experiment") <- per_experiment
  attr(all_genes, "n_per_experiment") <- vapply(per_experiment, length, integer(1))
  attr(all_genes, "source") <- source
  if (source == "all") {
    attr(all_genes, "parameters") <- list(
      l2fc_thresh = l2fc_thresh,
      p_thresh = p_thresh
    )
  }

  return(all_genes)
}


#' Load DEA results from export directories
#'
#' Loads previously exported DEA results from one or more experiment
#' subdirectories. Prioritizes RDS files (complete objects) and falls back
#' to CSV reconstruction. Optionally loads associated PART results.
#'
#' @param dea_dir Path to the parent directory containing DEA result
#'   subdirectories (one per experiment/comparison).
#' @param experiments Character vector. Specific subdirectory names to load.
#'   If NULL (default), all subdirectories are scanned.
#' @param prefix Character. Filename prefix used during export. Default: "dea".
#' @param what Character vector specifying what to load per experiment:
#'   - "dea": the DEA result (default)
#'   - "part": the PART clustering result
#'   - "all": both dea and part
#' @param strip_version Logical. If TRUE, strips Ensembl version suffixes
#'   from gene IDs. Default: FALSE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A named list (one element per experiment). Each element is a list
#'   containing the requested objects:
#'   \describe{
#'     \item{dea}{DEA result list with: results (data.frame), summary, comparison,
#'       norm_counts. If loaded from RDS, includes dds object.}
#'     \item{part}{artemis_part object (if requested and available)}
#'   }
#'   If only "dea" is requested, each element is the DEA result directly
#'   (not wrapped in a sub-list) for convenience.
#'
#' @examples
#' \dontrun{
#' # Load all DEA results
#' dea_list <- DEMETER_load_dea("results/DEA/")
#' # Returns: list(C5RO = <dea_result>, CS1AN = <dea_result>, ...)
#'
#' # Pass directly to circos prep
#' circos_data <- DEMETER_prepare_part_wgcna_circos(
#'   part_result, wgcna$modules, wgcna$trait_cor,
#'   dea_results = dea_list
#' )
#'
#' # Load specific experiments with PART results
#' results <- DEMETER_load_dea("results/DEA/",
#'                              experiments = c("C5RO", "CS1AN"),
#'                              what = "all")
#' # results$C5RO$dea, results$C5RO$part
#'
#' }
#' @export
DEMETER_load_dea <- function(dea_dir,
                              experiments = NULL,
                              prefix = "dea",
                              what = "dea",
                              strip_version = FALSE,
                              verbose = TRUE) {

  if (!dir.exists(dea_dir)) {
    stop("Directory not found: ", dea_dir)
  }

  # Parse 'what'
  if ("all" %in% what) {
    what <- c("dea", "part")
  }
  valid_what <- c("dea", "part")
  invalid <- what[!what %in% valid_what]
  if (length(invalid) > 0) {
    stop("Invalid 'what' values: ", paste(invalid, collapse = ", "),
         "\nValid options: dea, part, all")
  }

  load_dea <- "dea" %in% what
  load_part <- "part" %in% what
  dea_only <- load_dea && !load_part

  if (verbose) {
    cat("=== Loading DEA Results ===\n")
    cat("Directory:", dea_dir, "\n")
    cat("Loading:", paste(what, collapse = ", "), "\n")
  }

  # Identify subdirectories
  all_dirs <- list.dirs(dea_dir, full.names = FALSE, recursive = FALSE)

  if (!is.null(experiments)) {
    missing <- experiments[!experiments %in% all_dirs]
    if (length(missing) > 0) {
      warning("Experiment directories not found: ", paste(missing, collapse = ", "))
    }
    all_dirs <- all_dirs[all_dirs %in% experiments]
  }

  if (length(all_dirs) == 0) {
    stop("No experiment subdirectories found in: ", dea_dir)
  }

  results <- list()

  for (exp_name in all_dirs) {
    exp_dir <- file.path(dea_dir, exp_name)
    exp_result <- list()

    if (verbose) cat("\n", exp_name, ":\n", sep = "")

    # --- Load DEA ---
    if (load_dea) {
      dea <- .load_dea_result(exp_dir, prefix, verbose)
      if (!is.null(dea)) {
        exp_result$dea <- dea
      }
    }

    # --- Load PART ---
    if (load_part) {
      part <- .load_part_result(exp_dir, prefix, verbose)
      if (!is.null(part)) {
        exp_result$part <- part
      }
    }

    # --- Strip Ensembl version suffixes if requested ---
    if (strip_version && length(exp_result) > 0) {

      if (!is.null(exp_result$dea)) {
        d <- exp_result$dea
        # Strip from results data.frame
        if (!is.null(d$results)) {
          gene_col <- intersect(c("feature_id", "gene_id", "gene"), colnames(d$results))[1]
          if (!is.na(gene_col)) {
            d$results[[gene_col]] <- .strip_ensembl_version(d$results[[gene_col]])
          }
        }
        # Strip from norm_counts rownames
        if (!is.null(d$norm_counts)) {
          rownames(d$norm_counts) <- .strip_ensembl_version(rownames(d$norm_counts))
        }
        exp_result$dea <- d
      }

      if (!is.null(exp_result$part)) {
        p <- exp_result$part
        # Strip from cluster_map
        if (!is.null(p$cluster_map)) {
          p$cluster_map$gene <- .strip_ensembl_version(p$cluster_map$gene)
          rownames(p$cluster_map) <- p$cluster_map$gene
        }
        # Strip from clusters named vector
        if (!is.null(p$clusters)) {
          names(p$clusters) <- .strip_ensembl_version(names(p$clusters))
        }
        # Strip from data matrix rownames
        if (!is.null(p$data)) {
          rownames(p$data) <- .strip_ensembl_version(rownames(p$data))
        }
        exp_result$part <- p
      }
    }

    if (length(exp_result) > 0) {
      # If only DEA requested, return it directly (not nested)
      if (dea_only) {
        results[[exp_name]] <- exp_result$dea
      } else {
        results[[exp_name]] <- exp_result
      }
    }
  }

  if (verbose) {
    if (strip_version) cat("Ensembl version suffixes stripped\n")
    cat("\n--- Summary ---\n")
    cat("Experiments loaded:", length(results), "/", length(all_dirs), "\n")
  }

  return(results)
}


#' Load single DEA result (RDS or CSV)
#' @noRd
.load_dea_result <- function(exp_dir, prefix, verbose) {

  # Try RDS first
  rds_file <- file.path(exp_dir, paste0(prefix, "_dea_result.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("  DEA: loaded from RDS\n")
    return(readRDS(rds_file))
  }

  # Fall back to CSV reconstruction
  csv_file <- file.path(exp_dir, paste0(prefix, "_all_results.csv"))
  if (!file.exists(csv_file)) {
    if (verbose) cat("  DEA: not found (no RDS or CSV)\n")
    return(NULL)
  }

  if (verbose) cat("  DEA: reconstructed from CSV\n")

  results_df <- read.csv(csv_file, stringsAsFactors = FALSE)

  # Read metadata for comparison info
  comparison <- NA
  summary_stats <- list()
  meta_file <- file.path(exp_dir, paste0(prefix, "_metadata.txt"))
  if (file.exists(meta_file)) {
    meta_lines <- readLines(meta_file, warn = FALSE)

    # Extract comparison
    comp_line <- grep("^  .+ vs .+", meta_lines, value = TRUE)
    if (length(comp_line) > 0) {
      comparison <- trimws(comp_line[1])
    }

    # Extract summary counts
    total_line <- grep("Total features tested:", meta_lines, value = TRUE)
    up_line <- grep("Significant \\(up\\):", meta_lines, value = TRUE)
    down_line <- grep("Significant \\(down\\):", meta_lines, value = TRUE)
    sig_line <- grep("Total significant:", meta_lines, value = TRUE)

    if (length(total_line) > 0) {
      summary_stats$total <- as.integer(gsub("\\D", "", total_line[1]))
    }
    if (length(up_line) > 0) {
      summary_stats$up <- as.integer(gsub("\\D", "", up_line[1]))
    }
    if (length(down_line) > 0) {
      summary_stats$down <- as.integer(gsub("\\D", "", down_line[1]))
    }
    if (length(sig_line) > 0) {
      summary_stats$significant <- as.integer(gsub("\\D", "", sig_line[1]))
    }
  }

  list(
    results = results_df,
    dds = NULL,
    summary = summary_stats,
    comparison = comparison,
    norm_counts = NULL,
    .reconstructed = TRUE
  )
}


#' Load single PART result (RDS or CSV)
#' @noRd
.load_part_result <- function(exp_dir, prefix, verbose) {

  # Try RDS
  rds_file <- file.path(exp_dir, paste0(prefix, "_part_result.rds"))
  if (file.exists(rds_file)) {
    if (verbose) cat("  PART: loaded from RDS\n")
    return(readRDS(rds_file))
  }

  # CSV reconstruction (limited — no dendrogram or full data matrix)
  clusters_file <- file.path(exp_dir, paste0(prefix, "_part_clusters.csv"))
  if (!file.exists(clusters_file)) {
    if (verbose) cat("  PART: not found\n")
    return(NULL)
  }

  if (verbose) cat("  PART: reconstructed from CSV (limited)\n")

  cluster_df <- read.csv(clusters_file, stringsAsFactors = FALSE)

  # Detect gene and cluster columns
  gene_col <- intersect(c("gene", "feature_id", "gene_id"), colnames(cluster_df))[1]
  clust_col <- intersect(c("cluster"), colnames(cluster_df))[1]

  if (is.na(gene_col) || is.na(clust_col)) {
    if (verbose) cat("  PART: could not identify gene/cluster columns\n")
    return(NULL)
  }

  clusters <- setNames(cluster_df[[clust_col]], cluster_df[[gene_col]])
  cluster_names <- sort(unique(clusters))

  # Generate colors
  n_colors <- length(cluster_names)
  set3_base <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                 "#FDB462", "#B3DE69", "#FCCDE5")
  cols <- grDevices::colorRampPalette(set3_base)(n_colors)
  names(cols) <- cluster_names
  if ("C0" %in% cluster_names) cols["C0"] <- "#D9D9D9"

  # Build cluster_map
  cluster_map <- data.frame(
    gene = names(clusters),
    cluster = unname(clusters),
    cluster_color = cols[clusters],
    stringsAsFactors = FALSE,
    row.names = names(clusters)
  )

  cluster_sizes <- table(clusters)
  cluster_sizes <- cluster_sizes[order(names(cluster_sizes))]

  result <- list(
    clusters = clusters,
    cluster_map = cluster_map,
    n_clusters = sum(!cluster_names %in% "C0"),
    n_outliers = sum(clusters == "C0"),
    cluster_sizes = setNames(as.integer(cluster_sizes), names(cluster_sizes)),
    cluster_colors = cols,
    dendrogram = NULL,
    data = NULL,
    .reconstructed = TRUE
  )

  class(result) <- c("artemis_part", "list")
  return(result)
}


# ==============================================================================
# COLOR UTILITIES
# ==============================================================================

#' Generate distinct visible colors
#'
#' Generates a set of distinct colors suitable for plotting, filtering out
#' near-white colors that would be invisible on light backgrounds. Useful for
#' remapping WGCNA module colors or any categorical color assignment.
#'
#' @param n Integer. Number of colors needed.
#' @param palette Character. HCL color palette name. Default: "Dark 3".
#'   See \code{grDevices::hcl.pals("qualitative")} for options.
#' @param max_brightness Numeric (0-255). Maximum mean RGB brightness allowed.
#'   Colors brighter than this are excluded. Default: 200.
#'
#' @return Character vector of \code{n} hex color codes.
#'
#' @examples
#' \dontrun{
#' cols <- DEMETER_generate_colors(10)
#' cols <- DEMETER_generate_colors(5, palette = "Set 2")
#'
#' # Remap WGCNA module colors for circos
#' mod_names <- rownames(circos_data$module_df)
#' mod_names <- mod_names[mod_names != "dummy"]
#' new_colors <- setNames(DEMETER_generate_colors(length(mod_names)), mod_names)
#' identity_rename <- setNames(mod_names, mod_names)
#' circos_data <- DEMETER_rename_circos_modules(circos_data, identity_rename, new_colors)
#'
#' }
#' @export
DEMETER_generate_colors <- function(n, palette = "Dark 3", max_brightness = 200) {
  # Over-generate to have replacements after filtering
  pool_size <- max(n * 3, 50)
  pool <- grDevices::hcl.colors(pool_size, palette = palette)

  # Filter out near-white colors
  rgb_vals <- grDevices::col2rgb(pool)
  brightness <- colMeans(rgb_vals)
  pool <- pool[brightness <= max_brightness]

  if (length(pool) < n) {
    warning("Only ", length(pool), " colors passed brightness filter (requested ", n,
            "). Try a darker palette or increase max_brightness.")
    pool <- c(pool, rep(pool, length.out = n - length(pool)))
  }

  pool[seq_len(n)]
}


# ==============================================================================
# PART-WGCNA CIRCOS DATA PREPARATION
# ==============================================================================

#' Prepare data for PART-WGCNA circos plot
#'
#' Combines PART clustering results, WGCNA module assignments, trait correlations,
#' and DEA log2 fold-change values into a structured object for circos visualization.
#'
#' @param part_result An \code{artemis_part} object from \code{ARTEMIS_part_clustering()}.
#' @param wgcna_modules A \code{wgcna_modules} object from WGCNA pipeline or
#'   \code{DEMETER_load_wgcna()}.
#' @param wgcna_trait_cor A \code{wgcna_trait_cor} object from
#'   \code{ARTEMIS_wgcna_module_traits()} or \code{DEMETER_load_wgcna()}.
#'   NULL if no trait correlations available.
#' @param dea_results Named list of DEA results. Each element should be either:
#'   - A result from \code{ARTEMIS_differential_counts()} (has \code{$results})
#'   - A data.frame with feature_id and log2FoldChange columns
#'   Names become column labels in the circos (e.g., list(C5RO = dea_c5ro, CS1AN = dea_cs1an)).
#'   NULL if no DEA results to display.
#' @param strip_version Logical. Strip version suffix from Ensembl gene IDs
#'   (e.g., ENSG00000103888.18 -> ENSG00000103888) for matching between
#'   PART and WGCNA. Default: TRUE.
#' @param p_thresh Numeric. P-value threshold for trait significance in module_df.
#'   Default: 0.05.
#' @param sig_type Character. Type of significance to use for filtering:
#'   "padj" (adjusted p-value, default) or "pvalue" (raw p-value).
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{part_wgcna_circos} containing:
#'   \describe{
#'     \item{PART_df}{Data.frame with genes as rows, columns: cluster, module,
#'       and one L2FC column per DEA comparison}
#'     \item{module_df}{Data.frame with modules as rows, columns: module, size,
#'       trait correlations, and significance flags}
#'     \item{association_matrix}{Integer matrix (clusters x modules) of gene counts}
#'     \item{cluster_colors}{Named character vector of cluster colors}
#'     \item{module_colors}{Named character vector of module colors}
#'     \item{module_gene}{Data.frame with gene and module columns}
#'   }
#'
#' @examples
#' \dontrun{
#' # Prepare circos data from GAIA objects
#' circos_data <- DEMETER_prepare_part_wgcna_circos(
#'   part_result = part_res,
#'   wgcna_modules = wgcna$modules,
#'   wgcna_trait_cor = wgcna$trait_cor,
#'   dea_results = list(C5RO = dea_c5ro, CS1AN = dea_cs1an)
#' )
#'
#' # Filter and plot
#' circos_data <- DEMETER_filter_circos_modules(circos_data, c("blue", "brown", "turquoise"))
#' AETHER_plot_part_wgcna_circos(circos_data, output_file = "circos.pdf")
#'
#' }
#' @export
DEMETER_prepare_part_wgcna_circos <- function(part_result,
                                               wgcna_modules,
                                               wgcna_trait_cor = NULL,
                                               dea_results = NULL,
                                               strip_version = TRUE,
                                               p_thresh = 0.05,
                                               sig_type = c("padj", "pvalue"),
                                               verbose = TRUE) {

  # --- Validate inputs ---
  if (!inherits(part_result, "artemis_part")) {
    stop("'part_result' must be an artemis_part object")
  }
  if (!inherits(wgcna_modules, "wgcna_modules")) {
    stop("'wgcna_modules' must be a wgcna_modules object")
  }

  if (verbose) cat("=== Preparing PART-WGCNA Circos Data ===\n")

  # --- Extract PART info ---
  cluster_map <- part_result$cluster_map
  part_genes <- cluster_map$gene
  cluster_colors <- part_result$cluster_colors

  # --- Extract WGCNA module-gene mapping ---
  if (!is.null(wgcna_modules$gene_module_df)) {
    module_gene <- wgcna_modules$gene_module_df
    # Normalize column name to "module"
    if (!"module" %in% colnames(module_gene) && "module_name" %in% colnames(module_gene)) {
      module_gene$module <- module_gene$module_name
    }
  } else if (!is.null(wgcna_modules$module_names)) {
    module_gene <- data.frame(
      gene = names(wgcna_modules$module_names),
      module = unname(wgcna_modules$module_names),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Cannot extract gene-module mapping from wgcna_modules")
  }
  rownames(module_gene) <- module_gene$gene

  # --- Gene ID matching ---
  if (strip_version) {
    part_genes_stripped <- sub("\\.\\d+$", "", part_genes)

    # Check if WGCNA genes also have versions
    wgcna_genes <- module_gene$gene
    wgcna_has_version <- any(grepl("\\.\\d+$", wgcna_genes))
    if (wgcna_has_version) {
      wgcna_genes_stripped <- sub("\\.\\d+$", "", wgcna_genes)
    } else {
      wgcna_genes_stripped <- wgcna_genes
    }

    # Build lookup: stripped WGCNA gene -> module
    wgcna_lookup <- setNames(module_gene$module, wgcna_genes_stripped)

    # Map PART genes to modules via stripped IDs
    part_modules <- wgcna_lookup[part_genes_stripped]
    names(part_modules) <- part_genes
  } else {
    part_modules <- setNames(
      module_gene[part_genes, "module"],
      part_genes
    )
  }

  # --- Build PART_df ---
  PART_df <- data.frame(
    cluster = cluster_map$cluster,
    stringsAsFactors = FALSE,
    row.names = part_genes
  )

  # Add DEA L2FC columns
  if (!is.null(dea_results)) {
    if (!is.list(dea_results)) stop("'dea_results' must be a named list")

    for (name in names(dea_results)) {
      dea <- dea_results[[name]]

      # Extract results data.frame
      if (is.data.frame(dea)) {
        res_df <- dea
      } else if (!is.null(dea$results)) {
        res_df <- dea$results
      } else {
        warning("DEA result '", name, "' has no recognizable format, skipping")
        next
      }

      # Detect gene ID column
      gene_col <- intersect(c("feature_id", "gene_id", "gene", "ensembl_id"),
                            colnames(res_df))[1]
      if (is.na(gene_col)) {
        warning("No gene ID column found in DEA result '", name, "', skipping")
        next
      }

      # Build lookup: gene -> l2fc
      lfc_lookup <- setNames(res_df$log2FoldChange, res_df[[gene_col]])

      # Match to PART genes
      PART_df[[name]] <- lfc_lookup[part_genes]
    }
  }

  # Add module assignment (after DEA columns so it's last)
  PART_df$module <- part_modules

  if (verbose) {
    n_matched <- sum(!is.na(PART_df$module))
    cat("Gene matching:\n")
    cat("  PART genes:", length(part_genes), "\n")
    cat("  WGCNA genes:", nrow(module_gene), "\n")
    cat("  Matched:", n_matched, "(", round(100 * n_matched / length(part_genes), 1), "%)\n")
    cat("  Unmatched:", sum(is.na(PART_df$module)), "\n")
  }

  # --- Build module_df ---
  if (!is.null(wgcna_trait_cor)) {
    sig_type <- match.arg(sig_type)
    cor_mat <- wgcna_trait_cor$cor_matrix
    if (sig_type == "pvalue") {
      pval_mat <- wgcna_trait_cor$pvalue_matrix
    } else {
      pval_mat <- if (!is.null(wgcna_trait_cor$padj_matrix)) {
        wgcna_trait_cor$padj_matrix
      } else {
        wgcna_trait_cor$pvalue_matrix
      }
    }

    module_names <- rownames(cor_mat)

    # Strip ME prefix if present
    clean_names <- sub("^ME", "", module_names)
    rownames(cor_mat) <- clean_names
    if (!is.null(pval_mat)) rownames(pval_mat) <- clean_names
    module_names <- clean_names

    module_df <- as.data.frame(cor_mat)

    if (!is.null(pval_mat)) {
      sig_df <- as.data.frame(pval_mat < p_thresh)
      colnames(sig_df) <- paste0(colnames(sig_df), "_sig")
      module_df <- cbind(module_df, sig_df)
    }

    module_df$module <- module_names
    rownames(module_df) <- module_names

    # Add module sizes
    module_sizes <- table(module_gene$module)
    module_df$size <- as.integer(module_sizes[module_names])
    module_df$size[is.na(module_df$size)] <- 0
  } else {
    # No trait correlations - just sizes
    module_names <- sort(unique(module_gene$module))
    module_sizes <- table(module_gene$module)
    module_df <- data.frame(
      module = module_names,
      size = as.integer(module_sizes[module_names]),
      row.names = module_names,
      stringsAsFactors = FALSE
    )
  }

  # --- Build association matrix ---
  cluster_names <- sort(unique(PART_df$cluster))
  all_modules <- unique(na.omit(PART_df$module))

  association_matrix <- matrix(0,
                               nrow = length(cluster_names),
                               ncol = length(all_modules),
                               dimnames = list(cluster_names, all_modules))

  valid_rows <- !is.na(PART_df$module)
  for (i in which(valid_rows)) {
    cl <- PART_df$cluster[i]
    mod <- PART_df$module[i]
    association_matrix[cl, mod] <- association_matrix[cl, mod] + 1
  }

  # --- Module colors ---
  # Use color map from wgcna_modules if available (module_name -> hex)
  module_color_vec <- character(length(all_modules))
  names(module_color_vec) <- all_modules

  if (!is.null(wgcna_modules$module_colors) && !is.null(names(wgcna_modules$module_colors))) {
    # New format: module_colors is a named vector (module_name -> hex color)
    for (m in all_modules) {
      if (m %in% names(wgcna_modules$module_colors)) {
        module_color_vec[m] <- wgcna_modules$module_colors[m]
      } else {
        module_color_vec[m] <- NA
      }
    }
  }

  # Fill any missing with generated colors
  needs_color <- is.na(module_color_vec) | module_color_vec == ""
  if (any(needs_color)) {
    n_needed <- sum(needs_color)
    generated <- DEMETER_generate_colors(n_needed)
    module_color_vec[needs_color] <- generated
  }

  # --- Return ---
  result <- list(
    PART_df = PART_df,
    module_df = module_df,
    association_matrix = association_matrix,
    cluster_colors = cluster_colors,
    module_colors = module_color_vec,
    module_gene = module_gene
  )

  class(result) <- c("part_wgcna_circos", "list")

  if (verbose) {
    sample_cols <- setdiff(colnames(PART_df), c("cluster", "module"))
    cat("\nCircos data ready:\n")
    cat("  Clusters:", length(cluster_names), "\n")
    cat("  Modules:", length(all_modules), "\n")
    cat("  DEA comparisons:", length(sample_cols), "\n")
    if (!is.null(wgcna_trait_cor)) {
      trait_cols <- colnames(cor_mat)[!grepl("_sig$", colnames(cor_mat))]
      cat("  Traits:", length(trait_cols), "\n")
    }
  }

  return(result)
}


#' @method print part_wgcna_circos
#' @export
print.part_wgcna_circos <- function(x, ...) {
  sample_cols <- setdiff(colnames(x$PART_df), c("cluster", "module"))
  trait_cols <- colnames(x$module_df)
  trait_cols <- trait_cols[!trait_cols %in% c("module", "size") & !grepl("_sig$", trait_cols)]

  cat("PART-WGCNA Circos Data\n")
  cat("----------------------\n")
  cat("Genes:", nrow(x$PART_df), "\n")
  cat("Clusters:", length(unique(x$PART_df$cluster)), "\n")
  cat("Modules:", nrow(x$module_df), "\n")
  cat("DEA comparisons:", length(sample_cols), "\n")
  cat("Traits:", length(trait_cols), "\n")
  invisible(x)
}


#' Select WGCNA modules for circos display
#'
#' Identifies modules worth displaying in the circos plot based on multiple
#' criteria: trait correlation significance and strength, module size, overlap
#' with PART clusters, and (optionally) enrichment results. Returns module
#' names suitable for passing to \code{DEMETER_filter_circos_modules()}.
#'
#' All enabled criteria are combined with AND logic: a module must pass every
#' applicable filter. Criteria are skipped when their data is unavailable
#' (e.g., enrichment filtering is skipped if \code{enrichment = NULL}).
#'
#' @param circos_data A \code{part_wgcna_circos} object.
#' @param enrichment Optional. Enrichment results from \code{APOLLO_enrich_gost()}
#'   or similar. Must have a \code{module} column. Modules with at least one
#'   significant term pass the filter.
#' @param min_cor Numeric. Minimum absolute trait correlation with at least one
#'   trait. Default: 0.3. Set to 0 to disable.
#' @param require_sig_cor Logical. Require at least one significant trait
#'   correlation (based on the _sig columns in module_df). Default: TRUE.
#'   Ignored if no trait data is available.
#' @param min_size Integer. Minimum module size (gene count). Default: 30.
#' @param min_overlap Integer. Minimum number of genes overlapping with any
#'   PART cluster (from the association matrix). Default: 1.
#' @param exclude_grey Logical. Exclude the grey (unassigned) module. Default: TRUE.
#' @param verbose Logical. Print filtering summary. Default: TRUE.
#'
#' @return Character vector of module names passing all criteria.
#'
#' @examples
#' \dontrun{
#' selected <- DEMETER_select_circos_modules(circos_data, min_cor = 0.4)
#' circos_data <- DEMETER_filter_circos_modules(circos_data, selected)
#'
#' }
#' @export
DEMETER_select_circos_modules <- function(circos_data,
                                           enrichment = NULL,
                                           min_cor = 0.3,
                                           require_sig_cor = TRUE,
                                           min_size = 30,
                                           min_overlap = 1,
                                           exclude_grey = TRUE,
                                           verbose = TRUE) {

  if (!inherits(circos_data, "part_wgcna_circos")) {
    stop("'circos_data' must be a part_wgcna_circos object")
  }

  module_df <- circos_data$module_df
  assoc <- circos_data$association_matrix
  all_modules <- rownames(module_df)
  keep <- rep(TRUE, length(all_modules))
  names(keep) <- all_modules

  if (verbose) cat("=== Selecting Circos Modules ===\n")
  if (verbose) cat("Starting modules:", length(all_modules), "\n")

  # --- 1. Exclude grey / dummy ---
  if (exclude_grey) {
    grey_modules <- intersect(all_modules, c("grey", "dummy", "module_0"))
    keep[grey_modules] <- FALSE
    if (verbose && length(grey_modules) > 0) {
      cat("Excluded grey/dummy:", length(grey_modules), "->",
          sum(keep), "remaining\n")
    }
  }

  # --- 2. Module size ---
  if (min_size > 0 && "size" %in% colnames(module_df)) {
    too_small <- all_modules[module_df$size < min_size]
    keep[too_small] <- FALSE
    if (verbose) cat("Min size (", min_size, "):", sum(!keep[too_small]),
                     "removed ->", sum(keep), "remaining\n", sep = "")
  }

  # --- 3. PART overlap ---
  if (min_overlap > 0 && !is.null(assoc)) {
    # Total overlap per module = column sums of association matrix
    modules_in_assoc <- intersect(all_modules, colnames(assoc))
    module_overlap <- colSums(assoc[, modules_in_assoc, drop = FALSE])
    no_overlap <- names(module_overlap[module_overlap < min_overlap])
    # Also flag modules not in the association matrix at all
    not_in_assoc <- setdiff(all_modules, colnames(assoc))
    no_overlap <- union(no_overlap, not_in_assoc)
    keep[no_overlap] <- FALSE
    if (verbose) cat("Min overlap (", min_overlap, "): ",
                     length(no_overlap), " removed -> ",
                     sum(keep), " remaining\n", sep = "")
  }

  # --- 4. Trait correlation ---
  # Detect trait columns (not _sig, not module, not size)
  trait_cols <- colnames(module_df)[!colnames(module_df) %in% c("module", "size") &
                                     !grepl("_sig$", colnames(module_df))]
  sig_cols <- paste0(trait_cols, "_sig")
  sig_cols <- sig_cols[sig_cols %in% colnames(module_df)]
  has_traits <- length(trait_cols) > 0

  if (has_traits) {
    # Correlation strength
    if (min_cor > 0) {
      max_abs_cor <- apply(abs(module_df[, trait_cols, drop = FALSE]), 1, max, na.rm = TRUE)
      weak_cor <- all_modules[max_abs_cor < min_cor]
      keep[weak_cor] <- FALSE
      if (verbose) cat("Min |cor| (", min_cor, "): ",
                       length(weak_cor), " removed -> ",
                       sum(keep), " remaining\n", sep = "")
    }

    # Significance
    if (require_sig_cor && length(sig_cols) > 0) {
      any_sig <- apply(module_df[, sig_cols, drop = FALSE], 1, any, na.rm = TRUE)
      not_sig <- all_modules[!any_sig]
      keep[not_sig] <- FALSE
      if (verbose) cat("Significant correlation: ",
                       length(not_sig), " removed -> ",
                       sum(keep), " remaining\n", sep = "")
    }
  } else if (verbose) {
    cat("No trait data available, skipping correlation filters\n")
  }

  # --- 5. Enrichment ---
  if (!is.null(enrichment)) {
    if ("module" %in% colnames(enrichment)) {
      enriched_modules <- unique(enrichment$module)
      not_enriched <- setdiff(all_modules[keep], enriched_modules)
      keep[not_enriched] <- FALSE
      if (verbose) cat("Has enrichment: ",
                       length(not_enriched), " removed -> ",
                       sum(keep), " remaining\n", sep = "")
    } else {
      warning("enrichment data has no 'module' column, skipping enrichment filter")
    }
  }

  selected <- names(keep[keep])

  if (verbose) {
    cat("---\nSelected:", length(selected), "of", length(all_modules), "modules\n")
    if (length(selected) > 0) cat("Modules:", paste(selected, collapse = ", "), "\n")
  }

  selected
}


#' Filter circos data to specific WGCNA modules
#'
#' Keeps selected modules and assigns all other genes to a "dummy" module.
#' The dummy module appears in the circos but is visually muted, allowing
#' focus on modules of interest while preserving all PART cluster data.
#'
#' @param circos_data A \code{part_wgcna_circos} object.
#' @param selected_modules Character vector of module names to keep.
#' @param dummy_name Character. Name for the catch-all module. Default: "dummy".
#' @param dummy_color Character. Color for the dummy module. Default: "#CCCCCC".
#'
#' @return Updated \code{part_wgcna_circos} object.
#'
#' @export
DEMETER_filter_circos_modules <- function(circos_data,
                                           selected_modules,
                                           dummy_name = "dummy",
                                           dummy_color = "#CCCCCC") {

  if (!inherits(circos_data, "part_wgcna_circos")) {
    stop("'circos_data' must be a part_wgcna_circos object")
  }

  PART_df <- circos_data$PART_df
  module_df <- circos_data$module_df
  module_gene <- circos_data$module_gene

  # Reassign non-selected modules to dummy
  PART_df$module <- ifelse(
    PART_df$module %in% selected_modules,
    PART_df$module,
    dummy_name
  )
  # Also catch NAs
  PART_df$module[is.na(PART_df$module)] <- dummy_name

  # Filter module_df to selected + dummy
  module_df <- module_df[module_df$module %in% selected_modules, , drop = FALSE]

  # Create dummy row
  dummy_row <- module_df[1, , drop = FALSE]
  dummy_row$module <- dummy_name
  dummy_row$size <- sum(PART_df$module == dummy_name)

  # Zero out trait values for dummy
  trait_cols <- colnames(module_df)
  trait_cols <- trait_cols[!trait_cols %in% c("module", "size") & !grepl("_sig$", trait_cols)]
  sig_cols <- grep("_sig$", colnames(module_df), value = TRUE)
  for (col in trait_cols) dummy_row[[col]] <- 0
  for (col in sig_cols) dummy_row[[col]] <- FALSE

  rownames(dummy_row) <- dummy_name
  module_df <- rbind(module_df, dummy_row)

  # Filter module_gene to selected only
  module_gene <- module_gene[module_gene$module %in% selected_modules, , drop = FALSE]

  # Rebuild association matrix
  cluster_names <- sort(unique(PART_df$cluster))
  module_names <- c(selected_modules, dummy_name)

  association_matrix <- matrix(0,
                               nrow = length(cluster_names),
                               ncol = length(module_names),
                               dimnames = list(cluster_names, module_names))

  for (i in 1:nrow(PART_df)) {
    cl <- PART_df$cluster[i]
    mod <- PART_df$module[i]
    association_matrix[cl, mod] <- association_matrix[cl, mod] + 1
  }

  # Update colors
  module_colors <- circos_data$module_colors
  module_colors <- module_colors[names(module_colors) %in% selected_modules]
  module_colors[dummy_name] <- dummy_color

  circos_data$PART_df <- PART_df
  circos_data$module_df <- module_df
  circos_data$association_matrix <- association_matrix
  circos_data$module_colors <- module_colors
  circos_data$module_gene <- module_gene

  return(circos_data)
}


#' Filter circos data to specific PART clusters
#'
#' @param circos_data A \code{part_wgcna_circos} object.
#' @param selected_clusters Character vector of cluster names to keep.
#'
#' @return Updated \code{part_wgcna_circos} object.
#'
#' @export
DEMETER_filter_circos_clusters <- function(circos_data, selected_clusters) {

  if (!inherits(circos_data, "part_wgcna_circos")) {
    stop("'circos_data' must be a part_wgcna_circos object")
  }

  PART_df <- circos_data$PART_df
  PART_df <- PART_df[PART_df$cluster %in% selected_clusters, , drop = FALSE]

  # Rebuild association matrix
  module_names <- colnames(circos_data$association_matrix)
  association_matrix <- matrix(0,
                               nrow = length(selected_clusters),
                               ncol = length(module_names),
                               dimnames = list(selected_clusters, module_names))

  for (i in 1:nrow(PART_df)) {
    cl <- PART_df$cluster[i]
    mod <- PART_df$module[i]
    if (!is.na(mod) && mod %in% module_names) {
      association_matrix[cl, mod] <- association_matrix[cl, mod] + 1
    }
  }

  # Filter cluster colors
  cluster_colors <- circos_data$cluster_colors
  cluster_colors <- cluster_colors[names(cluster_colors) %in% selected_clusters]

  circos_data$PART_df <- PART_df
  circos_data$association_matrix <- association_matrix
  circos_data$cluster_colors <- cluster_colors

  return(circos_data)
}


#' Rename and recolor WGCNA modules in circos data
#'
#' @param circos_data A \code{part_wgcna_circos} object.
#' @param rename_vector Named character vector: old names as names, new names as values.
#' @param recolor_vector Named character vector: new names as names, hex colors as values.
#'   If NULL, colors are auto-generated for renamed modules.
#'
#' @return Updated \code{part_wgcna_circos} object.
#'
#' @examples
#' \dontrun{
#' # Rename WGCNA color-based names to generic labels
#' rename_vec <- c(blue = "module_1", brown = "module_2", turquoise = "module_3",
#'                 dummy = "dummy")
#' recolor_vec <- c(module_1 = "dodgerblue2", module_2 = "#E31A1C",
#'                  module_3 = "green4", dummy = "#CCCCCC")
#' circos_data <- DEMETER_rename_circos_modules(circos_data, rename_vec, recolor_vec)
#'
#' }
#' @export
DEMETER_rename_circos_modules <- function(circos_data,
                                           rename_vector,
                                           recolor_vector = NULL) {

  if (!inherits(circos_data, "part_wgcna_circos")) {
    stop("'circos_data' must be a part_wgcna_circos object")
  }

  # Rename in PART_df
  circos_data$PART_df$module <- rename_vector[circos_data$PART_df$module]

  # Rename in module_df
  rownames(circos_data$module_df) <- rename_vector[rownames(circos_data$module_df)]
  circos_data$module_df$module <- rename_vector[circos_data$module_df$module]

  # Rename in association matrix
  colnames(circos_data$association_matrix) <- rename_vector[colnames(circos_data$association_matrix)]

  # Rename in module_gene
  circos_data$module_gene$module <- rename_vector[circos_data$module_gene$module]

  # Update module colors
  if (!is.null(recolor_vector)) {
    circos_data$module_colors <- recolor_vector
  } else {
    new_names <- rename_vector[names(circos_data$module_colors)]
    names(circos_data$module_colors) <- new_names
  }

  return(circos_data)
}
