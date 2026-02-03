# GAIA/Apollo/gene_enrichment.R
# Gene list enrichment functions
#
# Functions for enrichment analysis on gene lists (e.g., from WGCNA modules).
# Complements peak_annotation.R which handles peak-based enrichment.
#
# Uses gprofiler2 for enrichment - accepts multiple ID types (ENSEMBL, SYMBOL, etc.)
# and includes GO, KEGG, Reactome, WikiPathways in one call.

library(gprofiler2)


#' Convert gene IDs with fallback
#'
#' Converts gene IDs from one type to another (e.g., ENSEMBL to SYMBOL),
#' preserving the original ID (optionally stripped of version numbers) when
#' conversion fails. Useful for converting count matrix rownames while
#' keeping all genes.
#'
#' @param genes Character vector of gene IDs (e.g., rownames of a count matrix).
#' @param org_db OrgDb object for ID conversion (e.g., org.Hs.eg.db).
#' @param from_type Type of input gene IDs. Default: "ENSEMBL".
#' @param to_type Type of output gene IDs. Default: "SYMBOL".
#' @param strip_version Logical. Remove version numbers from Ensembl-style IDs
#'   before conversion (e.g., ENSG00000141510.16 -> ENSG00000141510). Default: TRUE.
#' @param verbose Logical. Print conversion statistics. Default: TRUE.
#'
#' @return Character vector of converted IDs, same length as input.
#'   Genes that couldn't be converted retain their original ID (stripped if requested).
#'   Has attribute "conversion_stats" with success/failure counts.
#'
#' @details
#' Unlike APOLLO_prepare_genelist() which drops unmapped genes (suitable for
#' enrichment analysis), this function preserves all genes - useful when you
#' need to maintain the same number of rows in a count matrix.
#'
#' When duplicates occur after conversion (multiple ENSEMBL IDs mapping to the
#' same SYMBOL), a suffix is added to make names unique (e.g., "TP53", "TP53.1").
#'
#' @examples
#' library(org.Hs.eg.db)
#'
#' # Convert ENSEMBL to SYMBOL for count matrix rownames
#' ensembl_ids <- c("ENSG00000141510.16", "ENSG00000012048.23", "ENSG00000000003.15")
#' symbols <- APOLLO_convert_genes(ensembl_ids, org.Hs.eg.db)
#'
#' # Apply to count matrix
#' rownames(counts) <- APOLLO_convert_genes(rownames(counts), org.Hs.eg.db)
#'
#' @export
APOLLO_convert_genes <- function(genes,
                                  org_db,
                                  from_type = "ENSEMBL",
                                  to_type = "SYMBOL",
                                  strip_version = TRUE,
                                  verbose = TRUE) {

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required. Install from Bioconductor.")
  }

  genes_original <- genes
  n_input <- length(genes)

  # Strip version numbers if requested
  if (strip_version) {
    genes_stripped <- sub("\\.[0-9]+$", "", genes)
  } else {
    genes_stripped <- genes
  }

  # Convert IDs
  converted <- tryCatch({
    AnnotationDbi::mapIds(
      org_db,
      keys = genes_stripped,
      keytype = from_type,
      column = to_type,
      multiVals = "first"
    )
  }, error = function(e) {
    warning("Conversion failed: ", e$message)
    setNames(rep(NA, length(genes_stripped)), genes_stripped)
  })

  # Count successes
  n_converted <- sum(!is.na(converted))
  n_failed <- n_input - n_converted

  # Replace NA with original (stripped) ID
  result <- converted
  result[is.na(result)] <- genes_stripped[is.na(result)]

  # Handle duplicates by making names unique
  if (any(duplicated(result))) {
    n_dups <- sum(duplicated(result))
    result <- make.unique(result, sep = ".")
    if (verbose) {
      cat("Note:", n_dups, "duplicate names made unique with suffix\n")
    }
  }

  if (verbose) {
    cat("Converted:", n_converted, "/", n_input, "genes (",
        round(100 * n_converted / n_input, 1), "%)\n", sep = "")
    if (n_failed > 0) {
      cat("Preserved original ID for", n_failed, "unmapped genes\n")
    }
  }

  # Add stats as attribute
  attr(result, "conversion_stats") <- list(
    input = n_input,
    converted = n_converted,
    preserved = n_failed
  )

  return(as.character(result))
}


#' Prepare gene list for enrichment analysis
#'
#' Converts gene symbols to ENTREZID format required by most enrichment tools.
#' Can work with a single gene list or a named list of gene lists (e.g., modules).
#'
#' @param genes Character vector of gene symbols, or a named list of character
#'   vectors (e.g., from WGCNA modules).
#' @param org_db OrgDb object for ID conversion (e.g., org.Hs.eg.db).
#'   If NULL, skips conversion and just cleans IDs (useful with strip_version).
#' @param from_type Type of input gene IDs. Default: "SYMBOL".
#' @param to_type Type of output gene IDs. Default: "ENTREZID".
#' @param strip_version Logical. Remove version numbers from Ensembl-style IDs
#'   (e.g., ENSG00000141510.16 -> ENSG00000141510). Default: FALSE.
#' @param drop_na Logical. Remove genes that couldn't be converted. Default: TRUE.
#' @param verbose Logical. Print conversion statistics. Default: TRUE.
#'
#' @return If genes is a vector, returns converted vector.
#'   If genes is a list, returns named list with same structure.
#'   Attribute "conversion_stats" contains success/failure counts.
#'
#' @examples
#' library(org.Hs.eg.db)
#' genes <- c("TP53", "BRCA1", "EGFR")
#' entrez <- APOLLO_prepare_genelist(genes, org.Hs.eg.db)
#'
#' # From WGCNA modules
#' module_genes <- list(blue = c("TP53", "BRCA1"), red = c("EGFR", "MYC"))
#' module_entrez <- APOLLO_prepare_genelist(module_genes, org.Hs.eg.db)
#'
#' # With Ensembl IDs that have version numbers
#' ensembl_genes <- c("ENSG00000141510.16", "ENSG00000012048.23")
#' clean_genes <- APOLLO_prepare_genelist(ensembl_genes, org_db = NULL, strip_version = TRUE)
#'
#' @export
APOLLO_prepare_genelist <- function(genes,
                                     org_db = NULL,
                                     from_type = "SYMBOL",
                                     to_type = "ENTREZID",
                                     strip_version = FALSE,
                                     drop_na = TRUE,
                                     verbose = TRUE) {

  # Handle list input (e.g., module gene lists)
  if (is.list(genes) && !is.data.frame(genes)) {
    results <- lapply(names(genes), function(name) {
      if (verbose) cat("Converting:", name, "\n")
      APOLLO_prepare_genelist(
        genes[[name]],
        org_db = org_db,
        from_type = from_type,
        to_type = to_type,
        strip_version = strip_version,
        drop_na = drop_na,
        verbose = FALSE
      )
    })
    names(results) <- names(genes)

    if (verbose) {
      total_in <- sum(sapply(genes, length))
      total_out <- sum(sapply(results, length))
      cat("\nTotal: ", total_out, "/", total_in, " genes processed (",
          round(100 * total_out / total_in, 1), "%)\n", sep = "")
    }

    return(results)
  }

  # Single vector processing

genes <- unique(as.character(genes))
  n_input <- length(genes)

  # Strip version numbers if requested (e.g., ENSG00000141510.16 -> ENSG00000141510)
  if (strip_version) {
    genes <- sub("\\.[0-9]+$", "", genes)
    genes <- unique(genes)  # May have duplicates after stripping
    if (verbose && length(genes) < n_input) {
      cat("Stripped version numbers, ", n_input - length(genes),
          " duplicates removed\n", sep = "")
    }
    n_input <- length(genes)
  }

  # If no org_db, just return cleaned genes
  if (is.null(org_db)) {
    if (verbose) cat("No conversion (org_db = NULL), returning ", n_input, " genes\n", sep = "")
    return(genes)
  }

  # Convert IDs
  converted <- tryCatch({
    AnnotationDbi::mapIds(
      org_db,
      keys = genes,
      keytype = from_type,
      column = to_type,
      multiVals = "first"
    )
  }, error = function(e) {
    warning("Conversion failed: ", e$message)
    setNames(rep(NA, length(genes)), genes)
  })

  n_success <- sum(!is.na(converted))

  if (verbose) {
    cat("Converted: ", n_success, "/", n_input, " genes (",
        round(100 * n_success / n_input, 1), "%)\n", sep = "")
  }

  if (drop_na) {
    converted <- converted[!is.na(converted)]
  }

  attr(converted, "conversion_stats") <- list(
    input = n_input,
    converted = n_success,
    failed = n_input - n_success
  )

  return(as.character(converted))
}


#' Extract gene lists from WGCNA modules
#'
#' Convenience function to extract gene lists from a wgcna_modules object,
#' optionally converting to ENTREZID.
#'
#' @param modules A wgcna_modules object from ARTEMIS_wgcna_detect_modules().
#' @param modules_of_interest Character vector of module colors to extract.
#'   Default: NULL (all modules except grey).
#' @param exclude_grey Logical. Exclude the grey (unassigned) module. Default: TRUE.
#' @param strip_version Logical. Remove version numbers from Ensembl-style IDs
#'   (e.g., ENSG00000141510.16 -> ENSG00000141510). Default: FALSE.
#' @param org_db Optional OrgDb for ID conversion. If provided, converts to ENTREZID.
#' @param from_type Type of input gene IDs (for conversion). Default: "ENSEMBL".
#' @param to_type Type of output gene IDs (for conversion). Default: "ENTREZID".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Named list of gene vectors, one per module.
#'
#' @examples
#' # Get gene symbols (no conversion)
#' module_genes <- APOLLO_extract_module_genes(modules)
#'
#' # Strip Ensembl version numbers for gprofiler
#' module_genes <- APOLLO_extract_module_genes(modules, strip_version = TRUE)
#'
#' # Get ENTREZID from Ensembl IDs
#' library(org.Hs.eg.db)
#' module_entrez <- APOLLO_extract_module_genes(modules, strip_version = TRUE,
#'                                               org_db = org.Hs.eg.db,
#'                                               from_type = "ENSEMBL")
#'
#' @export
APOLLO_extract_module_genes <- function(modules,
                                         modules_of_interest = NULL,
                                         exclude_grey = TRUE,
                                         strip_version = FALSE,
                                         org_db = NULL,
                                         from_type = "ENSEMBL",
                                         to_type = "ENTREZID",
                                         verbose = TRUE) {

  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules object")
  }

  module_colors <- modules$module_colors

  # Determine modules to extract
  all_modules <- unique(module_colors)
  if (exclude_grey) {
    all_modules <- all_modules[all_modules != "grey"]
  }

  if (!is.null(modules_of_interest)) {
    # Validate
    invalid <- setdiff(modules_of_interest, all_modules)
    if (length(invalid) > 0) {
      warning("Modules not found: ", paste(invalid, collapse = ", "))
    }
    all_modules <- intersect(modules_of_interest, all_modules)
  }

  if (verbose) {
    cat("Extracting genes from", length(all_modules), "modules\n")
  }

  # Extract gene lists
  gene_lists <- lapply(all_modules, function(mod) {
    names(module_colors)[module_colors == mod]
  })
  names(gene_lists) <- all_modules

  if (verbose) {
    sizes <- sapply(gene_lists, length)
    cat("Module sizes:", paste(names(sizes), sizes, sep = "=", collapse = ", "), "\n")
  }

  # Strip version numbers if requested (before conversion)
  if (strip_version) {
    if (verbose) cat("Stripping version numbers from gene IDs...\n")
    gene_lists <- lapply(gene_lists, function(genes) {
      unique(sub("\\.[0-9]+$", "", genes))
    })
  }

  # Convert IDs if org_db provided
  if (!is.null(org_db)) {
    if (verbose) cat("\nConverting", from_type, "to", to_type, "...\n")
    gene_lists <- APOLLO_prepare_genelist(
      gene_lists,
      org_db = org_db,
      from_type = from_type,
      to_type = to_type,
      strip_version = FALSE,  # Already done above
      verbose = verbose
    )
  }

  return(gene_lists)
}


#' Run enrichment analysis using gprofiler2
#'
#' Performs functional enrichment analysis on gene lists using gprofiler2::gost().
#' Supports GO, KEGG, Reactome, WikiPathways, and more in a single call.
#' Accepts multiple ID types directly (ENSEMBL, SYMBOL, etc.) - no conversion needed.
#'
#' @param gene_lists Named list of gene vectors. Can be ENSEMBL IDs (with or without
#'   version numbers), gene symbols, or other identifiers.
#' @param organism gprofiler2 organism code. Common values:
#'   "hsapiens" (human), "mmusculus" (mouse), "rnorvegicus" (rat).
#'   Default: "hsapiens".
#' @param sources Character vector of data sources to query. Options include:
#'   "GO:BP", "GO:MF", "GO:CC" (Gene Ontology),
#'   "KEGG" (KEGG pathways),
#'   "REAC" (Reactome),
#'   "WP" (WikiPathways),
#'   "TF" (TRANSFAC transcription factors),
#'   "MIRNA" (miRTarBase miRNA targets),
#'   "CORUM" (protein complexes),
#'   "HP" (Human Phenotype Ontology).
#'   Default: c("GO:BP", "KEGG", "REAC").
#' @param user_threshold Significance threshold for term filtering. Default: 0.05.
#' @param correction_method Multiple testing correction method:
#'   "g_SCS" (gprofiler's default), "fdr", "bonferroni". Default: "g_SCS".
#' @param domain_scope Background for statistical test:
#'   "annotated" (genes with any annotation),
#'   "known" (all known genes),
#'   "custom" (provide custom_bg). Default: "annotated".
#' @param custom_bg Custom background gene set (when domain_scope = "custom").
#' @param min_term_size Minimum term/pathway size. Default: 10.
#' @param max_term_size Maximum term/pathway size. Default: 500.
#' @param significant Logical. Only return significant results. Default: TRUE.
#' @param exclude_iea Logical. Exclude GO terms inferred from electronic annotation. Default: FALSE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list with class "gost_enrichment" containing:
#'   \item{results}{Named list of gost result objects, one per gene list}
#'   \item{combined}{Combined data.frame of all results with 'query' column indicating source}
#'   \item{summary}{Summary data.frame with counts per gene list and source}
#'   \item{metadata}{List of analysis parameters}
#'
#' @details
#' gprofiler2 automatically detects the input ID type in most cases. For best results
#' with Ensembl IDs that have version numbers (e.g., ENSG00000141510.16), use
#' APOLLO_extract_module_genes() with strip_version = TRUE before calling this function.
#'
#' @examples
#' # Basic usage with gene symbols
#' gene_lists <- list(
#'   module1 = c("TP53", "BRCA1", "EGFR", "MYC"),
#'   module2 = c("IL6", "TNF", "IL1B", "CXCL8")
#' )
#' results <- APOLLO_enrich_gost(gene_lists)
#'
#' # From WGCNA modules (Ensembl IDs)
#' module_genes <- APOLLO_extract_module_genes(modules, strip_version = TRUE)
#' results <- APOLLO_enrich_gost(module_genes, sources = c("GO:BP", "GO:MF", "KEGG", "REAC"))
#'
#' # Mouse data
#' results <- APOLLO_enrich_gost(gene_lists, organism = "mmusculus")
#'
#' @export
APOLLO_enrich_gost <- function(gene_lists,
                                organism = "hsapiens",
                                sources = c("GO:BP", "KEGG", "REAC"),
                                user_threshold = 0.05,
                                correction_method = "g_SCS",
                                domain_scope = "annotated",
                                custom_bg = NULL,
                                min_term_size = 10,
                                max_term_size = 500,
                                significant = TRUE,
                                exclude_iea = FALSE,
                                verbose = TRUE) {

  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required. Install with: install.packages('gprofiler2')")
  }

  if (verbose) {
    cat("=== Functional Enrichment Analysis (gprofiler2) ===\n")
    cat("Organism:", organism, "\n")
    cat("Sources:", paste(sources, collapse = ", "), "\n")
    cat("Gene lists:", length(gene_lists), "\n")
    cat("Threshold:", user_threshold, "(", correction_method, ")\n\n")
  }

  results <- list()
  combined_list <- list()

  for (name in names(gene_lists)) {
    genes <- gene_lists[[name]]

    if (length(genes) < 3) {
      if (verbose) cat(name, ": Skipping (< 3 genes)\n")
      next
    }

    if (verbose) cat(name, " (", length(genes), " genes): ", sep = "")

    gost_result <- tryCatch({
      gprofiler2::gost(
        query = genes,
        organism = organism,
        sources = sources,
        user_threshold = user_threshold,
        correction_method = correction_method,
        domain_scope = domain_scope,
        custom_bg = custom_bg,
        evcodes = TRUE,  # Include gene IDs in results
        significant = significant,
        exclude_iea = exclude_iea
      )
    }, error = function(e) {
      if (verbose) cat("Error - ", e$message, "\n")
      NULL
    })

    if (!is.null(gost_result) && !is.null(gost_result$result) && nrow(gost_result$result) > 0) {
      # Filter by term size
      result_df <- gost_result$result
      result_df <- result_df[result_df$term_size >= min_term_size &
                             result_df$term_size <= max_term_size, ]

      if (nrow(result_df) > 0) {
        n_sig <- nrow(result_df)
        if (verbose) cat(n_sig, "significant terms\n")

        # Store full gost object
        results[[name]] <- gost_result

        # Add module identifier and combine
        result_df$module <- name
        combined_list[[name]] <- result_df
      } else {
        if (verbose) cat("No terms within size range\n")
      }
    } else {
      if (verbose) cat("No significant terms\n")
    }
  }

  # Combine results
  if (length(combined_list) > 0) {
    combined <- do.call(rbind, combined_list)
    rownames(combined) <- NULL
  } else {
    combined <- data.frame()
  }

  # Create summary by module and source
  summary_rows <- list()
  for (name in names(gene_lists)) {
    for (src in sources) {
      if (nrow(combined) > 0) {
        n_terms <- sum(combined$module == name & combined$source == src)
      } else {
        n_terms <- 0
      }
      summary_rows[[paste(name, src, sep = "_")]] <- data.frame(
        module = name,
        source = src,
        n_genes = length(gene_lists[[name]]),
        n_terms = n_terms,
        stringsAsFactors = FALSE
      )
    }
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  if (verbose) {
    cat("\n--- Summary ---\n")
    cat("Modules with results:", length(results), "/", length(gene_lists), "\n")
    cat("Total significant terms:", nrow(combined), "\n")

    # Breakdown by source
    if (nrow(combined) > 0) {
      cat("\nBy source:\n")
      for (src in sources) {
        cat("  ", src, ": ", sum(combined$source == src), "\n", sep = "")
      }
    }
  }

  result <- list(
    results = results,
    combined = combined,
    summary = summary_df,
    metadata = list(
      organism = organism,
      sources = sources,
      user_threshold = user_threshold,
      correction_method = correction_method,
      min_term_size = min_term_size,
      max_term_size = max_term_size
    )
  )
  class(result) <- c("gost_enrichment", "list")

  return(result)
}


#' Filter enrichment results by source or module
#'
#' Extract specific subsets from a gost_enrichment result object.
#'
#' @param enrich_result A gost_enrichment object from APOLLO_enrich_gost().
#' @param sources Character vector of sources to keep (e.g., c("GO:BP", "KEGG")).
#'   Default: NULL (all sources).
#' @param modules Character vector of modules to keep. Default: NULL (all modules).
#' @param top_n Integer. Keep only top N terms per module per source. Default: NULL (all).
#'
#' @return Filtered data.frame of enrichment results.
#'
#' @examples
#' # Get only KEGG results from blue module
#' kegg_blue <- APOLLO_filter_enrichment(results, sources = "KEGG", modules = "blue")
#'
#' # Get top 10 GO:BP terms per module
#' top_bp <- APOLLO_filter_enrichment(results, sources = "GO:BP", top_n = 10)
#'
#' @export
APOLLO_filter_enrichment <- function(enrich_result,
                                      sources = NULL,
                                      modules = NULL,
                                      top_n = NULL) {

  if (!inherits(enrich_result, "gost_enrichment")) {
    stop("enrich_result must be a gost_enrichment object from APOLLO_enrich_gost()")
  }

  df <- enrich_result$combined

  if (nrow(df) == 0) {
    return(df)
  }

  # Filter by source
  if (!is.null(sources)) {
    df <- df[df$source %in% sources, ]
  }

  # Filter by module
  if (!is.null(modules)) {
    df <- df[df$module %in% modules, ]
  }

  # Keep top N per module per source
  if (!is.null(top_n) && nrow(df) > 0) {
    df <- df[order(df$p_value), ]
    df <- do.call(rbind, lapply(split(df, list(df$module, df$source)), function(x) {
      if (nrow(x) > 0) head(x, top_n) else x
    }))
    rownames(df) <- NULL
  }

  return(df)
}
