# GAIA/Apollo/gene_enrichment.R
# Gene list enrichment functions
#
# Functions for enrichment analysis on gene lists (e.g., from WGCNA modules).
# Complements peak_annotation.R which handles peak-based enrichment.
#
# Uses gprofiler2 for enrichment - accepts multiple ID types (ENSEMBL, SYMBOL, etc.)
# and includes GO, KEGG, Reactome, WikiPathways in one call.



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
#' \dontrun{
#' library(org.Hs.eg.db)
#'
#' # Convert ENSEMBL to SYMBOL for count matrix rownames
#' ensembl_ids <- c("ENSG00000141510.16", "ENSG00000012048.23", "ENSG00000000003.15")
#' symbols <- APOLLO_convert_genes(ensembl_ids, org.Hs.eg.db)
#'
#' # Apply to count matrix
#' rownames(counts) <- APOLLO_convert_genes(rownames(counts), org.Hs.eg.db)
#'
#' }
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
      cat("[APOLLO] Note:", n_dups, "duplicate names made unique with suffix\n")
    }
  }

  if (verbose) {
    cat("[APOLLO] Converted:", n_converted, "/", n_input, "genes (",
        round(100 * n_converted / n_input, 1), "%)\n")
    if (n_failed > 0) {
      cat("   Preserved original ID for", n_failed, "unmapped genes\n")
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
#' \dontrun{
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
#' }
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
      if (verbose) cat("[APOLLO] Converting:", name, "\n")
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
      cat("    Total:", total_out, "/", total_in, "genes processed (",
          round(100 * total_out / total_in, 1), "%)\n")
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
      cat("[APOLLO] Stripped version numbers,", n_input - length(genes),
          "duplicates removed\n")
    }
    n_input <- length(genes)
  }

  # If no org_db, just return cleaned genes
  if (is.null(org_db)) {
    if (verbose) cat("[APOLLO] No conversion (org_db = NULL), returning ", n_input, " genes\n", sep = "")
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
    cat("[APOLLO] Converted: ", n_success, "/", n_input, " genes (",
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


#' Extract gene lists from clustering results
#'
#' Extracts gene lists from clustering results, returning a named list suitable
#' for enrichment analysis. Accepts \code{wgcna_modules}, \code{wgcna_hubs}
#' (from ARTEMIS_wgcna_*), or \code{artemis_part} objects. When passed a
#' wgcna_hubs object (from ARTEMIS_wgcna_hub_genes() with n_top = NULL and
#' a kME threshold), only hub genes passing the threshold are returned —
#' suitable for enrichment analysis on large modules.
#'
#' @param modules A \code{wgcna_modules} object from \code{ARTEMIS_wgcna_detect_modules()},
#'   a \code{wgcna_hubs} object from \code{ARTEMIS_wgcna_hub_genes()}, or an
#'   \code{artemis_part} object from \code{ARTEMIS_part()}.
#' @param modules_of_interest Character vector of module/cluster names to extract.
#'   Default: NULL (all clusters, excluding unassigned).
#' @param exclude_grey Logical. Exclude unassigned genes: the grey module for
#'   WGCNA (\code{module_0}), or the outlier cluster for PART (\code{C0}).
#'   Default: TRUE.
#' @param strip_version Logical. Remove version numbers from Ensembl-style IDs
#'   (e.g., ENSG00000141510.16 -> ENSG00000141510). Default: FALSE.
#' @param org_db Optional OrgDb for ID conversion. If provided, converts to ENTREZID.
#' @param from_type Type of input gene IDs (for conversion). Default: "ENSEMBL".
#' @param to_type Type of output gene IDs (for conversion). Default: "ENTREZID".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Named list of gene vectors, one per cluster/module.
#'
#' @examples
#' \dontrun{
#' # From WGCNA modules
#' cluster_genes <- APOLLO_extract_cluster_genes(modules)
#'
#' # From PART clustering (C0 outliers excluded by default)
#' part_result <- ARTEMIS_part(mat, seed = 42)
#' cluster_genes <- APOLLO_extract_cluster_genes(part_result)
#'
#' # Strip Ensembl version numbers for gprofiler
#' cluster_genes <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE)
#'
#' # Get ENTREZID from Ensembl IDs
#' library(org.Hs.eg.db)
#' cluster_entrez <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE,
#'                                                 org_db = org.Hs.eg.db,
#'                                                 from_type = "ENSEMBL")
#'
#' }
#' @export
APOLLO_extract_cluster_genes <- function(modules,
                                          modules_of_interest = NULL,
                                          exclude_grey = TRUE,
                                         strip_version = FALSE,
                                         org_db = NULL,
                                         from_type = "ENSEMBL",
                                         to_type = "ENTREZID",
                                         verbose = TRUE) {

  # --- wgcna_hubs dispatch ---
  if (inherits(modules, "wgcna_hubs")) {
    hub_df <- modules$hub_genes  # already filtered by is_hub / mm_threshold

    if (nrow(hub_df) == 0) {
      stop("wgcna_hubs object contains no hub genes (hub_genes is empty). ",
           "Consider lowering mm_threshold or setting n_top = NULL.")
    }

    gene_lists <- split(hub_df$gene, hub_df$module)

    if (!is.null(modules_of_interest)) {
      invalid <- setdiff(modules_of_interest, names(gene_lists))
      if (length(invalid) > 0) {
        warning("Modules not found in hub genes: ", paste(invalid, collapse = ", "))
      }
      gene_lists <- gene_lists[intersect(modules_of_interest, names(gene_lists))]
    }

    if (verbose) {
      sizes <- sapply(gene_lists, length)
      cat("[APOLLO] Hub genes per module (MM threshold =",
          modules$criteria$mm_threshold, "):\n")
      cat("[APOLLO]", paste(names(sizes), sizes, sep = "=", collapse = ", "), "\n")
    }

    if (strip_version) {
      gene_lists <- lapply(gene_lists, function(genes) unique(sub("\\.[0-9]+$", "", genes)))
    }

    if (!is.null(org_db)) {
      gene_lists <- APOLLO_prepare_genelist(
        gene_lists,
        org_db = org_db,
        from_type = from_type,
        to_type = to_type,
        strip_version = FALSE,
        verbose = verbose
      )
    }

    return(gene_lists)
  }

  # --- artemis_part dispatch ---
  if (inherits(modules, "artemis_part")) {
    cluster_vec <- modules$clusters
    gene_lists  <- split(names(cluster_vec), cluster_vec)

    if (exclude_grey) gene_lists[["C0"]] <- NULL

    if (!is.null(modules_of_interest)) {
      invalid <- setdiff(modules_of_interest, names(gene_lists))
      if (length(invalid) > 0) {
        warning("Clusters not found: ", paste(invalid, collapse = ", "))
      }
      gene_lists <- gene_lists[intersect(modules_of_interest, names(gene_lists))]
    }

    if (verbose) {
      sizes <- sapply(gene_lists, length)
      cat("[APOLLO] Genes per cluster:", paste(names(sizes), sizes, sep = "=", collapse = ", "), "\n")
      if (modules$n_outliers > 0 && exclude_grey) {
        cat("   (", modules$n_outliers, "outlier genes in C0 excluded)\n")
      }
    }

    if (strip_version) {
      gene_lists <- lapply(gene_lists, function(genes) unique(sub("\\.[0-9]+$", "", genes)))
    }

    if (!is.null(org_db)) {
      gene_lists <- APOLLO_prepare_genelist(
        gene_lists,
        org_db       = org_db,
        from_type    = from_type,
        to_type      = to_type,
        strip_version = FALSE,
        verbose      = verbose
      )
    }

    return(gene_lists)
  }

  if (!inherits(modules, "wgcna_modules")) {
    stop("modules must be a wgcna_modules, wgcna_hubs, or artemis_part object")
  }

  module_names_vec <- modules$module_names

  # Determine modules to extract
  all_modules <- unique(module_names_vec)
  if (exclude_grey) {
    all_modules <- all_modules[all_modules != "module_0"]
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
    cat("[APOLLO] Extracting genes from", length(all_modules), "modules\n")
  }

  # Extract gene lists
  gene_lists <- lapply(all_modules, function(mod) {
    names(module_names_vec)[module_names_vec == mod]
  })
  names(gene_lists) <- all_modules

  if (verbose) {
    sizes <- sapply(gene_lists, length)
    cat("[APOLLO] Module sizes:", paste(names(sizes), sizes, sep = "=", collapse = ", "), "\n")
  }

  # Strip version numbers if requested (before conversion)
  if (strip_version) {
    if (verbose) cat("[APOLLO] Stripping version numbers from gene IDs...\n")
    gene_lists <- lapply(gene_lists, function(genes) {
      unique(sub("\\.[0-9]+$", "", genes))
    })
  }

  # Convert IDs if org_db provided
  if (!is.null(org_db)) {
    if (verbose) cat("[APOLLO] Converting", from_type, "to", to_type, "...\n")
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
#' @param min_term_size Minimum term/pathway size. Default: 1 (no filtering,
#'   matches gprofiler online defaults).
#' @param max_term_size Maximum term/pathway size. Default: 1e6 (no effective
#'   upper limit, matches gprofiler online defaults).
#' @param max_query_size Maximum number of genes per query. Gene lists exceeding
#'   this are skipped. Prevents slow, uninformative enrichment on very large
#'   gene sets (e.g., WGCNA grey module). Default: 10000. Set to NULL to disable.
#' @param significant Logical. Only include significant results in
#'   \code{combined}/\code{summary}. All terms are always fetched from
#'   gprofiler2 in a single API call regardless of this setting (filtering is
#'   applied client-side), so toggling this does not incur extra requests.
#'   Default: TRUE.
#' @param exclude_iea Logical. Exclude GO terms inferred from electronic annotation. Default: FALSE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list with class "gost_enrichment" containing:
#'   \item{results}{Named list of FULL (unfiltered) gost result objects, one
#'     per gene list — includes every evaluated term regardless of
#'     significance, for use with \code{AETHER_plot_gost_full()}}
#'   \item{combined}{Combined data.frame of \code{significant}-filtered results with 'query' column indicating source}
#'   \item{summary}{Summary data.frame with counts per gene list and source}
#'   \item{metadata}{List of analysis parameters}
#'
#' @details
#' gprofiler2 automatically detects the input ID type in most cases. For best results
#' with Ensembl IDs that have version numbers (e.g., ENSG00000141510.16), use
#' APOLLO_extract_cluster_genes() with strip_version = TRUE before calling this function.
#'
#' @examples
#' \dontrun{
#' # Basic usage with gene symbols
#' gene_lists <- list(
#'   module1 = c("TP53", "BRCA1", "EGFR", "MYC"),
#'   module2 = c("IL6", "TNF", "IL1B", "CXCL8")
#' )
#' results <- APOLLO_enrich_gost(gene_lists)
#'
#' # From WGCNA modules (Ensembl IDs)
#' module_genes <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE)
#' results <- APOLLO_enrich_gost(module_genes, sources = c("GO:BP", "GO:MF", "KEGG", "REAC"))
#'
#' # Mouse data
#' results <- APOLLO_enrich_gost(gene_lists, organism = "mmusculus")
#'
#' }
#' @export
APOLLO_enrich_gost <- function(gene_lists,
                                organism = "hsapiens",
                                sources = c("GO:BP", "KEGG", "REAC"),
                                user_threshold = 0.05,
                                correction_method = "g_SCS",
                                domain_scope = "annotated",
                                custom_bg = NULL,
                                min_term_size = 1,
                                max_term_size = 1e6,
                                max_query_size = 10000,
                                significant = TRUE,
                                exclude_iea = FALSE,
                                verbose = TRUE) {

  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required. Install with: install.packages('gprofiler2')")
  }

  if (verbose) {
    cat("[APOLLO] Functional Enrichment Analysis (gprofiler2) \n")
    cat("    Organism:", organism, "\n")
    cat("    Sources:", paste(sources, collapse = ", "), "\n")
    cat("    Gene lists:", length(gene_lists), "\n")
    cat("    Threshold:", user_threshold, "(", correction_method, ")\n\n")
  }

  results <- list()
  combined_list <- list()

  for (name in names(gene_lists)) {
    genes <- gene_lists[[name]]

    if (length(genes) < 3) {
      if (verbose) cat("    ", name, ": Skipping (< 3 genes)\n")
      next
    }

    if (!is.null(max_query_size) && length(genes) > max_query_size) {
      if (verbose) cat("    ", name, " (", length(genes), " genes): Skipping (exceeds max_query_size = ",
                       max_query_size, ")\n", sep = "")
      next
    }

    if (verbose) cat("    ", name, " (", length(genes), " genes): ", sep = "")

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
        # Always fetch ALL evaluated terms (not just significant ones): the
        # `significant` column gprofiler2 computes is filtered on client-side
        # below, and the full object is kept in $results for
        # AETHER_plot_gost_full(). Avoids a second API call to get unfiltered
        # terms for that plot.
        significant = FALSE,
        exclude_iea = exclude_iea
      )
    }, error = function(e) {
      if (verbose) cat("     Error - ", e$message, "\n")
      NULL
    })

    if (!is.null(gost_result) && !is.null(gost_result$result) && nrow(gost_result$result) > 0) {
      # Store the FULL (unfiltered) gost object regardless of significance —
      # used by AETHER_plot_gost_full() to plot every evaluated term.
      results[[name]] <- gost_result

      # Filter by term size, then (optionally) by significance for combined/summary
      result_df <- gost_result$result
      result_df <- result_df[result_df$term_size >= min_term_size &
                             result_df$term_size <= max_term_size, ]

      if (significant) {
        if ("significant" %in% colnames(result_df)) {
          result_df <- result_df[result_df$significant, ]
        } else {
          result_df <- result_df[!is.na(result_df$p_value) &
                                   result_df$p_value <= user_threshold, ]
        }
      }

      if (nrow(result_df) > 0) {
        n_sig <- nrow(result_df)
        if (verbose) cat(n_sig, if (significant) "significant terms\n" else "terms (unfiltered)\n")

        # Add module identifier and combine
        result_df$module <- name
        combined_list[[name]] <- result_df
      } else {
        if (verbose) cat("No", if (significant) "significant terms\n" else "terms within size range\n")
      }
    } else {
      if (verbose) cat("No terms returned\n")
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
    cat("[APOLLO] --- Summary ---\n")
    cat("    Modules with results:", length(results), "/", length(gene_lists), "\n")
    cat("    Total significant terms:", nrow(combined), "\n")

    # Breakdown by source
    if (nrow(combined) > 0) {
      cat("    By source:\n")
      for (src in sources) {
        cat("    ", src, ": ", sum(combined$source == src), "\n", sep = "")
      }
    }
  }

  # Store original gene list sizes for annotation coverage calculation
  module_sizes <- sapply(gene_lists, length)

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
      max_term_size = max_term_size,
      module_sizes = module_sizes
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
#' \dontrun{
#' # Get only KEGG results from blue module
#' kegg_blue <- APOLLO_filter_enrichment(results, sources = "KEGG", modules = "blue")
#'
#' # Get top 10 GO:BP terms per module
#' top_bp <- APOLLO_filter_enrichment(results, sources = "GO:BP", top_n = 10)
#'
#' }
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


###############################################################################
########### Gene Set Enrichment Analysis (GSEA, fgsea) ###########
###############################################################################

#' Build a ranked gene vector for GSEA from a differential expression result
#'
#' Extracts a named numeric vector (gene -> ranking statistic) from an
#' \code{artemis_limma} or \code{artemis_ts_de} result, suitable for
#' \code{APOLLO_gsea()}.
#'
#' @param de_result An \code{artemis_limma} object (from
#'   \code{ARTEMIS_limma_de()}), an \code{artemis_ts_de} object (from
#'   \code{ARTEMIS_timeseries_conditional()}/\code{_temporal()} or their
#'   \code{_limma} counterparts), or a plain data.frame with the relevant
#'   columns.
#' @param rank_by Character. Ranking statistic:
#'   \describe{
#'     \item{\code{"t"}}{The test statistic - limma's \code{t} column or
#'       DESeq2's \code{stat} column (whichever is present is used
#'       automatically). Accounts for both effect size and precision; the
#'       recommended default.}
#'     \item{\code{"log2FoldChange"}}{Fold change only - ignores variance.}
#'     \item{\code{"signed_neg_log10p"}}{\code{sign(log2FoldChange) *
#'       -log10(pvalue)}. Useful when no test statistic column is available.}
#'   }
#'   Default: \code{"t"}.
#' @param gene_col Character. Column in the results data.frame to use as gene
#'   names. Default: \code{"feature_id"}. Use \code{"original_id"} for Olink
#'   data when the gene sets are expected to match the matrix's original IDs
#'   (e.g. OlinkIDs) rather than the substituted \code{feature_id_col} symbols.
#' @param comparison Character, integer, or \code{NULL}. For
#'   \code{artemis_ts_de} objects (which contain multiple comparisons): the
#'   comparison name or index to extract. \code{NULL} uses the first
#'   comparison and prints a note. Ignored for \code{artemis_limma} and
#'   data.frame inputs.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Named numeric vector, sorted descending, NA values removed,
#'   duplicate gene names removed (first occurrence kept).
#'
#' @examples
#' \dontrun{
#' de <- ARTEMIS_limma_de(ol$wide, ol$sample_meta, group_col = "Group",
#'                         reference = "Control", experiment = "Sepsis")
#' ranked <- APOLLO_rank_from_de(de, rank_by = "t")
#'
#' # From a timeseries result, specific comparison
#' ranked_tp2 <- APOLLO_rank_from_de(temp_de, rank_by = "t",
#'                                    comparison = "2_vs_1")
#' }
#' @export
APOLLO_rank_from_de <- function(de_result,
                                 rank_by    = c("t", "log2FoldChange", "signed_neg_log10p"),
                                 gene_col   = "feature_id",
                                 comparison = NULL,
                                 verbose    = TRUE) {

  rank_by <- match.arg(rank_by)

  # ---------------------------------------------------------------------------
  # Resolve the results data.frame
  # ---------------------------------------------------------------------------
  if (inherits(de_result, "artemis_limma")) {
    df <- de_result$results

  } else if (inherits(de_result, "artemis_ts_de")) {
    comp_names <- names(de_result$results)
    if (is.null(comparison)) {
      comparison <- comp_names[1]
      if (verbose) {
        cat("[APOLLO] No comparison specified; using '", comparison, "'\n", sep = "")
      }
    } else if (is.numeric(comparison)) {
      comparison <- comp_names[comparison]
    }
    if (!comparison %in% comp_names) {
      stop("comparison '", comparison, "' not found. Available: ",
           paste(comp_names, collapse = ", "))
    }
    df <- de_result$results[[comparison]]$results

  } else if (is.data.frame(de_result)) {
    df <- de_result

  } else {
    stop("'de_result' must be an artemis_limma object, an artemis_ts_de ",
         "object, or a data.frame.")
  }

  if (!gene_col %in% colnames(df)) {
    stop("gene_col '", gene_col, "' not found in results. Available columns: ",
         paste(colnames(df), collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Build the ranking statistic
  # ---------------------------------------------------------------------------
  if (rank_by == "t") {
    stat_col <- if ("t" %in% colnames(df)) "t" else if ("stat" %in% colnames(df)) "stat" else NULL
    if (is.null(stat_col)) {
      stop("rank_by='t' requires a 't' (limma) or 'stat' (DESeq2) column ",
           "in the results. Use rank_by='log2FoldChange' or ",
           "'signed_neg_log10p' instead.")
    }
    stat_vec <- df[[stat_col]]
    if (verbose && stat_col == "stat") {
      cat("[APOLLO] Using DESeq2 'stat' column as the ranking statistic ",
          "(equivalent role to limma 't')\n", sep = "")
    }

  } else if (rank_by == "log2FoldChange") {
    if (!"log2FoldChange" %in% colnames(df)) {
      stop("rank_by='log2FoldChange' requires a 'log2FoldChange' column.")
    }
    stat_vec <- df$log2FoldChange

  } else {
    if (!all(c("log2FoldChange", "pvalue") %in% colnames(df))) {
      stop("rank_by='signed_neg_log10p' requires 'log2FoldChange' and ",
           "'pvalue' columns.")
    }
    stat_vec <- sign(df$log2FoldChange) * -log10(df$pvalue)
  }

  ranked <- stats::setNames(stat_vec, as.character(df[[gene_col]]))

  # ---------------------------------------------------------------------------
  # Clean: drop NA, deduplicate, sort descending
  # ---------------------------------------------------------------------------
  ranked <- ranked[!is.na(ranked)]
  if (anyDuplicated(names(ranked))) {
    n_dup  <- sum(duplicated(names(ranked)))
    ranked <- ranked[!duplicated(names(ranked))]
    if (verbose) cat("[APOLLO] Removed", n_dup, "duplicate gene names (kept first occurrence)\n")
  }
  ranked <- sort(ranked, decreasing = TRUE)

  if (verbose) {
    cat("[APOLLO] Ranked", length(ranked), "genes by '", rank_by, "'\n", sep = "")
  }

  return(ranked)
}


#' Gene Set Enrichment Analysis (fgsea)
#'
#' Runs \code{fgsea::fgseaMultilevel()} on a pre-ranked gene list against a
#' collection of gene sets. Gene sets can be supplied as a named list or
#' fetched automatically from MSigDB via \code{msigdbr}. Use
#' \code{APOLLO_rank_from_de()} to build \code{ranked_genes} from an
#' \code{artemis_limma} or \code{artemis_ts_de} result.
#'
#' @param ranked_genes Named numeric vector. Names are gene symbols (or IDs
#'   matching the gene sets); values are the ranking statistic (e.g. limma's
#'   \code{t}, DESeq2's \code{stat}, or signed -log10(p)). Must have names.
#' @param gene_sets Named list of character vectors (gene set name -> member
#'   genes), or \code{NULL} to fetch from MSigDB via \code{msigdbr}.
#' @param species Character. Species name passed to \code{msigdbr::msigdbr()}
#'   when \code{gene_sets = NULL}. Default: \code{"Homo sapiens"}.
#' @param collection Character vector. MSigDB collection codes to fetch when
#'   \code{gene_sets = NULL}. Default: \code{c("H", "C2", "C5")}. A
#'   subcategory can be appended with a colon (e.g. \code{"C2:CP:REACTOME"},
#'   \code{"C5:GO:BP"}).
#' @param min_size Integer. Minimum gene set size (after overlap with ranked
#'   genes). Default: 15.
#' @param max_size Integer. Maximum gene set size. Default: 500.
#' @param fdr_threshold Numeric. FDR threshold for \code{$significant}.
#'   Default: 0.05.
#' @param verbose Logical. Print progress and summary. Default: TRUE.
#'
#' @return An S3 object of class \code{"apollo_gsea"} containing:
#'   \describe{
#'     \item{results}{Full fgsea result data.frame (ordered by padj):
#'       pathway, pval, padj, ES, NES, size, leadingEdge (collapsed string),
#'       collection (if gene sets came from msigdbr).}
#'     \item{significant}{\code{results} filtered to \code{fdr_threshold}.}
#'     \item{ranked_genes}{The input vector, NA-removed and sorted descending.}
#'     \item{gene_sets}{Named list of only the gene sets actually tested
#'       (passed min/max size) - used by \code{AETHER_plot_gsea_enrichment()}.}
#'     \item{params}{List of parameters used.}
#'   }
#'
#' @details
#' Requires the \code{fgsea} package (Suggests). When \code{gene_sets = NULL},
#' also requires \code{msigdbr} (Suggests).
#'
#' @examples
#' \dontrun{
#' ranked <- APOLLO_rank_from_de(de, rank_by = "t")
#' gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))
#'
#' # Custom gene sets
#' my_sets <- list(my_pathway = c("TP53", "BRCA1", "EGFR"))
#' gsea_result <- APOLLO_gsea(ranked, gene_sets = my_sets)
#' }
#' @export
APOLLO_gsea <- function(ranked_genes,
                         gene_sets     = NULL,
                         species       = "Homo sapiens",
                         collection    = c("H", "C2", "C5"),
                         min_size      = 15L,
                         max_size      = 500L,
                         fdr_threshold = 0.05,
                         verbose       = TRUE) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
  }

  if (is.null(names(ranked_genes)) || length(ranked_genes) == 0L) {
    stop("'ranked_genes' must be a named numeric vector (names = gene IDs).")
  }

  # Remove NA, deduplicate names (keep first occurrence), sort descending
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  if (anyDuplicated(names(ranked_genes))) {
    n_dup <- sum(duplicated(names(ranked_genes)))
    ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
    if (verbose) cat("[APOLLO] Removed", n_dup, "duplicate gene names (kept first occurrence)\n")
  }
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  if (verbose) {
    cat("[APOLLO] Gene Set Enrichment Analysis (fgsea)\n")
    cat("    Ranked genes:", length(ranked_genes), "\n")
  }

  # ---------------------------------------------------------------------------
  # Resolve gene sets
  # ---------------------------------------------------------------------------
  pathway_to_col <- NULL   # populated below when gene_sets come from msigdbr
  if (is.null(gene_sets)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Package 'msigdbr' is required when gene_sets = NULL. ",
           "Install with: install.packages('msigdbr')")
    }

    if (verbose) {
      cat("    Gene sets   : MSigDB (", paste(collection, collapse = " + "),
          ") for ", species, "\n", sep = "")
    }

    gs_dfs <- lapply(collection, function(col) {
      # msigdbr subcollection codes can themselves contain a colon (e.g.
      # "GO:BP", "CP:REACTOME"), so everything after the first ":" must be
      # rejoined rather than truncated to the second token only -- taking
      # just cat_sub[2] silently fetched the wrong (broader) gene set for
      # "C2:CP:REACTOME" (all of C2:CP, not just REACTOME) and an invalid
      # subcategory for "C5:GO:BP" ("GO" instead of "GO:BP").
      cat_sub  <- strsplit(col, ":", fixed = TRUE)[[1]]
      cat_code <- cat_sub[1]
      subcat   <- if (length(cat_sub) > 1) paste(cat_sub[-1], collapse = ":") else NULL

      tryCatch(
        msigdbr::msigdbr(species = species, category = cat_code, subcategory = subcat),
        error = function(e) {
          warning("msigdbr failed for collection '", col, "': ", conditionMessage(e))
          NULL
        }
      )
    })
    gs_dfs <- Filter(Negate(is.null), gs_dfs)

    if (length(gs_dfs) == 0L) {
      stop("No gene sets retrieved from MSigDB. Check species name and collection codes.")
    }

    gs_df     <- do.call(rbind, gs_dfs)
    gene_sets <- split(gs_df$gene_symbol, gs_df$gs_name)

    unique_gs      <- gs_df[!duplicated(gs_df$gs_name), ]
    pathway_to_col <- setNames(unique_gs$gs_collection, unique_gs$gs_name)
  }

  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("'gene_sets' must be a named list of character vectors.")
  }

  if (verbose) {
    cat("    Testing", length(gene_sets), "gene sets against",
        format(length(ranked_genes), big.mark = ","), "ranked genes...\n\n")
  }

  res <- fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats    = ranked_genes,
    minSize  = as.integer(min_size),
    maxSize  = as.integer(max_size),
    eps      = 0
  )

  res_df <- as.data.frame(res)
  res_df$leadingEdge <- vapply(
    res$leadingEdge,
    function(x) paste(x, collapse = ", "),
    character(1)
  )
  res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
  rownames(res_df) <- NULL

  if (!is.null(pathway_to_col)) {
    res_df$collection <- pathway_to_col[res_df$pathway]
  }

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < fdr_threshold, , drop = FALSE]

  if (verbose) {
    n_enriched <- sum(sig$NES > 0, na.rm = TRUE)
    n_depleted <- sum(sig$NES < 0, na.rm = TRUE)
    cat("[APOLLO] --- Summary ---\n")
    cat("    Gene sets tested:", nrow(res_df), "\n")
    cat("    FDR <", fdr_threshold, ":", nrow(sig), "significant (",
        n_enriched, "enriched,", n_depleted, "depleted)\n")
  }

  # Keep only tested gene sets (passed min/max size) for enrichment plots
  tested_sets <- gene_sets[names(gene_sets) %in% res_df$pathway]

  result <- list(
    results      = res_df,
    significant  = sig,
    ranked_genes = ranked_genes,
    gene_sets    = tested_sets,
    params       = list(
      species       = species,
      collection    = collection,
      min_size      = min_size,
      max_size      = max_size,
      n_genes       = length(ranked_genes),
      n_sets        = nrow(res_df),
      fdr_threshold = fdr_threshold
    )
  )
  class(result) <- c("apollo_gsea", "list")

  return(result)
}


#' @method print apollo_gsea
#' @export
print.apollo_gsea <- function(x, ...) {
  cat("GSEA Result (fgsea)\n")
  cat("------------------------------\n")
  cat("Ranked genes : ", x$params$n_genes, "\n", sep = "")
  cat("Gene sets    : ", x$params$n_sets, " tested\n", sep = "")
  n_enr <- sum(x$significant$NES > 0, na.rm = TRUE)
  n_dep <- sum(x$significant$NES < 0, na.rm = TRUE)
  cat("Significant (FDR < ", x$params$fdr_threshold, "): ", nrow(x$significant),
      " (enriched: ", n_enr, ", depleted: ", n_dep, ")\n", sep = "")
  cat("\nSlots: $results, $significant, $ranked_genes, $gene_sets, $params\n")
  invisible(x)
}
