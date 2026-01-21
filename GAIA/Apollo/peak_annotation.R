#' Apollo - Peak Annotation & Enrichment Functions
#'
#' @description Functions for annotating genomic peaks with gene information
#' and performing pathway/GO enrichment analysis.


#' Create or Load TxDb from GTF/GFF Annotation
#'
#' @description Creates a TxDb object from a GTF/GFF annotation file, with
#' automatic caching to avoid recreating it on subsequent runs.
#'
#' @param gtf_path Character string. Path to GTF or GFF annotation file.
#' @param cache_dir Character string. Directory to store cached TxDb SQLite file.
#'   Default is "data/" relative to ZERO_DAWN root.
#' @param cache_name Character string or NULL. Name for the cached SQLite file
#'   (without .sqlite extension). If NULL, derives from GTF filename.
#' @param chr_mapping Character string or named vector. Chromosome name mapping
#'   to apply. Can be:
#'   \itemize{
#'     \item "T2T" - Use built-in T2T-CHM13 NCBI->UCSC mapping
#'     \item A named vector: c("NC_060925.1" = "chr1", ...)
#'     \item NULL - No renaming (default)
#'   }
#' @param organism Character string. Organism name for TxDb metadata
#'   (default = "Homo sapiens").
#' @param force Logical. If TRUE, recreate TxDb even if cache exists
#'   (default = FALSE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A TxDb object.
#'
#' @details
#' The function checks for a cached TxDb in cache_dir:
#' - If found, loads and returns it (fast)
#' - If not found, creates TxDb from GTF and saves to cache (slow, but one-time)
#'
#' TxDb objects are saved using AnnotationDbi::saveDb() as SQLite databases,
#' which is the proper serialization method for these objects.
#'
#' This is particularly useful for non-standard genomes like T2T-CHM13 that
#' don't have pre-built TxDb packages on Bioconductor.
#'
#' @export
#'
#' @examples
#' # First run: creates TxDb from GTF (slow)
#' txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T")
#'
#' # Subsequent runs: loads from cache (fast)
#' txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T")
#'
#' # Force recreation
#' txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T", force = TRUE)
#'
APOLLO_make_txdb <- function(gtf_path,
                              cache_dir = NULL,
                              cache_name = NULL,
                              chr_mapping = NULL,
                              organism = "Homo sapiens",
                              force = FALSE,
                              verbose = TRUE) {

  # Check for required packages
  if (!requireNamespace("txdbmaker", quietly = TRUE)) {
    stop("Package 'txdbmaker' is required. Install from Bioconductor:\n",
         "  BiocManager::install('txdbmaker')")
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required. Install from Bioconductor:\n",
         "  BiocManager::install('AnnotationDbi')")
  }

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required. Install from Bioconductor:\n",
         "  BiocManager::install('GenomicFeatures')")
  }

  # Validate GTF path
  if (!file.exists(gtf_path)) {
    stop("GTF file not found: ", gtf_path)
  }

  # Set default cache directory (ZERO_DAWN/data/)
  if (is.null(cache_dir)) {
    # Try to find ZERO_DAWN root by looking for GAIA directory
    script_dir <- getwd()
    if (dir.exists(file.path(script_dir, "GAIA"))) {
      cache_dir <- file.path(script_dir, "data")
    } else if (dir.exists(file.path(dirname(script_dir), "GAIA"))) {
      cache_dir <- file.path(dirname(script_dir), "data")
    } else if (dir.exists(file.path(dirname(dirname(script_dir)), "GAIA"))) {
      cache_dir <- file.path(dirname(dirname(script_dir)), "data")
    } else {
      # Fall back to current directory
      cache_dir <- file.path(script_dir, "data")
    }
  }

  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    if (verbose) cat("Creating cache directory:", cache_dir, "\n")
    dir.create(cache_dir, recursive = TRUE)
  }

  # Determine cache filename
  if (is.null(cache_name)) {
    # Derive from GTF filename
    cache_name <- gsub("\\.(gtf|gff|gff3)(\\.gz)?$", "", basename(gtf_path),
                       ignore.case = TRUE)
    cache_name <- paste0(cache_name, "_TxDb")
  }

  # TxDb uses SQLite format
  cache_path <- file.path(cache_dir, paste0(cache_name, ".sqlite"))

  # Check for existing cache
  if (file.exists(cache_path) && !force) {
    if (verbose) {
      cat("Loading cached TxDb from:", cache_path, "\n")
    }

    txdb <- tryCatch({
      AnnotationDbi::loadDb(cache_path)
    }, error = function(e) {
      warning("Failed to load cached TxDb: ", e$message, "\nRecreating...")
      NULL
    })

    # Verify it's a valid TxDb
    if (!is.null(txdb) && inherits(txdb, "TxDb")) {
      if (verbose) cat("TxDb loaded successfully.\n")
      return(txdb)
    }
  }

  # Create TxDb from GTF
  if (verbose) {
    cat("Creating TxDb from GTF annotation...\n")
    cat("  Source:", gtf_path, "\n")
    cat("  This may take several minutes for large annotation files.\n")
  }

  # Determine format
  is_gff <- grepl("\\.(gff|gff3)(\\.gz)?$", gtf_path, ignore.case = TRUE)
  format <- if (is_gff) "gff3" else "gtf"

  # Use txdbmaker (the modern package for this)
  txdb <- txdbmaker::makeTxDbFromGFF(
    file = gtf_path,
    format = format,
    organism = organism,
    dataSource = basename(gtf_path)
  )

  # Apply chromosome name mapping if requested
  if (!is.null(chr_mapping)) {
    if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
      stop("Package 'GenomeInfoDb' is required for chromosome renaming.\n",
           "  BiocManager::install('GenomeInfoDb')")
    }

    # Get mapping vector
    if (is.character(chr_mapping) && length(chr_mapping) == 1) {
      # Use built-in mapping
      if (verbose) cat("  Applying", chr_mapping, "chromosome name mapping...\n")
      mapping_vec <- APOLLO_get_chr_mapping(chr_mapping)
    } else if (is.character(chr_mapping) && !is.null(names(chr_mapping))) {
      mapping_vec <- chr_mapping
    } else {
      stop("chr_mapping must be a genome name (e.g., 'T2T') or a named character vector")
    }

    # Get current seqlevels
    current_seqlevels <- GenomeInfoDb::seqlevels(txdb)

    # Build rename vector (only for chromosomes that exist in the TxDb)
    to_rename <- intersect(current_seqlevels, names(mapping_vec))

    if (length(to_rename) > 0) {
      rename_vec <- mapping_vec[to_rename]
      names(rename_vec) <- to_rename

      txdb <- GenomeInfoDb::renameSeqlevels(txdb, rename_vec)

      if (verbose) {
        cat("  Renamed", length(to_rename), "chromosomes to UCSC style\n")
      }
    } else {
      warning("No chromosome names matched the mapping. Seqlevels unchanged.")
    }
  }

  # Save to cache using AnnotationDbi::saveDb (proper method for TxDb)
  if (verbose) {
    cat("Saving TxDb to cache:", cache_path, "\n")
  }
  AnnotationDbi::saveDb(txdb, file = cache_path)

  if (verbose) {
    cat("TxDb created successfully.\n")
    # Print some basic stats
    cat("  Genes:", length(GenomicFeatures::genes(txdb)), "\n")
    cat("  Transcripts:", length(GenomicFeatures::transcripts(txdb)), "\n")
    cat("  Chromosomes:", length(GenomeInfoDb::seqlevels(txdb)), "\n")
  }

  return(txdb)
}


#' Annotate Peaks with Genomic Features and Nearest Genes
#'
#' @description Annotates genomic peaks with their genomic context (promoter,
#' exon, intron, intergenic, etc.) and nearest gene using ChIPseeker.
#'
#' @param regions Data.frame with at minimum: chr, start, end. Can also have
#'   peak_id or name column for identification.
#' @param txdb A TxDb object (from APOLLO_make_txdb or Bioconductor).
#' @param tss_region Numeric vector of length 2. Region around TSS to define
#'   as promoter, e.g., c(-3000, 3000) means 3kb upstream to 3kb downstream.
#' @param level Character. Annotation level - "transcript" or "gene"
#'   (default = "transcript").
#' @param verbose Logical. Print summary (default = TRUE).
#'
#' @return A data.frame with original region info plus annotation columns:
#' \describe{
#'   \item{annotation}{Genomic feature (Promoter, Exon, Intron, etc.)}
#'   \item{annotation_simple}{Simplified category}
#'   \item{gene_id}{Nearest gene ID}
#'   \item{gene_name}{Gene symbol (if available in TxDb)}
#'   \item{distance_to_tss}{Distance to nearest TSS}
#'   \item{transcript_id}{Associated transcript ID}
#' }
#'
#' @details
#' Uses ChIPseeker::annotatePeak() for annotation. The TSS region defines
#' what counts as a "promoter" - regions within this window of any TSS are
#' labeled as promoter regions.
#'
#' @export
#'
#' @examples
#' txdb <- APOLLO_make_txdb("annotation.gtf")
#' annotated <- APOLLO_annotate_peaks(significant_peaks, txdb)
#'
#' # See distribution of genomic features
#' table(annotated$annotation_simple)
#'
#' # Get genes for pathway analysis
#' genes <- unique(annotated$gene_id)
#'
APOLLO_annotate_peaks <- function(regions,
                                   txdb,
                                   tss_region = c(-3000, 3000),
                                   level = "transcript",
                                   verbose = TRUE) {

  # Check for required packages
  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    stop("Package 'ChIPseeker' is required. Install from Bioconductor:\n",
         "  BiocManager::install('ChIPseeker')")
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required. Install from Bioconductor:\n",
         "  BiocManager::install('GenomicRanges')")
  }

  # Validate inputs
  if (!is.data.frame(regions)) {
    stop("regions must be a data.frame")
  }

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop("regions missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(regions) == 0) {
    stop("regions data.frame is empty")
  }

  if (!inherits(txdb, "TxDb")) {
    stop("txdb must be a TxDb object")
  }

  # Determine peak ID column
  id_cols <- c("peak_id", "region_id", "name", "id")
  id_col <- intersect(id_cols, colnames(regions))[1]
  if (is.na(id_col)) {
    regions$`.peak_id` <- paste0("peak_", seq_len(nrow(regions)))
    id_col <- ".peak_id"
  }

  if (verbose) {
    cat("Annotating", nrow(regions), "peaks with genomic features...\n")
  }

  # Convert to GRanges
  peaks_gr <- GenomicRanges::GRanges(
    seqnames = regions$chr,
    ranges = IRanges::IRanges(
      start = regions$start + 1,  # Convert 0-based to 1-based
      end = regions$end
    ),
    peak_id = regions[[id_col]]
  )

  # Annotate with ChIPseeker
  # Suppress messages from ChIPseeker
  anno <- suppressMessages(
    ChIPseeker::annotatePeak(
      peaks_gr,
      TxDb = txdb,
      tssRegion = tss_region,
      level = level,
      verbose = FALSE
    )
  )

  # Convert to data.frame
  anno_df <- as.data.frame(anno)

  # Create simplified annotation categories
  anno_df$annotation_simple <- sapply(anno_df$annotation, function(x) {
    if (grepl("Promoter", x)) return("Promoter")
    if (grepl("5' UTR", x)) return("5' UTR")
    if (grepl("3' UTR", x)) return("3' UTR")
    if (grepl("Exon", x)) return("Exon")
    if (grepl("Intron", x)) return("Intron")
    if (grepl("Downstream", x)) return("Downstream")
    if (grepl("Intergenic", x)) return("Intergenic")
    return("Other")
  })

  # Build result data.frame with clean column names
  result <- data.frame(
    peak_id = anno_df$peak_id,
    chr = as.character(anno_df$seqnames),
    start = anno_df$start - 1,  # Convert back to 0-based
    end = anno_df$end,
    width = anno_df$width,
    annotation = anno_df$annotation,
    annotation_simple = anno_df$annotation_simple,
    gene_id = anno_df$geneId,
    distance_to_tss = anno_df$distanceToTSS,
    stringsAsFactors = FALSE
  )

  # Add transcript_id if available
  if ("transcriptId" %in% colnames(anno_df)) {
    result$transcript_id <- anno_df$transcriptId
  }

  # Add SYMBOL if available (depends on TxDb having this info)
  if ("SYMBOL" %in% colnames(anno_df)) {
    result$gene_name <- anno_df$SYMBOL
  } else {
    # Gene name not available - user may need to add via org.Db
    result$gene_name <- NA_character_
  }

  # Merge back any additional columns from original regions
  extra_cols <- setdiff(colnames(regions), c("chr", "start", "end", id_col))
  if (length(extra_cols) > 0) {
    # Match by peak_id
    match_idx <- match(result$peak_id, regions[[id_col]])
    for (col in extra_cols) {
      result[[col]] <- regions[[col]][match_idx]
    }
  }

  if (verbose) {
    cat("\nAnnotation summary:\n")
    cat("  Total peaks:", nrow(result), "\n")
    cat("  Unique genes:", length(unique(result$gene_id[!is.na(result$gene_id)])), "\n")

    # Feature distribution
    cat("\nGenomic feature distribution:\n")
    feature_table <- sort(table(result$annotation_simple), decreasing = TRUE)
    for (feat in names(feature_table)) {
      pct <- round(100 * feature_table[feat] / nrow(result), 1)
      cat(sprintf("  %-12s: %4d (%5.1f%%)\n", feat, feature_table[feat], pct))
    }
  }

  return(result)
}


#' GO Enrichment Analysis for Annotated Peaks
#'
#' @description Performs Gene Ontology enrichment analysis on genes associated
#' with annotated peaks using clusterProfiler.
#'
#' @param annotated_peaks Data.frame from APOLLO_annotate_peaks() with gene_id column.
#'   Alternatively, a character vector of gene IDs.
#' @param org_db Character string or OrgDb object. The organism annotation database.
#'   Default is "org.Hs.eg.db" for human.
#' @param ont Character. GO ontology: "BP" (Biological Process), "MF" (Molecular
#'   Function), "CC" (Cellular Component), or "ALL" (default = "BP").
#' @param gene_id_type Character. Type of gene IDs provided. One of "ENTREZID",
#'   "ENSEMBL", "SYMBOL", "REFSEQ" (default = "ENTREZID").
#' @param pval_cutoff Numeric. P-value cutoff for enrichment (default = 0.05).
#' @param qval_cutoff Numeric. Adjusted p-value cutoff (default = 0.1).
#' @param min_gs_size Integer. Minimum gene set size (default = 10).
#' @param max_gs_size Integer. Maximum gene set size (default = 500).
#' @param verbose Logical. Print summary (default = TRUE).
#'
#' @return An enrichResult object from clusterProfiler. Key methods:
#' \describe{
#'   \item{as.data.frame()}{Convert to data.frame of results}
#'   \item{dotplot()}{Create dot plot of top terms}
#'   \item{barplot()}{Create bar plot of top terms}
#'   \item{cnetplot()}{Gene-concept network plot}
#' }
#'
#' @details
#' If gene IDs are not ENTREZID, the function attempts to convert them using
#' the org.Db. This is necessary because GO enrichment requires ENTREZ IDs.
#'
#' For T2T-CHM13 or other genomes, gene symbols should still work with
#' org.Hs.eg.db since the gene names are the same.
#'
#' @export
#'
#' @examples
#' # From annotated peaks
#' annotated <- APOLLO_annotate_peaks(sig_peaks, txdb)
#' go_results <- APOLLO_enrich_go(annotated)
#'
#' # View results
#' head(as.data.frame(go_results))
#'
#' # Plot
#' dotplot(go_results, showCategory = 20)
#'
#' # From gene vector
#' genes <- c("TP53", "BRCA1", "EGFR")
#' go_results <- APOLLO_enrich_go(genes, gene_id_type = "SYMBOL")
#'
APOLLO_enrich_go <- function(annotated_peaks,
                              org_db = "org.Hs.eg.db",
                              ont = "BP",
                              gene_id_type = "ENTREZID",
                              pval_cutoff = 0.05,
                              qval_cutoff = 0.1,
                              min_gs_size = 10,
                              max_gs_size = 500,
                              verbose = TRUE) {

  # Check for required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required. Install from Bioconductor:\n",
         "  BiocManager::install('clusterProfiler')")
  }

  # Load org.Db if string provided
  if (is.character(org_db)) {
    if (!requireNamespace(org_db, quietly = TRUE)) {
      stop("Package '", org_db, "' is required. Install from Bioconductor:\n",
           "  BiocManager::install('", org_db, "')")
    }
    org_db <- get(org_db, envir = loadNamespace(org_db))
  }

  # Extract gene IDs
  if (is.data.frame(annotated_peaks)) {
    if (!"gene_id" %in% colnames(annotated_peaks)) {
      stop("annotated_peaks must have a 'gene_id' column. ",
           "Run APOLLO_annotate_peaks() first.")
    }
    gene_ids <- unique(annotated_peaks$gene_id)
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  } else if (is.character(annotated_peaks)) {
    gene_ids <- unique(annotated_peaks)
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  } else {
    stop("annotated_peaks must be a data.frame or character vector of gene IDs")
  }

  if (length(gene_ids) == 0) {
    stop("No valid gene IDs found")
  }

  if (verbose) {
    cat("GO Enrichment Analysis\n")
    cat("  Input genes:", length(gene_ids), "\n")
    cat("  Ontology:", ont, "\n")
    cat("  Gene ID type:", gene_id_type, "\n")
  }

  # Convert to ENTREZID if needed
  if (toupper(gene_id_type) != "ENTREZID") {
    if (verbose) cat("  Converting", gene_id_type, "to ENTREZID...\n")

    converted <- tryCatch({
      clusterProfiler::bitr(
        gene_ids,
        fromType = toupper(gene_id_type),
        toType = "ENTREZID",
        OrgDb = org_db
      )
    }, error = function(e) {
      stop("Failed to convert gene IDs: ", e$message)
    })

    if (nrow(converted) == 0) {
      stop("No gene IDs could be converted to ENTREZID")
    }

    if (verbose) {
      cat("  Successfully converted:", nrow(converted), "of", length(gene_ids), "\n")
    }

    gene_ids <- unique(converted$ENTREZID)
  }

  if (verbose) {
    cat("  Running enrichment...\n")
  }

  # Run GO enrichment
  ego <- clusterProfiler::enrichGO(
    gene = gene_ids,
    OrgDb = org_db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE  # Convert IDs to symbols in output
  )

  if (verbose) {
    n_sig <- sum(ego@result$p.adjust < qval_cutoff)
    cat("\nResults:\n")
    cat("  Significant terms (q <", qval_cutoff, "):", n_sig, "\n")

    if (n_sig > 0) {
      cat("\nTop 10 enriched terms:\n")
      top_terms <- head(ego@result[ego@result$p.adjust < qval_cutoff, ], 10)
      for (i in seq_len(nrow(top_terms))) {
        cat(sprintf("  %2d. %s (q=%.2e)\n",
                    i,
                    substr(top_terms$Description[i], 1, 50),
                    top_terms$p.adjust[i]))
      }
    }
  }

  return(ego)
}


#' KEGG Pathway Enrichment Analysis
#'
#' @description Performs KEGG pathway enrichment analysis on genes associated
#' with annotated peaks using clusterProfiler.
#'
#' @param annotated_peaks Data.frame from APOLLO_annotate_peaks() with gene_id column.
#'   Alternatively, a character vector of gene IDs.
#' @param org_db Character string or OrgDb object. The organism annotation database
#'   for ID conversion. Default is "org.Hs.eg.db" for human.
#' @param organism Character. KEGG organism code (default = "hsa" for human).
#' @param gene_id_type Character. Type of gene IDs provided. One of "ENTREZID",
#'   "ENSEMBL", "SYMBOL", "REFSEQ" (default = "ENTREZID").
#' @param pval_cutoff Numeric. P-value cutoff for enrichment (default = 0.05).
#' @param qval_cutoff Numeric. Adjusted p-value cutoff (default = 0.1).
#' @param verbose Logical. Print summary (default = TRUE).
#'
#' @return An enrichResult object from clusterProfiler.
#'
#' @export
#'
#' @examples
#' annotated <- APOLLO_annotate_peaks(sig_peaks, txdb)
#' kegg_results <- APOLLO_enrich_kegg(annotated)
#' dotplot(kegg_results)
#'
APOLLO_enrich_kegg <- function(annotated_peaks,
                                org_db = "org.Hs.eg.db",
                                organism = "hsa",
                                gene_id_type = "ENTREZID",
                                pval_cutoff = 0.05,
                                qval_cutoff = 0.1,
                                verbose = TRUE) {

  # Check for required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required. Install from Bioconductor:\n",
         "  BiocManager::install('clusterProfiler')")
  }

  # Load org.Db if string provided
  if (is.character(org_db)) {
    if (!requireNamespace(org_db, quietly = TRUE)) {
      stop("Package '", org_db, "' is required. Install from Bioconductor:\n",
           "  BiocManager::install('", org_db, "')")
    }
    org_db <- get(org_db, envir = loadNamespace(org_db))
  }

  # Extract gene IDs
  if (is.data.frame(annotated_peaks)) {
    if (!"gene_id" %in% colnames(annotated_peaks)) {
      stop("annotated_peaks must have a 'gene_id' column. ",
           "Run APOLLO_annotate_peaks() first.")
    }
    gene_ids <- unique(annotated_peaks$gene_id)
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  } else if (is.character(annotated_peaks)) {
    gene_ids <- unique(annotated_peaks)
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  } else {
    stop("annotated_peaks must be a data.frame or character vector of gene IDs")
  }

  if (length(gene_ids) == 0) {
    stop("No valid gene IDs found")
  }

  if (verbose) {
    cat("KEGG Pathway Enrichment Analysis\n")
    cat("  Input genes:", length(gene_ids), "\n")
    cat("  Organism:", organism, "\n")
    cat("  Gene ID type:", gene_id_type, "\n")
  }

  # Convert to ENTREZID if needed
  if (toupper(gene_id_type) != "ENTREZID") {
    if (verbose) cat("  Converting", gene_id_type, "to ENTREZID...\n")

    converted <- tryCatch({
      clusterProfiler::bitr(
        gene_ids,
        fromType = toupper(gene_id_type),
        toType = "ENTREZID",
        OrgDb = org_db
      )
    }, error = function(e) {
      stop("Failed to convert gene IDs: ", e$message)
    })

    if (nrow(converted) == 0) {
      stop("No gene IDs could be converted to ENTREZID")
    }

    if (verbose) {
      cat("  Successfully converted:", nrow(converted), "of", length(gene_ids), "\n")
    }

    gene_ids <- unique(converted$ENTREZID)
  }

  if (verbose) {
    cat("  Running enrichment...\n")
  }

  # Run KEGG enrichment
  ekegg <- clusterProfiler::enrichKEGG(
    gene = gene_ids,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff
  )

  if (verbose) {
    n_sig <- sum(ekegg@result$p.adjust < qval_cutoff)
    cat("\nResults:\n")
    cat("  Significant pathways (q <", qval_cutoff, "):", n_sig, "\n")

    if (n_sig > 0) {
      cat("\nTop 10 enriched pathways:\n")
      top_paths <- head(ekegg@result[ekegg@result$p.adjust < qval_cutoff, ], 10)
      for (i in seq_len(nrow(top_paths))) {
        cat(sprintf("  %2d. %s (q=%.2e)\n",
                    i,
                    substr(top_paths$Description[i], 1, 50),
                    top_paths$p.adjust[i]))
      }
    }
  }

  return(ekegg)
}
