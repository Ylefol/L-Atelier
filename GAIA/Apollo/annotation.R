###############################################################################
########### Apollo - Annotation & Enrichment Functions ###########
###############################################################################


#' Built-in Chromosome Name Mappings
#'
#' @description Returns a named vector mapping NCBI accession numbers to
#' UCSC-style chromosome names for supported genomes.
#'
#' @param genome Character string. Genome identifier. Currently supported:
#'   "T2T" or "T2T-CHM13v2.0" for the Telomere-to-Telomere human genome.
#'
#' @return A named character vector where names are NCBI accessions and
#'   values are UCSC-style chromosome names.
#'
#' @export
#'
#' @examples
#' mapping <- APOLLO_get_chr_mapping("T2T")
#' # Use with APOLLO_get_chromosome_sizes:
#' chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", name_mapping = "T2T")
#'
APOLLO_get_chr_mapping <- function(genome) {

  genome <- toupper(genome)

  if (genome %in% c("T2T", "T2T-CHM13V2.0", "T2T-CHM13V2", "CHM13")) {
    # T2T-CHM13v2.0 (GCF_009914755.1) NCBI accession to UCSC mapping
    mapping <- c(
      "NC_060925.1" = "chr1",
      "NC_060926.1" = "chr2",
      "NC_060927.1" = "chr3",
      "NC_060928.1" = "chr4",
      "NC_060929.1" = "chr5",
      "NC_060930.1" = "chr6",
      "NC_060931.1" = "chr7",
      "NC_060932.1" = "chr8",
      "NC_060933.1" = "chr9",
      "NC_060934.1" = "chr10",
      "NC_060935.1" = "chr11",
      "NC_060936.1" = "chr12",
      "NC_060937.1" = "chr13",
      "NC_060938.1" = "chr14",
      "NC_060939.1" = "chr15",
      "NC_060940.1" = "chr16",
      "NC_060941.1" = "chr17",
      "NC_060942.1" = "chr18",
      "NC_060943.1" = "chr19",
      "NC_060944.1" = "chr20",
      "NC_060945.1" = "chr21",
      "NC_060946.1" = "chr22",
      "NC_060947.1" = "chrX",
      "NC_060948.1" = "chrY",
      "NC_060949.1" = "chrM"
    )
    return(mapping)
  }

  stop("Unknown genome: ", genome, "\n",
       "Currently supported: T2T (T2T-CHM13v2.0)\n",
       "For other genomes, provide a custom mapping vector.")
}


#' Get Chromosome Sizes from GTF/GFF Annotation File
#'
#' @description Extracts chromosome names and sizes from a GTF or GFF annotation
#' file. Useful for circos plots and other visualizations that need chromosome
#' dimensions, especially for non-standard genomes (e.g., T2T instead of hg38).
#'
#' @param annotation_path Character string. Path to GTF or GFF file.
#'   Supports gzipped files (.gz extension).
#' @param chromosomes Optional character vector of chromosome names to include.
#'   If NULL (default), includes all chromosomes found. Use this to filter to
#'   standard chromosomes, e.g., c(paste0("chr", 1:22), "chrX", "chrY").
#'   Note: filtering happens AFTER name_mapping is applied.
#' @param name_mapping Optional. Rename chromosomes from NCBI accessions to
#'   UCSC-style names. Can be:
#'   \itemize{
#'     \item A character string for built-in mappings: "T2T" for T2T-CHM13v2.0
#'     \item A named character vector: c("NC_060925.1" = "chr1", ...)
#'   }
#' @param add_chr_prefix Logical. If TRUE and chromosome names don't start with
#'   "chr", adds the prefix. Applied AFTER name_mapping (default = FALSE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{chr}{Chromosome name}
#'   \item{size}{Chromosome size (max coordinate found in annotation)}
#' }
#'
#' @details
#' The function determines chromosome sizes by finding the maximum coordinate
#' (end position) for each chromosome in the annotation file. This is reliable
#' for well-annotated genomes where annotations extend near chromosome ends.
#'
#' For GFF3 files with ##sequence-region pragmas, those values are used directly
#' as they provide exact chromosome lengths.
#'
#' @export
#'
#' @examples
#' # Get all chromosome sizes
#' chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf.gz")
#'
#' # Filter to standard human chromosomes
#' standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
#' chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", chromosomes = standard_chrs)
#'
#' # For T2T genome with built-in mapping
#' chr_sizes <- APOLLO_get_chromosome_sizes(
#'   "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf",
#'   name_mapping = "T2T",
#'   chromosomes = c(paste0("chr", 1:22), "chrX", "chrY")
#' )
#'
#' # Custom mapping
#' my_mapping <- c("NC_000001.11" = "chr1", "NC_000002.12" = "chr2")
#' chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", name_mapping = my_mapping)
#'
APOLLO_get_chromosome_sizes <- function(annotation_path,
                                         chromosomes = NULL,
                                         name_mapping = NULL,
                                         add_chr_prefix = FALSE,
                                         verbose = TRUE) {

  if (!file.exists(annotation_path)) {
    stop("Annotation file not found: ", annotation_path)
  }

  if (verbose) {
    cat("Extracting chromosome sizes from:", basename(annotation_path), "\n")
  }

  # Determine if file is gzipped
  is_gzipped <- grepl("\\.gz$", annotation_path)

  # Determine file type (GTF or GFF)
  is_gtf <- grepl("\\.(gtf|gtf\\.gz)$", annotation_path, ignore.case = TRUE)
  is_gff <- grepl("\\.(gff|gff3|gff\\.gz|gff3\\.gz)$", annotation_path, ignore.case = TRUE)

  if (!is_gtf && !is_gff) {
    warning("Could not determine file type from extension. Assuming GTF format.")
    is_gtf <- TRUE
  }

  # Open connection
  if (is_gzipped) {
    con <- gzfile(annotation_path, "r")
  } else {
    con <- file(annotation_path, "r")
  }
  on.exit(close(con))

  # For GFF3, try to extract ##sequence-region pragmas first
  sequence_regions <- list()
  chr_max_coords <- list()

  if (verbose) {
    cat("  Scanning annotation file...\n")
  }

  # Read and process line by line to handle large files
  line_count <- 0
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break

    line_count <- line_count + 1

    # Progress indicator for large files
    if (verbose && line_count %% 500000 == 0) {
      cat("    Processed", format(line_count, big.mark = ","), "lines...\n")
    }

    # Skip empty lines
    if (nchar(trimws(line)) == 0) next

    # Check for GFF3 sequence-region pragma
    if (startsWith(line, "##sequence-region")) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) >= 4) {
        chr_name <- parts[2]
        chr_end <- as.numeric(parts[4])
        if (!is.na(chr_end)) {
          sequence_regions[[chr_name]] <- chr_end
        }
      }
      next
    }

    # Skip other comment lines
    if (startsWith(line, "#")) next

    # Parse data line (tab-separated)
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 5) next

    chr_name <- parts[1]
    end_pos <- suppressWarnings(as.numeric(parts[5]))

    if (is.na(end_pos)) next

    # Track max coordinate per chromosome
    if (is.null(chr_max_coords[[chr_name]])) {
      chr_max_coords[[chr_name]] <- end_pos
    } else {
      chr_max_coords[[chr_name]] <- max(chr_max_coords[[chr_name]], end_pos)
    }
  }

  if (verbose) {
    cat("  Scanned", format(line_count, big.mark = ","), "lines total.\n")
  }

  # Prefer sequence-region pragmas if available, otherwise use max coordinates
  if (length(sequence_regions) > 0) {
    if (verbose) {
      cat("  Using ##sequence-region pragmas for chromosome sizes.\n")
    }
    chr_sizes <- sequence_regions
  } else {
    if (verbose) {
      cat("  Using max coordinates for chromosome sizes.\n")
    }
    chr_sizes <- chr_max_coords
  }

  if (length(chr_sizes) == 0) {
    stop("No chromosome information found in annotation file.")
  }

  # Convert to data.frame
  result <- data.frame(
    chr = names(chr_sizes),
    size = as.numeric(unlist(chr_sizes)),
    stringsAsFactors = FALSE
  )

  # Apply name mapping if provided
  if (!is.null(name_mapping)) {
    # If string, get built-in mapping
    if (is.character(name_mapping) && length(name_mapping) == 1 &&
        !any(names(name_mapping) != "")) {
      if (verbose) {
        cat("  Applying built-in chromosome mapping for:", name_mapping, "\n")
      }
      name_mapping <- APOLLO_get_chr_mapping(name_mapping)
    } else if (verbose) {
      cat("  Applying custom chromosome name mapping...\n")
    }

    # Apply mapping - only rename chromosomes that are in the mapping
    mapped_idx <- result$chr %in% names(name_mapping)
    if (sum(mapped_idx) > 0) {
      result$chr[mapped_idx] <- name_mapping[result$chr[mapped_idx]]
      if (verbose) {
        cat("  Renamed", sum(mapped_idx), "chromosomes.\n")
      }
    } else {
      warning("No chromosomes matched the provided name_mapping.")
    }
  }

  # Add chr prefix if requested
  if (add_chr_prefix) {
    needs_prefix <- !grepl("^chr", result$chr, ignore.case = TRUE)
    result$chr[needs_prefix] <- paste0("chr", result$chr[needs_prefix])
  }

  # Filter to requested chromosomes
  if (!is.null(chromosomes)) {
    result <- result[result$chr %in% chromosomes, , drop = FALSE]

    if (nrow(result) == 0) {
      stop("None of the requested chromosomes found in annotation.\n",
           "  Requested: ", paste(head(chromosomes, 5), collapse = ", "),
           if (length(chromosomes) > 5) "..." else "", "\n",
           "  Found: ", paste(head(names(chr_sizes), 5), collapse = ", "),
           if (length(chr_sizes) > 5) "..." else "")
    }

    # Preserve requested order
    result <- result[match(chromosomes[chromosomes %in% result$chr], result$chr), ]
  }

  # Sort by chromosome (natural sort if possible)
  # Simple approach: extract numeric part for sorting
  chr_order <- order(
    !grepl("^chr[0-9]+$", result$chr),  # Standard numbered chrs first
    suppressWarnings(as.numeric(gsub("chr", "", result$chr))),  # By number
    result$chr  # Then alphabetically
  )
  result <- result[chr_order, ]
  rownames(result) <- NULL

  if (verbose) {
    cat("  Found", nrow(result), "chromosomes.\n")
    cat("  Total genome size:", format(sum(result$size), big.mark = ","), "bp\n")
  }

  return(result)
}


#' Convert Chromosome Sizes to Seqinfo Object
#'
#' @description Converts a chromosome sizes data.frame to a Seqinfo object
#' for use with GenomicRanges and related packages.
#'
#' @param chr_sizes A data.frame with columns 'chr' and 'size', typically
#'   from APOLLO_get_chromosome_sizes().
#' @param genome Optional character string for genome name (e.g., "hg38", "T2T").
#'
#' @return A Seqinfo object.
#'
#' @export
#'
APOLLO_chr_sizes_to_seqinfo <- function(chr_sizes, genome = NA_character_) {

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package 'GenomeInfoDb' is required. ",
         "Install with: BiocManager::install('GenomeInfoDb')")
  }

  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = chr_sizes$chr,
    seqlengths = chr_sizes$size,
    isCircular = rep(FALSE, nrow(chr_sizes)),
    genome = genome
  )

  return(seqinfo)
}
