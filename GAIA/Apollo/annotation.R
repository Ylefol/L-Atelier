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
#' Supports caching to avoid re-parsing large annotation files.
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
#' @param cache_dir Character string. Directory to store cached chromosome sizes.
#'   Default is "data/" relative to ZERO_DAWN root. Set to NULL to disable caching.
#' @param cache_name Character string or NULL. Name for the cached RDS file
#'   (without .rds extension). If NULL, derives from annotation filename.
#' @param force Logical. If TRUE, re-extract sizes even if cache exists
#'   (default = FALSE).
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
#' **Caching**: Results are cached as RDS files in cache_dir. The cache filename
#' includes the name_mapping used, so different mapping configurations get
#' separate cache files. Use force=TRUE to regenerate the cache.
#'
#' @export
#'
#' @examples
#' # Get all chromosome sizes (first run slow, subsequent runs fast)
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
#' # Force regeneration of cache
#' chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", force = TRUE)
#'
APOLLO_get_chromosome_sizes <- function(annotation_path,
                                         chromosomes = NULL,
                                         name_mapping = NULL,
                                         add_chr_prefix = FALSE,
                                         cache_dir = NULL,
                                         cache_name = NULL,
                                         force = FALSE,
                                         verbose = TRUE) {

  if (!file.exists(annotation_path)) {
    stop("Annotation file not found: ", annotation_path)
  }

  # ---------------------------------------------------------------------------
  # Set up caching
  # ---------------------------------------------------------------------------
  use_cache <- !is.null(cache_dir) || is.null(cache_dir)  # Default to using cache

  if (use_cache) {
    # Set default cache directory (ZERO_DAWN/data/)
    if (is.null(cache_dir)) {
      script_dir <- getwd()
      if (dir.exists(file.path(script_dir, "GAIA"))) {
        cache_dir <- file.path(script_dir, "data")
      } else if (dir.exists(file.path(dirname(script_dir), "GAIA"))) {
        cache_dir <- file.path(dirname(script_dir), "data")
      } else if (dir.exists(file.path(dirname(dirname(script_dir)), "GAIA"))) {
        cache_dir <- file.path(dirname(dirname(script_dir)), "data")
      } else {
        cache_dir <- file.path(script_dir, "data")
      }
    }

    # Create cache directory if needed
    if (!dir.exists(cache_dir)) {
      if (verbose) cat("Creating cache directory:", cache_dir, "\n")
      dir.create(cache_dir, recursive = TRUE)
    }

    # Determine cache filename
    if (is.null(cache_name)) {
      cache_name <- gsub("\\.(gtf|gff|gff3)(\\.gz)?$", "", basename(annotation_path),
                         ignore.case = TRUE)
      cache_name <- paste0(cache_name, "_chrom_sizes")
    }

    # Include name_mapping in cache filename
    if (!is.null(name_mapping)) {
      if (is.character(name_mapping) && length(name_mapping) == 1 &&
          (is.null(names(name_mapping)) || names(name_mapping)[1] == "")) {
        mapping_suffix <- paste0("_", toupper(name_mapping))
      } else {
        mapping_suffix <- "_chrMapped"
      }
      cache_name <- paste0(cache_name, mapping_suffix)
    }

    cache_path <- file.path(cache_dir, paste0(cache_name, ".rds"))

    # Check for existing cache
    if (file.exists(cache_path) && !force) {
      if (verbose) {
        cat("Loading cached chromosome sizes from:", basename(cache_path), "\n")
      }

      result <- tryCatch({
        readRDS(cache_path)
      }, error = function(e) {
        warning("Failed to load cached sizes: ", e$message, "\nRecreating...")
        NULL
      })

      if (!is.null(result) && is.data.frame(result)) {
        # Apply chromosome filter (not cached since user may want different subsets)
        if (!is.null(chromosomes)) {
          result <- result[result$chr %in% chromosomes, , drop = FALSE]
          if (nrow(result) > 0) {
            result <- result[match(chromosomes[chromosomes %in% result$chr], result$chr), ]
          }
        }

        if (verbose) {
          cat("  Loaded", nrow(result), "chromosomes.\n")
        }
        return(result)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Extract chromosome sizes from annotation file
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Save to cache (before chromosome filtering, so full set is cached)
  # ---------------------------------------------------------------------------
  if (use_cache) {
    if (verbose) {
      cat("  Saving chromosome sizes to cache:", basename(cache_path), "\n")
    }
    saveRDS(result, cache_path)
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


#' Calculate Sequence Composition for Genomic Regions
#'
#' @description Extracts DNA sequences for genomic regions and calculates
#' composition metrics including GC/AT content and repeat patterns.
#'
#' @param regions A data.frame with columns: chr, start, end. Additional columns
#'   are preserved. Coordinates should be 0-based (BED format).
#' @param fasta_path Path to a FASTA file. Must have an accompanying .fai index
#'   (create with `samtools faidx`).
#' @param chr_mapping Either:
#'   \itemize{
#'     \item NULL (default) - no chromosome name translation
#'     \item A genome name (e.g., "T2T") - uses built-in mapping
#'     \item A named character vector mapping region chr names to FASTA chr names
#'   }
#' @param extend Integer. Extend regions by this many bp on each side before
#'   extracting sequence (default = 0). Useful for standardizing region sizes.
#' @param min_width Integer. Minimum region width (after extension) to calculate
#'   composition (default = 100). Regions below this return NA values.
#' @param include_repeats Logical. Calculate homopolymer and dinucleotide repeat
#'   metrics (default = TRUE).
#' @param include_sequence Logical. Include the extracted sequence in output
#'   (default = TRUE). Set to FALSE to reduce memory usage for large datasets.
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with one row per input region containing:
#' \describe{
#'   \item{peak_id}{Region identifier (from input or generated)}
#'   \item{width}{Final region width after extension}
#'   \item{gc_percent}{Percentage of G+C bases}
#'   \item{at_percent}{Percentage of A+T bases}
#'   \item{n_Nbases}{Number of N (ambiguous) bases}
#'   \item{longest_homopolymer}{Length of longest single-base repeat (if include_repeats)}
#'   \item{homopolymer_base}{Base forming the longest homopolymer (if include_repeats)}
#'   \item{longest_dinucleotide}{Length of longest dinucleotide repeat in bp (if include_repeats)}
#'   \item{dinucleotide_motif}{Motif of longest dinucleotide repeat (if include_repeats)}
#'   \item{sequence}{The actual nucleotide sequence extracted (if include_sequence)}
#' }
#'
#' @details
#' **Chromosome mapping**: When using genomes like T2T-CHM13 where your regions
#' use UCSC-style names (chr1, chr2) but the FASTA uses NCBI accessions
#' (NC_060925.1, NC_060926.1), use chr_mapping = "T2T" to automatically translate.
#' The mapping is applied internally - your output will retain the original
#' region chromosome names.
#'
#' **Extension behavior**: Consistent with ELEUTHIA_expand_regions(), the extend
#' parameter adds the specified bp to EACH side. So extend = 50 adds 100bp total.
#' Start coordinates are clamped at 0.
#'
#' **Minimum width threshold**: Regions smaller than min_width (after extension)
#' return NA for composition values. This avoids unreliable percentages from very
#' short sequences (e.g., 10bp regions where each base is 10%).
#'
#' **Repeat detection**:
#' \itemize{
#'   \item Homopolymers: Consecutive identical bases (e.g., AAAA, TTTTTT)
#'   \item Dinucleotide repeats: Alternating two-base patterns (e.g., ATATAT, CGCGCG)
#' }
#'
#' @export
#'
#' @examples
#' # Basic usage with genomic regions
#' comp <- APOLLO_sequence_composition(regions, "reference.fa")
#'
#' # With T2T genome (regions have chr1, FASTA has NC_060925.1)
#' comp <- APOLLO_sequence_composition(regions, "T2T.fna", chr_mapping = "T2T")
#'
#' # With extension to standardize region size
#' comp <- APOLLO_sequence_composition(regions, "reference.fa", extend = 100)
#'
#' # Without sequence column (saves memory)
#' comp <- APOLLO_sequence_composition(regions, "ref.fa", include_sequence = FALSE)
#'
#' # Merge back to original data
#' regions_with_comp <- cbind(regions, comp[, c("gc_percent", "at_percent")])
#'
APOLLO_sequence_composition <- function(regions,
                                          fasta_path,
                                          chr_mapping = NULL,
                                          extend = 0,
                                          min_width = 100,
                                          include_repeats = TRUE,
                                          include_sequence = TRUE,
                                          verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check for required packages
  # ---------------------------------------------------------------------------
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required. Install from Bioconductor:\n",
         "  BiocManager::install('Rsamtools')")
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install from Bioconductor:\n",
         "  BiocManager::install('Biostrings')")
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required. Install from Bioconductor:\n",
         "  BiocManager::install('GenomicRanges')")
  }

  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(regions)) {
    stop("regions must be a data.frame")
  }

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop("regions missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  n_regions <- nrow(regions)
  if (n_regions == 0) {
    stop("regions data.frame is empty")
  }

  if (!file.exists(fasta_path)) {
    stop("FASTA file not found: ", fasta_path)
  }

  # Check for index file
  fai_path <- paste0(fasta_path, ".fai")
  if (!file.exists(fai_path)) {
    stop("FASTA index (.fai) not found: ", fai_path, "\n",
         "  Create with: samtools faidx ", fasta_path)
  }

  if (!is.numeric(extend) || length(extend) != 1 || extend < 0) {
    stop("extend must be a single non-negative number")
  }
  extend <- as.integer(extend)

  if (!is.numeric(min_width) || length(min_width) != 1 || min_width < 1) {
    stop("min_width must be a positive integer")
  }
  min_width <- as.integer(min_width)

  # ---------------------------------------------------------------------------
  # Determine peak ID column
  # ---------------------------------------------------------------------------
  id_cols <- c("peak_id", "region_id", "name", "id")
  id_col <- intersect(id_cols, colnames(regions))[1]
  if (is.na(id_col)) {
    peak_ids <- paste0("region_", seq_len(n_regions))
  } else {
    peak_ids <- regions[[id_col]]
  }

  # ---------------------------------------------------------------------------
  # Process chromosome mapping
  # ---------------------------------------------------------------------------
  # chr_mapping translates region chromosome names to FASTA chromosome names
  # e.g., for T2T: chr1 -> NC_060925.1
  chr_map_vec <- NULL

  if (!is.null(chr_mapping)) {
    if (is.character(chr_mapping) && length(chr_mapping) == 1 &&
        (is.null(names(chr_mapping)) || names(chr_mapping)[1] == "")) {
      # Built-in mapping name (e.g., "T2T")
      # Get the NCBI -> UCSC mapping and reverse it to UCSC -> NCBI
      ncbi_to_ucsc <- APOLLO_get_chr_mapping(chr_mapping)
      chr_map_vec <- setNames(names(ncbi_to_ucsc), as.character(ncbi_to_ucsc))
    } else if (is.character(chr_mapping) && !is.null(names(chr_mapping))) {
      # Custom mapping provided directly (region_name -> fasta_name)
      chr_map_vec <- chr_mapping
    } else {
      stop("chr_mapping must be a genome name (e.g., 'T2T') or a named character vector")
    }
  }

  if (verbose) {
    cat("Sequence composition analysis\n")
    cat("  Regions:", n_regions, "\n")
    cat("  FASTA:", basename(fasta_path), "\n")
    if (!is.null(chr_map_vec)) {
      cat("  Chromosome mapping: enabled (", length(chr_map_vec), " mappings)\n", sep = "")
    }
    if (extend > 0) {
      cat("  Extension: +/-", extend, "bp (", extend * 2, "bp total)\n")
    }
    cat("  Min width threshold:", min_width, "bp\n")
  }

  # ---------------------------------------------------------------------------
  # Apply extension (consistent with ELEUTHIA_expand_regions)
  # ---------------------------------------------------------------------------
  start_ext <- pmax(0L, as.integer(regions$start) - extend)
  end_ext <- as.integer(regions$end) + extend
  widths <- end_ext - start_ext

  # Identify regions below min_width threshold
  below_threshold <- widths < min_width
  n_below <- sum(below_threshold)

  if (verbose && n_below > 0) {
    cat("  Regions below min_width:", n_below, "(will return NA)\n")
  }

  # ---------------------------------------------------------------------------
  # Open FASTA file
  # ---------------------------------------------------------------------------
  fa <- Rsamtools::FaFile(fasta_path)
  open(fa)
  on.exit(close(fa), add = TRUE)

  # Get available chromosomes in FASTA
  fa_seqinfo <- Rsamtools::seqinfo(fa)
  fa_chroms <- GenomeInfoDb::seqnames(fa_seqinfo)

  # ---------------------------------------------------------------------------
  # Translate chromosome names if mapping provided
  # ---------------------------------------------------------------------------
  # Create vector of FASTA-compatible chromosome names
  if (!is.null(chr_map_vec)) {
    # Translate region chr names to FASTA names where mapping exists
    chr_for_fasta <- ifelse(
      regions$chr %in% names(chr_map_vec),
      chr_map_vec[regions$chr],
      regions$chr  # Keep original if not in mapping
    )
  } else {
    chr_for_fasta <- regions$chr
  }

  # ---------------------------------------------------------------------------
  # Create GRanges for sequence extraction
  # ---------------------------------------------------------------------------
  # Only include regions that are above threshold and have valid chromosomes
  valid_idx <- which(!below_threshold & chr_for_fasta %in% fa_chroms)
  n_valid <- length(valid_idx)

  if (verbose) {
    n_missing_chr <- sum(!chr_for_fasta %in% fa_chroms & !below_threshold)
    if (n_missing_chr > 0) {
      cat("  Regions with missing chromosomes:", n_missing_chr, "\n")
    }
    cat("  Valid regions for extraction:", n_valid, "\n")
  }

  # ---------------------------------------------------------------------------
  # Initialize result vectors
  # ---------------------------------------------------------------------------
  gc_percent <- rep(NA_real_, n_regions)
  at_percent <- rep(NA_real_, n_regions)
  n_Nbases <- rep(NA_integer_, n_regions)

  if (include_sequence) {
    sequences <- rep(NA_character_, n_regions)
  }

  if (include_repeats) {
    longest_homo <- rep(NA_integer_, n_regions)
    homo_base <- rep(NA_character_, n_regions)
    longest_di <- rep(NA_integer_, n_regions)
    di_motif <- rep(NA_character_, n_regions)
  }

  # ---------------------------------------------------------------------------
  # Extract sequences and calculate composition
  # ---------------------------------------------------------------------------
  if (n_valid > 0) {
    if (verbose) cat("  Extracting sequences...\n")

    # Create GRanges for valid regions (1-based for Bioconductor)
    # Use translated chromosome names (chr_for_fasta) for FASTA lookup
    gr <- GenomicRanges::GRanges(
      seqnames = chr_for_fasta[valid_idx],
      ranges = IRanges::IRanges(
        start = start_ext[valid_idx] + 1,  # Convert to 1-based
        end = end_ext[valid_idx]
      )
    )

    # Extract sequences
    seqs <- Biostrings::getSeq(fa, gr)

    if (verbose) cat("  Calculating composition...\n")

    # Process each sequence
    for (i in seq_along(valid_idx)) {
      idx <- valid_idx[i]
      seq_str <- as.character(seqs[[i]])
      seq_upper <- toupper(seq_str)
      seq_len <- nchar(seq_upper)

      # Store the sequence if requested
      if (include_sequence) {
        sequences[idx] <- seq_upper
      }

      # Count bases
      a_count <- lengths(regmatches(seq_upper, gregexpr("A", seq_upper)))
      t_count <- lengths(regmatches(seq_upper, gregexpr("T", seq_upper)))
      g_count <- lengths(regmatches(seq_upper, gregexpr("G", seq_upper)))
      c_count <- lengths(regmatches(seq_upper, gregexpr("C", seq_upper)))
      n_Nbases[idx] <- lengths(regmatches(seq_upper, gregexpr("N", seq_upper)))

      # Calculate percentages (excluding Ns from denominator)
      effective_len <- seq_len - n_Nbases[idx]
      if (effective_len > 0) {
        gc_percent[idx] <- 100 * (g_count + c_count) / effective_len
        at_percent[idx] <- 100 * (a_count + t_count) / effective_len
      }

      # Repeat detection
      if (include_repeats) {
        # Homopolymer detection
        homo_result <- .find_longest_homopolymer(seq_upper)
        longest_homo[idx] <- homo_result$length
        homo_base[idx] <- homo_result$base

        # Dinucleotide repeat detection
        di_result <- .find_longest_dinucleotide(seq_upper)
        longest_di[idx] <- di_result$length
        di_motif[idx] <- di_result$motif
      }

      # Progress indicator for large datasets
      if (verbose && i %% 1000 == 0) {
        cat("    Processed", i, "of", n_valid, "regions\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Build result data.frame
  # ---------------------------------------------------------------------------
  result <- data.frame(
    peak_id = peak_ids,
    width = widths,
    gc_percent = round(gc_percent, 2),
    at_percent = round(at_percent, 2),
    n_Nbases = n_Nbases,
    stringsAsFactors = FALSE
  )

  if (include_repeats) {
    result$longest_homopolymer <- longest_homo
    result$homopolymer_base <- homo_base
    result$longest_dinucleotide <- longest_di
    result$dinucleotide_motif <- di_motif
  }

  # Add sequence column if requested (at the end since it can be long)
  if (include_sequence) {
    result$sequence <- sequences
  }

  if (verbose) {
    valid_gc <- gc_percent[!is.na(gc_percent)]
    if (length(valid_gc) > 0) {
      cat("\nComposition summary (", length(valid_gc), " regions):\n", sep = "")
      cat("  GC%: mean =", round(mean(valid_gc), 1),
          ", median =", round(median(valid_gc), 1),
          ", range =", round(min(valid_gc), 1), "-", round(max(valid_gc), 1), "\n")
      cat("  AT%: mean =", round(mean(at_percent[!is.na(at_percent)]), 1), "\n")

      if (include_repeats) {
        valid_homo <- longest_homo[!is.na(longest_homo)]
        if (length(valid_homo) > 0) {
          cat("  Longest homopolymer: max =", max(valid_homo),
              ", mean =", round(mean(valid_homo), 1), "\n")
        }
        valid_di <- longest_di[!is.na(longest_di)]
        if (length(valid_di) > 0) {
          cat("  Longest dinucleotide repeat: max =", max(valid_di),
              "bp, mean =", round(mean(valid_di), 1), "bp\n")
        }
      }
    }
  }

  return(result)
}


#' Find Longest Homopolymer Run
#'
#' @description Internal helper function to find the longest run of consecutive
#' identical bases in a sequence.
#'
#' @param seq Character string. DNA sequence (uppercase).
#'
#' @return List with components:
#' \describe{
#'   \item{length}{Length of longest homopolymer}
#'   \item{base}{The base forming the longest run (A, T, G, C, or N)}
#' }
#'
#' @keywords internal
#'
.find_longest_homopolymer <- function(seq) {
  # Find runs of each base
  bases <- c("A", "T", "G", "C")
  max_len <- 0
  max_base <- NA_character_

  for (base in bases) {
    # Pattern for runs of this base
    pattern <- paste0(base, "+")
    matches <- regmatches(seq, gregexpr(pattern, seq))[[1]]

    if (length(matches) > 0) {
      longest <- max(nchar(matches))
      if (longest > max_len) {
        max_len <- longest
        max_base <- base
      }
    }
  }

  return(list(length = as.integer(max_len), base = max_base))
}


#' Find Longest Dinucleotide Repeat
#'
#' @description Internal helper function to find the longest dinucleotide
#' repeat in a sequence (e.g., ATATAT, CGCGCG).
#'
#' @param seq Character string. DNA sequence (uppercase).
#'
#' @return List with components:
#' \describe{
#'   \item{length}{Length in base pairs of the longest dinucleotide repeat}
#'   \item{motif}{The two-base motif (e.g., "AT", "CG")}
#' }
#'
#' @keywords internal
#'
.find_longest_dinucleotide <- function(seq) {
  # All possible dinucleotide motifs (excluding same-base like AA, TT)
  bases <- c("A", "T", "G", "C")
  motifs <- character()
  for (b1 in bases) {
    for (b2 in bases) {
      if (b1 != b2) {
        motifs <- c(motifs, paste0(b1, b2))
      }
    }
  }

  max_len <- 0
  max_motif <- NA_character_

  for (motif in motifs) {
    # Pattern: motif repeated 2+ times
    # e.g., for "AT": (AT){2,}
    pattern <- paste0("(", motif, "){2,}")
    matches <- regmatches(seq, gregexpr(pattern, seq))[[1]]

    if (length(matches) > 0) {
      longest <- max(nchar(matches))
      if (longest > max_len) {
        max_len <- longest
        max_motif <- motif
      }
    }
  }

  # Return 0 if no dinucleotide repeats found (need at least 2 repeats = 4bp)
  if (max_len < 4) {
    return(list(length = 0L, motif = NA_character_))
  }

  return(list(length = as.integer(max_len), motif = max_motif))
}
