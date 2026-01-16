#' Eleuthia - Peak Loading Functions
#'
#' @description Functions for loading and processing peak files from
#' ATAC-seq and ChIP-seq analyses (narrowPeak, broadPeak, BED formats).


#' Load a narrowPeak File
#'
#' @description Reads a MACS2 narrowPeak file into a data.frame with
#' standardized column names.
#'
#' @param file_path Character string. Path to the narrowPeak file.
#' @param min_score Numeric. Minimum peak score to retain (default = 0, keep all).
#' @param min_qvalue Numeric. Minimum -log10(qvalue) to retain (default = 0).
#'
#' @return A data.frame with columns: chr, start, end, name, score, strand,
#'   signal_value, pvalue, qvalue, peak_summit, and summit_pos (absolute position).
#'
#' @details
#' narrowPeak format (BED6+4):
#' \enumerate{
#'   \item chr - Chromosome
#'   \item start - Start position (0-based)
#'   \item end - End position
#'   \item name - Peak name
#'   \item score - Peak score (0-1000)
#'   \item strand - Strand (usually ".")
#'   \item signalValue - Fold enrichment
#'   \item pValue - -log10(p-value)
#'   \item qValue - -log10(q-value)
#'   \item peak - Offset from start to peak summit
#' }
#'
#' @export
#'
#' @examples
#' peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak")
#' peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak", min_qvalue = 2)
#'
ELEUTHIA_load_narrowpeak <- function(file_path,
                                      min_score = 0,
                                      min_qvalue = 0) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Read narrowPeak file (no header)
  peaks <- read.table(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
  )

  # Assign column names
  colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand",
                       "signal_value", "pvalue", "qvalue", "peak_summit")

  # Calculate absolute summit position
  peaks$summit_pos <- peaks$start + peaks$peak_summit

  # Apply filters
  if (min_score > 0) {
    peaks <- peaks[peaks$score >= min_score, ]
  }

  if (min_qvalue > 0) {
    peaks <- peaks[peaks$qvalue >= min_qvalue, ]
  }

  # Reset row names
  rownames(peaks) <- NULL

  return(peaks)
}


#' Load Multiple Peak Files from Sample Sheet
#'
#' @description Loads all peak files for a given omics type from a validated
#' sample sheet into a named list.
#'
#' @param sample_sheet A validated sample sheet data.frame with sample_id column.
#' @param omics Character string. Omics type to load ("ATACseq" or "CHIPseq").
#' @param min_score Numeric. Minimum peak score to retain (default = 0).
#' @param min_qvalue Numeric. Minimum -log10(qvalue) to retain (default = 0).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A named list of peak data.frames, with names corresponding to sample_id.
#'
#' @export
#'
#' @examples
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' atac_peaks <- ELEUTHIA_load_peaks_from_sheet(result$sample_sheet, "ATACseq")
#'
ELEUTHIA_load_peaks_from_sheet <- function(sample_sheet,
                                            omics,
                                            min_score = 0,
                                            min_qvalue = 0,
                                            verbose = TRUE) {

  # Get subset for this omics type
  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type: ", omics)
  }

  # Check format is appropriate
  valid_formats <- c("peaks", "bed")
  if (!all(subset_df$format %in% valid_formats)) {
    stop("All samples must have format 'peaks' or 'bed' for peak loading")
  }

  if (verbose) {
    cat("Loading", nrow(subset_df), omics, "peak files...\n")
  }

  # Load each file

  peak_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    file_path <- file.path(subset_df$file_loc[i], subset_df$file_name[i])

    if (verbose) {
      cat("  Loading:", sample_id, "\n")
    }

    # Determine file type from format or extension
    if (subset_df$format[i] == "peaks" || grepl("\\.narrowPeak$", file_path)) {
      peak_list[[sample_id]] <- ELEUTHIA_load_narrowpeak(
        file_path,
        min_score = min_score,
        min_qvalue = min_qvalue
      )
    } else {
      # Generic BED loading
      peak_list[[sample_id]] <- ELEUTHIA_load_bed(file_path)
    }
  }

  if (verbose) {
    cat("Loaded", length(peak_list), "peak files.\n")
    total_peaks <- sum(sapply(peak_list, nrow))
    cat("Total peaks across all samples:", total_peaks, "\n")
  }

  return(peak_list)
}


#' Load a BED File
#'
#' @description Reads a BED file into a data.frame with standardized column names.
#'
#' @param file_path Character string. Path to the BED file.
#'
#' @return A data.frame with at minimum columns: chr, start, end.
#'   Additional columns depend on BED format (BED3, BED6, etc.).
#'
#' @export
#'
ELEUTHIA_load_bed <- function(file_path) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Read BED file
  bed <- read.table(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
  )

  # Assign column names based on number of columns
  ncols <- ncol(bed)

  if (ncols >= 3) {
    colnames(bed)[1:3] <- c("chr", "start", "end")
  }
  if (ncols >= 4) {
    colnames(bed)[4] <- "name"
  }
  if (ncols >= 5) {
    colnames(bed)[5] <- "score"
  }
  if (ncols >= 6) {
    colnames(bed)[6] <- "strand"
  }

  return(bed)
}


#' Create Consensus Peak Set
#'
#' @description Merges overlapping peaks across multiple samples to create
#' a consensus peak set for quantification.
#'
#' @param peak_list Named list of peak data.frames (output from
#'   ELEUTHIA_load_peaks_from_sheet).
#' @param min_overlap Integer. Minimum number of samples a peak must appear
#'   in to be included in consensus (default = 1, include all).
#' @param merge_distance Integer. Maximum distance (bp) between peaks to merge
#'   them into a single region (default = 0, only merge overlapping peaks).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, peak_id, n_samples.
#'   Each row represents a consensus peak region.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Combines all peaks from all samples
#'   \item Sorts by chromosome and position
#'   \item Merges overlapping/adjacent peaks (within merge_distance)
#'   \item Counts how many samples contributed to each merged region
#'   \item Filters by min_overlap threshold
#'   \item Assigns unique peak_id to each consensus region
#' }
#'
#' @export
#'
#' @examples
#' peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' consensus <- ELEUTHIA_create_consensus_peaks(peak_list, min_overlap = 2)
#'
ELEUTHIA_create_consensus_peaks <- function(peak_list,
                                             min_overlap = 1,
                                             merge_distance = 0,
                                             verbose = TRUE) {

  if (length(peak_list) == 0) {
    stop("peak_list is empty")
  }

  if (verbose) {
    cat("Creating consensus peak set from", length(peak_list), "samples...\n")
  }

  # Combine all peaks with sample source
  all_peaks <- do.call(rbind, lapply(names(peak_list), function(sample) {
    df <- peak_list[[sample]][, c("chr", "start", "end")]
    df$sample <- sample
    return(df)
  }))

  if (verbose) {
    cat("  Total peaks before merging:", nrow(all_peaks), "\n")
  }

  # Sort by chromosome and position
  all_peaks <- all_peaks[order(all_peaks$chr, all_peaks$start), ]

  # Merge overlapping peaks using a simple iterative approach
  # Group by chromosome first
  chromosomes <- unique(all_peaks$chr)
  consensus_list <- list()

  for (chrom in chromosomes) {
    chr_peaks <- all_peaks[all_peaks$chr == chrom, ]

    if (nrow(chr_peaks) == 0) next

    # Initialize merged regions
    merged <- data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      samples = character(),
      stringsAsFactors = FALSE
    )

    current_start <- chr_peaks$start[1]
    current_end <- chr_peaks$end[1]
    current_samples <- chr_peaks$sample[1]

    for (i in seq_len(nrow(chr_peaks))[-1]) {
      peak_start <- chr_peaks$start[i]
      peak_end <- chr_peaks$end[i]
      peak_sample <- chr_peaks$sample[i]

      # Check if this peak overlaps with current region (within merge_distance)
      if (peak_start <= current_end + merge_distance) {
        # Extend current region
        current_end <- max(current_end, peak_end)
        current_samples <- paste(current_samples, peak_sample, sep = ",")
      } else {
        # Save current region and start new one
        merged <- rbind(merged, data.frame(
          chr = chrom,
          start = current_start,
          end = current_end,
          samples = current_samples,
          stringsAsFactors = FALSE
        ))

        current_start <- peak_start
        current_end <- peak_end
        current_samples <- peak_sample
      }
    }

    # Don't forget the last region
    merged <- rbind(merged, data.frame(
      chr = chrom,
      start = current_start,
      end = current_end,
      samples = current_samples,
      stringsAsFactors = FALSE
    ))

    consensus_list[[chrom]] <- merged
  }

  # Combine all chromosomes
  consensus <- do.call(rbind, consensus_list)
  rownames(consensus) <- NULL

  # Count unique samples per region
  consensus$n_samples <- sapply(strsplit(consensus$samples, ","), function(x) {
    length(unique(x))
  })

  # Filter by min_overlap
  if (min_overlap > 1) {
    n_before <- nrow(consensus)
    consensus <- consensus[consensus$n_samples >= min_overlap, ]
    if (verbose) {
      cat("  Filtered from", n_before, "to", nrow(consensus),
          "peaks (min_overlap =", min_overlap, ")\n")
    }
  }

  # Create peak IDs
  consensus$peak_id <- paste0("peak_", seq_len(nrow(consensus)))

  # Remove samples column (was just for counting)
  consensus$samples <- NULL

  # Reorder columns
  consensus <- consensus[, c("chr", "start", "end", "peak_id", "n_samples")]

  if (verbose) {
    cat("  Final consensus peaks:", nrow(consensus), "\n")
    cat("  Peaks per chromosome:\n")
    chr_counts <- table(consensus$chr)
    # Show top 5 chromosomes
    top_chrs <- head(sort(chr_counts, decreasing = TRUE), 5)
    for (ch in names(top_chrs)) {
      cat("    ", ch, ":", top_chrs[ch], "\n")
    }
  }

  return(consensus)
}


#' Call Regions from Fragment BED Files (Sparse Data)
#'
#' @description Identifies enriched regions from fragment BED files when
#' traditional peak callers fail due to low binding event counts.
#' Designed for factors with sparse binding (e.g., ~100-500 expected sites).
#'
#' @param bed_files Character vector of paths to BED files, or a named list
#'   of loaded BED data.frames.
#' @param merge_distance Integer. Maximum distance (bp) between fragments to
#'   merge into a single region (default = 150).
#' @param min_fragments Integer. Minimum number of fragments required in a
#'   region to call it as a "peak" (default = 5).
#' @param extend Integer. Extend each fragment by this many bp on each side
#'   before merging (default = 75). Helps connect nearby fragments.
#' @param pool_samples Logical. If TRUE, pool all samples before calling
#'   regions (recommended for sparse data). If FALSE, call per-sample and
#'   take union (default = TRUE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, peak_id, n_fragments,
#'   width. Each row represents a called region.
#'
#' @details
#' This function provides a simple alternative to statistical peak callers
#' for scenarios where:
#' \itemize{
#'   \item The factor has few binding sites (< 500 expected)
#'   \item Traditional peak callers cannot build reliable background models
#'   \item Visual inspection confirms signal at expected locations
#' }
#'
#' The algorithm:
#' \enumerate{
#'   \item Optionally pools fragments from all samples
#'   \item Extends each fragment by \code{extend} bp on each side
#'   \item Merges overlapping/nearby extended fragments (within \code{merge_distance})
#'   \item Counts original fragments falling within each merged region
#'   \item Filters regions by \code{min_fragments} threshold
#' }
#'
#' For CUT&TAG/ChIP-seq data with sparse binding, pooling samples is
#' recommended as it increases power to detect true binding sites.
#'
#' @export
#'
#' @examples
#' # Call regions from multiple BED files
#' bed_files <- c("sample1.bed", "sample2.bed", "sample3.bed")
#' regions <- ELEUTHIA_call_regions_from_fragments(bed_files, min_fragments = 5)
#'
#' # More stringent calling
#' regions <- ELEUTHIA_call_regions_from_fragments(bed_files,
#'                                                  merge_distance = 100,
#'                                                  min_fragments = 10)
#'
ELEUTHIA_call_regions_from_fragments <- function(bed_files,
                                                   merge_distance = 150,
                                                   min_fragments = 5,
                                                   extend = 75,
                                                   pool_samples = TRUE,
                                                   verbose = TRUE) {

  # Load BED files if paths provided
  if (is.character(bed_files)) {
    if (verbose) cat("Loading", length(bed_files), "BED files...\n")

    bed_list <- lapply(bed_files, function(f) {
      if (!file.exists(f)) stop("File not found: ", f)
      ELEUTHIA_load_bed(f)
    })
    names(bed_list) <- basename(bed_files)
  } else if (is.list(bed_files)) {
    bed_list <- bed_files
  } else {
    stop("bed_files must be a character vector of paths or a list of data.frames")
  }

  if (verbose) {
    total_frags <- sum(sapply(bed_list, nrow))
    cat("  Total fragments:", format(total_frags, big.mark = ","), "\n")
  }

  # Pool or process separately
  if (pool_samples) {
    if (verbose) cat("Pooling fragments from all samples...\n")

    all_frags <- do.call(rbind, lapply(bed_list, function(df) {
      df[, c("chr", "start", "end")]
    }))
  } else {
    # For non-pooled, we'd call per sample and merge - implement if needed
    stop("Non-pooled calling not yet implemented. Use pool_samples = TRUE")
  }

  if (verbose) {
    cat("Calling regions (merge_distance =", merge_distance,
        ", min_fragments =", min_fragments, ")...\n")
  }

  # Extend fragments
  all_frags$start_ext <- pmax(0, all_frags$start - extend)
  all_frags$end_ext <- all_frags$end + extend

  # Sort by chromosome and position
  all_frags <- all_frags[order(all_frags$chr, all_frags$start_ext), ]

  # Merge overlapping extended regions per chromosome
  chromosomes <- unique(all_frags$chr)
  regions_list <- list()

  for (chrom in chromosomes) {
    chr_frags <- all_frags[all_frags$chr == chrom, ]

    if (nrow(chr_frags) == 0) next

    # Initialize
    merged <- data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      n_fragments = integer(),
      stringsAsFactors = FALSE
    )

    current_start <- chr_frags$start_ext[1]
    current_end <- chr_frags$end_ext[1]
    current_count <- 1

    for (i in seq_len(nrow(chr_frags))[-1]) {
      frag_start <- chr_frags$start_ext[i]
      frag_end <- chr_frags$end_ext[i]

      # Check if this fragment overlaps/is near current region
      if (frag_start <= current_end + merge_distance) {
        # Extend region
        current_end <- max(current_end, frag_end)
        current_count <- current_count + 1
      } else {
        # Save current region if it meets threshold
        if (current_count >= min_fragments) {
          merged <- rbind(merged, data.frame(
            chr = chrom,
            start = current_start,
            end = current_end,
            n_fragments = current_count,
            stringsAsFactors = FALSE
          ))
        }

        # Start new region
        current_start <- frag_start
        current_end <- frag_end
        current_count <- 1
      }
    }

    # Don't forget last region
    if (current_count >= min_fragments) {
      merged <- rbind(merged, data.frame(
        chr = chrom,
        start = current_start,
        end = current_end,
        n_fragments = current_count,
        stringsAsFactors = FALSE
      ))
    }

    regions_list[[chrom]] <- merged
  }

  # Combine all chromosomes
  regions <- do.call(rbind, regions_list)
  rownames(regions) <- NULL

  if (nrow(regions) == 0) {
    warning("No regions found meeting the criteria. ",
            "Try lowering min_fragments or increasing extend/merge_distance.")
    return(data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      peak_id = character(),
      n_fragments = integer(),
      width = integer()
    ))
  }

  # Add peak IDs and width
  regions$peak_id <- paste0("region_", seq_len(nrow(regions)))
  regions$width <- regions$end - regions$start

  # Reorder columns
  regions <- regions[, c("chr", "start", "end", "peak_id", "n_fragments", "width")]

  if (verbose) {
    cat("\nRegion calling complete:\n")
    cat("  Regions found:", nrow(regions), "\n")
    cat("  Mean width:", round(mean(regions$width)), "bp\n")
    cat("  Median fragments per region:", median(regions$n_fragments), "\n")
    cat("  Total fragments in regions:",
        format(sum(regions$n_fragments), big.mark = ","), "\n")
  }

  return(regions)
}


#' Load ChIP/CUT&TAG BED Files from Sample Sheet
#'
#' @description Loads fragment BED files for ChIP-seq or CUT&TAG data
#' from a validated sample sheet.
#'
#' @param sample_sheet A validated sample sheet data.frame.
#' @param omics Character string. Omics type to load (default = "CHIPseq").
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A named list of BED data.frames, with names corresponding to sample_id.
#'
#' @export
#'
ELEUTHIA_load_bed_from_sheet <- function(sample_sheet,
                                          omics = "CHIPseq",
                                          verbose = TRUE) {

  # Get subset for this omics type
  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type: ", omics)
  }

  if (verbose) {
    cat("Loading", nrow(subset_df), omics, "BED files...\n")
  }

  bed_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    file_path <- file.path(subset_df$file_loc[i], subset_df$file_name[i])

    if (verbose) {
      cat("  Loading:", sample_id, "\n")
    }

    bed_list[[sample_id]] <- ELEUTHIA_load_bed(file_path)
  }

  if (verbose) {
    total_frags <- sum(sapply(bed_list, nrow))
    cat("Loaded", length(bed_list), "BED files.\n")
    cat("Total fragments:", format(total_frags, big.mark = ","), "\n")
  }

  return(bed_list)
}


#' Convert Consensus Peaks to SAF Format
#'
#' @description Converts a consensus peak data.frame to SAF (Simplified
#' Annotation Format) for use with Rsubread::featureCounts().
#'
#' @param consensus_peaks A consensus peak data.frame from
#'   ELEUTHIA_create_consensus_peaks().
#'
#' @return A data.frame in SAF format with columns: GeneID, Chr, Start, End, Strand.
#'
#' @details
#' SAF format is required by Rsubread::featureCounts() for counting reads
#' in custom genomic regions. The "GeneID" column contains the peak_id.
#'
#' @export
#'
#' @examples
#' consensus <- ELEUTHIA_create_consensus_peaks(peak_list)
#' saf <- ELEUTHIA_peaks_to_saf(consensus)
#' counts <- Rsubread::featureCounts(bam_files, annot.ext = saf, isGTFAnnotationFile = FALSE)
#'
ELEUTHIA_peaks_to_saf <- function(consensus_peaks) {

  saf <- data.frame(
    GeneID = consensus_peaks$peak_id,
    Chr = consensus_peaks$chr,
    Start = consensus_peaks$start + 1,  # SAF is 1-based
    End = consensus_peaks$end,
    Strand = "*",
    stringsAsFactors = FALSE
  )

  return(saf)
}
