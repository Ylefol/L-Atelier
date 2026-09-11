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
#' \dontrun{
#' peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak")
#' peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak", min_qvalue = 2)
#'
#' }
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
#' \dontrun{
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' atac_peaks <- ELEUTHIA_load_peaks_from_sheet(result$sample_sheet, "ATACseq")
#'
#' }
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
  valid_formats <- c("peaks", "bed", "seacr")
  if (!all(subset_df$format %in% valid_formats)) {
    stop("All samples must have format 'peaks', 'bed', or 'seacr' for peak loading")
  }

  if (verbose) {
    cat("[ELEUTHIA] Loading", nrow(subset_df), omics, "peak files...\n")
  }

  # Load each file

  peak_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    file_path <- file.path(subset_df$file_loc[i], subset_df$file_name[i])

    if (verbose) {
      cat("    Loading:", sample_id, "\n")
    }

    # Determine file type from format or extension
    if (subset_df$format[i] == "seacr" || grepl("_seacr_peaks\\.bed$", file_path)) {
      peak_list[[sample_id]] <- .eleuthia_load_seacr_peaks(file_path)
    } else if (subset_df$format[i] == "peaks" || grepl("\\.narrowPeak$", file_path)) {
      peak_list[[sample_id]] <- ELEUTHIA_load_narrowpeak(
        file_path,
        min_score = min_score,
        min_qvalue = min_qvalue
      )
    } else {
      peak_list[[sample_id]] <- ELEUTHIA_load_bed(file_path)
    }
  }

  if (verbose) {
    cat("[ELEUTHIA] Loaded", length(peak_list), "peak files.\n")
    total_peaks <- sum(sapply(peak_list, nrow))
    cat("[ELEUTHIA] Total peaks across all samples:", total_peaks, "\n")
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
.eleuthia_load_seacr_peaks <- function(file_path) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  peaks <- read.table(file_path, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE,
                      col.names = c("peak_id", "chr", "start", "end",
                                    "total_signal", "max_signal"))
  return(peaks)
}


#' Load a BED File
#'
#' @description Reads a BED format file (BED3, BED6, narrowPeak, broadPeak, etc.)
#'   into a data.frame with standardized column names. The first three columns are
#'   always named \code{chr}, \code{start}, \code{end}; additional columns are named
#'   \code{name}, \code{score}, \code{strand} where present.
#'
#' @param file_path Character string. Path to the BED file.
#'
#' @return A data.frame with at minimum columns: \code{chr}, \code{start}, \code{end}.
#'   Coordinates are returned as-read (BED 0-based half-open intervals).
#'
#' @export
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
#' @param sample_groups Named character vector or NULL. Maps sample IDs (names)
#'   to group labels (values). When provided, a \code{groups} column is added to
#'   the output containing sorted, comma-separated unique group labels that
#'   contributed to each consensus peak (default = NULL, no group column added).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, peak_id, n_samples, and
#'   optionally groups (when \code{sample_groups} is supplied).
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
#' \dontrun{
#' peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' consensus <- ELEUTHIA_create_consensus_peaks(peak_list, min_overlap = 2)
#'
#' }
ELEUTHIA_create_consensus_peaks <- function(peak_list,
                                             min_overlap = 1,
                                             merge_distance = 0,
                                             sample_groups = NULL,
                                             verbose = TRUE) {

  if (length(peak_list) == 0) {
    stop("peak_list is empty")
  }

  if (verbose) {
    cat("[ELEUTHIA] Creating consensus peak set from", length(peak_list), "samples...\n")
  }

  # Combine all peaks with sample source
  all_peaks <- do.call(rbind, lapply(names(peak_list), function(sample) {
    df <- peak_list[[sample]][, c("chr", "start", "end")]
    df$sample <- sample
    return(df)
  }))

  if (verbose) {
    cat("    Total peaks before merging:", nrow(all_peaks), "\n")
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
      cat("    Filtered from", n_before, "to", nrow(consensus),
          "peaks (min_overlap =", min_overlap, ")\n")
    }
  }

  # Create peak IDs
  consensus$peak_id <- paste0("peak_", seq_len(nrow(consensus)))

  # Optionally map samples → group labels before dropping the samples column
  if (!is.null(sample_groups)) {
    consensus$groups <- sapply(strsplit(consensus$samples, ","), function(samps) {
      grps <- unique(sample_groups[unique(samps)])
      grps <- grps[!is.na(grps)]
      paste(sort(grps), collapse = ",")
    })
  }

  consensus$samples <- NULL

  col_order <- c("chr", "start", "end", "peak_id", "n_samples")
  if (!is.null(sample_groups)) col_order <- c(col_order, "groups")
  consensus <- consensus[, col_order]

  if (verbose) {
    cat("    Final consensus peaks:", nrow(consensus), "\n")
    cat("    Chromosomes with peaks:", length(unique(consensus$chr)), "\n")
  }

  return(consensus)
}


#' Merge Fragments into Candidate Regions
#'
#' @description Pools and merges fragments from BED files into candidate regions
#' with a baseline noise filter. This is the first step of a two-step workflow
#' where you can then examine the distribution and choose a final threshold.
#'
#' @param bed_files Character vector of paths to BED files, or a named list
#'   of loaded BED data.frames.
#' @param merge_distance Integer. Maximum distance (bp) between fragments to
#'   merge into a single region (default = 150).
#' @param min_fragments Integer. Minimum fragments to retain a candidate region
#'   (default = 10). This is a noise filter - regions below this are discarded
#'   during merging. Set to 1 to keep all regions.
#' @param extend Integer. Extend each fragment by this many bp on each side
#'   before merging (default = 75). Helps connect nearby fragments.
#' @param chrom_sizes Optional. Chromosome sizes for clamping extended coordinates.
#'   Prevents regions from extending beyond chromosome boundaries. Can be:
#'   \itemize{
#'     \item A named numeric vector where names are chromosome names and values
#'       are chromosome lengths
#'     \item A data.frame with columns 'chr' and 'size'
#'     \item NULL (default) - no clamping, extended coordinates may exceed
#'       chromosome length
#'   }
#' @param pool_samples Logical. If TRUE, pool all samples before merging
#'   (recommended for sparse data). If FALSE, not yet implemented.
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, n_fragments, width.
#'   Contains candidate regions passing the min_fragments noise filter.
#'   Use ELEUTHIA_select_regions() for further quantile-based filtering.
#'
#' @details
#' This function performs the merging step with a baseline noise filter:
#' \enumerate{
#'   \item Pools fragments from all samples (if pool_samples = TRUE)
#'   \item Extends each fragment by \code{extend} bp on each side
#'   \item Merges overlapping/nearby extended fragments (within \code{merge_distance})
#'   \item Counts fragments in each merged region
#'   \item Discards regions with fewer than \code{min_fragments} (noise filter)
#' }
#'
#' The min_fragments parameter removes obvious noise (regions with very few
#' fragments) to make the quantile distribution more meaningful. This is
#' different from the final selection threshold in ELEUTHIA_select_regions().
#'
#' @seealso \code{\link{ELEUTHIA_select_regions}} for quantile-based filtering,
#'   \code{\link{ELEUTHIA_plot_region_distribution}} for visualization,
#'   \code{\link{ELEUTHIA_call_regions_from_fragments}} for one-step workflow
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Two-step workflow:
#' # 1. Merge fragments into candidate regions (with noise filter)
#' candidates <- ELEUTHIA_merge_fragments(bed_files, min_fragments = 10)
#'
#' # 2. Examine distribution of candidates
#' ELEUTHIA_plot_region_distribution(candidates)
#'
#' # 3. Select top regions based on distribution
#' regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)
#'
#' }
ELEUTHIA_merge_fragments <- function(bed_files,
                                      merge_distance = 150,
                                      min_fragments = 10,
                                      extend = 75,
                                      chrom_sizes = NULL,
                                      pool_samples = TRUE,
                                      verbose = TRUE) {

  # Load BED files if paths provided
  if (is.character(bed_files)) {
    if (verbose) cat("[ELEUTHIA] Loading", length(bed_files), "BED files...\n")

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

  total_frags <- sum(sapply(bed_list, nrow))
  if (verbose) {
    cat("[ELEUTHIA]   Total fragments:", format(total_frags, big.mark = ","), "\n")
  }

  # Pool or process separately
  if (pool_samples) {
    if (verbose) cat("[ELEUTHIA] Pooling fragments from all samples...\n")

    all_frags <- do.call(rbind, lapply(seq_along(bed_list), function(i) {
      df <- bed_list[[i]][, c("chr", "start", "end")]
      df$sample_idx <- i
      df
    }))
  } else {
    stop("Non-pooled calling not yet implemented. Use pool_samples = TRUE")
  }

  if (verbose) {
    cat("[ELEUTHIA] Merging regions (merge_distance =", merge_distance,
        ", extend =", extend, ", min_fragments =", min_fragments, ")...\n")
  }


  # Process chromosome sizes if provided
  chrom_size_vec <- NULL
  if (!is.null(chrom_sizes)) {
    if (is.data.frame(chrom_sizes)) {
      if (!all(c("chr", "size") %in% colnames(chrom_sizes))) {
        stop("chrom_sizes data.frame must have 'chr' and 'size' columns")
      }
      chrom_size_vec <- setNames(chrom_sizes$size, chrom_sizes$chr)
    } else if (is.numeric(chrom_sizes) && !is.null(names(chrom_sizes))) {
      chrom_size_vec <- chrom_sizes
    } else {
      stop("chrom_sizes must be a named numeric vector or data.frame with 'chr' and 'size' columns")
    }
    if (verbose) {
      cat("[ELEUTHIA]   Chromosome sizes provided for", length(chrom_size_vec), "chromosomes\n")
    }
  }

  # Extend fragments
  all_frags$start_ext <- pmax(0, all_frags$start - extend)
  all_frags$end_ext <- all_frags$end + extend

  # Clamp end coordinates to chromosome sizes if provided
  if (!is.null(chrom_size_vec)) {
    for (chr_name in unique(all_frags$chr)) {
      if (chr_name %in% names(chrom_size_vec)) {
        chr_idx <- all_frags$chr == chr_name
        chr_len <- chrom_size_vec[chr_name]
        all_frags$end_ext[chr_idx] <- pmin(all_frags$end_ext[chr_idx], chr_len)
      }
    }
    n_clamped <- sum(all_frags$end_ext < (all_frags$end + extend))
    if (verbose && n_clamped > 0) {
      cat("[ELEUTHIA]   Fragments clamped at chromosome boundaries:", n_clamped, "\n")
    }
  }

  # Sort by chromosome and position
  all_frags <- all_frags[order(all_frags$chr, all_frags$start_ext), ]

  # Merge overlapping extended regions per chromosome
  # Use list accumulation instead of rbind in loop for performance
  n_samp        <- length(bed_list)
  chromosomes   <- unique(all_frags$chr)
  regions_list  <- list()

  for (chrom in chromosomes) {
    chr_frags <- all_frags[all_frags$chr == chrom, ]

    if (nrow(chr_frags) == 0) next

    # Pre-allocate vectors (upper bound is number of fragments)
    n_max          <- nrow(chr_frags)
    starts         <- integer(n_max)
    ends           <- integer(n_max)
    counts         <- integer(n_max)
    sample_counts  <- integer(n_max)
    region_idx     <- 0

    current_start      <- chr_frags$start_ext[1]
    current_end        <- chr_frags$end_ext[1]
    current_count      <- 1
    current_sample_set <- logical(n_samp)
    current_sample_set[chr_frags$sample_idx[1]] <- TRUE

    for (i in seq_len(nrow(chr_frags))[-1]) {
      frag_start <- chr_frags$start_ext[i]
      frag_end   <- chr_frags$end_ext[i]

      # Check if this fragment overlaps/is near current region
      if (frag_start <= current_end + merge_distance) {
        # Extend region
        current_end <- max(current_end, frag_end)
        current_count <- current_count + 1
        current_sample_set[chr_frags$sample_idx[i]] <- TRUE
      } else {
        # Save current region if it passes noise filter
        if (current_count >= min_fragments) {
          region_idx <- region_idx + 1
          starts[region_idx]        <- current_start
          ends[region_idx]          <- current_end
          counts[region_idx]        <- current_count
          sample_counts[region_idx] <- sum(current_sample_set)
        }

        # Start new region
        current_start      <- frag_start
        current_end        <- frag_end
        current_count      <- 1
        current_sample_set <- logical(n_samp)
        current_sample_set[chr_frags$sample_idx[i]] <- TRUE
      }
    }

    # Don't forget last region (with filter)
    if (current_count >= min_fragments) {
      region_idx <- region_idx + 1
      starts[region_idx]        <- current_start
      ends[region_idx]          <- current_end
      counts[region_idx]        <- current_count
      sample_counts[region_idx] <- sum(current_sample_set)
    }

    # Create data.frame from pre-allocated vectors (single allocation)
    if (region_idx > 0) {
      regions_list[[chrom]] <- data.frame(
        chr        = rep(chrom, region_idx),
        start      = starts[1:region_idx],
        end        = ends[1:region_idx],
        n_fragments = counts[1:region_idx],
        n_samples  = sample_counts[1:region_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all chromosomes
  regions <- do.call(rbind, regions_list)

  if (is.null(regions) || nrow(regions) == 0) {
    warning("No regions created. Check input data.")
    return(data.frame(
      chr        = character(),
      start      = integer(),
      end        = integer(),
      n_fragments = integer(),
      n_samples  = integer(),
      width      = integer()
    ))
  }
  rownames(regions) <- NULL

  # Add width
  regions$width <- regions$end - regions$start

  if (verbose) {
    cat("[ELEUTHIA] Merging complete:\n")
    cat("    Candidate regions:", nrow(regions), "\n")
    cat("    Mean width:", round(mean(regions$width)), "bp\n")
    cat("    Fragment count range:", min(regions$n_fragments), "-",
        max(regions$n_fragments), "\n")
    cat("    Median fragments:", median(regions$n_fragments), "\n")
    cat("    Sample coverage: regions in all", n_samp, "samples:",
        sum(regions$n_samples == n_samp), "\n")
    cat("    Total fragments in regions:",
        format(sum(regions$n_fragments), big.mark = ","),
        sprintf("(%.1f%% of input)\n", 100 * sum(regions$n_fragments) / total_frags))
  }

  return(regions)
}


#' Select Regions Based on Fragment Count or Density Threshold
#'
#' @description Filters candidate regions from ELEUTHIA_merge_fragments() based
#' on fragment count or fragment density, with optional maximum width and
#' minimum sample support filters.
#'
#' @param candidate_regions Data.frame from ELEUTHIA_merge_fragments() with
#'   columns: chr, start, end, n_fragments, n_samples.
#' @param min_fragments Integer or NULL. Absolute minimum fragment count to keep
#'   a region. Ignored when \code{by_density = TRUE}; use \code{min_density}
#'   instead. If NULL, only quantile_threshold is used.
#' @param quantile_threshold Numeric between 0 and 1, or NULL. Keep regions
#'   above this quantile. When \code{by_density = FALSE} (default) the quantile
#'   is computed on raw fragment count; when \code{by_density = TRUE} it is
#'   computed on fragment density (fragments per bp). E.g., 0.90 keeps the top
#'   10% of regions. If NULL, only the absolute threshold is used.
#' @param max_width Integer or NULL. Maximum region width in bp. Regions wider
#'   than this are removed before any other filter. Useful for discarding broad,
#'   diffuse accumulations that inflate raw fragment counts without genuine
#'   local enrichment. Default NULL (no cap).
#' @param min_samples Integer or NULL. Minimum number of samples that must
#'   have contributed at least one fragment to a region for it to be retained.
#'   Requires the \code{n_samples} column produced by
#'   \code{ELEUTHIA_merge_fragments()}. This filter is applied after
#'   \code{max_width} and before the fragment count/density threshold, so it
#'   can substantially reduce the candidate pool before quantile computation.
#'   Default NULL (no sample support requirement).
#' @param by_density Logical. When TRUE, threshold filtering is performed on
#'   fragment density (n_fragments / width in bp) rather than raw fragment count.
#'   This normalises for region width so that narrow, concentrated peaks are
#'   preferred over broad, diffuse ones with equivalent total counts.
#'   Default FALSE (backward-compatible behaviour).
#' @param min_density Numeric or NULL. Absolute minimum fragment density
#'   (fragments per bp) to keep a region. Only used when
#'   \code{by_density = TRUE}. If NULL, only quantile_threshold is used.
#' @param verbose Logical. Print selection summary (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, peak_id, n_fragments,
#'   n_samples, fragment_density, width. Only regions passing all filters are
#'   included. \code{fragment_density} is always returned for inspection.
#'
#' @details
#' Filters are applied in the following order:
#' \enumerate{
#'   \item \code{max_width}: remove regions wider than the cap
#'   \item \code{min_samples}: require fragment contribution from at least N samples
#'   \item Fragment count or density threshold (whichever of the absolute and
#'     quantile thresholds is more stringent)
#' }
#'
#' Applying \code{min_samples} before the quantile threshold is intentional:
#' the quantile is then computed only on the biologically supported subset,
#' preventing low-support regions from diluting the distribution.
#'
#' \strong{Choosing between count and density mode:}
#' Raw fragment count (\code{by_density = FALSE}) reflects total signal
#' accumulation and is appropriate when region widths are comparable. Fragment
#' density (\code{by_density = TRUE}) normalises for region width and better
#' matches RPKM-normalised BigWig tracks, making it preferable when the merging
#' step produces regions of highly variable size.
#'
#' @seealso \code{\link{ELEUTHIA_merge_fragments}} for the merging step,
#'   \code{\link{ELEUTHIA_plot_region_distribution}} to visualize before selecting
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Merge first
#' candidates <- ELEUTHIA_merge_fragments(bed_files)
#'
#' # Select by absolute fragment count (original behaviour)
#' regions <- ELEUTHIA_select_regions(candidates, min_fragments = 100)
#'
#' # Require signal in at least 3 of N samples + top 10% by density
#' regions <- ELEUTHIA_select_regions(candidates,
#'                                     quantile_threshold = 0.90,
#'                                     min_samples        = 3,
#'                                     by_density         = TRUE)
#'
#' # Full combination: width cap + sample support + density
#' regions <- ELEUTHIA_select_regions(candidates,
#'                                     quantile_threshold = 0.90,
#'                                     max_width          = 2000,
#'                                     min_samples        = 3,
#'                                     by_density         = TRUE)
#' }
ELEUTHIA_select_regions <- function(candidate_regions,
                                     min_fragments      = NULL,
                                     quantile_threshold = NULL,
                                     max_width          = NULL,
                                     min_samples        = NULL,
                                     by_density         = FALSE,
                                     min_density        = NULL,
                                     verbose            = TRUE) {

  # Validate inputs
  if (!is.data.frame(candidate_regions)) {
    stop("candidate_regions must be a data.frame")
  }

  if (!"n_fragments" %in% colnames(candidate_regions)) {
    stop("candidate_regions must have 'n_fragments' column")
  }

  if (!is.null(min_samples)) {
    if (!"n_samples" %in% colnames(candidate_regions)) {
      stop("min_samples requires an 'n_samples' column in candidate_regions. ",
           "Ensure ELEUTHIA_merge_fragments() was used to produce the input.")
    }
    if (!is.numeric(min_samples) || min_samples < 1) {
      stop("min_samples must be a positive integer")
    }
  }

  if (by_density) {
    if (is.null(min_density) && is.null(quantile_threshold)) {
      stop("When by_density = TRUE, at least one of min_density or ",
           "quantile_threshold must be provided")
    }
  } else {
    if (is.null(min_fragments) && is.null(quantile_threshold)) {
      stop("At least one of min_fragments or quantile_threshold must be provided")
    }
  }

  if (!is.null(quantile_threshold) &&
      (quantile_threshold < 0 || quantile_threshold > 1)) {
    stop("quantile_threshold must be between 0 and 1")
  }

  if (!is.null(max_width) && (!is.numeric(max_width) || max_width <= 0)) {
    stop("max_width must be a positive number")
  }

  n_candidates <- nrow(candidate_regions)

  if (n_candidates == 0) {
    warning("No candidate regions provided")
    return(data.frame(
      chr              = character(),
      start            = integer(),
      end              = integer(),
      peak_id          = character(),
      n_fragments      = integer(),
      n_samples        = integer(),
      fragment_density = numeric(),
      width            = integer()
    ))
  }

  # Ensure width and density are present
  if (!"width" %in% colnames(candidate_regions)) {
    candidate_regions$width <- candidate_regions$end - candidate_regions$start
  }
  candidate_regions$fragment_density <- candidate_regions$n_fragments /
    candidate_regions$width

  # Ensure n_samples column exists (may be absent if input is not from
  # ELEUTHIA_merge_fragments, e.g. legacy data)
  has_n_samples <- "n_samples" %in% colnames(candidate_regions)

  # --- Step 1: max_width filter ---
  if (!is.null(max_width)) {
    n_before <- nrow(candidate_regions)
    candidate_regions <- candidate_regions[candidate_regions$width <= max_width, ]
    n_removed <- n_before - nrow(candidate_regions)
    if (verbose) {
      cat("[ELEUTHIA] Width filter (max_width =", max_width, "bp):\n")
      cat("    Removed:", n_removed, "regions",
          sprintf("(%.1f%% of candidates)\n", 100 * n_removed / n_candidates))
      cat("    Remaining:", nrow(candidate_regions), "\n\n")
    }
    if (nrow(candidate_regions) == 0) {
      warning("No regions remain after max_width filter. ",
              "Consider increasing max_width.")
      return(data.frame(
        chr              = character(),
        start            = integer(),
        end              = integer(),
        peak_id          = character(),
        n_fragments      = integer(),
        n_samples        = integer(),
        fragment_density = numeric(),
        width            = integer()
      ))
    }
  }

  # --- Step 2: min_samples filter ---
  if (!is.null(min_samples) && has_n_samples) {
    n_before <- nrow(candidate_regions)
    candidate_regions <- candidate_regions[
      candidate_regions$n_samples >= min_samples, ]
    n_removed <- n_before - nrow(candidate_regions)
    if (verbose) {
      cat("[ELEUTHIA] Sample support filter (min_samples =", min_samples, "):\n")
      cat("    Removed:", n_removed, "regions",
          sprintf("(%.1f%% of candidates)\n", 100 * n_removed / n_candidates))
      cat("    Remaining:", nrow(candidate_regions), "\n\n")
    }
    if (nrow(candidate_regions) == 0) {
      warning("No regions remain after min_samples filter. ",
              "Consider lowering min_samples.")
      return(data.frame(
        chr              = character(),
        start            = integer(),
        end              = integer(),
        peak_id          = character(),
        n_fragments      = integer(),
        n_samples        = integer(),
        fragment_density = numeric(),
        width            = integer()
      ))
    }
  }

  # --- Step 3: fragment count or density threshold ---
  if (by_density) {
    metric          <- candidate_regions$fragment_density
    abs_threshold   <- if (!is.null(min_density)) min_density else 0
    quant_threshold <- if (!is.null(quantile_threshold)) {
      quantile(metric, probs = quantile_threshold)
    } else {
      0
    }
    effective_threshold <- max(abs_threshold, quant_threshold)

    if (verbose) {
      cat("[ELEUTHIA] Region selection (by density):\n")
      cat("    Candidate regions:", nrow(candidate_regions), "\n")
      if (!is.null(min_density)) {
        cat("    Absolute density threshold:", min_density, "fragments/bp\n")
      }
      if (!is.null(quantile_threshold)) {
        cat("    Quantile threshold:", quantile_threshold,
            "(=", round(quant_threshold, 4), "fragments/bp)\n")
      }
      cat("    Effective threshold:", round(effective_threshold, 4),
          "fragments/bp\n")
    }

    regions <- candidate_regions[metric >= effective_threshold, ]

  } else {
    metric          <- candidate_regions$n_fragments
    abs_threshold   <- if (!is.null(min_fragments)) min_fragments else 0
    quant_threshold <- if (!is.null(quantile_threshold)) {
      quantile(metric, probs = quantile_threshold)
    } else {
      0
    }
    effective_threshold <- max(abs_threshold, quant_threshold)

    if (verbose) {
      cat("[ELEUTHIA] Region selection (by fragment count):\n")
      cat("    Candidate regions:", nrow(candidate_regions), "\n")
      if (!is.null(min_fragments)) {
        cat("    Absolute threshold:", min_fragments, "fragments\n")
      }
      if (!is.null(quantile_threshold)) {
        cat("    Quantile threshold:", quantile_threshold,
            "(=", round(quant_threshold, 1), "fragments)\n")
      }
      cat("    Effective threshold:", round(effective_threshold, 1), "fragments\n")
    }

    regions <- candidate_regions[metric >= effective_threshold, ]
  }

  rownames(regions) <- NULL

  if (nrow(regions) == 0) {
    warning("No regions passed the threshold. Consider lowering threshold.")
    return(data.frame(
      chr              = character(),
      start            = integer(),
      end              = integer(),
      peak_id          = character(),
      n_fragments      = integer(),
      n_samples        = integer(),
      fragment_density = numeric(),
      width            = integer()
    ))
  }

  # Add peak IDs
  regions$peak_id <- paste0("region_", seq_len(nrow(regions)))

  # Reorder columns — include n_samples only if present
  base_cols <- c("chr", "start", "end", "peak_id", "n_fragments")
  if (has_n_samples) base_cols <- c(base_cols, "n_samples")
  regions <- regions[, c(base_cols, "fragment_density", "width")]

  if (verbose) {
    cat("[ELEUTHIA] Selection complete:\n")
    cat("    Regions selected:", nrow(regions),
        sprintf("(%.1f%% of candidates)\n", 100 * nrow(regions) / n_candidates))
    cat("    Mean width:", round(mean(regions$width)), "bp\n")
    cat("    Fragment count range:", min(regions$n_fragments), "-",
        max(regions$n_fragments), "\n")
    cat("    Density range:", round(min(regions$fragment_density), 4), "-",
        round(max(regions$fragment_density), 4), "fragments/bp\n")
    if (has_n_samples) {
      cat("    Sample support range:", min(regions$n_samples), "-",
          max(regions$n_samples), "samples\n")
    }
  }

  return(regions)
}


#' Plot Distribution of Fragment Counts in Candidate Regions
#'
#' @description Visualizes the distribution of fragment counts across candidate
#' regions to help choose an appropriate selection threshold.
#'
#' @param candidate_regions Data.frame from ELEUTHIA_merge_fragments() with
#'   n_fragments column.
#' @param show_quantiles Numeric vector of quantiles to mark on the plot
#'   (default = c(0.50, 0.75, 0.90, 0.95, 0.99)).
#' @param log_scale Logical. Use log10 scale for x-axis (default = TRUE).
#'   Useful when fragment counts span several orders of magnitude.
#' @param title Character. Plot title.
#'
#' @return A ggplot object. Also prints quantile summary to console.
#'
#' @details
#' The plot shows:
#' \itemize{
#'   \item Histogram of fragment counts across all candidate regions
#'   \item Vertical lines marking key quantiles
#'   \item A table showing how many regions would be kept at each quantile
#' }
#'
#' Use this to decide on a threshold before calling ELEUTHIA_select_regions().
#'
#' @export
#'
#' @examples
#' \dontrun{
#' candidates <- ELEUTHIA_merge_fragments(bed_files)
#' p <- ELEUTHIA_plot_region_distribution(candidates)
#' print(p)
#'
#' # Then select based on what you see
#' regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)
#'
#' }
ELEUTHIA_plot_region_distribution <- function(candidate_regions,
                                               show_quantiles = c(0.50, 0.75, 0.90, 0.95, 0.99),
                                               log_scale = TRUE,
                                               title = "Distribution of Fragment Counts") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (!"n_fragments" %in% colnames(candidate_regions)) {
    stop("candidate_regions must have 'n_fragments' column")
  }

  n_total <- nrow(candidate_regions)
  frags <- candidate_regions$n_fragments

  # Calculate quantiles
  quant_values <- quantile(frags, probs = show_quantiles)
  quant_df <- data.frame(
    quantile = show_quantiles,
    threshold = quant_values,
    n_above = sapply(quant_values, function(q) sum(frags >= q)),
    stringsAsFactors = FALSE
  )
  quant_df$pct_above <- round(100 * quant_df$n_above / n_total, 1)

  # Print summary
  cat("[ELEUTHIA] Fragment count distribution:\n")
  cat("    Total candidate regions:", n_total, "\n")
  cat("    Range:", min(frags), "-", max(frags), "\n")
  cat("    Mean:", round(mean(frags), 1), "\n")
  cat("    Median:", median(frags), "\n\n")

  cat("[ELEUTHIA] Quantile thresholds:\n")
  for (i in seq_len(nrow(quant_df))) {
    cat(sprintf("    %5.1f%% quantile: >= %6.0f fragments -> %5d regions (%5.1f%%)\n",
                100 * quant_df$quantile[i],
                quant_df$threshold[i],
                quant_df$n_above[i],
                quant_df$pct_above[i]))
  }

  # Build density plot (cleaner than histogram, naturally bounded to data range)
  df <- data.frame(n_fragments = frags)

  # Use log-transformed data for density if log_scale requested
  if (log_scale && max(frags) > 10 * min(frags[frags > 0])) {
    df$x_plot <- log10(frags)
    quant_x <- log10(quant_values)
    x_label <- "Fragment count per region (log10 scale)"
    use_log <- TRUE
  } else {
    df$x_plot <- frags
    quant_x <- quant_values
    x_label <- "Fragment count per region"
    use_log <- FALSE
  }

  # Calculate density to get y-values for label positioning and x-axis cutoff
  dens <- density(df$x_plot)
  max_density <- max(dens$y)

  # Find x-axis cutoff where density drops below 1% of peak
  density_threshold <- max_density * 0.01
  cutoff_idx <- max(which(dens$y >= density_threshold))
  x_max <- dens$x[cutoff_idx]

  # Ensure we include at least the 99th quantile
  x_max <- max(x_max, max(quant_x) * 1.05)
  x_min <- min(df$x_plot) * 0.95

  # Set up x-axis breaks and labels
  if (use_log) {
    possible_breaks <- log10(c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000))
    x_breaks <- possible_breaks[possible_breaks >= x_min & possible_breaks <= x_max]
    x_labels <- as.character(round(10^x_breaks))
  } else {
    x_breaks <- ggplot2::waiver()
    x_labels <- ggplot2::waiver()
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.5, color = "steelblue4") +
    ggplot2::geom_vline(xintercept = quant_x,
                        linetype = "dashed", color = "red", linewidth = 0.6) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max)) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("n = ", format(n_total, big.mark = ","), " candidate regions"),
      x = x_label,
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )

  # Add custom x-axis breaks/labels
  if (use_log) {
    p <- p + ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels)
  }

  # Add quantile labels at top of plot
  for (i in seq_len(nrow(quant_df))) {
    p <- p + ggplot2::annotate(
      "text",
      x = quant_x[i],
      y = max_density * (1.0 - 0.07 * (i - 1)),
      label = paste0("Q", quant_df$quantile[i] * 100),
      color = "red",
      size = 3,
      hjust = -0.1
    )
  }

  return(p)
}


#' Call Regions from Fragment BED Files (Sparse Data)
#'
#' @description Convenience wrapper that combines ELEUTHIA_merge_fragments()
#' and ELEUTHIA_select_regions() into a single call. For more control, use
#' the two-step workflow with those functions directly.
#'
#' @param bed_files Character vector of paths to BED files, or a named list
#'   of loaded BED data.frames.
#' @param merge_distance Integer. Maximum distance (bp) between fragments to
#'   merge into a single region (default = 150).
#' @param min_fragments Integer. Minimum number of fragments required in a
#'   region to call it as a "peak" (default = 5).
#' @param extend Integer. Extend each fragment by this many bp on each side
#'   before merging (default = 75). Helps connect nearby fragments.
#' @param chrom_sizes Optional. Chromosome sizes for clamping extended coordinates.
#' @param pool_samples Logical. If TRUE, pool all samples before calling
#'   regions (recommended for sparse data). (default = TRUE).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A data.frame with columns: chr, start, end, peak_id, n_fragments,
#'   width. Each row represents a called region.
#'
#' @details
#' This is a convenience function that wraps the two-step workflow:
#' \enumerate{
#'   \item ELEUTHIA_merge_fragments() - merge fragments into candidate regions
#'   \item ELEUTHIA_select_regions() - filter by min_fragments threshold
#' }
#'
#' For more control (e.g., quantile-based selection, examining distribution),
#' use the two-step workflow directly.
#'
#' @seealso \code{\link{ELEUTHIA_merge_fragments}}, \code{\link{ELEUTHIA_select_regions}},
#'   \code{\link{ELEUTHIA_plot_region_distribution}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # One-step workflow (quick)
#' regions <- ELEUTHIA_call_regions_from_fragments(bed_files, min_fragments = 100)
#'
#' # Two-step workflow (more control)
#' candidates <- ELEUTHIA_merge_fragments(bed_files)
#' ELEUTHIA_plot_region_distribution(candidates)  # examine distribution
#' regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)
#'
#' }
ELEUTHIA_call_regions_from_fragments <- function(bed_files,
                                                   merge_distance = 150,
                                                   min_fragments = 5,
                                                   extend = 75,
                                                   chrom_sizes = NULL,
                                                   pool_samples = TRUE,
                                                   verbose = TRUE) {

  # Step 1: Merge fragments
  candidates <- ELEUTHIA_merge_fragments(
    bed_files = bed_files,
    merge_distance = merge_distance,
    extend = extend,
    chrom_sizes = chrom_sizes,
    pool_samples = pool_samples,
    verbose = verbose
  )

  # Step 2: Select regions
  regions <- ELEUTHIA_select_regions(
    candidate_regions = candidates,
    min_fragments = min_fragments,
    verbose = verbose
  )

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
    cat("[ELEUTHIA] Loading", nrow(subset_df), omics, "BED files...\n")
  }

  bed_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    file_path <- file.path(subset_df$file_loc[i], subset_df$file_name[i])

    if (verbose) {
      cat("    Loading:", sample_id, "\n")
    }

    bed_list[[sample_id]] <- ELEUTHIA_load_bed(file_path)
  }

  if (verbose) {
    total_frags <- sum(sapply(bed_list, nrow))
    cat("[ELEUTHIA] Loaded", length(bed_list), "BED files.\n")
    cat("[ELEUTHIA] Total fragments:", format(total_frags, big.mark = ","), "\n")
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
#' \dontrun{
#' consensus <- ELEUTHIA_create_consensus_peaks(peak_list)
#' saf <- ELEUTHIA_peaks_to_saf(consensus)
#' counts <- Rsubread::featureCounts(bam_files, annot.ext = saf, isGTFAnnotationFile = FALSE)
#'
#' }
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


#' Expand Genomic Regions by a Window Size
#'
#' @description Symmetrically expands genomic regions by a specified window size.
#' Handles chromosome boundaries by clamping start to 0 and optionally clamping
#' end to chromosome length if sizes are provided.
#'
#' @param regions A data.frame with at minimum columns: chr, start, end.
#'   Typically from ELEUTHIA_create_consensus_peaks() or
#'   ELEUTHIA_call_regions_from_fragments().
#' @param window_size Integer. Number of base pairs to extend on each side
#'   (default = 0, no expansion). Must be >= 0.
#' @param chrom_sizes Optional. Either:
#'   \itemize{
#'     \item A named numeric vector where names are chromosome names and values
#'       are chromosome lengths
#'     \item A data.frame with columns 'chr' and 'size'
#'     \item NULL (default) - end coordinates are not clamped to chromosome length
#'   }
#' @param verbose Logical. Print summary messages (default = TRUE).
#'
#' @return A data.frame with the same structure as input, but with start and end
#'   coordinates expanded by window_size (respecting chromosome boundaries).
#'   An additional column 'original_width' is added showing the pre-expansion width,
#'   and 'window_size' records the expansion used.
#'
#' @details
#' For each region:
#' \itemize{
#'   \item new_start = max(0, start - window_size)
#'   \item new_end = end + window_size (or min(end + window_size, chrom_length) if sizes provided)
#' }
#'
#' If window_size = 0, the function returns the regions unchanged (with added
#' metadata columns). This allows consistent usage in pipelines that test
#' multiple window sizes including zero.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # No expansion (window_size = 0)
#' regions_0 <- ELEUTHIA_expand_regions(chip_regions, window_size = 0)
#'
#' # Expand by 500bp on each side
#' regions_500 <- ELEUTHIA_expand_regions(chip_regions, window_size = 500)
#'
#' # With chromosome size clamping
#' chrom_sizes <- c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559)
#' regions_500 <- ELEUTHIA_expand_regions(chip_regions, window_size = 500,
#'                                         chrom_sizes = chrom_sizes)
#'
#' # Use in anchored analysis workflow
#' expanded <- ELEUTHIA_expand_regions(chip_regions, window_size = 250)
#' atac_at_anchors <- ELEUTHIA_quantify_bed(sample_sheet, expanded, "ATACseq")
#'
#' }
ELEUTHIA_expand_regions <- function(regions,
                                     window_size = 0,
                                     chrom_sizes = NULL,
                                     verbose = TRUE) {

  # Validate inputs
  if (!is.data.frame(regions)) {
    stop("regions must be a data.frame")
  }

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop("regions missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(window_size) || length(window_size) != 1 || window_size < 0) {
    stop("window_size must be a single non-negative number")
  }

  window_size <- as.integer(window_size)
  n_regions <- nrow(regions)

  if (n_regions == 0) {
    warning("regions data.frame is empty, returning as-is")
    return(regions)
  }

  # Process chromosome sizes if provided
  chrom_size_vec <- NULL
  if (!is.null(chrom_sizes)) {
    if (is.data.frame(chrom_sizes)) {
      if (!all(c("chr", "size") %in% colnames(chrom_sizes))) {
        stop("chrom_sizes data.frame must have 'chr' and 'size' columns")
      }
      chrom_size_vec <- setNames(chrom_sizes$size, chrom_sizes$chr)
    } else if (is.numeric(chrom_sizes) && !is.null(names(chrom_sizes))) {
      chrom_size_vec <- chrom_sizes
    } else {
      stop("chrom_sizes must be a named numeric vector or data.frame with 'chr' and 'size' columns")
    }
  }

  # Make a copy to avoid modifying input
  expanded <- regions

  # Store original width
  expanded$original_width <- expanded$end - expanded$start

  # Expand coordinates
  expanded$start <- pmax(0L, as.integer(expanded$start) - window_size)
  expanded$end <- as.integer(expanded$end) + window_size

  # Clamp end to chromosome sizes if provided
  if (!is.null(chrom_size_vec)) {
    for (i in seq_len(n_regions)) {
      chr <- expanded$chr[i]
      if (chr %in% names(chrom_size_vec)) {
        expanded$end[i] <- min(expanded$end[i], chrom_size_vec[chr])
      }
    }
  }

  # Record window size used
  expanded$window_size <- window_size

  if (verbose) {
    new_widths <- expanded$end - expanded$start
    cat("[ELEUTHIA] Region expansion complete:\n")
    cat("    Regions:", n_regions, "\n")
    cat("    Window size: +/-", window_size, "bp\n")
    cat("    Original mean width:", round(mean(expanded$original_width)), "bp\n")
    cat("    Expanded mean width:", round(mean(new_widths)), "bp\n")

    # Report any regions that hit chromosome start boundary
    n_clamped_start <- sum(regions$start - window_size < 0)
    if (n_clamped_start > 0) {
      cat("    Regions clamped at chromosome start:", n_clamped_start, "\n")
    }

    if (!is.null(chrom_size_vec)) {
      # Count end-clamped regions
      n_clamped_end <- sum(mapply(function(chr, end_orig, end_new) {
        if (chr %in% names(chrom_size_vec)) {
          return(end_orig + window_size > chrom_size_vec[chr])
        }
        return(FALSE)
      }, regions$chr, regions$end, expanded$end))

      if (n_clamped_end > 0) {
        cat("    Regions clamped at chromosome end:", n_clamped_end, "\n")
      }
    }
  }

  return(expanded)
}


#' Quantify BED Fragments Against Regions
#'
#' @description Counts the number of BED file fragments overlapping each region.
#' Memory-efficient: loads one BED file at a time from the sample sheet.
#' Works for any omics type (ATACseq, CHIPseq, etc.).
#'
#' @param sample_sheet A validated sample sheet data.frame with bed_loc column.
#' @param regions A regions data.frame with chr, start, end, and peak_id columns
#'   (from ELEUTHIA_call_regions_from_fragments or ELEUTHIA_create_consensus_peaks).
#' @param omics Character string. Omics type to quantify ("ATACseq", "CHIPseq", etc.).
#' @param min_overlap Integer. Minimum bp overlap required to count a fragment
#'   (default = 1, any overlap counts).
#' @param stranded Logical. If TRUE, only count a fragment toward a region
#'   when the fragment's strand (BED column 6) matches the region's strand
#'   (requires a \code{strand} column in \code{regions}, values \code{"+"}/
#'   \code{"-"}). Use this for strand-specific peak sets (e.g. DRIPc-seq
#'   R-loop peaks, which are called separately per strand) where a region at
#'   the same coordinates on the opposite strand is a distinct feature.
#'   Default \code{FALSE} (strand-blind, matches prior behavior).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{counts}{Matrix of fragment counts (regions x samples)}
#'   \item{annotation}{Data.frame with region annotations}
#'   \item{targets}{Data.frame with sample metadata}
#' }
#'
#' @details
#' This function quantifies fragment overlaps with genomic regions using BED
#' files specified in the sample sheet's bed_loc column. It is memory-efficient
#' because it loads and processes one BED file at a time, discarding each before
#' loading the next.
#'
#' A fragment is counted if it overlaps the region by at least \code{min_overlap}
#' base pairs. Each fragment is counted at most once per region.
#'
#' This function works for any omics type - it simply filters the sample sheet
#' by the specified omics value and loads BED files from bed_loc.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ATAC-seq: quantify against consensus peaks
#' atac_peaks <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
#' consensus <- ELEUTHIA_create_consensus_peaks(atac_peaks, min_overlap = 2)
#' atac_counts <- ELEUTHIA_quantify_bed(sample_sheet, consensus, "ATACseq")
#'
#' # ChIP-seq: quantify against called regions
#' chip_beds <- ELEUTHIA_load_bed_from_sheet(sample_sheet, "CHIPseq")
#' chip_regions <- ELEUTHIA_call_regions_from_fragments(chip_beds, min_fragments = 10)
#' chip_counts <- ELEUTHIA_quantify_bed(sample_sheet, chip_regions, "CHIPseq")
#'
#' }
ELEUTHIA_quantify_bed <- function(sample_sheet,
                                   regions,
                                   omics,
                                   min_overlap = 1,
                                   stranded = FALSE,
                                   verbose = TRUE) {

  # Check for data.table
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for fast interval overlaps. ",
         "Install with: install.packages('data.table')")
  }

  # Validate inputs
  if (!"bed_loc" %in% colnames(sample_sheet)) {
    stop("sample_sheet must have a 'bed_loc' column")
  }

  if (nrow(regions) == 0) {
    stop("regions data.frame is empty")
  }

  # Ensure regions have required columns
  required_cols <- c("chr", "start", "end", "peak_id")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop("regions missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (isTRUE(stranded) && !"strand" %in% colnames(regions)) {
    stop("regions must have a 'strand' column when stranded = TRUE")
  }

  # Get subset for this omics type
  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type: ", omics)
  }

  # Filter to rows with valid bed_loc
  valid_bed <- !is.na(subset_df$bed_loc) & subset_df$bed_loc != "NA" &
               trimws(subset_df$bed_loc) != ""

  if (sum(valid_bed) == 0) {
    stop("No valid bed_loc paths found for omics type: ", omics)
  }

  subset_df <- subset_df[valid_bed, , drop = FALSE]

  sample_names <- subset_df$sample_id
  n_samples <- length(sample_names)
  n_regions <- nrow(regions)

  if (verbose) {
    cat("[ELEUTHIA] Quantifying", n_samples, omics, "samples against", n_regions, "regions...\n")
    cat("[ELEUTHIA] (Memory-efficient mode with fast interval overlaps)\n\n")
  }

  # Initialize count matrix
  counts <- matrix(0L, nrow = n_regions, ncol = n_samples)
  rownames(counts) <- regions$peak_id
  colnames(counts) <- sample_names

  # Prepare regions as data.table for foverlaps
  # foverlaps requires: key columns, and start <= end
  regions_dt <- data.table::data.table(
    chr = regions$chr,
    start = as.integer(regions$start),
    end = as.integer(regions$end),
    region_idx = seq_len(n_regions)
  )
  if (isTRUE(stranded)) {
    data.table::set(regions_dt, j = "region_strand", value = regions$strand)
  }
  data.table::setkey(regions_dt, chr, start, end)

  # Process each sample one at a time
  for (s in seq_len(n_samples)) {
    sample_id <- sample_names[s]
    bed_path <- subset_df$bed_loc[s]

    if (verbose) {
      cat("    [", s, "/", n_samples, "] ", sample_id, ": ", sep = "")
    }

    # Check file exists
    if (!file.exists(bed_path)) {
      stop("BED file not found: ", bed_path)
    }

    # Load BED file using data.table::fread (much faster than read.table)
    bed_dt <- data.table::fread(
      bed_path,
      header = FALSE,
      sep = "\t",
      select = if (isTRUE(stranded)) c(1, 2, 3, 6) else 1:3,  # +strand (col 6) when stranded
      col.names = if (isTRUE(stranded)) c("chr", "start", "end", "frag_strand")
                  else c("chr", "start", "end"),
      showProgress = FALSE
    )

    if (verbose) {
      cat(format(nrow(bed_dt), big.mark = ","), "fragments ... ")
    }

    # Ensure integer types for overlap calculation
    # Use data.table::set() to avoid CEDTA check (data.table is in Suggests)
    data.table::set(bed_dt, j = "start", value = as.integer(bed_dt$start))
    data.table::set(bed_dt, j = "end",   value = as.integer(bed_dt$end))

    # Set key for foverlaps (requires start, end naming)
    data.table::setkey(bed_dt, chr, start, end)

    # Find all overlaps using foverlaps (very fast interval join)
    # type="any" finds any overlap, nomatch=NULL excludes non-matches
    overlaps <- data.table::foverlaps(
      bed_dt,
      regions_dt,
      type = "any",
      nomatch = NULL
    )

    if (nrow(overlaps) > 0) {
      # Strand match required first when stranded=TRUE -- a fragment on the
      # opposite strand from the region doesn't count as an overlap at all
      # (e.g. DRIPc-seq peaks called separately per strand)
      if (isTRUE(stranded)) {
        overlaps <- overlaps[overlaps$frag_strand == overlaps$region_strand, ]
      }

      # Apply minimum overlap filter if needed
      if (nrow(overlaps) > 0 && min_overlap > 1) {
        # Calculate actual overlap bp
        data.table::set(overlaps, j = "overlap_bp",
                        value = pmin(overlaps$end, overlaps$i.end) -
                                pmax(overlaps$start, overlaps$i.start))
        # Use $ to avoid [.data.table dispatch (data.table is in Suggests)
        overlaps <- overlaps[overlaps$overlap_bp >= min_overlap, ]
      }

      # Count fragments per region using tabulate() — avoids [.data.table dispatch
      # tabulate(x, nbins) returns a vector of length nbins where [i] = count of i in x
      if (nrow(overlaps) > 0) {
        counts[, s] <- tabulate(overlaps$region_idx, nbins = nrow(regions_dt))
      }
    }

    if (verbose) {
      cat(" done\n")
    }

    # Clean up
    rm(bed_dt, overlaps)
  }

  # Create targets data.frame — include replicate/batch columns only if present
  meta_cols <- intersect(c("sample_id", "group", "bio_rep", "tech_rep", "batch"),
                         colnames(subset_df))
  targets <- subset_df[, meta_cols, drop = FALSE]
  rownames(targets) <- targets$sample_id

  if (verbose) {
    cat("[ELEUTHIA] Quantification complete:\n")
    cat("    Count matrix:", n_regions, "regions x", n_samples, "samples\n")
    cat("    Total counts:", format(sum(counts), big.mark = ","), "\n")
    cat("    Mean counts per region:", round(mean(rowSums(counts)), 1), "\n")
    cat("    Regions with zero counts:", sum(rowSums(counts) == 0), "\n")
  }

  return(list(
    counts = counts,
    annotation = regions,
    targets = targets
  ))
}
