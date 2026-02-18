###############################################################################
########### Circos Plot Functions ###########
###############################################################################

#' Prepare Data for Circos Plots
#'
#' @description Bins fragment counts from BED files into genomic tiles for
#' circos plot visualization. Processes files one at a time for memory
#' efficiency, caps outliers using mean + 3*SD threshold.
#'
#' @param sample_sheet Data frame with sample metadata. Must contain columns:
#'   sample_id, bed_loc, omics, group.
#' @param chr_sizes Data frame with chromosome sizes. Must contain columns:
#'   chr, size. Typically from APOLLO_get_chromosome_sizes().
#' @param omics Character string or vector. Which omics types to include
#'   (e.g., "ATACseq", "CHIPseq", or c("ATACseq", "CHIPseq")).
#' @param bin_size Numeric. Bin size in base pairs (default = 1e6 = 1Mb).
#' @param aggregate Logical. If TRUE, aggregate counts by group (one value per
#'   group). If FALSE (default), keep individual sample values.
#' @param cap_outliers Logical. If TRUE (default), cap outliers at mean + 3*SD
#'   per sample/group.
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{bins}{Data frame with chr, start, end for each genomic bin}
#'   \item{signal}{Matrix of binned signal values (bins x samples/groups)}
#'   \item{sample_info}{Data frame with sample/group metadata and colors}
#'   \item{chr_sizes}{The input chr_sizes (passed through for plotting)}
#'   \item{bin_size}{The bin size used}
#'   \item{omics}{The omics type(s) processed}
#' }
#'
#' @details
#' The function operates directly from the sample_sheet without loading all
#' data into memory. BED files are processed one at a time:
#' 1. Load BED file
#' 2. Count fragments in each genomic bin
#' 3. Discard BED data
#' 4. Move to next sample
#'
#' This approach minimizes memory usage for large datasets.
#'
#' For aggregate = TRUE, samples within each group are summed, then normalized
#' by the number of samples in that group.
#'
#' @note If more than 6 samples are processed in non-aggregated mode, a message
#' suggests using aggregate = TRUE for cleaner visualization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load chromosome sizes
#' chr_sizes <- APOLLO_get_chromosome_sizes(
#'   "annotation.gtf",
#'   name_mapping = "T2T",
#'   chromosomes = c(paste0("chr", 1:22), "chrX", "chrY")
#' )
#'
#' # Prepare circos data for ATAC-seq
#' circos_data <- AETHER_prepare_circos_data(
#'   sample_sheet = my_sample_sheet,
#'   chr_sizes = chr_sizes,
#'   omics = "ATACseq",
#'   bin_size = 1e6,
#'   aggregate = FALSE
#' )
#'
#' # Aggregated by group
#' circos_data_agg <- AETHER_prepare_circos_data(
#'   sample_sheet = my_sample_sheet,
#'   chr_sizes = chr_sizes,
#'   omics = c("ATACseq", "CHIPseq"),
#'   aggregate = TRUE
#' )
#' }
#'
AETHER_prepare_circos_data <- function(sample_sheet,
                                        chr_sizes,
                                        omics,
                                        bin_size = 1e6,
                                        aggregate = FALSE,
                                        cap_outliers = TRUE,
                                        verbose = TRUE) {

  # Check for data.table
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for fast interval overlaps. ",
         "Install with: install.packages('data.table')")
  }

  # Validate inputs
  required_ss_cols <- c("sample_id", "bed_loc", "omics", "group")
  missing_ss <- setdiff(required_ss_cols, colnames(sample_sheet))
  if (length(missing_ss) > 0) {
    stop("sample_sheet missing required columns: ", paste(missing_ss, collapse = ", "))
  }

  required_chr_cols <- c("chr", "size")
  missing_chr <- setdiff(required_chr_cols, colnames(chr_sizes))
  if (length(missing_chr) > 0) {
    stop("chr_sizes missing required columns: ", paste(missing_chr, collapse = ", "))
  }

  if (nrow(chr_sizes) == 0) {
    stop("chr_sizes is empty")
  }

  # Filter sample_sheet for requested omics
  subset_df <- sample_sheet[sample_sheet$omics %in% omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type(s): ", paste(omics, collapse = ", "))
  }

  # Filter to rows with valid bed_loc
  valid_bed <- !is.na(subset_df$bed_loc) & subset_df$bed_loc != "NA" &
    trimws(subset_df$bed_loc) != ""

  if (sum(valid_bed) == 0) {
    stop("No valid bed_loc paths found for omics type(s): ", paste(omics, collapse = ", "))
  }

  subset_df <- subset_df[valid_bed, , drop = FALSE]

  sample_names <- subset_df$sample_id
  n_samples <- length(sample_names)

  if (verbose) {
    cat("Preparing circos data for", n_samples, "samples\n")
    cat("  Omics:", paste(unique(subset_df$omics), collapse = ", "), "\n")
    cat("  Groups:", paste(unique(subset_df$group), collapse = ", "), "\n")
    cat("  Bin size:", format(bin_size, big.mark = ","), "bp\n")
    cat("  Mode:", if (aggregate) "Aggregated by group" else "Individual samples", "\n\n")
  }

  # Suggest aggregation for many samples
  if (!aggregate && n_samples > 6) {
    message("Note: ", n_samples, " samples in individual mode. ",
            "Consider using aggregate = TRUE for cleaner visualization.")
  }

  # Create genomic bins
  bins_list <- list()

  for (i in seq_len(nrow(chr_sizes))) {
    chr_name <- chr_sizes$chr[i]
    chr_len <- chr_sizes$size[i]

    # Create bin starts (0-based)
    bin_starts <- seq(0, chr_len - 1, by = bin_size)
    bin_ends <- pmin(bin_starts + bin_size, chr_len)

    bins_list[[i]] <- data.frame(
      chr = chr_name,
      start = as.integer(bin_starts),
      end = as.integer(bin_ends),
      stringsAsFactors = FALSE
    )
  }

  bins <- do.call(rbind, bins_list)
  bins$bin_id <- seq_len(nrow(bins))
  n_bins <- nrow(bins)

  if (verbose) {
    cat("Created", format(n_bins, big.mark = ","), "genomic bins across",
        nrow(chr_sizes), "chromosomes\n\n")
  }

  # Prepare bins as data.table for foverlaps
  bins_dt <- data.table::data.table(
    chr = bins$chr,
    start = bins$start,
    end = bins$end,
    bin_id = bins$bin_id
  )
  data.table::setkey(bins_dt, chr, start, end)

  # Initialize count matrix
  counts <- matrix(0L, nrow = n_bins, ncol = n_samples)
  colnames(counts) <- sample_names

  # Process each sample one at a time
  for (s in seq_len(n_samples)) {
    sample_id <- sample_names[s]
    bed_path <- subset_df$bed_loc[s]

    if (verbose) {
      cat("  [", s, "/", n_samples, "] ", sample_id, ": ", sep = "")
    }

    # Check file exists
    if (!file.exists(bed_path)) {
      stop("BED file not found: ", bed_path)
    }

    # Load BED file
    bed_dt <- data.table::fread(
      bed_path,
      header = FALSE,
      sep = "\t",
      select = 1:3,
      col.names = c("chr", "start", "end"),
      showProgress = FALSE
    )

    if (verbose) {
      cat(format(nrow(bed_dt), big.mark = ","), "fragments ... ")
    }

    # Filter to chromosomes in chr_sizes
    bed_dt <- bed_dt[bed_dt$chr %in% chr_sizes$chr, ]

    if (nrow(bed_dt) == 0) {
      if (verbose) cat("no fragments in target chromosomes\n")
      next
    }

    # Ensure integer types
    bed_dt[, `:=`(start = as.integer(start), end = as.integer(end))]
    data.table::setkey(bed_dt, chr, start, end)

    # Find overlaps with bins
    overlaps <- data.table::foverlaps(
      bed_dt,
      bins_dt,
      type = "any",
      nomatch = NULL
    )

    if (nrow(overlaps) > 0) {
      # Count fragments per bin
      bin_counts <- overlaps[, .N, by = bin_id]
      counts[bin_counts$bin_id, s] <- bin_counts$N
    }

    if (verbose) {
      cat("done\n")
    }

    # Clean up
    rm(bed_dt, overlaps)
  }

  if (verbose) {
    cat("\nBinning complete. Total counts:", format(sum(counts), big.mark = ","), "\n")
  }

  # Aggregate by group if requested
  if (aggregate) {
    groups <- unique(subset_df$group)
    agg_counts <- matrix(0, nrow = n_bins, ncol = length(groups))
    colnames(agg_counts) <- groups

    for (g in seq_along(groups)) {
      group_samples <- subset_df$sample_id[subset_df$group == groups[g]]
      group_cols <- which(colnames(counts) %in% group_samples)

      if (length(group_cols) == 1) {
        agg_counts[, g] <- counts[, group_cols]
      } else {
        # Sum and normalize by number of samples
        agg_counts[, g] <- rowSums(counts[, group_cols]) / length(group_cols)
      }
    }

    signal <- agg_counts

    # Parse batch from sample_ids (last element after splitting on '_')
    all_batches <- sapply(subset_df$sample_id, function(sid) {
      parts <- strsplit(sid, "_")[[1]]
      parts[length(parts)]
    })

    # Create sample_info for groups
    sample_info <- data.frame(
      id = groups,
      group = groups,
      omics = sapply(groups, function(g) {
        paste(unique(subset_df$omics[subset_df$group == g]), collapse = "+")
      }),
      n_samples = sapply(groups, function(g) sum(subset_df$group == g)),
      batch = sapply(groups, function(g) {
        group_batches <- all_batches[subset_df$group == g]
        paste(unique(group_batches), collapse = ",")
      }),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat("Aggregated to", length(groups), "groups\n")
    }

  } else {
    signal <- counts

    # Parse batch from sample_id (last element after splitting on '_', e.g., "b1")
    parsed_batch <- sapply(subset_df$sample_id, function(sid) {
      parts <- strsplit(sid, "_")[[1]]
      parts[length(parts)]
    })

    # Create sample_info
    sample_info <- data.frame(
      id = subset_df$sample_id,
      group = subset_df$group,
      omics = subset_df$omics,
      batch = parsed_batch,
      stringsAsFactors = FALSE
    )
  }

  # Cap outliers (mean + 3*SD per column)
  if (cap_outliers) {
    for (col in seq_len(ncol(signal))) {
      col_vals <- signal[, col]
      col_mean <- mean(col_vals, na.rm = TRUE)
      col_sd <- sd(col_vals, na.rm = TRUE)

      if (!is.na(col_sd) && col_sd > 0) {
        cap_value <- col_mean + 3 * col_sd
        signal[signal[, col] > cap_value, col] <- cap_value
      }
    }

    if (verbose) {
      cat("Outliers capped at mean + 3*SD per sample/group\n")
    }
  }

  # Remove bin_id column from bins before returning
  bins <- bins[, c("chr", "start", "end")]

  return(list(
    bins = bins,
    signal = signal,
    sample_info = sample_info,
    chr_sizes = chr_sizes,
    bin_size = bin_size,
    omics = omics
  ))
}


#' Create Circos Plot
#'
#' @description Creates a circos plot from prepared data, with signal tracks
#' showing genomic coverage across chromosomes. Supports multiple samples/groups
#' with customizable colors and track organization.
#'
#' @param circos_data List from AETHER_prepare_circos_data() containing bins,
#'   signal matrix, sample_info, and chr_sizes.
#' @param output Character string. Output file path. Format is inferred from
#'   extension: ".png" for PNG, ".pdf" for PDF. If NULL, plot is displayed
#'   but not saved.
#' @param split_by Character vector specifying how to split tracks. Options:
#'   \itemize{
#'     \item "omics" - one track per omics type (always applied)
#'     \item "group" - additionally split by group
#'     \item "batch" - additionally split by batch
#'   }
#'   Can combine: c("omics", "group") creates one track per omics+group.
#'   Default is "omics" only.
#' @param color_by Character string specifying color granularity. Options:
#'   \itemize{
#'     \item "track" (default) - one base color per track, gradients for samples
#'     \item "group" - colors keyed by group name (e.g., "WT", "KO")
#'     \item "omics_group" - colors keyed by omics_group (e.g., "ATACseq_WT")
#'   }
#'   Allows finer color control than track splitting. For example, split_by="omics"
#'   with color_by="omics_group" puts WT and KO in same track but different colors.
#'   When multiple samples share the same color key, gradients are generated.
#' @param invert_gradient Logical. Controls gradient direction relative to signal:
#'   \itemize{
#'     \item FALSE (default) - highest signal = darkest color (plotted on top)
#'     \item TRUE - lowest signal = darkest color (plotted on top)
#'   }
#' @param colors Named character vector of colors. Names should match the
#'   color_by scheme:
#'   \itemize{
#'     \item color_by="track": names like "ATACseq", "CHIPseq"
#'     \item color_by="group": names like "WT", "KO"
#'     \item color_by="omics_group": names like "ATACseq_WT", "ATACseq_KO"
#'   }
#'   If NULL, assigns colors automatically.
#' @param track_height Numeric. Height of each signal track as fraction of
#'   total radius (default = 0.12).
#' @param genome Character. Genome assembly name for cytoband data (default = "hg38").
#' @param track_margin Numeric. Margin between tracks (default = 0.01).
#' @param line_width Numeric. Width of signal lines (default = 1.0).
#' @param ideogram_height Numeric. Height of chromosome ideogram (default = 0.04).
#' @param title Character string. Plot title displayed in center (default = NULL).
#' @param width Numeric. Output width in inches (default = 10).
#' @param height Numeric. Output height in inches (default = 10).
#' @param res Numeric. Resolution for PNG output in DPI (default = 300).
#' @param highlights Named list of highlight definitions. Each element should be
#'   a list with:
#'   \itemize{
#'     \item regions: data.frame with chr, start, end columns
#'     \item tracks: character - track name(s) to highlight, or "all" for all signal tracks
#'     \item color: hex color string (should include alpha for transparency, e.g., "#FF000033")
#'   }
#'   Example: list("peaks" = list(regions = peaks_df, tracks = "ATACseq", color = "#FF000033"))
#' @param highlight_padding Numeric. Padding for highlights as fraction of chromosome
#'   length (default = 0.005, i.e., 0.5% on each side). Increase to make small regions
#'   more visible.
#' @param annotation_regions Data frame with annotated regions for chord connections.
#'   Must contain columns: chr, start, end, annotation. The annotation column should
#'   contain category labels (e.g., "Promoter", "Intron", "Exon"). Compatible with
#'   APOLLO_annotate_peaks() output. User can simplify annotation labels before passing.
#' @param annotation_colors Named character vector of colors for annotation categories.
#'   Names should match values in annotation_regions$annotation. If NULL, default
#'   colors are used.
#' @param show_chords Logical. If TRUE (default), draw chord connections from regions
#'   to annotation sector. If FALSE, only show annotation sector without chords.
#' @param chord_alpha Numeric. Transparency for chord fill (0-1, default = 0.5).
#' @param annotation_sector_size Numeric. Size of the annotation pseudo-chromosome
#'   in base pairs (default = 2e8). Adjust for visual balance.
#' @param start_degree Numeric. Starting angle in degrees for the first sector
#'   (default = 90, which is top/12 o'clock). Use smaller values to rotate clockwise.
#'   For annotation sector at center-right (3 o'clock), try values around -50 to -70
#'   depending on the number of chromosomes.
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return Invisibly returns NULL. Side effect is creating the plot.
#'
#' @details
#' The function creates a circos plot with:
#' \itemize{
#'   \item Outer ring: Chromosome ideogram with labels and alternating colors
#'   \item Inner rings: Signal tracks, each independently scaled
#' }
#'
#' Track organization examples:
#' \itemize{
#'   \item split_by = "omics": ATACseq track, CHIPseq track
#'   \item split_by = c("omics", "group"): ATACseq_WT, ATACseq_KO, CHIPseq_WT, CHIPseq_KO
#'   \item split_by = c("omics", "batch"): ATACseq_b1, ATACseq_b2, etc.
#' }
#'
#' Each track is independently scaled (y-axis) to prevent squishing when
#' combining different omics types with different signal ranges.
#'
#' Color assignment depends on the color_by parameter:
#' \itemize{
#'   \item color_by="track": One base color per track, with gradient shades for
#'     multiple samples (darkest = highest signal, plotted on top)
#'   \item color_by="group": Colors assigned by group name, allowing different
#'     colors within the same track
#'   \item color_by="omics_group": Colors assigned by omics+group combination,
#'     e.g., "ATACseq_WT" and "ATACseq_KO" can have different colors even in
#'     the same track when split_by="omics"
#' }
#'
#' Samples are always plotted in order of mean signal (lowest first, highest
#' last on top) to preserve data visibility.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' circos_data <- AETHER_prepare_circos_data(
#'   sample_sheet, chr_sizes, omics = c("ATACseq", "CHIPseq")
#' )
#'
#' # One track per omics (default)
#' AETHER_create_circos(circos_data, output = "circos.png")
#'
#' # Split by omics and group
#' AETHER_create_circos(circos_data, output = "circos.png",
#'                       split_by = c("omics", "group"))
#'
#' # With centered title
#' AETHER_create_circos(circos_data, output = "circos.png",
#'                       title = "Multi-omics Overview")
#'
#' # With region highlights
#' my_highlights <- list(
#'   "ATAC_peaks" = list(
#'     regions = atac_peaks,  # data.frame with chr, start, end
#'     tracks = "ATACseq",
#'     color = "#FF000033"    # red with transparency
#'   ),
#'   "shared_peaks" = list(
#'     regions = shared_peaks,
#'     tracks = "all",        # highlight on all signal tracks
#'     color = "#0000FF33"    # blue with transparency
#'   )
#' )
#' AETHER_create_circos(circos_data, output = "circos.png",
#'                       highlights = my_highlights)
#'
#' # With annotation chords (Phase 2B)
#' # annotation_regions should have chr, start, end, annotation columns
#' # (compatible with APOLLO_annotate_peaks() output)
#' annotated_peaks <- APOLLO_annotate_peaks(my_peaks, txdb)
#'
#' # Optional: simplify annotation labels
#' annotated_peaks$annotation <- gsub("Distal Intergenic", "Intergenic",
#'                                     annotated_peaks$annotation)
#'
#' AETHER_create_circos(circos_data, output = "circos_with_chords.png",
#'                       annotation_regions = annotated_peaks,
#'                       show_chords = TRUE,
#'                       chord_alpha = 0.4)
#'
#' # Custom annotation colors
#' my_anno_colors <- c(
#'   "Promoter" = "#E78AC3",
#'   "Intron" = "#66C2A5",
#'   "Exon" = "#8DA0CB",
#'   "Intergenic" = "#FC8D62"
#' )
#' AETHER_create_circos(circos_data, output = "circos_custom.png",
#'                       annotation_regions = annotated_peaks,
#'                       annotation_colors = my_anno_colors)
#' }
#'
AETHER_create_circos <- function(circos_data,
                                  genome='hg38',
                                  output = NULL,
                                  split_by = "omics",
                                  color_by = "track",
                                  invert_gradient = FALSE,
                                  colors = NULL,
                                  track_height = 0.12,
                                  track_margin = 0.01,
                                  line_width = 1.0,
                                  ideogram_height = 0.04,
                                  title = NULL,
                                  width = 10,
                                  height = 10,
                                  res = 300,
                                  highlights = NULL,
                                  highlight_padding = 0.005,
                                  annotation_regions = NULL,
                                  annotation_colors = NULL,
                                  show_chords = TRUE,
                                  chord_alpha = 0.5,
                                  annotation_sector_size = 2e8,
                                  start_degree = 90,
                                  verbose = TRUE) {

  # Check for circlize
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for circos plots. ",
         "Install with: install.packages('circlize')")
  }

  # Validate circos_data structure
  required_elements <- c("bins", "signal", "sample_info", "chr_sizes")
  missing_elements <- setdiff(required_elements, names(circos_data))
  if (length(missing_elements) > 0) {
    stop("circos_data missing required elements: ",
         paste(missing_elements, collapse = ", "))
  }

  bins <- circos_data$bins
  signal <- circos_data$signal
  sample_info <- circos_data$sample_info
  chr_sizes <- circos_data$chr_sizes

  n_samples <- ncol(signal)
  sample_ids <- colnames(signal)

  # Validate split_by
  valid_splits <- c("omics", "group", "batch")
  invalid_splits <- setdiff(split_by, valid_splits)
  if (length(invalid_splits) > 0) {
    stop("Invalid split_by values: ", paste(invalid_splits, collapse = ", "),
         "\nValid options: ", paste(valid_splits, collapse = ", "))
  }

  # Ensure omics is always included
  if (!"omics" %in% split_by) {
    split_by <- c("omics", split_by)
  }

  # Check if required columns exist for splitting
  if ("group" %in% split_by && !"group" %in% colnames(sample_info)) {
    stop("split_by includes 'group' but sample_info has no 'group' column")
  }
  if ("batch" %in% split_by && !"batch" %in% colnames(sample_info)) {
    stop("split_by includes 'batch' but sample_info has no 'batch' column. ",
         "Batch is parsed from sample_id (last element after '_').")
  }

  if (verbose) {
    cat("Creating circos plot\n")
    cat("  Samples/groups:", n_samples, "\n")
    cat("  Chromosomes:", nrow(chr_sizes), "\n")
    cat("  Bins:", nrow(bins), "\n")
    cat("  Split by:", paste(split_by, collapse = " + "), "\n")
  }

  # Create track groupings based on split_by
  track_ids <- apply(sample_info[, split_by, drop = FALSE], 1, function(x) {
    paste(x, collapse = "_")
  })
  sample_info$track_id <- track_ids

  # Get unique tracks in order
  unique_tracks <- unique(track_ids)
  n_tracks <- length(unique_tracks)

  if (verbose) {
    cat("  Number of tracks:", n_tracks, "\n")
    cat("  Tracks:", paste(unique_tracks, collapse = ", "), "\n")
  }

  # Validate color_by
  valid_color_by <- c("track", "group", "omics_group")
  if (!color_by %in% valid_color_by) {
    stop("Invalid color_by value: ", color_by,
         "\nValid options: ", paste(valid_color_by, collapse = ", "))
  }

  # Build sample-level color keys based on color_by mode
  # These keys are used to look up colors for each sample
  if (color_by == "track") {
    sample_color_keys <- sample_info$track_id
  } else if (color_by == "group") {
    sample_color_keys <- sample_info$group
  } else if (color_by == "omics_group") {
    sample_color_keys <- paste(sample_info$omics, sample_info$group, sep = "_")
  }
  names(sample_color_keys) <- sample_ids

  # Get unique color keys for legend
  unique_color_keys <- unique(sample_color_keys)

  # Assign colors based on color_by mode
  omics_palette <- AETHER_default_palette("omics")
  group_palette <- AETHER_default_palette("groups")

  if (is.null(colors)) {
    # Auto-assign colors based on color_by mode
    key_colors <- character(length(unique_color_keys))
    names(key_colors) <- unique_color_keys

    for (i in seq_along(unique_color_keys)) {
      key <- unique_color_keys[i]
      key_parts <- strsplit(key, "_")[[1]]

      # Try omics first, then group, then rainbow
      if (key_parts[1] %in% names(omics_palette)) {
        key_colors[i] <- omics_palette[key_parts[1]]
      } else if (key %in% names(group_palette)) {
        key_colors[i] <- group_palette[key]
      } else if (length(key_parts) > 1 && key_parts[2] %in% names(group_palette)) {
        key_colors[i] <- group_palette[key_parts[2]]
      } else {
        key_colors[i] <- grDevices::rainbow(length(unique_color_keys))[i]
      }
    }
  } else {
    # User provided colors
    key_colors <- colors
  }

  # Build sample-level colors by looking up each sample's key
  sample_colors <- key_colors[sample_color_keys]
  names(sample_colors) <- sample_ids

  if (verbose) {
    cat("  Color by:", color_by, "\n")
    cat("  Colors assigned for", length(unique_color_keys), "groups\n")
  }

  # Set up alternating chromosome colors for ideogram
  n_chr <- nrow(chr_sizes)
  chr_colors <- rep(c("#E8E8E8", "#C8C8C8"), length.out = n_chr)
  names(chr_colors) <- chr_sizes$chr

  # ---------------------------------------------------------------------------
  # Process annotation regions for chords (Phase 2B)
  # ---------------------------------------------------------------------------
  use_annotations <- !is.null(annotation_regions)
  annotation_segments <- NULL
  chord_data <- NULL
  custom_cytoband <- NULL

  if (use_annotations) {
    # Validate annotation_regions
    required_anno_cols <- c("chr", "start", "end", "annotation")
    missing_anno <- setdiff(required_anno_cols, colnames(annotation_regions))
    if (length(missing_anno) > 0) {
      stop("annotation_regions missing required columns: ",
           paste(missing_anno, collapse = ", "))
    }

    if (verbose) {
      cat("  Processing annotation regions for chords...\n")
      cat("    Regions:", nrow(annotation_regions), "\n")
      cat("    Categories:", paste(unique(annotation_regions$annotation), collapse = ", "), "\n")
    }

    # Calculate proportional segment allocation
    anno_counts <- table(annotation_regions$annotation)
    total_regions <- sum(anno_counts)
    anno_proportions <- anno_counts / total_regions
    anno_ranges <- anno_proportions * annotation_sector_size

    # Build segment boundaries
    anno_categories <- names(anno_counts)
    segment_starts <- numeric(length(anno_categories))
    segment_ends <- numeric(length(anno_categories))
    current_pos <- 0

    for (i in seq_along(anno_categories)) {
      segment_starts[i] <- current_pos
      segment_ends[i] <- current_pos + anno_ranges[i]
      current_pos <- segment_ends[i]
    }

    annotation_segments <- data.frame(
      category = anno_categories,
      count = as.numeric(anno_counts),
      proportion = as.numeric(anno_proportions),
      start = segment_starts,
      end = segment_ends,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat("    Segment allocation:\n")
      for (i in seq_len(nrow(annotation_segments))) {
        cat("      ", annotation_segments$category[i], ": ",
            round(annotation_segments$proportion[i] * 100, 1), "% (",
            annotation_segments$count[i], " regions)\n", sep = "")
      }
    }

    # Set default annotation colors if not provided
    if (is.null(annotation_colors)) {
      # Default colors from original implementation
      default_anno_colors <- c(
        "Promoter" = rgb(231/255, 138/255, 195/255, 0.7),    # pink
        "Exon" = rgb(141/255, 160/255, 203/255, 0.7),        # blue
        "Intron" = rgb(102/255, 194/255, 165/255, 0.7),      # green
        "Distal Intergenic" = rgb(252/255, 141/255, 98/255, 0.7),  # orange
        "Intergenic" = rgb(252/255, 141/255, 98/255, 0.7),   # orange (alias)
        "3' UTR" = rgb(166/255, 206/255, 227/255, 0.7),      # light blue
        "5' UTR" = rgb(253/255, 191/255, 111/255, 0.7),      # light orange
        "Downstream" = rgb(178/255, 178/255, 178/255, 0.7)   # gray
      )
      annotation_colors <- default_anno_colors
    }

    # Assign colors to segments (use rainbow for unknown categories)
    segment_colors <- character(length(anno_categories))
    for (i in seq_along(anno_categories)) {
      cat_name <- anno_categories[i]
      if (cat_name %in% names(annotation_colors)) {
        segment_colors[i] <- annotation_colors[cat_name]
      } else {
        # Try partial matching
        matched <- FALSE
        for (known in names(annotation_colors)) {
          if (grepl(known, cat_name, ignore.case = TRUE) ||
              grepl(cat_name, known, ignore.case = TRUE)) {
            segment_colors[i] <- annotation_colors[known]
            matched <- TRUE
            break
          }
        }
        if (!matched) {
          segment_colors[i] <- grDevices::rainbow(length(anno_categories))[i]
        }
      }
    }
    annotation_segments$color <- segment_colors

    # Build chord data: map each region to its annotation segment
    chord_data <- data.frame(
      chr = annotation_regions$chr,
      start = annotation_regions$start,
      end = annotation_regions$end,
      annotation = annotation_regions$annotation,
      stringsAsFactors = FALSE
    )

    # Add target coordinates in annotation sector
    chord_data$to_chr <- "annotations"
    chord_data$to_start <- NA_real_
    chord_data$to_end <- NA_real_
    chord_data$color <- NA_character_

    for (i in seq_len(nrow(chord_data))) {
      anno <- chord_data$annotation[i]
      seg_idx <- which(annotation_segments$category == anno)
      if (length(seg_idx) > 0) {
        chord_data$to_start[i] <- annotation_segments$start[seg_idx]
        chord_data$to_end[i] <- annotation_segments$end[seg_idx]
        chord_data$color[i] <- annotation_segments$color[seg_idx]
      }
    }

    # Filter out any regions with unknown annotations
    chord_data <- chord_data[!is.na(chord_data$to_start), ]

    # Create custom cytoband by reading default and adding annotations
    tryCatch({
      default_cytoband <- circlize::read.cytoband(species = genome)$df
      custom_cytoband <- rbind(
        default_cytoband,
        data.frame(V1 = "annotations", V2 = 0, V3 = annotation_sector_size,
                   V4 = "", V5 = "", stringsAsFactors = FALSE)
      )
    }, error = function(e) {
      warning("Could not read default cytoband for genome '", genome, "'. ",
              "Creating minimal cytoband from chr_sizes.")
      # Create minimal cytoband from chr_sizes
      custom_cytoband <<- data.frame(
        V1 = c(chr_sizes$chr, "annotations"),
        V2 = 0,
        V3 = c(chr_sizes$size, annotation_sector_size),
        V4 = "",
        V5 = "",
        stringsAsFactors = FALSE
      )
    })

    if (verbose) {
      cat("    Custom cytoband created with 'annotations' sector\n")
    }
  }

  # Determine output format
  if (!is.null(output)) {
    ext <- tolower(tools::file_ext(output))
    if (ext == "png") {
      grDevices::png(output, width = width, height = height, units = "in", res = res)
    } else if (ext == "pdf") {
      grDevices::pdf(output, width = width, height = height)
    } else {
      warning("Unknown file extension '", ext, "'. Defaulting to PNG.")
      grDevices::png(output, width = width, height = height, units = "in", res = res)
    }
    on.exit(grDevices::dev.off())
  }

  # Initialize circos
  circlize::circos.clear()

  # Set circos parameters
  circlize::circos.par(
    start.degree = start_degree,
    gap.degree = 2,
    track.margin = c(track_margin, track_margin),
    cell.padding = c(0, 0, 0, 0)
  )

  # Initialize with chromosome labels only (no axis ticks)
  # Use custom cytoband if annotations are present
  if (use_annotations && !is.null(custom_cytoband)) {
    chr_order <- c(chr_sizes$chr, "annotations")
    circlize::circos.initializeWithIdeogram(custom_cytoband, chromosome.index = chr_order)
  } else {
    circlize::circos.initializeWithIdeogram(species = genome, chromosome.index = chr_sizes$chr)
  }

  # Create signal tracks - one per track grouping
  # Each track has independent y-axis scaling
  # Store y-limits for later use in highlights
  track_ylim_store <- list()

  for (t in seq_along(unique_tracks)) {
    track_name <- unique_tracks[t]

    # Get samples in this track
    track_sample_idx <- which(sample_info$track_id == track_name)
    track_sample_ids <- sample_ids[track_sample_idx]

    # Calculate y_max for this track independently
    track_y_max <- max(signal[, track_sample_idx, drop = FALSE], na.rm = TRUE)

    # Store y-limits for this track (track index = t + 1 because ideogram is track 1)
    track_ylim_store[[t + 1]] <- c(-0.03 * track_y_max, track_y_max)

    # Calculate mean signal for each sample in this track
    track_means <- colMeans(signal[, track_sample_idx, drop = FALSE], na.rm = TRUE)

    # Order samples by mean signal (ascending: lowest first, highest last)
    # Highest signal plotted last (on top)
    plot_order <- track_sample_idx[order(track_means)]

    # Assign colors based on color_by mode
    n_in_track <- length(plot_order)

    if (color_by == "track") {
      # Track mode: generate gradient from track's base color
      track_col <- key_colors[track_name]

      if (n_in_track > 1) {
        # Generate gradient
        track_sample_cols <- AETHER_generate_gradient(track_col, n = n_in_track,
                                                       direction = "to_light")
        # Control gradient direction: invert_gradient=FALSE means darkest=highest signal
        if (invert_gradient) {
          # Darkest first (lowest signal), lightest last (highest signal)
          # No reversal needed - gradient already goes dark to light
        } else {
          # Darkest last (highest signal) - reverse so dark is at end
          track_sample_cols <- rev(track_sample_cols)
        }
      } else {
        track_sample_cols <- track_col
      }
    } else {
      # Group or omics_group mode: generate gradients per color key within track
      # Get color keys for samples in this track (in plot_order)
      track_sample_ids <- sample_ids[plot_order]
      track_color_keys <- sample_color_keys[track_sample_ids]

      # Initialize color vector for this track
      track_sample_cols <- character(n_in_track)

      # For each unique color key in this track, generate gradient if needed
      unique_keys_in_track <- unique(track_color_keys)

      for (key in unique_keys_in_track) {
        key_indices <- which(track_color_keys == key)
        n_with_key <- length(key_indices)
        base_col <- key_colors[key]

        if (n_with_key > 1) {
          # Generate gradient for samples sharing this color key
          key_gradient <- AETHER_generate_gradient(base_col, n = n_with_key,
                                                    direction = "to_light")
          # Control gradient direction
          if (invert_gradient) {
            # Darkest first (lowest signal)
          } else {
            # Darkest last (highest signal)
            key_gradient <- rev(key_gradient)
          }
          track_sample_cols[key_indices] <- key_gradient
        } else {
          track_sample_cols[key_indices] <- base_col
        }
      }
    }

    # Create the track using local() to capture variables
    local({
      t_name <- track_name
      t_y_max <- track_y_max
      t_plot_order <- plot_order
      t_sample_cols <- track_sample_cols
      t_bins <- bins
      t_signal <- signal
      t_line_width <- line_width

      circlize::circos.track(
        ylim = c(-0.03 * t_y_max, t_y_max),  # Small negative offset to raise data above bottom border
        track.height = track_height,
        bg.border = "gray70",
        bg.col = NA,
        panel.fun = function(x, y) {
          chr <- circlize::CELL_META$sector.index
          xlim <- circlize::CELL_META$xlim

          # Get bins for this chromosome
          chr_bins <- t_bins[t_bins$chr == chr, ]
          if (nrow(chr_bins) == 0) return()

          # Calculate bin midpoints
          midpoints <- (chr_bins$start + chr_bins$end) / 2

          # Get row indices for this chromosome
          chr_rows <- which(t_bins$chr == chr)

          # Filter to bins within sector boundaries
          in_bounds <- midpoints >= xlim[1] & midpoints <= xlim[2]
          if (sum(in_bounds) == 0) return()

          midpoints <- midpoints[in_bounds]
          chr_rows <- chr_rows[in_bounds]

          # Plot each sample in order (lowest signal first, highest last)
          for (i in seq_along(t_plot_order)) {
            s <- t_plot_order[i]
            values <- t_signal[chr_rows, s]

            circlize::circos.lines(
              midpoints,
              values,
              col = t_sample_cols[i],
              lwd = t_line_width,
              area = FALSE
            )
          }
        }
      )
    })
  }

  # ---------------------------------------------------------------------------
  # Add highlights (Phase 2A)
  # ---------------------------------------------------------------------------
  if (!is.null(highlights) && length(highlights) > 0) {
    if (verbose) {
      cat("  Adding", length(highlights), "highlight layer(s)...\n")
    }

    # Build mapping from track names to track indices
    # Track 1 is ideogram, signal tracks start at 2
    track_index_map <- setNames(seq_along(unique_tracks) + 2, unique_tracks)

    for (hl_name in names(highlights)) {
      hl <- highlights[[hl_name]]

      # Validate highlight structure
      if (!is.list(hl) || is.null(hl$regions) || is.null(hl$tracks) || is.null(hl$color)) {
        warning("Highlight '", hl_name, "' missing required elements (regions, tracks, color). Skipping.")
        next
      }

      regions <- hl$regions
      target_tracks <- hl$tracks
      hl_color <- hl$color

      # Validate regions data.frame
      if (!is.data.frame(regions) || !all(c("chr", "start", "end") %in% colnames(regions))) {
        warning("Highlight '", hl_name, "' regions must be a data.frame with chr, start, end. Skipping.")
        next
      }

      # Determine target track indices
      if (length(target_tracks) == 1 && target_tracks == "all") {
        # "all" means all signal tracks (offset by 2: ideogram is 1, first signal is 2)
        target_indices <- seq_along(unique_tracks) + 2
      } else {
        # Map track names to indices
        target_indices <- track_index_map[target_tracks]
        target_indices <- target_indices[!is.na(target_indices)]

        if (length(target_indices) == 0) {
          warning("Highlight '", hl_name, "' has no valid target tracks. ",
                  "Available: ", paste(unique_tracks, collapse = ", "))
          next
        }
      }

      # Draw rectangles for each region on each target track
      # Snap to bin coordinates for visual alignment with signal data
      n_drawn <- 0
      for (i in seq_len(nrow(regions))) {
        region_chr <- regions$chr[i]
        region_start <- regions$start[i]
        region_end <- regions$end[i]

        # Skip if chromosome not in plot
        if (!region_chr %in% chr_sizes$chr) next

        # Find overlapping bins for visual alignment
        chr_bins <- bins[bins$chr == region_chr, ]
        overlapping_bins <- chr_bins[chr_bins$end > region_start & chr_bins$start < region_end, ]

        if (nrow(overlapping_bins) > 0) {
          # Snap to bin boundaries for visual alignment
          region_start_snapped <- min(overlapping_bins$start)
          region_end_snapped <- max(overlapping_bins$end)
        } else {
          # No overlapping bins - find nearest bin
          bin_midpoints <- (chr_bins$start + chr_bins$end) / 2
          region_midpoint <- (region_start + region_end) / 2
          nearest_idx <- which.min(abs(bin_midpoints - region_midpoint))
          region_start_snapped <- chr_bins$start[nearest_idx]
          region_end_snapped <- chr_bins$end[nearest_idx]
        }

        # Apply padding for visibility (as fraction of chromosome length)
        chr_len <- chr_sizes$size[chr_sizes$chr == region_chr]
        padding <- chr_len * highlight_padding
        region_start_final <- max(0, region_start_snapped - padding)
        region_end_final <- min(chr_len, region_end_snapped + padding)

        n_drawn <- n_drawn + 1

        for (track_idx in target_indices) {
          tryCatch({
            # Enter the track/sector context before drawing
            circlize::set.current.cell(sector.index = region_chr, track.index = track_idx)

            # Get y-limits from current cell context (like original implementation)
            ylim <- circlize::get.cell.meta.data('cell.ylim')

            # Debug first rectangle
            if (verbose && n_drawn == 1 && track_idx == target_indices[1]) {
              cat("      DEBUG: First rect - chr=", region_chr,
                  " original=", region_start, "-", region_end,
                  " snapped=", region_start_snapped, "-", region_end_snapped,
                  " final=", region_start_final, "-", region_end_final,
                  " track=", track_idx, "\n", sep="")
            }

            circlize::circos.rect(
              xleft = region_start_final,
              xright = region_end_final,
              ybottom = ylim[1],
              ytop = ylim[2],
              col = hl_color,
              border = NA
            )
          }, error = function(e) {
            if (verbose) {
              warning("Highlight error on ", region_chr, " track ", track_idx, ": ", e$message)
            }
          })
        }
      }

      if (verbose) {
        cat("    '", hl_name, "': ", nrow(regions), " regions on ",
            length(target_indices), " track(s)\n", sep = "")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Draw annotation sector and chords (Phase 2B)
  # ---------------------------------------------------------------------------
  if (use_annotations && !is.null(annotation_segments) && !is.null(chord_data)) {
    if (verbose) {
      cat("  Adding annotation sector and chords...\n")
    }

    # Determine track indices to mask (ideogram + all signal tracks)
    # Track 1 = ideogram, signal tracks use +2 offset (tracks 3 to n_tracks+2)
    tracks_to_mask <- seq_len(n_tracks + 2)

    # Step 1: Mask signal tracks in annotation sector with white
    # First update each track's plot region
    for (track_idx in tracks_to_mask) {
      tryCatch({
        circlize::circos.updatePlotRegion(
          sector.index = "annotations",
          track.index = track_idx,
          bg.border = "white",
          bg.col = "white"
        )
      }, error = function(e) {
        # Silently ignore if track doesn't exist in annotation sector
      })
    }

    # Use highlight.sector to fully mask with white
    tryCatch({
      circlize::highlight.sector(
        sector.index = "annotations",
        track.index = tracks_to_mask,
        col = "white",
        padding = c(0.01, 0.01, 0.01, 0.01)
      )
    }, error = function(e) {
      if (verbose) {
        warning("Could not fully mask annotation sector: ", e$message)
      }
    })

    # Step 2: Define innermost track for annotation labels
    # Use the innermost track (n_tracks + 2) - same offset as highlights
    # Note: We no longer draw colored segments here, just use for label positioning
    innermost_track <- n_tracks + 2

    # Step 3: Draw chords if enabled
    if (show_chords && nrow(chord_data) > 0) {
      if (verbose) {
        cat("    Drawing ", nrow(chord_data), " chords...\n", sep = "")
      }

      # Apply alpha to chord colors
      for (i in seq_len(nrow(chord_data))) {
        region_chr <- chord_data$chr[i]
        region_start <- chord_data$start[i]
        region_end <- chord_data$end[i]
        to_start <- chord_data$to_start[i]
        to_end <- chord_data$to_end[i]
        chord_col <- chord_data$color[i]

        # Skip if chromosome not in plot
        if (!region_chr %in% chr_sizes$chr) next

        # Apply alpha to color
        if (!is.na(chord_col)) {
          # Extract RGB and apply alpha
          rgb_vals <- grDevices::col2rgb(chord_col)
          chord_col_alpha <- grDevices::rgb(
            rgb_vals[1, 1], rgb_vals[2, 1], rgb_vals[3, 1],
            alpha = chord_alpha * 255,
            maxColorValue = 255
          )
        } else {
          chord_col_alpha <- grDevices::rgb(0.5, 0.5, 0.5, chord_alpha)
        }

        tryCatch({
          circlize::circos.link(
            sector.index1 = region_chr,
            point1 = c(region_start, region_end),
            sector.index2 = "annotations",
            point2 = c(to_start, to_end),
            col = chord_col_alpha,
            border = chord_col_alpha
          )
        }, error = function(e) {
          # Silently skip failed chords
        })
      }

      if (verbose) {
        cat("    Chords drawn\n")
      }
    }

    # Step 4: Add category labels to annotation sector (radial orientation, at chord centers)
    if (verbose) {
      cat("    Adding annotation labels...\n")
    }

    # Determine text orientation based on annotation sector position
    # Calculate where the annotation sector's midpoint falls angularly
    total_size <- sum(chr_sizes$size) + annotation_sector_size
    chr_proportion <- sum(chr_sizes$size) / total_size
    anno_proportion <- annotation_sector_size / total_size

    # Account for gaps: n_sectors gaps of gap.degree each
    n_sectors <- nrow(chr_sizes) + 1  # chromosomes + annotations
    total_gap <- n_sectors * 2  # gap.degree = 2
    available_degrees <- 360 - total_gap

    # Annotation sector starts after all chromosomes (counter-clockwise from start_degree)
    chr_span <- chr_proportion * available_degrees
    anno_span <- anno_proportion * available_degrees

    # Midpoint of annotation sector (in degrees from start)
    # Sectors go counter-clockwise, so we subtract from start_degree
    anno_start_angle <- start_degree - chr_span - (nrow(chr_sizes) * 2)  # subtract chr span and gaps
    anno_mid_angle <- anno_start_angle - (anno_span / 2)

    # Normalize to 0-360 range
    anno_mid_angle <- anno_mid_angle %% 360
    if (anno_mid_angle < 0) anno_mid_angle <- anno_mid_angle + 360

    # Determine facing: right half (315-360, 0-45) or (270-90 via 0) = clockwise
    # left half (90-270) = reverse.clockwise
    # Right side: 270 < angle <= 360 OR 0 <= angle < 90
    if ((anno_mid_angle > 270 && anno_mid_angle <= 360) ||
        (anno_mid_angle >= 0 && anno_mid_angle < 90)) {
      anno_text_facing <- "clockwise"
    } else {
      anno_text_facing <- "reverse.clockwise"
    }

    if (verbose) {
      cat("    Annotation sector midpoint:", round(anno_mid_angle, 1),
          "deg -> text facing:", anno_text_facing, "\n")
    }

    for (i in seq_len(nrow(annotation_segments))) {
      seg_start <- annotation_segments$start[i]
      seg_end <- annotation_segments$end[i]
      label <- annotation_segments$category[i]

      # Calculate midpoint for label placement (center of chord connection area)
      label_pos <- (seg_start + seg_end) / 2

      tryCatch({
        circlize::set.current.cell(sector.index = "annotations", track.index = innermost_track)
        ylim <- circlize::get.cell.meta.data("cell.ylim")

        # Only add label if segment is large enough
        segment_width <- seg_end - seg_start
        if (segment_width > annotation_sector_size * 0.05) {  # At least 5% of sector
          circlize::circos.text(
            x = label_pos,
            y = mean(ylim),  # Center of track
            labels = label,
            facing = anno_text_facing,  # Auto-determined based on position
            niceFacing = FALSE,  # Don't auto-adjust orientation
            cex = 0.6,
            adj = c(0.5, 0.5)  # Center the text
          )
        }
      }, error = function(e) {
        # Silently skip failed labels
      })
    }

    if (verbose) {
      cat("  Annotation sector complete\n")
    }
  }

  # Add centered title if provided
  if (!is.null(title)) {
    text(0, 0, title, cex = 1.2, font = 2)
  }

  # Add legend showing color keys (based on color_by mode)
  graphics::legend(
    "bottomright",
    legend = unique_color_keys,
    col = key_colors[unique_color_keys],
    lwd = 3,
    cex = 0.7,
    bg = "white",
    box.lwd = 0.5,
    title = if (color_by == "track") "Tracks" else if (color_by == "group") "Groups" else "Omics/Group"
  )

  circlize::circos.clear()

  if (verbose) {
    if (!is.null(output)) {
      cat("Plot saved to:", output, "\n")
    } else {
      cat("Plot displayed\n")
    }
  }

  invisible(NULL)
}


#' Assign Colors for Circos Plot (Internal)
#'
#' @description Internal helper function to assign colors to samples for
#' circos plots.
#'
#' @param sample_info Data frame with sample metadata
#' @param sample_ids Character vector of sample IDs
#'
#' @return Named character vector of colors
#'
#' @keywords internal
#'
.AETHER_assign_circos_colors <- function(sample_info, sample_ids) {

  n_samples <- length(sample_ids)
  is_aggregated <- "n_samples" %in% colnames(sample_info)

  colors <- character(n_samples)
  names(colors) <- sample_ids

  # Get default palettes
  group_colors <- AETHER_default_palette("groups")
  omics_colors <- AETHER_default_palette("omics")

  if (is_aggregated) {
    # Use group colors directly
    for (i in seq_len(n_samples)) {
      group <- sample_info$group[i]
      if (group %in% names(group_colors)) {
        colors[i] <- group_colors[group]
      } else {
        omics_type <- sample_info$omics[i]
        if (omics_type %in% names(omics_colors)) {
          colors[i] <- omics_colors[omics_type]
        } else {
          colors[i] <- grDevices::rainbow(n_samples)[i]
        }
      }
    }
  } else {
    # Individual samples: generate gradients per group
    groups <- unique(sample_info$group)

    for (g in groups) {
      group_idx <- which(sample_info$group == g)
      n_in_group <- length(group_idx)

      # Get base color for group
      if (g %in% names(group_colors)) {
        base_color <- group_colors[g]
      } else {
        omics_type <- sample_info$omics[group_idx[1]]
        if (omics_type %in% names(omics_colors)) {
          base_color <- omics_colors[omics_type]
        } else {
          base_color <- grDevices::rainbow(length(groups))[which(groups == g)]
        }
      }

      # Generate gradient
      if (n_in_group == 1) {
        colors[group_idx] <- base_color
      } else {
        gradient <- AETHER_generate_gradient(base_color, n = n_in_group,
                                              direction = "to_light")
        colors[group_idx] <- gradient
      }
    }
  }

  return(colors)
}