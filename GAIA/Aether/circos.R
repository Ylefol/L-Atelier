library(ggplot2)
source('~/A_Projects/ZERO_DAWN/GAIA/Aether/colors.R')

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
#' @param track_margin Numeric. Margin between tracks (default = 0.01).
#' @param line_width Numeric. Width of signal lines (default = 1.0).
#' @param ideogram_height Numeric. Height of chromosome ideogram (default = 0.04).
#' @param title Character string. Plot title displayed in center (default = NULL).
#' @param width Numeric. Output width in inches (default = 10).
#' @param height Numeric. Output height in inches (default = 10).
#' @param res Numeric. Resolution for PNG output in DPI (default = 300).
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
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(track_margin, track_margin),
    cell.padding = c(0, 0, 0, 0)
  )

  # Initialize with chromosome labels only (no axis ticks)
  circlize::circos.initializeWithIdeogram(species = genome, chromosome.index = chr_sizes$chr)

  # Create signal tracks - one per track grouping
  # Each track has independent y-axis scaling
  for (t in seq_along(unique_tracks)) {
    track_name <- unique_tracks[t]

    # Get samples in this track
    track_sample_idx <- which(sample_info$track_id == track_name)
    track_sample_ids <- sample_ids[track_sample_idx]

    # Calculate y_max for this track independently
    track_y_max <- max(signal[, track_sample_idx, drop = FALSE], na.rm = TRUE)

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