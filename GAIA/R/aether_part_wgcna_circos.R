# GAIA/Aether/part_wgcna_circos.R
# Circos plot comparing PART clustering with WGCNA modules
#
# Requires: circlize, ComplexHeatmap (grid)


#' Plot PART-WGCNA comparison circos
#'
#' Creates a circos plot showing the relationship between PART clustering
#' (DEA-based) and WGCNA modules. PART clusters occupy 270 degrees showing
#' per-gene L2FC values, WGCNA modules occupy 90 degrees showing trait
#' correlations. Chords connect clusters to modules based on shared genes.
#'
#' @param circos_data A \code{part_wgcna_circos} object from
#'   \code{DEMETER_prepare_part_wgcna_circos()}.
#' @param output_file Character. Output file path. Extension determines format
#'   (.pdf or .png). Default: "PART_WGCNA_circos.pdf".
#' @param plot_width Numeric. Plot width in inches. Default: 15.
#' @param plot_height Numeric. Plot height in inches. Default: 15.
#' @param data_track_height Numeric. Height of each data track (expression/trait).
#'   Default: 0.05.
#' @param label_track_height Numeric. Height of the innermost label track.
#'   Default: 0.03.
#' @param expression_limit Numeric or NULL. Symmetric cap for the L2FC color
#'   scale. If numeric (e.g., 3), the scale runs from -3 to 3 and values beyond
#'   are clamped to the extreme colors. If NULL (default), uses the full data
#'   range.
#' @param expression_colors Character vector of 2-3 colors for L2FC values
#'   (low, center, high). Default: c("blue", "white", "red").
#' @param trait_colors Character vector of 2-3 colors for trait correlations.
#'   Default: c("#276419", "white", "#8e0152").
#' @param na_color Character. Color for NA expression values. Default: "white".
#' @param show_trait_values Logical. Show correlation values as text in module
#'   sectors. Default: FALSE.
#' @param show_trait_significance Logical. Show significance markers in module
#'   sectors. Default: TRUE.
#' @param show_sector_borders Logical. Show borders on label track sectors.
#'   Default: TRUE.
#' @param show_legend Logical. Add a legend page to the output. Default: TRUE.
#' @param show_track_labels Logical. Add text labels in the gap identifying each
#'   track (sample names for PART, trait names for WGCNA). Experimental — may
#'   need adjustment depending on track count and label length. Default: FALSE.
#' @param track_label_cex Numeric. Text size for track labels. Default: 0.4.
#' @param position_grouping_threshold Numeric (0-1). Fraction of module sector
#'   within which same-cluster genes are grouped into a single chord.
#'   Default: 0.25.
#' @param min_chord_width_pct Numeric (0-1). Minimum chord width as fraction
#'   of module sector size. Default: 0.1.
#' @param min_module_visual_size Integer. Minimum visual size for small modules
#'   so they remain visible. Default: 1000.
#' @param part_degrees Numeric. Degrees allocated to PART section. Default: 270.
#' @param wgcna_degrees Numeric. Degrees allocated to WGCNA section. Default: 90.
#' @param gap_degrees Numeric. Gap between PART and WGCNA sections. Default: 15.
#' @param gap_after_track Integer or NULL. Insert an empty spacer track after
#'   this data track number (1-indexed). Affects all sectors equally. Useful for
#'   visually grouping tracks. Default: NULL (no gap).
#' @param gap_track_height Numeric. Height of the spacer track. Default: 0.02.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisible NULL. Plot is saved to \code{output_file}.
#'
#' @examples
#' \dontrun{
#' AETHER_plot_part_wgcna_circos(circos_data, output_file = "circos.pdf")
#'
#' }
#' @export
AETHER_plot_part_wgcna_circos <- function(circos_data,
                                           output_file = "PART_WGCNA_circos.pdf",
                                           plot_width = 15,
                                           plot_height = 15,
                                           data_track_height = 0.05,
                                           label_track_height = 0.03,
                                           expression_limit = NULL,
                                           expression_colors = c("blue", "white", "red"),
                                           trait_colors = c("#276419", "white", "#8e0152"),
                                           na_color = "white",
                                           show_trait_values = FALSE,
                                           show_trait_significance = TRUE,
                                           show_sector_borders = TRUE,
                                           show_legend = TRUE,
                                           show_track_labels = FALSE,
                                           track_label_cex = 0.4,
                                           position_grouping_threshold = 0.25,
                                           min_chord_width_pct = 0.1,
                                           min_module_visual_size = 1000,
                                           part_degrees = 270,
                                           wgcna_degrees = 90,
                                           gap_degrees = 15,
                                           gap_after_track = NULL,
                                           gap_track_height = 0.02,
                                           verbose = TRUE) {

  if (!inherits(circos_data, "part_wgcna_circos")) {
    stop("'circos_data' must be a part_wgcna_circos object")
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install with: install.packages('circlize')")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required.")
  }

  # =========================================================================
  # Data extraction
  # =========================================================================

  PART_df <- circos_data$PART_df
  modules_df <- circos_data$module_df
  association_matrix <- circos_data$association_matrix
  module_gene_df <- circos_data$module_gene
  cluster_colors <- circos_data$cluster_colors
  module_colors <- circos_data$module_colors

  genes <- rownames(PART_df)
  cluster_names <- sort(unique(PART_df$cluster))
  module_names <- rownames(modules_df)

  # Hide dummy module from display entirely
  if ("dummy" %in% module_names) {
    if (verbose) cat("[AETHER] Hiding dummy module from display\n")
    module_names <- module_names[module_names != "dummy"]
    modules_df <- modules_df[module_names, , drop = FALSE]
    module_colors <- module_colors[names(module_colors) != "dummy"]
  }

  # Detect sample columns (L2FC data) and trait columns
  sample_cols <- setdiff(colnames(PART_df), c("cluster", "module"))
  n_samples <- length(sample_cols)

  all_mod_cols <- colnames(modules_df)
  trait_cols <- all_mod_cols[!all_mod_cols %in% c("module", "size") &
                             !grepl("_sig$", all_mod_cols)]
  n_traits <- length(trait_cols)
  total_tracks <- max(n_samples, n_traits)

  # Validate gap_after_track
  has_gap <- !is.null(gap_after_track)
  if (has_gap) {
    gap_after_track <- as.integer(gap_after_track)
    if (gap_after_track < 1 || gap_after_track >= total_tracks) {
      warning("gap_after_track must be between 1 and total_tracks - 1. Ignoring.")
      has_gap <- FALSE
    }
  }

  cluster_sizes <- table(PART_df$cluster)[cluster_names]

  if (verbose) {
    cat("[AETHER] Plotting PART-WGCNA Circos\n")
    cat("Clusters:", length(cluster_names), "| Modules:", length(module_names), "\n")
    cat("Sample tracks:", n_samples, "| Trait tracks:", n_traits, "\n")
  }

  # =========================================================================
  # Sector layout: 270/90 degree split
  # =========================================================================

  sectors <- c(cluster_names, module_names)

  # Apply minimum visual size to modules (except dummy)
  modules_df_adj <- modules_df
  for (module in module_names) {
    if (module != "dummy" && modules_df_adj[module, "size"] < min_module_visual_size) {
      modules_df_adj[module, "size"] <- min_module_visual_size
    }
  }

  total_part_genes <- sum(cluster_sizes)
  total_wgcna_genes <- sum(modules_df_adj$size)

  # Scaling to achieve the degree split
  wgcna_scale <- 1
  part_scale <- (part_degrees / wgcna_degrees) *
    (total_wgcna_genes / total_part_genes) * wgcna_scale

  # Build xlim matrix
  xlim_mat <- matrix(0, ncol = 2, nrow = length(sectors),
                     dimnames = list(sectors, NULL))
  for (cluster in cluster_names) {
    xlim_mat[cluster, ] <- c(0, cluster_sizes[cluster] * part_scale)
  }
  for (module in module_names) {
    xlim_mat[module, ] <- c(0, modules_df_adj[module, "size"] * wgcna_scale)
  }

  # Remove zero-size modules
  zero_modules <- module_names[modules_df_adj[module_names, "size"] == 0]
  if (length(zero_modules) > 0) {
    if (verbose) cat("[AETHER] Removing zero-size modules:", paste(zero_modules, collapse = ", "), "\n")
    sectors <- sectors[!sectors %in% zero_modules]
    module_names <- module_names[!module_names %in% zero_modules]
    xlim_mat <- xlim_mat[sectors, , drop = FALSE]
  }

  # Gap degrees: gaps only between the two sections
  gap_deg <- c(
    rep(0, length(cluster_names) - 1),
    gap_degrees,
    rep(0, length(module_names) - 1),
    gap_degrees
  )

  # =========================================================================
  # Open device
  # =========================================================================

  ext <- tolower(tools::file_ext(output_file))
  if (ext == "pdf") {
    grDevices::pdf(file = output_file, width = plot_width, height = plot_height)
  } else if (ext == "png") {
    grDevices::png(filename = output_file,
                   width = plot_width, height = plot_height,
                   units = "in", res = 300)
  } else {
    grDevices::pdf(file = output_file, width = plot_width, height = plot_height)
  }

  # =========================================================================
  # Initialize circos
  # =========================================================================

  circlize::circos.clear()
  circlize::circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.001, 0.001),
    start.degree = 90,
    gap.degree = gap_deg
  )
  circlize::circos.initialize(factors = sectors, xlim = xlim_mat)

  # =========================================================================
  # Color functions
  # =========================================================================

  # Expression color function (symmetric around 0)
  all_exp_vals <- unlist(lapply(sample_cols, function(col) PART_df[[col]]))
  all_exp_vals <- all_exp_vals[!is.na(all_exp_vals)]
  max_abs_exp <- if (!is.null(expression_limit)) expression_limit
                 else max(abs(all_exp_vals), na.rm = TRUE)

  if (length(expression_colors) == 3) {
    exp_col_fun <- circlize::colorRamp2(
      c(-max_abs_exp, 0, max_abs_exp), expression_colors
    )
  } else if (length(expression_colors) == 2) {
    exp_col_fun <- circlize::colorRamp2(
      c(-max_abs_exp, max_abs_exp), expression_colors
    )
  } else {
    stop("expression_colors must be 2 or 3 colors")
  }

  # Trait color function (correlation: -1 to 1)
  if (length(trait_colors) == 3) {
    trait_col_fun <- circlize::colorRamp2(c(-1, 0, 1), trait_colors)
  } else if (length(trait_colors) == 2) {
    trait_col_fun <- circlize::colorRamp2(c(-1, 1), trait_colors)
  } else {
    stop("trait_colors must be 2 or 3 colors")
  }

  # =========================================================================
  # Data tracks
  # =========================================================================

  # Alignment: samples go outer, traits go inner (or vice versa if one has more)
  sample_offset <- max(0, n_traits - n_samples)
  trait_offset <- max(0, n_samples - n_traits)

  circos_track_idx <- 0
  # Map from data track number (1..total_tracks) to circos track index
  data_track_to_circos <- integer(total_tracks)

  for (track_idx in seq_len(total_tracks)) {
    circos_track_idx <- circos_track_idx + 1
    data_track_to_circos[track_idx] <- circos_track_idx

    sample_idx <- track_idx - sample_offset
    trait_idx <- track_idx - trait_offset

    has_sample <- sample_idx >= 1 && sample_idx <= n_samples
    has_trait <- trait_idx >= 1 && trait_idx <= n_traits

    if (has_sample) current_sample_col <- sample_cols[sample_idx]
    if (has_trait) {
      current_trait_col <- trait_cols[trait_idx]
      current_sig_col <- paste0(current_trait_col, "_sig")
    }

    circlize::circos.trackPlotRegion(
      track.index = circos_track_idx,
      ylim = c(0, 1),
      bg.border = NA,
      track.height = data_track_height,
      panel.fun = function(x, y) {
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")

        if (sector %in% cluster_names) {
          # --- PART cluster: per-gene expression ---
          if (has_sample) {
            cluster_genes <- rownames(PART_df)[PART_df$cluster == sector]
            pos <- 0
            for (gene in cluster_genes) {
              gene_exp <- PART_df[gene, current_sample_col]
              col <- if (is.na(gene_exp)) na_color else exp_col_fun(gene_exp)
              circlize::circos.rect(pos, 0, pos + part_scale, 1,
                                    col = col, border = NA)
              pos <- pos + part_scale
            }
            circlize::circos.rect(xlim[1], 0, xlim[2], 1,
                                  col = NA, border = "black", lwd = 0.5)
          }

        } else if (sector %in% module_names) {
          # --- WGCNA module: trait correlation ---
          if (has_trait) {
            trait_val <- modules_df[sector, current_trait_col]
            col <- trait_col_fun(trait_val)
            border_col <- if (show_sector_borders) "black" else NA
            circlize::circos.rect(xlim[1], 0, xlim[2], 1,
                                  col = col, border = border_col, lwd = 0.5)

            # Text label
            trait_text <- ""
            if (show_trait_values) {
              trait_text <- as.character(round(trait_val, 2))
            }

            has_sig_col <- current_sig_col %in% colnames(modules_df)
            if (show_trait_significance && has_sig_col) {
              is_sig <- modules_df[sector, current_sig_col]
              if (isTRUE(is_sig)) {
                trait_text <- if (nchar(trait_text) > 0) {
                  paste0(trait_text, " *")
                } else {
                  "*"
                }
              }
            }

            if (nchar(trait_text) > 0) {
              circlize::circos.text(mean(xlim), 0.5, labels = trait_text,
                                    facing = "bending.inside", cex = 0.6,
                                    col = "black", font = 1)
            }
          }
        }
      }
    )

    # Insert spacer track after the specified data track
    if (has_gap && track_idx == gap_after_track) {
      circos_track_idx <- circos_track_idx + 1
      circlize::circos.trackPlotRegion(
        track.index = circos_track_idx,
        ylim = c(0, 1),
        bg.border = NA,
        track.height = gap_track_height,
        panel.fun = function(x, y) { invisible(NULL) }
      )
      if (verbose) cat("[AETHER] Spacer track inserted after data track", gap_after_track, "\n")
    }
  }

  if (verbose) cat("[AETHER] Data tracks drawn\n")

  # =========================================================================
  # Track labels (experimental — set show_track_labels = TRUE to enable)
  # Places text in the gap between PART and WGCNA sections to identify
  # which sample/trait each concentric track represents.
  # =========================================================================

  if (show_track_labels) {
    first_cluster <- cluster_names[1]
    first_module <- module_names[1]

    for (track_idx in seq_len(total_tracks)) {
      ci <- data_track_to_circos[track_idx]
      sample_idx <- track_idx - sample_offset
      trait_idx <- track_idx - trait_offset

      has_sample <- sample_idx >= 1 && sample_idx <= n_samples
      has_trait <- trait_idx >= 1 && trait_idx <= n_traits

      # PART track label (at gap before first cluster)
      # facing = "inside" reads radially (perpendicular to tracks)
      if (has_sample) {
        circlize::circos.text(
          x = 0, y = 0.5,
          labels = sample_cols[sample_idx],
          sector.index = first_cluster,
          track.index = ci,
          facing = "inside",
          niceFacing = TRUE,
          cex = track_label_cex,
          adj = c(1.1, 0.5)
        )
      }

      # WGCNA track label (at gap before first module)
      # facing = "outside" so text reads outward (away from plot) on this side
      if (has_trait) {
        circlize::circos.text(
          x = 0, y = 0.5,
          labels = trait_cols[trait_idx],
          sector.index = first_module,
          track.index = ci,
          facing = "outside",
          niceFacing = TRUE,
          cex = track_label_cex,
          adj = c(-0.1, 0.5)
        )
      }
    }

    if (verbose) cat("[AETHER] Track labels drawn\n")
  }

  # =========================================================================
  # Label track (innermost)
  # =========================================================================

  label_track_idx <- circos_track_idx + 1
  circlize::circos.trackPlotRegion(
    track.index = label_track_idx,
    ylim = c(0, 1),
    bg.border = NA,
    track.height = label_track_height,
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      border_col <- if (show_sector_borders) "black" else NA

      if (sector %in% cluster_names) {
        col <- cluster_colors[sector]
      } else {
        col <- module_colors[sector]
      }
      circlize::circos.rect(xlim[1], 0, xlim[2], 1,
                             col = col, border = border_col)
    }
  )

  if (verbose) cat("[AETHER] Label track drawn\n")

  # =========================================================================
  # Chords: position-based connections between clusters and modules
  # =========================================================================

  cluster_cumul_pos <- setNames(rep(0, length(cluster_names)), cluster_names)

  for (module in module_names) {
    # Get genes in this module (ordered as in module_gene_df)
    module_genes_ordered <- module_gene_df$gene[module_gene_df$module == module]
    n_genes_mod <- length(module_genes_ordered)
    if (n_genes_mod == 0) next

    # Map each gene to its PART cluster
    gene_positions <- data.frame(
      gene = module_genes_ordered,
      position = 0:(n_genes_mod - 1),
      position_frac = (0:(n_genes_mod - 1)) / n_genes_mod,
      cluster = NA_character_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(gene_positions))) {
      g <- gene_positions$gene[i]
      if (g %in% rownames(PART_df)) {
        gene_positions$cluster[i] <- as.character(PART_df[g, "cluster"])
      }
    }

    gene_positions <- gene_positions[!is.na(gene_positions$cluster), , drop = FALSE]
    if (nrow(gene_positions) == 0) next

    # --- Group consecutive same-cluster genes ---
    groups <- list()
    current_group <- list(
      cluster = gene_positions$cluster[1],
      start_frac = gene_positions$position_frac[1],
      end_frac = gene_positions$position_frac[1],
      n_genes = 1
    )

    if (nrow(gene_positions) > 1) {
      for (i in 2:nrow(gene_positions)) {
        same_cluster <- gene_positions$cluster[i] == current_group$cluster
        close_enough <- (gene_positions$position_frac[i] - current_group$end_frac) <=
          position_grouping_threshold

        if (same_cluster && close_enough) {
          current_group$end_frac <- gene_positions$position_frac[i]
          current_group$n_genes <- current_group$n_genes + 1
        } else {
          groups[[length(groups) + 1]] <- current_group
          current_group <- list(
            cluster = gene_positions$cluster[i],
            start_frac = gene_positions$position_frac[i],
            end_frac = gene_positions$position_frac[i],
            n_genes = 1
          )
        }
      }
    }
    groups[[length(groups) + 1]] <- current_group

    # --- Calculate chord positions with minimum width ---
    module_sector_size <- modules_df_adj[module, "size"] * wgcna_scale
    min_chord_width <- min_chord_width_pct * module_sector_size

    # Pass 1: initial positions
    group_positions <- lapply(groups, function(g) {
      start <- max(0, min(module_sector_size, g$start_frac * module_sector_size))
      end <- max(0, min(module_sector_size, g$end_frac * module_sector_size))
      list(group = g, initial_start = start, initial_end = end,
           final_start = start, final_end = end)
    })

    # Pass 2: expand to minimum width
    for (i in seq_along(group_positions)) {
      gp <- group_positions[[i]]
      current_width <- gp$initial_end - gp$initial_start

      if (current_width < min_chord_width) {
        needed <- min_chord_width - current_width

        space_before <- if (i > 1) {
          gp$initial_start - group_positions[[i - 1]]$initial_end
        } else {
          gp$initial_start
        }
        space_after <- if (i < length(group_positions)) {
          group_positions[[i + 1]]$initial_start - gp$initial_end
        } else {
          module_sector_size - gp$initial_end
        }

        exp_fwd <- min(space_after, needed / 2)
        exp_bwd <- min(space_before, needed / 2)

        # Asymmetric expansion if symmetric isn't enough
        total_exp <- exp_fwd + exp_bwd
        if (total_exp < needed) {
          remaining <- needed - total_exp
          add_fwd <- min(space_after - exp_fwd, remaining)
          exp_fwd <- exp_fwd + add_fwd
          remaining <- remaining - add_fwd
          if (remaining > 0) {
            exp_bwd <- exp_bwd + min(space_before - exp_bwd, remaining)
          }
        }

        gp$final_start <- max(0, gp$initial_start - exp_bwd)
        gp$final_end <- min(module_sector_size, gp$initial_end + exp_fwd)
        group_positions[[i]] <- gp
      }
    }

    # Pass 3: draw chords
    for (i in seq_along(group_positions)) {
      gp <- group_positions[[i]]
      cluster <- gp$group$cluster

      cluster_start <- cluster_cumul_pos[cluster]
      cluster_end <- cluster_cumul_pos[cluster] + (gp$group$n_genes * part_scale)

      chord_col <- circlize::add_transparency(module_colors[module], 0.5)

      circlize::circos.link(
        sector.index1 = cluster,
        point1 = c(cluster_start, cluster_end),
        sector.index2 = module,
        point2 = c(gp$final_start, gp$final_end),
        col = chord_col,
        border = NA
      )

      cluster_cumul_pos[cluster] <- cluster_end
    }
  }

  if (verbose) cat("[AETHER] Chords drawn\n")

  circlize::circos.clear()

  # =========================================================================
  # Legend (separate page)
  # =========================================================================

  if (show_legend) {
    grid::grid.newpage()

    legend_list <- list()

    # Expression color scale
    exp_at <- pretty(c(-max_abs_exp, max_abs_exp))
    legend_list[[1]] <- ComplexHeatmap::Legend(
      col_fun = exp_col_fun,
      at = exp_at,
      title = if (!is.null(expression_limit)) "L2FC (capped)" else "L2FC",
      direction = "horizontal",
      title_position = "topcenter",
      legend_width = grid::unit(4, "cm")
    )

    # Trait color scale
    if (n_traits > 0) {
      legend_list[[length(legend_list) + 1]] <- ComplexHeatmap::Legend(
        col_fun = trait_col_fun,
        at = c(-1, -0.5, 0, 0.5, 1),
        title = "Trait Correlation",
        direction = "horizontal",
        title_position = "topcenter",
        legend_width = grid::unit(4, "cm")
      )
    }

    # Module colors — match sector order (modules appear after clusters)
    ordered_modules <- intersect(sectors, names(module_colors))
    legend_list[[length(legend_list) + 1]] <- ComplexHeatmap::Legend(
      labels = ordered_modules,
      legend_gp = grid::gpar(fill = module_colors[ordered_modules]),
      title = "WGCNA Modules",
      ncol = 2
    )

    # Cluster colors — match sector order (clusters appear first)
    ordered_clusters <- intersect(sectors, names(cluster_colors))
    legend_list[[length(legend_list) + 1]] <- ComplexHeatmap::Legend(
      labels = ordered_clusters,
      legend_gp = grid::gpar(fill = cluster_colors[ordered_clusters]),
      title = "PART Clusters",
      ncol = 2
    )

    packed <- ComplexHeatmap::packLegend(
      list = legend_list,
      direction = "vertical",
      gap = grid::unit(0.5, "cm")
    )
    ComplexHeatmap::draw(packed,
                         x = grid::unit(0.5, "npc"),
                         y = grid::unit(0.5, "npc"),
                         just = c("center", "center"))
  }

  grDevices::dev.off()

  if (verbose) cat("[AETHER] Saved to:", output_file, "\n")

  invisible(NULL)
}
