# ==============================================================================
# KHALKOS - Shared Utilities
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Bronze (khalkos) — the fundamental working material of Hephaestus's forge.
# Every creation depends on it. KHALKOS provides the shared helpers,
# object manipulation tools, and configuration that underpin all other
# modules. Direct parallel to GAIA's DEMETER module.
#
# Core responsibilities:
#   - SingleCellExperiment object manipulation helpers
#   - Common data structure utilities
#   - Configuration and parameter validation
#   - Logging and verbose output helpers
#   - Colour palettes and theme defaults
#   - Miscellaneous helpers shared across modules
#
# All functions prefixed: KHALKOS_
# ==============================================================================


# ── Internal constants ────────────────────────────────────────────────────────

# 20-colour hand-curated categorical palette (Tableau + extensions, distinct).
# Must stay identical to ASPIS's .aspis_discrete_palette() (aspis_embedding.R)
# -- that function delegates to KHALKOS_default_palette() specifically so
# there is only ever one canonical categorical palette in CAULDRON. Do not
# fork this list; every ASPIS plotting function's un-styled default (no
# palette= supplied) and every KHALKOS_assign_metadata_colours() map must
# keep producing the same colours for the same position, or a category's
# colour silently disagrees between plots that pass an explicit named
# palette and plots that don't.
.khalkos_cat_palette_20 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#79706E", "#D4A6C8", "#86BCB6", "#FFBE7D", "#8CD17D",
  "#499894", "#E6D16A", "#D37295", "#FABFD2", "#B6992D"
)


# ==============================================================================
# Object validation & safe access
# ==============================================================================

#' Validate a SingleCellExperiment object
#'
#' Checks that an object is a \code{SingleCellExperiment} and optionally that
#' required assays, reduced dimensions, \code{colData} columns, and
#' \code{metadata} keys are present.  All failures are collected and reported
#' together in a single informative error.
#'
#' @param sce A \code{SingleCellExperiment} to validate.
#' @param require_assays Character vector of assay names that must be present.
#' @param require_reduceddims Character vector of \code{reducedDim} names that
#'   must be present.
#' @param require_coldata Character vector of \code{colData} column names that
#'   must be present.
#' @param require_metadata Character vector of \code{metadata} keys that must
#'   be present.
#' @param verbose Logical. Print a short confirmation on success.
#'   Default \code{TRUE}.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass; throws an
#'   informative error otherwise.
#' @export
KHALKOS_validate_sce <- function(sce,
                                  require_assays      = NULL,
                                  require_reduceddims = NULL,
                                  require_coldata     = NULL,
                                  require_metadata    = NULL,
                                  verbose             = TRUE) {

  if (!is(sce, "SingleCellExperiment"))
    stop("Object is not a SingleCellExperiment (class: ",
         paste(class(sce), collapse = ", "), ").", call. = FALSE)

  errors <- character(0)

  if (!is.null(require_assays)) {
    miss <- setdiff(require_assays, assayNames(sce))
    if (length(miss))
      errors <- c(errors, paste0(
        "Missing assay(s): ", paste(miss, collapse = ", "),
        ".  Available: ", paste(assayNames(sce), collapse = ", ")
      ))
  }

  if (!is.null(require_reduceddims)) {
    miss <- setdiff(require_reduceddims, reducedDimNames(sce))
    if (length(miss))
      errors <- c(errors, paste0(
        "Missing reducedDim(s): ", paste(miss, collapse = ", "),
        ".  Available: ", paste(reducedDimNames(sce), collapse = ", ")
      ))
  }

  if (!is.null(require_coldata)) {
    miss <- setdiff(require_coldata, names(colData(sce)))
    if (length(miss))
      errors <- c(errors, paste0(
        "Missing colData column(s): ", paste(miss, collapse = ", ")
      ))
  }

  if (!is.null(require_metadata)) {
    miss <- setdiff(require_metadata, names(metadata(sce)))
    if (length(miss))
      errors <- c(errors, paste0(
        "Missing metadata key(s): ", paste(miss, collapse = ", ")
      ))
  }

  if (length(errors))
    stop("SCE validation failed:\n",
         paste0("  - ", errors, collapse = "\n"), call. = FALSE)

  if (isTRUE(verbose))
    message("SCE validated: ", format(ncol(sce), big.mark = ","),
            " cells \u00d7 ", format(nrow(sce), big.mark = ","), " genes.")

  invisible(TRUE)
}


#' Safe assay retrieval
#'
#' Returns \code{assay(sce, assay_name)} with a clear error message listing
#' available assays when the requested name is not found.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param assay_name Character scalar. Name of the assay to retrieve.
#' @param error_if_missing Logical. If \code{TRUE} (default), throw an error
#'   when the assay is absent.  If \code{FALSE}, return \code{NULL} silently.
#'
#' @return The assay matrix, or \code{NULL} if missing and
#'   \code{error_if_missing = FALSE}.
#' @export
KHALKOS_get_assay <- function(sce, assay_name, error_if_missing = TRUE) {
  if (assay_name %in% assayNames(sce))
    return(assay(sce, assay_name))

  if (isTRUE(error_if_missing))
    stop("Assay '", assay_name, "' not found in SCE.  ",
         "Available: ", paste(assayNames(sce), collapse = ", "), call. = FALSE)

  NULL
}


#' Safe assay assignment
#'
#' A thin wrapper around \code{assay(sce, assay_name) <- value} that returns
#' the modified SCE, making it pipe-friendly.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param assay_name Character scalar. Name of the assay to set (created if
#'   absent, replaced if present).
#' @param value Matrix to store.
#'
#' @return SCE with the assay updated.
#' @export
KHALKOS_set_assay <- function(sce, assay_name, value) {
  assay(sce, assay_name) <- value
  sce
}


# ==============================================================================
# Colour palettes
# ==============================================================================

#' Default CAULDRON colour palettes
#'
#' Returns a vector of \code{n} colours from the CAULDRON palette suite.
#' Categorical palettes are based on a hand-curated 20-colour set extended
#' via \code{colorRampPalette} for \code{n > 20}.  Sequential and diverging
#' palettes are interpolated between key anchor colours.
#'
#' @param n Integer. Number of colours to return.
#' @param type Character. Palette type: \code{"categorical"} (default),
#'   \code{"sequential"}, or \code{"diverging"}.
#'
#' @return Character vector of \code{n} hex colour codes.
#' @export
KHALKOS_default_palette <- function(n, type = c("categorical", "sequential",
                                                  "diverging")) {
  type <- match.arg(type)
  n    <- as.integer(n)

  if (type == "categorical") {
    base <- .khalkos_cat_palette_20
    if (n <= length(base)) return(base[seq_len(n)])
    return(colorRampPalette(base)(n))
  }

  if (type == "sequential") {
    # Purple → teal → yellow (viridis-inspired)
    return(colorRampPalette(
      c("#440154", "#31688E", "#35B779", "#FDE725")
    )(n))
  }

  # diverging: blue → white → red
  colorRampPalette(c("#313695", "#F7F7F7", "#D73027"))(n)
}


#' Assign a stable colour to every level of qualifying colData columns
#'
#' CAULDRON's plotting functions otherwise assign discrete colours
#' positionally, fresh on every call (see \code{\link{KHALKOS_default_palette}}
#' and how e.g. \code{ASPIS_plot_umap} uses it) -- so the same category (e.g.
#' cluster \code{"4"}, genotype \code{"KO"}) can end up a different colour in
#' two different plots if the set/order of levels present differs between
#' calls. This builds one colour map per qualifying \code{colData(sce)}
#' column instead, meant to be computed once and threaded through to every
#' plotting call in a report (e.g. \code{ASPIS_plot_umap(palette = ...)},
#' \code{KERAUNOS_propeller_proportions(cluster_colors = ..., group_colors =
#' ...)}) so the same category renders identically everywhere.
#'
#' Only character/factor columns with at most \code{max_levels} unique
#' values qualify -- this excludes continuous metadata (e.g. a UMI count
#' column) by type, and high-cardinality identifier columns (e.g. cell
#' barcodes) by level count, without needing either named explicitly.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param palette Character vector of hex colours to draw from, recycled if
#'   a column has more levels than colours supplied. Default \code{NULL}
#'   uses \code{\link{KHALKOS_default_palette}}.
#' @param max_levels Integer. Only \code{colData} columns with at most this
#'   many unique values are assigned colours. Default \code{30}.
#' @param overrides Named list, one entry per \code{colData} column to
#'   override (need not cover every column). Each entry is itself a named
#'   character vector (category label -> hex colour); only the categories
#'   named are overwritten -- the rest of that column keeps its
#'   palette-assigned colour. Default: empty list (no overrides).
#' @param verbose Logical. Print which columns were included/skipped.
#'   Default \code{TRUE}.
#'
#' @return Named list, one entry per qualifying \code{colData} column, each
#'   a named character vector (category label -> hex colour). A column
#'   absent from the result (excluded by type or \code{max_levels}) simply
#'   isn't a name in this list -- indexing it (e.g. \code{colours$gene_count})
#'   returns \code{NULL}, which every consumer here already treats as "use
#'   the default palette instead".
#' @export
KHALKOS_assign_metadata_colours <- function(sce, palette = NULL, max_levels = 30,
                                             overrides = list(), verbose = TRUE) {

  if (!inherits(sce, "SingleCellExperiment"))
    stop("'sce' must be a SingleCellExperiment.", call. = FALSE)

  col_data <- as.data.frame(colData(sce))

  is_categorical <- function(x) {
    (is.character(x) || is.factor(x)) && length(unique(x)) <= max_levels
  }
  candidate_cols <- names(col_data)[vapply(col_data, is_categorical, logical(1))]
  skipped_cols   <- setdiff(names(col_data), candidate_cols)

  if (isTRUE(verbose)) {
    cat("[KHALKOS] Assigning metadata colours\n")
    cat("    Included (character/factor, <=", max_levels, "levels):",
        paste(candidate_cols, collapse = ", "), "\n")
    cat("    Skipped:", paste(skipped_cols, collapse = ", "), "\n\n")
  }

  colours <- lapply(candidate_cols, function(col) {
    values <- as.character(col_data[[col]])
    lv <- unique(values)
    # Numeric-looking levels (e.g. cluster IDs "1".."11") sort numerically;
    # everything else sorts alphabetically -- matches the level-ordering
    # convention already used by TALARIA_export_gene_expr_widget(), so
    # colour assignment is stable/reproducible across a report.
    lv_num <- suppressWarnings(as.numeric(lv))
    lv <- if (!anyNA(lv_num)) as.character(sort(lv_num)) else sort(lv)

    base_pal <- if (!is.null(palette)) palette else KHALKOS_default_palette(length(lv))
    base_pal <- rep_len(base_pal, length(lv))
    map <- stats::setNames(base_pal, lv)

    if (!is.null(overrides[[col]]))
      map[names(overrides[[col]])] <- overrides[[col]]

    map
  })
  names(colours) <- candidate_cols

  colours
}


# ==============================================================================
# ggplot2 theme
# ==============================================================================

#' CAULDRON ggplot2 theme
#'
#' A clean, publication-ready ggplot2 theme: white background, subtle grey
#' gridlines and borders, no minor grid, and bold plot titles.  Use as a
#' drop-in replacement for \code{theme_bw()}.
#'
#' @param base_size Numeric. Base font size in points. Default \code{11}.
#' @param base_family Character. Base font family. Default \code{""}.
#'
#' @return A ggplot2 \code{theme} object.
#' @export
KHALKOS_theme_cauldron <- function(base_size = 11, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
      panel.border      = element_rect(colour = "grey70", fill = NA,
                                       linewidth = 0.5),
      axis.ticks        = element_line(colour = "grey70", linewidth = 0.3),
      strip.background  = element_rect(fill = "grey95", colour = "grey70"),
      strip.text        = element_text(size = rel(0.9)),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.title        = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle     = element_text(colour = "grey40", size = rel(0.9))
    )
}


# ==============================================================================
# Logging / messaging
# ==============================================================================

#' Verbose-aware logging helper
#'
#' Emits a message, warning, or error depending on \code{level}.  At the
#' \code{"info"} level the message is only printed when \code{verbose = TRUE}.
#'
#' @param msg Character. Message text.
#' @param verbose Logical. Print info-level messages. Default \code{TRUE}.
#' @param level Character. One of \code{"info"} (default), \code{"warning"},
#'   or \code{"error"}.  \code{"error"} always throws; \code{"warning"} always
#'   warns.
#'
#' @return Invisibly \code{NULL}.
#' @export
KHALKOS_log_message <- function(msg, verbose = TRUE,
                                 level = c("info", "warning", "error")) {
  level <- match.arg(level)
  if (level == "error")  stop(msg,    call. = FALSE)
  if (level == "warning") warning(msg, call. = FALSE)
  if (isTRUE(verbose) && level == "info") message(msg)
  invisible(NULL)
}


# ==============================================================================
# Package availability check
# ==============================================================================

#' Check optional package availability
#'
#' Tests whether each package in \code{packages} can be loaded with
#' \code{requireNamespace()}.  Missing packages are reported together with an
#' installation hint.
#'
#' @param packages Character vector of package names to check.
#' @param install_hint Character scalar.  Custom install instructions printed
#'   in the error/warning.  When \code{NULL} (default), a generic
#'   \code{BiocManager::install()} call is suggested.
#' @param error Logical. If \code{TRUE} (default), throw an error for missing
#'   packages.  If \code{FALSE}, issue a warning instead.
#'
#' @return Invisibly \code{TRUE} if all packages are available; invisibly
#'   \code{FALSE} otherwise (only when \code{error = FALSE}).
#' @export
KHALKOS_check_packages <- function(packages, install_hint = NULL,
                                    error = TRUE) {
  missing_pkgs <- packages[!vapply(packages, requireNamespace,
                                   logical(1L), quietly = TRUE)]

  if (length(missing_pkgs) == 0L) return(invisible(TRUE))

  hint <- if (is.null(install_hint)) {
    paste0('BiocManager::install(c(',
           paste0('"', missing_pkgs, '"', collapse = ", "), '))')
  } else {
    install_hint
  }

  msg <- paste0(
    "Required package(s) not available: ",
    paste(missing_pkgs, collapse = ", "),
    "\n  Install with: ", hint
  )

  if (isTRUE(error)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  invisible(FALSE)
}


# ==============================================================================
# SCE merging
# ==============================================================================

#' Merge multiple SingleCellExperiment objects
#'
#' Column-binds a list of \code{SingleCellExperiment} objects into a single
#' combined SCE.  Optionally adds a \code{colData} column identifying which
#' input object each cell came from.
#'
#' All SCEs must share the same genes (rows) in the same order.  If assay
#' names differ across objects, \code{cbind} will fill missing assays with
#' \code{NA} — ensure all objects have been processed to the same point before
#' merging.
#'
#' @param sce_list A named or unnamed list of \code{SingleCellExperiment}
#'   objects.  Must contain at least 2 elements.
#' @param sample_col Character scalar or \code{NULL}.  If non-\code{NULL}, a
#'   new \code{colData} column with this name is added to each SCE before
#'   merging, recording which sample it came from.  Default \code{"sample_id"}.
#' @param sample_names Character vector of length \code{length(sce_list)}.
#'   Labels for each SCE.  Defaults to \code{names(sce_list)} or
#'   \code{"sample_1"}, \code{"sample_2"}, … if unnamed.
#' @param verbose Logical. Print a summary. Default \code{TRUE}.
#'
#' @return A merged \code{SingleCellExperiment}.
#' @export
KHALKOS_merge_sce <- function(sce_list,
                               sample_col   = "sample_id",
                               sample_names = NULL,
                               verbose      = TRUE) {

  if (!is.list(sce_list) || length(sce_list) < 2L)
    stop("'sce_list' must be a list of at least 2 SingleCellExperiment objects.",
         call. = FALSE)

  for (i in seq_along(sce_list)) {
    if (!is(sce_list[[i]], "SingleCellExperiment"))
      stop("Element ", i, " of sce_list is not a SingleCellExperiment.",
           call. = FALSE)
  }

  if (is.null(sample_names)) {
    sample_names <- if (!is.null(names(sce_list)))
      names(sce_list)
    else
      paste0("sample_", seq_along(sce_list))
  }

  if (length(sample_names) != length(sce_list))
    stop("'sample_names' must have the same length as 'sce_list'.",
         call. = FALSE)

  # Tag each SCE with its sample name before merging
  if (!is.null(sample_col)) {
    sce_list <- lapply(seq_along(sce_list), function(i) {
      sce <- sce_list[[i]]
      colData(sce)[[sample_col]] <- sample_names[i]
      sce
    })
  }

  merged <- do.call(BiocGenerics::cbind, sce_list)

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\u2500\u2500 KHALKOS: Merged %d SCE objects %s\n  Cells   : %s\n  Genes   : %s\n  Samples : %s\n%s\n",
      length(sce_list),
      strrep("\u2500", 30),
      format(ncol(merged), big.mark = ","),
      format(nrow(merged), big.mark = ","),
      paste(sample_names, collapse = ", "),
      strrep("\u2500", 56)
    ))
  }

  merged
}
