#' Load Olink NPX Data
#'
#' @description Reads an Olink Explore NPX file (.parquet or .csv), optionally
#' joins plate layout metadata, and returns a structured \code{olink_data}
#' object containing all samples and proteins as-loaded.
#'
#' This function performs no QC filtering. To remove failed samples or
#' non-biological controls, pass the result to \code{HADES_filter_olink()}.
#'
#' @param npx_file Character. Path to the Olink NPX file (.parquet or .csv).
#' @param metadata_file Character or NULL. Path to plate layout metadata
#'   (.xlsx, .csv, or .tsv). Joined onto data by \code{sample_col}. Default: NULL.
#' @param metadata_sheet Integer or character. Sheet to read from .xlsx files.
#'   Default: 1.
#' @param sample_col Character. Column name present in both the NPX data and
#'   metadata used as the join key. Default: "SampleID".
#' @param npx_col Character. Which NPX representation to use as the primary
#'   value. One of \code{"NPX"}, \code{"ExtNPX"}, or \code{"PCNormalizedNPX"}.
#'   Default: \code{"NPX"}. NPX is log2-scale and approximately normally
#'   distributed — do NOT log-transform again before analysis.
#' @param restrict_to_metadata Logical. If \code{TRUE} and a
#'   \code{metadata_file} is provided, removes any samples from the NPX data
#'   that have no corresponding entry in the metadata. Useful when the NPX file
#'   contains additional samples (e.g. validation runs, reference samples) that
#'   are intentionally absent from the study metadata. Default: \code{FALSE}.
#' @param verbose Logical. Print loading progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"olink_data"} containing:
#' \describe{
#'   \item{data}{Long-format data.frame: SampleID, OlinkID, Assay, UniProt,
#'     Panel, NPX (the selected NPX column, renamed to "NPX"), AssayQC.
#'     Metadata columns are not replicated here — join from \code{$sample_meta}
#'     when needed.}
#'   \item{wide}{Numeric matrix: proteins (rows) x samples (cols).
#'     Rownames = OlinkID. Colnames = SampleID. NPX NAs propagated as-is.}
#'   \item{sample_meta}{Per-sample data.frame for all samples. Contains
#'     SampleID, SampleType, PlateID, WellID, SampleQC, plus any joined
#'     metadata columns.}
#'   \item{assay_meta}{Per-protein data.frame: OlinkID, Assay, UniProt, Panel,
#'     warn_fraction (proportion of samples with AssayQC == "WARN").}
#'   \item{params}{List of loader parameters used (for provenance).}
#' }
#'
#' @details
#' Olink process control probes (AssayType != "assay") are excluded during
#' parsing — these are internal instrument rows (ext_ctrl, inc_ctrl, amp_ctrl),
#' not proteins, and are not meaningful for analysis.
#'
#' All biological samples and QC control samples are retained in the returned
#' object. Use \code{HADES_filter_olink()} to remove SampleQC failures or
#' restrict to biological samples only.
#'
#' AssayQC == "WARN" proteins are flagged but retained. \code{$assay_meta}
#' carries a \code{warn_fraction} column (proportion of samples with WARN).
#'
#' Requires \code{arrow} for .parquet input; requires \code{readxl} for .xlsx
#' metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load only
#' ol_raw <- ELEUTHIA_load_olink("data.parquet")
#'
#' # Load with plate layout metadata, then filter
#' ol_raw <- ELEUTHIA_load_olink(
#'   npx_file      = "data.parquet",
#'   metadata_file = "plate_layout.xlsx"
#' )
#' ol <- HADES_filter_olink(ol_raw, keep_controls = FALSE, filter_qc = TRUE)
#' }
ELEUTHIA_load_olink <- function(npx_file,
                                 metadata_file        = NULL,
                                 metadata_sheet       = 1,
                                 sample_col           = "SampleID",
                                 npx_col              = "NPX",
                                 restrict_to_metadata = FALSE,
                                 verbose              = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (!file.exists(npx_file)) stop("npx_file not found: ", npx_file)

  npx_col <- match.arg(npx_col, c("NPX", "ExtNPX", "PCNormalizedNPX"))

  ext <- tolower(sub(".*\\.", "", basename(npx_file)))

  # ---------------------------------------------------------------------------
  # Read NPX file
  # ---------------------------------------------------------------------------

  if (verbose) cat("[ELEUTHIA] Reading NPX file:", npx_file, "\n")

  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required to read .parquet files.\n",
           "Install with: install.packages('arrow')")
    }
    raw <- as.data.frame(arrow::read_parquet(npx_file))
  } else if (ext == "csv") {
    raw <- utils::read.csv(npx_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported NPX file format: '.", ext, "'. Use .parquet or .csv.")
  }

  if (verbose) cat("[ELEUTHIA]   Rows read:", nrow(raw), " | Columns:", ncol(raw), "\n")

  # ---------------------------------------------------------------------------
  # Validate required columns
  # ---------------------------------------------------------------------------

  required_cols <- c(sample_col, "OlinkID", "Assay", "UniProt", "Panel",
                     npx_col, "AssayQC", "SampleQC", "AssayType", "SampleType")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("Required columns missing from NPX file: ",
         paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Remove Olink process control rows (structural, not QC filtering)
  # AssayType == "assay": actual proteins
  # ext_ctrl, inc_ctrl, amp_ctrl: internal instrument rows, not proteins
  # ---------------------------------------------------------------------------

  n_before <- nrow(raw)
  raw      <- raw[raw$AssayType == "assay", ]
  if (verbose && nrow(raw) < n_before) {
    cat("[ELEUTHIA]   Parsed", n_before - nrow(raw),
        "instrument control rows (AssayType != 'assay') as non-protein\n")
  }

  # Drop rows with missing or blank SampleID — can occur in some parquet exports
  # as internal calibration rows that survive the AssayType filter.
  valid_sid <- !is.na(raw[[sample_col]]) &
               nzchar(trimws(as.character(raw[[sample_col]])))
  n_invalid <- sum(!valid_sid)
  if (n_invalid > 0) {
    raw <- raw[valid_sid, ]
    if (verbose)
      cat("[ELEUTHIA]   Dropped", n_invalid, "rows with missing/blank", sample_col, "\n")
  }

  # ---------------------------------------------------------------------------
  # Read metadata file — kept separate, NOT merged into raw
  # Merging into the long-format raw would replicate every metadata value
  # once per protein per sample, causing severe memory amplification.
  # Metadata is joined at the sample level onto $sample_meta only.
  # ---------------------------------------------------------------------------

  meta <- NULL
  if (!is.null(metadata_file)) {
    if (!file.exists(metadata_file))
      stop("metadata_file not found: ", metadata_file)

    if (verbose) cat("[ELEUTHIA] Reading metadata:", metadata_file, "\n")

    meta_ext <- tolower(sub(".*\\.", "", basename(metadata_file)))

    if (meta_ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("Package 'readxl' is required to read .xlsx files.\n",
             "Install with: install.packages('readxl')")
      meta <- as.data.frame(readxl::read_excel(metadata_file, sheet = metadata_sheet))
    } else if (meta_ext == "csv") {
      meta <- utils::read.csv(metadata_file, stringsAsFactors = FALSE,
                               check.names = FALSE)
    } else if (meta_ext %in% c("tsv", "txt")) {
      meta <- utils::read.delim(metadata_file, stringsAsFactors = FALSE,
                                 check.names = FALSE)
    } else {
      stop("Unsupported metadata format: '.", meta_ext,
           "'. Use .xlsx, .csv, or .tsv.")
    }

    if (!sample_col %in% colnames(meta))
      stop("sample_col '", sample_col, "' not found in metadata file.")

    # Drop blank/NA rows common in CSV exports from Excel
    meta <- meta[!is.na(meta[[sample_col]]) &
                   nzchar(trimws(as.character(meta[[sample_col]]))), ]

    if (verbose)
      cat("[ELEUTHIA]   Metadata:", ncol(meta) - 1L, "non-key columns,",
          nrow(meta), "rows\n")
  }

  # ---------------------------------------------------------------------------
  # Optionally restrict NPX data to samples present in metadata
  # ---------------------------------------------------------------------------

  if (isTRUE(restrict_to_metadata)) {
    if (is.null(meta))
      stop("restrict_to_metadata = TRUE requires a metadata_file to be provided.")
    meta_ids    <- as.character(meta[[sample_col]])
    npx_ids     <- as.character(raw[[sample_col]])
    not_in_meta <- unique(npx_ids[!npx_ids %in% meta_ids])
    if (length(not_in_meta) > 0) {
      raw <- raw[npx_ids %in% meta_ids, ]
      if (verbose) {
        cat("[ELEUTHIA]   restrict_to_metadata: removed", length(not_in_meta),
            "samples not in metadata:\n")
        cat("[ELEUTHIA]    ", paste(not_in_meta, collapse = ", "), "\n")
      }
    } else if (verbose) {
      cat("[ELEUTHIA]   restrict_to_metadata: all NPX samples present in metadata\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Build $sample_meta (one row per sample) then join metadata onto it
  # ---------------------------------------------------------------------------

  smeta_cols  <- intersect(c(sample_col, "SampleType", "PlateID",
                              "WellID", "SampleQC"),
                           colnames(raw))
  sample_meta <- unique(raw[, smeta_cols, drop = FALSE])
  rownames(sample_meta) <- NULL

  if (!is.null(meta)) {
    new_cols    <- setdiff(colnames(meta), colnames(sample_meta))
    meta_join   <- meta[, c(sample_col, new_cols), drop = FALSE]
    sample_meta <- merge(sample_meta, meta_join, by = sample_col,
                         all.x = TRUE, sort = FALSE)
    if (verbose)
      cat("[ELEUTHIA]   Joined", length(new_cols), "metadata columns onto sample_meta\n")
  }
  rownames(sample_meta) <- sample_meta[[sample_col]]

  # ---------------------------------------------------------------------------
  # Build $data (long format — measurement columns only, no metadata)
  # ---------------------------------------------------------------------------

  core_cols <- c(sample_col, "OlinkID", "Assay", "UniProt",
                 "Panel", npx_col, "AssayQC")
  data_df   <- raw[, intersect(core_cols, colnames(raw)), drop = FALSE]

  if (npx_col != "NPX")
    colnames(data_df)[colnames(data_df) == npx_col] <- "NPX"
  rownames(data_df) <- NULL

  # ---------------------------------------------------------------------------
  # Build $wide (proteins x samples matrix)
  # ---------------------------------------------------------------------------

  samples  <- unique(data_df[[sample_col]])
  proteins <- unique(data_df$OlinkID)

  wide_mat <- matrix(NA_real_,
                     nrow     = length(proteins),
                     ncol     = length(samples),
                     dimnames = list(proteins, samples))

  prot_idx   <- match(data_df$OlinkID,       proteins)
  sample_idx <- match(data_df[[sample_col]], samples)
  wide_mat[cbind(prot_idx, sample_idx)] <- data_df$NPX

  # ---------------------------------------------------------------------------
  # Build $assay_meta
  # ---------------------------------------------------------------------------

  assay_cols <- intersect(c("OlinkID", "Assay", "UniProt", "Panel"),
                          colnames(raw))
  assay_meta <- unique(raw[, assay_cols, drop = FALSE])
  rownames(assay_meta) <- NULL

  warn_list  <- tapply(raw$AssayQC, raw$OlinkID,
                       function(x) mean(x == "WARN", na.rm = TRUE))
  warn_df    <- data.frame(OlinkID       = names(warn_list),
                            warn_fraction = as.numeric(warn_list),
                            stringsAsFactors = FALSE)
  assay_meta <- merge(assay_meta, warn_df, by = "OlinkID",
                      all.x = TRUE, sort = FALSE)
  assay_meta <- assay_meta[match(proteins, assay_meta$OlinkID), ]
  rownames(assay_meta) <- NULL

  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------

  params <- list(
    npx_file      = npx_file,
    metadata_file = metadata_file,
    sample_col    = sample_col,
    npx_col       = npx_col
  )

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = sample_meta,
    assay_meta  = assay_meta,
    params      = params
  )
  class(result) <- c("olink_data", "list")

  if (verbose) {
    cat("[ELEUTHIA] olink_data object ready: ",
        nrow(wide_mat), " proteins x ",
        ncol(wide_mat), " samples\n", sep = "")
    cat("[ELEUTHIA]   Use HADES_filter_olink() to remove controls / QC failures\n")
  }

  return(result)
}


#' Export Olink QC Results
#'
#' @description Writes all outputs from the Olink QC/exploration pipeline to a
#' structured directory: filtered data as RDS and CSVs, QC plots, a
#' PC-metadata association heatmap, per-variable PCA plots (coloured by the two
#' PCs most significantly associated with each variable), and a plain-text QC
#' summary.
#'
#' @param olink_data An \code{olink_data} object — the filtered, QC-cleaned
#'   object to export.
#' @param assoc An \code{artemis_pc_assoc} object from
#'   \code{ARTEMIS_pc_metadata_association()}.
#' @param qc_plots A named list from \code{AETHER_plot_olink_qc()}.
#' @param output_dir Character. Root directory for all outputs. Created if it
#'   does not exist.
#' @param top_n_vars Integer. Number of top variables per PC (by p-value) to
#'   include in the PCA plot loop. Default: 5.
#' @param p_threshold Numeric. Significance threshold used in the association
#'   heatmap border annotation and in the summary file. Default: 0.05.
#' @param heatmap_top_n Integer or NULL. Passed to
#'   \code{AETHER_plot_pc_association(top_n)} to restrict the heatmap to the
#'   top N variables per PC. NULL shows all tested variables. Default: NULL.
#' @param ntop_pca Integer or NULL. Number of most variable proteins used for
#'   PCA. NULL uses all proteins. Default: NULL.
#' @param batch_plots Named list of ggplot objects for batch-check PCA plots.
#'   Each element is saved as \code{qc/pca_batch_{name}.png}. Supply one plot
#'   coloured by PlateID for a standard batch check, or a before/after pair
#'   when batch correction has been applied. Default: NULL (no plots saved).
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisibly returns NULL. All outputs are written to disk.
#'
#' @details
#' \strong{Directory layout:}
#' \preformatted{
#' output_dir/
#'   data/
#'     olink_filtered.rds   -- filtered olink_data object
#'     sample_meta.csv      -- per-sample metadata (with outlier annotations)
#'     assay_meta.csv       -- per-protein metadata (warn_fraction etc.)
#'   qc/
#'     npx_distributions.png
#'     sample_qc.csv        -- SampleQC pass/fail table
#'     warn_proteins.png    -- only if present in qc_plots
#'     summary.txt          -- key QC numbers and parameter record
#'   association/
#'     pc_metadata_heatmap.png
#'   pca/
#'     <variable>.png       -- one per selected variable
#' }
#'
#' Each PCA plot uses the two PCs most significantly associated with that
#' variable (lowest p-value in \code{assoc$pvalues}). The PC pair is shown in
#' the plot title.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol)
#' ol <- HADES_detect_outliers_olink(ol, filter = TRUE)
#'
#' qc_plots <- AETHER_plot_olink_qc(ol, plot_warn_proteins = TRUE)
#' assoc    <- ARTEMIS_pc_metadata_association(ol$wide, ol$sample_meta,
#'                                              sample_col = "SampleID")
#'
#' ELEUTHIA_export_olink_qc(ol, assoc, qc_plots, output_dir = "results/Olink_QC")
#' }
ELEUTHIA_export_olink_qc <- function(olink_data,
                                      assoc,
                                      qc_plots,
                                      output_dir,
                                      top_n_vars    = 5L,
                                      p_threshold   = 0.05,
                                      heatmap_top_n = NULL,
                                      ntop_pca      = NULL,
                                      batch_plots   = NULL,
                                      qc_steps      = NULL,
                                      verbose       = TRUE) {

  if (!inherits(olink_data, "olink_data"))
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  if (!inherits(assoc, "artemis_pc_assoc"))
    stop("assoc must be an artemis_pc_assoc object from ARTEMIS_pc_metadata_association().")
  if (!is.list(qc_plots))
    stop("qc_plots must be a list from AETHER_plot_olink_qc().")

  sid_col    <- olink_data$params$sample_col
  top_n_vars <- as.integer(top_n_vars)

  # ---------------------------------------------------------------------------
  # Create directory structure
  # ---------------------------------------------------------------------------

  dir_data  <- file.path(output_dir, "data")
  dir_qc    <- file.path(output_dir, "qc")
  dir_assoc <- file.path(output_dir, "association")
  dir_pca   <- file.path(output_dir, "pca")

  for (d in c(dir_data, dir_qc, dir_assoc, dir_pca))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  if (verbose) cat("[ELEUTHIA] ELEUTHIA_export_olink_qc -> ", output_dir, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # 1. Data exports
  # ---------------------------------------------------------------------------

  saveRDS(olink_data, file.path(dir_data, "olink_filtered.rds"))
  utils::write.csv(olink_data$sample_meta,
                   file.path(dir_data, "sample_meta.csv"),
                   row.names = FALSE)
  utils::write.csv(olink_data$assay_meta,
                   file.path(dir_data, "assay_meta.csv"),
                   row.names = FALSE)

  if (verbose) cat("[ELEUTHIA]   [data] olink_filtered.rds, sample_meta.csv, assay_meta.csv\n")

  # ---------------------------------------------------------------------------
  # 2. QC plots
  # Accepts either a single AETHER_plot_olink_qc() list (backward compatible)
  # or a named list of such lists for multiple passes (e.g. pre/post filtering).
  # A single set is detected by the presence of $npx_distributions at the top
  # level; it is wrapped internally so the saving loop is uniform.
  # ---------------------------------------------------------------------------

  if (!is.null(qc_plots)) {
    qc_sets <- if (!is.null(qc_plots$npx_distributions)) {
      setNames(list(qc_plots), "")  # single set — no filename prefix
    } else {
      qc_plots                     # named list of sets
    }

    for (set_nm in names(qc_sets)) {
      pfx   <- if (nzchar(set_nm)) paste0(set_nm, "_") else ""
      plots <- qc_sets[[set_nm]]
      saved <- character(0L)

      if (!is.null(plots$npx_distributions)) {
        fname <- paste0(pfx, "npx_distributions.png")
        ggplot2::ggsave(file.path(dir_qc, fname),
                        plots$npx_distributions, width = 9, height = 5, dpi = 150)
        saved <- c(saved, fname)
      }

      if (!is.null(plots$sample_qc_table)) {
        fname <- paste0(pfx, "sample_qc.csv")
        utils::write.csv(plots$sample_qc_table,
                         file.path(dir_qc, fname), row.names = FALSE)
        saved <- c(saved, fname)
      }

      if (!is.null(plots$warn_proteins)) {
        fname <- paste0(pfx, "warn_proteins.png")
        ggplot2::ggsave(file.path(dir_qc, fname),
                        plots$warn_proteins, width = 7, height = 4, dpi = 150)
        saved <- c(saved, fname)
      }

      if (verbose && length(saved) > 0L)
        cat("[ELEUTHIA]   [qc]  ", paste(saved, collapse = ", "), "\n")
    }
  }

  # ---------------------------------------------------------------------------
  # 2b. Batch-check PCA plots (user-supplied)
  # ---------------------------------------------------------------------------

  if (!is.null(batch_plots)) {
    if (!is.list(batch_plots) || is.null(names(batch_plots)))
      stop("batch_plots must be a named list of ggplot objects.")

    for (nm in names(batch_plots)) {
      out_path <- file.path(dir_qc, paste0("pca_batch_", nm, ".png"))
      ggplot2::ggsave(out_path, batch_plots[[nm]], width = 7, height = 5, dpi = 150)
    }

    if (verbose)
      cat("[ELEUTHIA]   [qc]  batch PCA:", length(batch_plots), "plot(s) —",
          paste0("pca_batch_", names(batch_plots), ".png", collapse = ", "), "\n")
  }

  # ---------------------------------------------------------------------------
  # 3. Summary text file
  # ---------------------------------------------------------------------------

  smeta <- olink_data$sample_meta
  ameta <- olink_data$assay_meta
  pmat  <- assoc$pvalues

  bar <- strrep("=", 72)
  sep <- strrep("-", 72)

  lines <- c(
    bar,
    "OLINK QC SUMMARY",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    bar,
    ""
  )

  # --- Filter chain (populated when qc_steps is provided) ---
  if (!is.null(qc_steps) && is.list(qc_steps) && length(qc_steps) > 0L) {

    raw_entry   <- qc_steps[["raw"]]
    filter_keys <- setdiff(names(qc_steps), "raw")

    if (!is.null(raw_entry)) {
      lines <- c(lines,
                 "RAW DATA (post-load)",
                 paste0("  Proteins : ", raw_entry$n_proteins),
                 paste0("  Samples  : ", raw_entry$n_samples,
                        "  (includes controls and QC failures)"),
                 "", sep, "FILTER CHAIN", sep, "")
    }

    fmt_dim <- function(before, after) {
      if (is.na(before)) return(as.character(after))
      d <- after - before
      if (d == 0) paste0(before, " -> ", after, "  (unchanged)")
      else        paste0(before, " -> ", after, "  (removed ", abs(d), ")")
    }

    prev_s <- if (!is.null(raw_entry)) raw_entry$n_samples  else NA_integer_
    prev_p <- if (!is.null(raw_entry)) raw_entry$n_proteins else NA_integer_

    for (i in seq_along(filter_keys)) {
      e   <- qc_steps[[filter_keys[[i]]]]
      thr <- if (!is.na(e$threshold) && nzchar(e$threshold)) e$threshold else "-"

      lines <- c(lines,
                 paste0("[", i, "] ", e$step),
                 paste0("    ", e$description),
                 paste0("    Threshold  : ", thr),
                 paste0("    Samples    : ", fmt_dim(prev_s, e$n_samples)),
                 paste0("    Proteins   : ", fmt_dim(prev_p, e$n_proteins)),
                 "")

      prev_s <- e$n_samples
      prev_p <- e$n_proteins
    }

    raw_s <- if (!is.null(raw_entry)) raw_entry$n_samples  else NA_integer_
    raw_p <- if (!is.null(raw_entry)) raw_entry$n_proteins else NA_integer_
    fin_s <- ncol(olink_data$wide)
    fin_p <- nrow(olink_data$wide)

    pct <- function(raw, fin) {
      if (is.na(raw) || raw == 0L) return("")
      paste0("  (removed ", raw - fin, " of ", raw, ", ",
             round(100 * (raw - fin) / raw, 1), "%)")
    }

    lines <- c(lines,
               sep,
               "FINAL DATA",
               paste0("  Proteins : ", fin_p, pct(raw_p, fin_p)),
               paste0("  Samples  : ", fin_s, pct(raw_s, fin_s)),
               "")

  } else {
    lines <- c(lines,
               "FINAL DATA",
               paste0("  Proteins : ", nrow(olink_data$wide)),
               paste0("  Samples  : ", ncol(olink_data$wide)),
               "")
  }

  # --- AssayQC warn (final set) ---
  if (!is.null(ameta) && "warn_fraction" %in% colnames(ameta)) {
    lines <- c(lines,
               sep,
               "ASSAYQC WARN (final set)",
               paste0("  Proteins with any WARN    : ",
                      sum(ameta$warn_fraction > 0,    na.rm = TRUE)),
               paste0("  Proteins with WARN >= 10% : ",
                      sum(ameta$warn_fraction >= 0.1, na.rm = TRUE)),
               "")
  }

  # --- SampleQC table (final set) ---
  if ("SampleQC" %in% colnames(smeta)) {
    qc_tbl <- table(smeta$SampleQC)
    lines  <- c(lines,
                sep,
                "SAMPLEQC (final set)",
                paste0("  ", names(qc_tbl), " : ", as.integer(qc_tbl)),
                "")
  }

  # --- PC-metadata associations ---
  n_sig_per_var <- colSums(pmat < p_threshold, na.rm = TRUE)
  sig_vars      <- names(n_sig_per_var)[n_sig_per_var > 0L]
  lines <- c(lines,
             sep,
             paste0("PC-METADATA ASSOCIATION (p < ", p_threshold, ")"),
             paste0("  Variables significant in >=1 PC : ", length(sig_vars)),
             if (length(sig_vars) > 0L)
               paste0("  Variables : ", paste(sig_vars, collapse = ", ")),
             "")

  # --- Export parameters ---
  lines <- c(lines,
             sep,
             "EXPORT PARAMETERS",
             paste0("  p_threshold   : ", p_threshold),
             paste0("  top_n_vars    : ", top_n_vars),
             paste0("  heatmap_top_n : ",
                    if (is.null(heatmap_top_n)) "NULL" else heatmap_top_n),
             paste0("  ntop_pca      : ",
                    if (is.null(ntop_pca)) "NULL" else ntop_pca),
             bar)

  writeLines(lines, file.path(dir_qc, "summary.txt"))
  if (verbose) cat("[ELEUTHIA]   [qc]  summary.txt\n")

  # ---------------------------------------------------------------------------
  # 4. Association heatmap
  # ---------------------------------------------------------------------------

  p_assoc <- AETHER_plot_pc_association(
    assoc_result = assoc,
    p_threshold  = p_threshold,
    top_n        = heatmap_top_n,
    order_vars   = TRUE
  )

  ggplot2::ggsave(file.path(dir_assoc, "pc_metadata_heatmap.png"),
                  p_assoc, width = 10, height = 5, dpi = 150)

  if (verbose) cat("[ELEUTHIA]   [association] pc_metadata_heatmap.png\n")

  # ---------------------------------------------------------------------------
  # 5. PCA plots — best 2 PCs per variable
  # ---------------------------------------------------------------------------

  selected_vars <- unique(unlist(lapply(seq_len(nrow(pmat)), function(i) {
    p_row   <- pmat[i, ]
    valid   <- which(!is.na(p_row))
    if (length(valid) == 0L) return(character(0L))
    top_idx <- valid[order(p_row[valid])][seq_len(min(top_n_vars, length(valid)))]
    colnames(pmat)[top_idx]
  })))

  for (var in selected_vars) {
    vals     <- olink_data$sample_meta[[var]]
    n_unique <- length(unique(stats::na.omit(vals)))
    var_type <- if (is.numeric(vals) && n_unique >= 10) "continuous" else "categorical"

    var_pvals <- pmat[, var]
    valid_pcs <- which(!is.na(var_pvals))
    top2_pcs  <- valid_pcs[order(var_pvals[valid_pcs])][seq_len(min(2L, length(valid_pcs)))]
    best_dims <- as.integer(sub("^PC", "", rownames(pmat)[top2_pcs]))
    if (length(best_dims) < 2L) best_dims <- c(best_dims, best_dims[1L] + 1L)

    p <- AETHER_plot_pca(
      counts        = olink_data$wide,
      sample_info   = olink_data$sample_meta,
      sample_col    = sid_col,
      group_col     = var,
      var_type      = var_type,
      dims          = best_dims,
      log_transform = FALSE,
      ntop          = ntop_pca,
      label_samples = FALSE,
      title         = paste0("PCA \u2014 ", var,
                             " (PC", best_dims[1L], " vs PC", best_dims[2L], ")"),
      verbose       = FALSE
    )

    # Dynamic width for categorical variables: the right-side legend can overflow
    # a fixed width when labels are long or many categories are present.
    # Width = 5in (plot area) + estimated legend width, floored at 7 and capped at 14.
    if (var_type == "categorical") {
      cats      <- unique(stats::na.omit(as.character(vals)))
      max_chars <- max(nchar(cats), 0L)
      n_cats    <- length(cats)
      # 0.07in per character + 0.02in per category (title/key overhead) + 0.5in fixed
      legend_w  <- max_chars * 0.07 + n_cats * 0.02 + 0.5
      png_width <- min(max(7, 5 + legend_w), 14)
    } else {
      png_width <- 7
    }

    ggplot2::ggsave(file.path(dir_pca, paste0(var, ".png")),
                    p, width = png_width, height = 5, dpi = 150)
  }

  if (verbose)
    cat("[ELEUTHIA]   [pca] ", length(selected_vars), "plots\n")

  invisible(NULL)
}


#' @method print olink_data
#' @export
print.olink_data <- function(x, ...) {
  cat("Olink Data Object\n")
  cat("-----------------\n")
  cat("Proteins : ", nrow(x$wide), "\n", sep = "")
  cat("Samples  : ", ncol(x$wide), "\n", sep = "")

  if ("PlateID" %in% colnames(x$sample_meta)) {
    n_plates <- length(unique(x$sample_meta$PlateID))
    cat("Plates   : ", n_plates, "\n", sep = "")
  }

  if ("SampleType" %in% colnames(x$sample_meta)) {
    st_tbl <- table(x$sample_meta$SampleType)
    cat("SampleTypes: ", paste(names(st_tbl), st_tbl,
                                sep = "=", collapse = ", "), "\n", sep = "")
  }

  if ("SampleQC" %in% colnames(x$sample_meta)) {
    qc_tbl <- table(x$sample_meta$SampleQC)
    cat("SampleQC : ", paste(names(qc_tbl), qc_tbl,
                              sep = "=", collapse = ", "), "\n", sep = "")
  }

  n_warn <- sum(x$assay_meta$warn_fraction > 0, na.rm = TRUE)
  cat("WARN proteins (>=1 sample): ", n_warn, "\n", sep = "")

  na_pct <- round(100 * mean(is.na(x$wide)), 2)
  cat("NPX NAs  : ", na_pct, "%\n", sep = "")

  npx_r <- range(x$data$NPX, na.rm = TRUE)
  cat("NPX range: [", round(npx_r[1], 2), ", ", round(npx_r[2], 2), "]\n",
      sep = "")

  cat("\nSlots: $data (long), $wide (matrix), $sample_meta, $assay_meta, $params\n")
  invisible(x)
}
