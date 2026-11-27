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
#'   \item{qc_log}{The QC record. Starts with one entry for the load itself
#'     (arguments, instrument rows and blank-ID rows dropped, samples removed by
#'     \code{restrict_to_metadata}); each \code{HADES_*_olink*} step appends one
#'     more. Read it with \code{HADES_get_olink_qc_log()}.}
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
  n_instrument_rows <- n_before - nrow(raw)
  if (verbose && n_instrument_rows > 0) {
    cat("[ELEUTHIA]   Parsed", n_instrument_rows,
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
  n_samples_read  <- length(unique(as.character(raw[[sample_col]])))
  not_in_meta_ids <- character(0)

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
      not_in_meta_ids <- as.character(not_in_meta)
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

  # First entry of the QC log: what the loader itself dropped. Later HADES_*_olink
  # steps append to this log (see HADES_get_olink_qc_log()).
  load_entry <- .olink_log_entry(
    step            = "ELEUTHIA_load_olink",
    action          = "loaded",
    params          = list(npx_file = basename(npx_file),
                           metadata_file = if (is.null(metadata_file)) NULL else basename(metadata_file),
                           sample_col = sample_col, npx_col = npx_col,
                           restrict_to_metadata = isTRUE(restrict_to_metadata)),
    samples_before  = n_samples_read, samples_after  = ncol(wide_mat),
    proteins_before = nrow(wide_mat), proteins_after = nrow(wide_mat),
    removed_samples = if (length(not_in_meta_ids) > 0L) {
      data.frame(sample_id = not_in_meta_ids, reason = "not in metadata file",
                 value = NA_real_, stringsAsFactors = FALSE)
    } else NULL,
    details         = list(instrument_rows_dropped = n_instrument_rows,
                           blank_sample_id_rows_dropped = n_invalid)
  )

  result <- list(
    data        = data_df,
    wide        = wide_mat,
    sample_meta = sample_meta,
    assay_meta  = assay_meta,
    params      = params,
    qc_log      = list(load_entry)
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
#' @description Writes the outputs of the Olink QC / exploration pipeline to a
#' structured directory: the filtered data (RDS and CSVs), the QC log of every
#' step that was run (with the IDs removed and why), QC plots, a
#' PC-metadata association heatmap, per-variable PCA plots (coloured by the two
#' PCs most significantly associated with each variable), and a plain-text QC
#' summary built from the log.
#'
#' The record of what was filtered comes from \code{olink_data$qc_log}, which
#' \code{ELEUTHIA_load_olink()} and the \code{HADES_*_olink*} functions maintain
#' as they run; nothing has to be assembled by hand.
#'
#' @param output_dir Character. Root directory for all outputs. Created if it
#'   does not exist.
#' @param olink_data An \code{olink_data} object - the filtered, QC-cleaned
#'   object to export.
#' @param assoc An \code{artemis_pc_assoc} object from
#'   \code{ARTEMIS_pc_metadata_association()}, or \code{NULL} to skip the
#'   association heatmap, per-variable PCA plots and the association section of
#'   the summary.
#' @param qc_plots A named list from \code{AETHER_plot_olink_qc()} (single pass),
#'   or a named list of such lists for several passes (e.g. \code{pre} and
#'   \code{post}), or \code{NULL}.
#' @param batch_plots Named list of ggplot objects for batch-check PCA plots.
#'   Each element is saved as \code{qc/<prefix>_pca_batch_<name>}. Default:
#'   \code{NULL}.
#' @param top_n_vars Integer. Number of top variables per PC (by p-value) to
#'   include in the PCA plot loop. Default: 5.
#' @param p_threshold Numeric. Significance threshold used in the association
#'   heatmap border annotation and in the summary file. Default: 0.05.
#' @param heatmap_top_n Integer or NULL. Passed to
#'   \code{AETHER_plot_pc_association(top_n)} to restrict the heatmap to the
#'   top N variables per PC. NULL shows all tested variables. Default: NULL.
#' @param ntop_pca Integer or NULL. Number of most variable proteins used for
#'   PCA. NULL uses all proteins. Default: NULL.
#' @param save_plots Logical. Save the plots. Default: TRUE.
#' @param plot_format Character. Plot file format (file extension), e.g.
#'   \code{"png"} or \code{"pdf"}. Default: \code{"png"}.
#' @param prefix Character. Prefix for the output file names. Default:
#'   \code{"olink"}.
#' @param save_rds Logical. Save the object as \code{data/<prefix>_filtered.rds}.
#'   Default: TRUE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisibly, a character vector of the files created.
#'
#' @details
#' \strong{Directory layout} (\code{prefix = "olink"}):
#' \preformatted{
#' output_dir/
#'   data/
#'     olink_filtered.rds          -- the olink_data object (incl. $qc_log)
#'     olink_sample_meta.csv       -- per-sample metadata (with outlier columns)
#'     olink_assay_meta.csv        -- per-protein metadata (warn_fraction etc.)
#'   qc/
#'     olink_qc_log.csv            -- one row per step, with its arguments
#'     olink_removed_samples.csv   -- sample, step, reason, value
#'     olink_removed_proteins.csv  -- OlinkID, Assay, step, reason, value
#'     olink_flagged_samples.csv   -- flagged but kept (outlier detection)
#'     olink_summary.txt           -- the QC record in readable form
#'     olink_<set>_npx_distributions.png, olink_<set>_sample_qc.csv,
#'     olink_<set>_warn_proteins.png   -- from qc_plots
#'     olink_pca_batch_<name>.png      -- from batch_plots
#'   association/
#'     olink_pc_metadata_heatmap.png
#'   pca/
#'     olink_pca_<variable>.png    -- one per selected variable
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
#' ol_raw <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.csv")
#' ol <- HADES_filter_olink(ol_raw)
#' ol <- HADES_filter_olink_proteins(ol)
#' ol <- HADES_detect_outliers_olink(ol, filter = TRUE)
#'
#' qc_plots <- AETHER_plot_olink_qc(ol, plot_warn_proteins = TRUE)
#' assoc    <- ARTEMIS_pc_metadata_association(ol$wide, ol$sample_meta,
#'                                              sample_col = "SampleID")
#'
#' ELEUTHIA_export_olink_qc("results/Olink_QC", ol, assoc = assoc,
#'                          qc_plots = qc_plots)
#' }
ELEUTHIA_export_olink_qc <- function(output_dir,
                                      olink_data,
                                      assoc         = NULL,
                                      qc_plots      = NULL,
                                      batch_plots   = NULL,
                                      top_n_vars    = 5L,
                                      p_threshold   = 0.05,
                                      heatmap_top_n = NULL,
                                      ntop_pca      = NULL,
                                      save_plots    = TRUE,
                                      plot_format   = "png",
                                      prefix        = "olink",
                                      save_rds      = TRUE,
                                      verbose       = TRUE) {

  if (!inherits(olink_data, "olink_data"))
    stop("olink_data must be an olink_data object from ELEUTHIA_load_olink().")
  if (!is.null(assoc) && !inherits(assoc, "artemis_pc_assoc"))
    stop("assoc must be an artemis_pc_assoc object from ARTEMIS_pc_metadata_association(), or NULL.")
  if (!is.null(qc_plots) && !is.list(qc_plots))
    stop("qc_plots must be a list from AETHER_plot_olink_qc(), or NULL.")
  if (!is.null(batch_plots) && (!is.list(batch_plots) || is.null(names(batch_plots))))
    stop("batch_plots must be a named list of ggplot objects.")

  sid_col    <- olink_data$params$sample_col
  top_n_vars <- as.integer(top_n_vars)
  files_created <- character(0)
  add_file <- function(f) files_created <<- c(files_created, f)
  pfile <- function(dir, name, ext) file.path(dir, paste0(prefix, "_", name, ".", ext))

  # ---------------------------------------------------------------------------
  # Directory structure
  # ---------------------------------------------------------------------------

  dir_data  <- file.path(output_dir, "data")
  dir_qc    <- file.path(output_dir, "qc")
  dir_assoc <- file.path(output_dir, "association")
  dir_pca   <- file.path(output_dir, "pca")

  dirs <- c(dir_data, dir_qc)
  if (!is.null(assoc)) dirs <- c(dirs, dir_assoc, dir_pca)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    cat("[ELEUTHIA] Exporting Olink QC results\n")
    cat("    Output directory:", output_dir, "\n\n")
  }

  # ---------------------------------------------------------------------------
  # 1. Data exports
  # ---------------------------------------------------------------------------

  if (save_rds) {
    f <- pfile(dir_data, "filtered", "rds")
    saveRDS(olink_data, f); add_file(f)
  }
  f <- pfile(dir_data, "sample_meta", "csv")
  utils::write.csv(olink_data$sample_meta, f, row.names = FALSE); add_file(f)
  f <- pfile(dir_data, "assay_meta", "csv")
  utils::write.csv(olink_data$assay_meta, f, row.names = FALSE); add_file(f)
  if (verbose) cat("    Data: filtered object, sample_meta, assay_meta\n")

  # ---------------------------------------------------------------------------
  # 2. QC log tables
  # ---------------------------------------------------------------------------

  has_log <- !is.null(olink_data$qc_log) && length(olink_data$qc_log) > 0L
  if (has_log) {
    for (what in c("steps", "removed_samples", "removed_proteins", "flagged_samples")) {
      nm <- if (what == "steps") "qc_log" else what
      f  <- pfile(dir_qc, nm, "csv")
      utils::write.csv(HADES_get_olink_qc_log(olink_data, what), f, row.names = FALSE)
      add_file(f)
    }
    if (verbose) cat("    QC log: qc_log, removed_samples, removed_proteins, flagged_samples\n")
  } else if (verbose) {
    cat("    QC log: none recorded in this object (no $qc_log)\n")
  }

  # ---------------------------------------------------------------------------
  # 3. QC plots (a single AETHER_plot_olink_qc() list, or a named list of them)
  # ---------------------------------------------------------------------------

  if (!is.null(qc_plots)) {
    qc_sets <- if (!is.null(qc_plots$npx_distributions)) {
      stats::setNames(list(qc_plots), "")   # single set - no name in the file
    } else {
      qc_plots
    }

    for (set_nm in names(qc_sets)) {
      pname <- function(x) if (nzchar(set_nm)) paste0(set_nm, "_", x) else x
      plots <- qc_sets[[set_nm]]
      saved <- character(0L)

      if (save_plots && !is.null(plots$npx_distributions)) {
        f <- pfile(dir_qc, pname("npx_distributions"), plot_format)
        ggplot2::ggsave(f, plots$npx_distributions, width = 9, height = 5, dpi = 150)
        add_file(f); saved <- c(saved, basename(f))
      }
      if (!is.null(plots$sample_qc_table)) {
        f <- pfile(dir_qc, pname("sample_qc"), "csv")
        utils::write.csv(plots$sample_qc_table, f, row.names = FALSE)
        add_file(f); saved <- c(saved, basename(f))
      }
      if (save_plots && !is.null(plots$warn_proteins)) {
        f <- pfile(dir_qc, pname("warn_proteins"), plot_format)
        ggplot2::ggsave(f, plots$warn_proteins, width = 7, height = 4, dpi = 150)
        add_file(f); saved <- c(saved, basename(f))
      }
      if (verbose && length(saved) > 0L)
        cat("    QC plots:", paste(saved, collapse = ", "), "\n")
    }
  }

  # Batch-check PCA plots (user-supplied)
  if (save_plots && !is.null(batch_plots)) {
    for (nm in names(batch_plots)) {
      f <- pfile(dir_qc, paste0("pca_batch_", nm), plot_format)
      ggplot2::ggsave(f, batch_plots[[nm]], width = 7, height = 5, dpi = 150)
      add_file(f)
    }
    if (verbose)
      cat("    Batch PCA:", length(batch_plots), "plot(s)\n")
  }

  # ---------------------------------------------------------------------------
  # 4. Summary text, built from the QC log
  # ---------------------------------------------------------------------------

  lines <- .olink_qc_summary_lines(
    olink_data, assoc, p_threshold = p_threshold, top_n_vars = top_n_vars,
    heatmap_top_n = heatmap_top_n, ntop_pca = ntop_pca, prefix = prefix)
  f <- pfile(dir_qc, "summary", "txt")
  writeLines(lines, f); add_file(f)
  if (verbose) cat("    Summary:", basename(f), "\n")

  # ---------------------------------------------------------------------------
  # 5. Association heatmap and per-variable PCA plots
  # ---------------------------------------------------------------------------

  if (!is.null(assoc) && save_plots) {
    pmat <- assoc$pvalues

    p_assoc <- AETHER_plot_pc_association(
      assoc_result = assoc,
      p_threshold  = p_threshold,
      top_n        = heatmap_top_n,
      order_vars   = TRUE
    )
    f <- pfile(dir_assoc, "pc_metadata_heatmap", plot_format)
    ggplot2::ggsave(f, p_assoc, width = 10, height = 5, dpi = 150)
    add_file(f)
    if (verbose) cat("    Association: pc_metadata_heatmap\n")

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
        title         = paste0("PCA - ", var,
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

      f <- pfile(dir_pca, paste0("pca_", var), plot_format)
      ggplot2::ggsave(f, p, width = png_width, height = 5, dpi = 150)
      add_file(f)
    }

    if (verbose)
      cat("    PCA plots:", length(selected_vars), "\n")
  }

  if (verbose)
    cat("\n[ELEUTHIA] Olink QC export complete:", length(files_created), "files\n")

  invisible(files_created)
}


# Internal: the text of <prefix>_summary.txt, built from olink_data$qc_log.
#' @noRd
.olink_qc_summary_lines <- function(olink_data, assoc, p_threshold, top_n_vars,
                                    heatmap_top_n, ntop_pca, prefix) {

  smeta <- olink_data$sample_meta
  ameta <- olink_data$assay_meta
  log   <- olink_data$qc_log

  bar <- strrep("=", 72)
  sep <- strrep("-", 72)

  # "n -> m  (removed k)", used for the before/after of each step
  fmt_dim <- function(before, after) {
    if (is.na(before)) return(as.character(after))
    d <- after - before
    if (d == 0) paste0(before, " -> ", after, "  (unchanged)")
    else        paste0(before, " -> ", after, "  (removed ", abs(d), ")")
  }

  # up to 20 IDs inline, the rest summarised (the CSVs hold the full lists)
  fmt_ids <- function(ids) {
    ids <- as.character(ids)
    if (length(ids) <= 20L) return(paste(ids, collapse = ", "))
    paste0(paste(ids[1:20], collapse = ", "), ", ... (+", length(ids) - 20L, " more)")
  }

  lines <- c(bar, "OLINK QC SUMMARY",
             paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             bar, "")

  if (!is.null(log) && length(log) > 0L) {

    first <- log[[1L]]
    steps <- if (first$action == "loaded") log[-1L] else log

    if (first$action == "loaded") {
      lines <- c(lines,
                 "LOAD",
                 paste0("  Parameters : ", .olink_fmt_params(first$params)),
                 paste0("  Proteins   : ", first$proteins_after),
                 paste0("  Samples    : ", first$samples_after,
                        "  (includes controls and QC failures)"))
      d <- first$details
      if (!is.null(d$instrument_rows_dropped))
        lines <- c(lines, paste0("  Instrument control rows dropped : ", d$instrument_rows_dropped))
      if (!is.null(d$blank_sample_id_rows_dropped))
        lines <- c(lines, paste0("  Rows with blank sample ID dropped : ", d$blank_sample_id_rows_dropped))
      if (nrow(first$removed_samples) > 0L)
        lines <- c(lines, paste0("  Samples not in metadata, removed (", nrow(first$removed_samples),
                                 ") : ", fmt_ids(first$removed_samples$sample_id)))
      lines <- c(lines, "")
    }

    lines <- c(lines, sep, "FILTER CHAIN", sep, "")
    if (length(steps) == 0L) lines <- c(lines, "  (no filter steps recorded)", "")

    for (i in seq_along(steps)) {
      e <- steps[[i]]
      lines <- c(lines,
                 paste0("[", i, "] ", e$step, "  (", e$action, ")"),
                 paste0("    Parameters : ", .olink_fmt_params(e$params)))
      if (!is.null(e$details$cutoff))
        lines <- c(lines, paste0("    Cutoff     : ", round(e$details$cutoff, 3),
                                 "  (distance in PC space)"))
      lines <- c(lines,
                 paste0("    Samples    : ", fmt_dim(e$samples_before,  e$samples_after)),
                 paste0("    Proteins   : ", fmt_dim(e$proteins_before, e$proteins_after)))
      if (nrow(e$removed_samples) > 0L)
        lines <- c(lines, paste0("    Samples removed (", nrow(e$removed_samples), "): ",
                                 fmt_ids(e$removed_samples$sample_id)))
      if (nrow(e$flagged_samples) > 0L)
        lines <- c(lines, paste0("    Samples flagged, kept in data (", nrow(e$flagged_samples), "): ",
                                 fmt_ids(e$flagged_samples$sample_id)))
      if (nrow(e$removed_proteins) > 0L)
        lines <- c(lines, paste0("    Proteins removed (", nrow(e$removed_proteins), "): ",
                                 fmt_ids(ifelse(is.na(e$removed_proteins$Assay),
                                                e$removed_proteins$OlinkID,
                                                e$removed_proteins$Assay))))
      lines <- c(lines, "")
    }

    raw_s <- first$samples_after
    raw_p <- first$proteins_after
    fin_s <- ncol(olink_data$wide)
    fin_p <- nrow(olink_data$wide)
    pct <- function(raw, fin) {
      if (is.na(raw) || raw == 0L) return("")
      paste0("  (removed ", raw - fin, " of ", raw, ", ",
             round(100 * (raw - fin) / raw, 1), "%)")
    }
    lines <- c(lines, sep, "FINAL DATA",
               paste0("  Proteins : ", fin_p, pct(raw_p, fin_p)),
               paste0("  Samples  : ", fin_s, pct(raw_s, fin_s)),
               "",
               paste0("  Full lists of removed / flagged IDs: ", prefix,
                      "_removed_samples.csv, ", prefix, "_removed_proteins.csv, ",
                      prefix, "_flagged_samples.csv"),
               "")

  } else {
    lines <- c(lines,
               "No QC log is recorded in this object ($qc_log is empty).",
               "Create the object with ELEUTHIA_load_olink() and filter it with the",
               "HADES_*_olink* functions to record each step.",
               "", sep, "FINAL DATA",
               paste0("  Proteins : ", nrow(olink_data$wide)),
               paste0("  Samples  : ", ncol(olink_data$wide)),
               "")
  }

  # --- AssayQC warn (final set) ---
  if (!is.null(ameta) && "warn_fraction" %in% colnames(ameta)) {
    lines <- c(lines, sep, "ASSAYQC WARN (final set)",
               paste0("  Proteins with any WARN    : ",
                      sum(ameta$warn_fraction > 0,    na.rm = TRUE)),
               paste0("  Proteins with WARN >= 10% : ",
                      sum(ameta$warn_fraction >= 0.1, na.rm = TRUE)),
               "")
  }

  # --- SampleQC table (final set) ---
  if ("SampleQC" %in% colnames(smeta)) {
    qc_tbl <- table(smeta$SampleQC)
    lines  <- c(lines, sep, "SAMPLEQC (final set)",
                paste0("  ", names(qc_tbl), " : ", as.integer(qc_tbl)),
                "")
  }

  # --- PC-metadata associations ---
  if (!is.null(assoc)) {
    pmat          <- assoc$pvalues
    n_sig_per_var <- colSums(pmat < p_threshold, na.rm = TRUE)
    sig_vars      <- names(n_sig_per_var)[n_sig_per_var > 0L]
    lines <- c(lines, sep,
               paste0("PC-METADATA ASSOCIATION (p < ", p_threshold, ")"),
               paste0("  Variables significant in >=1 PC : ", length(sig_vars)),
               if (length(sig_vars) > 0L)
                 paste0("  Variables : ", paste(sig_vars, collapse = ", ")),
               "")
  }

  # --- Export parameters ---
  c(lines, sep, "EXPORT PARAMETERS",
    paste0("  p_threshold   : ", p_threshold),
    paste0("  top_n_vars    : ", top_n_vars),
    paste0("  heatmap_top_n : ", if (is.null(heatmap_top_n)) "NULL" else heatmap_top_n),
    paste0("  ntop_pca      : ", if (is.null(ntop_pca)) "NULL" else ntop_pca),
    bar)
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

  if (!is.null(x$qc_log))
    cat("QC log   : ", length(x$qc_log), " step(s) recorded\n", sep = "")

  cat("\nSlots: $data (long), $wide (matrix), $sample_meta, $assay_meta, $params, $qc_log\n")
  invisible(x)
}
