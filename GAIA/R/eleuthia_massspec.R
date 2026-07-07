#' Load DIA Mass Spectrometry Data (DIA-NN report.pg_matrix)
#'
#' @description Reads a DIA-NN \code{report.pg_matrix} file (.xlsx or tab-delimited
#' text), parses sample identifiers from the column names (which are Windows or
#' POSIX file paths in DIA-NN output), log2-transforms intensities, and returns
#' a structured \code{massspec_data} object.
#'
#' This function performs no QC filtering. To remove controls, pool samples, and
#' high-missingness samples/proteins, pass the result to
#' \code{HADES_filter_massspec()} and \code{HADES_filter_massspec_proteins()}.
#'
#' @param ms_file Character. Path to the DIA-NN \code{report.pg_matrix} file.
#'   Accepted formats: \code{.xlsx}, \code{.tsv}, \code{.txt}, \code{.csv}.
#' @param sheet Integer or character. Sheet to read if \code{ms_file} is an
#'   xlsx workbook. Default: 1.
#' @param metadata_file Character or NULL. Path to a sample metadata file
#'   (.csv, .tsv, or .xlsx). Joined onto \code{$sample_meta} by
#'   \code{sample_col}. Default: NULL.
#' @param metadata_sheet Integer or character. Sheet for xlsx metadata files.
#'   Default: 1.
#' @param sample_col Character. Column in \code{metadata_file} used as the
#'   join key. Parsed sample IDs (numeric strings, SEP###-# labels, etc.) are
#'   matched against this column. Default: \code{"SampleID"}.
#' @param n_annotation_cols Integer or NULL. Number of leading columns that are
#'   protein annotations rather than sample intensities. \code{NULL} (default)
#'   uses auto-detection: annotation columns are those whose names do not contain
#'   a file-path separator (\code{/} or \code{\\}) or a \code{.d} extension.
#' @param verbose Logical. Print loading progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"massspec_data"} containing:
#' \describe{
#'   \item{wide}{Numeric matrix: proteins (rows) x samples (cols).
#'     Rownames = Protein.Group. Colnames = parsed SampleID.
#'     Values are log2-transformed intensities; zeros and negatives become NA.}
#'   \item{sample_meta}{Per-sample data.frame: SampleID, SampleType, PlateID,
#'     plus any columns joined from \code{metadata_file}.}
#'   \item{protein_meta}{Per-protein data.frame with annotation columns from
#'     the input file (e.g. Protein.Group, Protein.Names, Genes,
#'     First.Protein.Description). Rownames match \code{rownames($wide)}.}
#'   \item{params}{List of loader parameters used (for provenance).}
#' }
#'
#' @details
#' \strong{Expected input format:} Column names must follow the standardised
#' convention produced by the project's \code{prep_MS_data.R} script:
#' \code{{SampleID}_Plate {N}} (e.g. \code{"1027_Plate 1"},
#' \code{"SEP001-3_Plate 7"}, \code{"pool1_Plate 12"}). Raw DIA-NN file paths
#' are not supported; data must be prepared through \code{prep_MS_data.R} first.
#'
#' \strong{Sample ID parsing:} \code{PlateID} is extracted from the
#' \code{_Plate N} suffix; \code{SampleID} is everything before it. Column
#' names that do not match the \code{_Plate N} suffix pattern yield
#' \code{PlateID = NA} and \code{SampleID = original column name}.
#'
#' \strong{Sample type classification:}
#' \itemize{
#'   \item \code{"SAMPLE"} — any biological sample with a recognised SampleID
#'     pattern: pure numeric IDs (e.g. \code{1027}), \code{SEP###-#}
#'     longitudinal samples (e.g. \code{SEP001-3}), and \code{SEP_CTR_*}
#'     control samples. All three map to \code{"SAMPLE"} because their
#'     biological distinction (group, timepoint, cohort) is captured in the
#'     metadata via the \code{Group} column — SampleType does not need to
#'     re-encode it.
#'   \item \code{"POOL"} — pooled QC injections (\code{pool*}). Present on
#'     every plate for batch drift monitoring. No metadata entry.
#'   \item \code{"OTHER"} — any column whose parsed SampleID does not match
#'     any of the above patterns. Most commonly caused by samples that were
#'     not renamed via \code{prep_MS_data.R}, ad-hoc sub-cohort naming, or
#'     samples intentionally absent from the metadata file. These columns
#'     cannot be joined to metadata and are dropped by
#'     \code{HADES_filter_massspec()} by default.
#' }
#'
#' \strong{Log2 transformation:} Raw DIA intensities are on a linear scale.
#' Intensities \eqn{\le 0} (zeros = not detected; negatives = artefact) are
#' set to NA before \code{log2()} is applied. Do NOT transform again before
#' passing to \code{ARTEMIS_limma_de()}.
#'
#' \strong{Normalization:} Not performed at load time. Pass \code{$wide} to
#' \code{POSEIDON_normalize_massspec()} after filtering.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load only
#' ms_raw <- ELEUTHIA_load_massspec("report.pg_matrix.xlsx")
#'
#' # Load with Olink-matched metadata, then filter
#' ms_raw <- ELEUTHIA_load_massspec(
#'   ms_file       = "report.pg_matrix.xlsx",
#'   metadata_file = "olink_sample_meta.csv",
#'   sample_col    = "SampleID"
#' )
#' ms <- HADES_filter_massspec(ms_raw, keep_types = "SAMPLE")
#' ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
#' ms$wide <- POSEIDON_normalize_massspec(ms$wide)
#' }
ELEUTHIA_load_massspec <- function(ms_file,
                                    sheet             = 1,
                                    metadata_file     = NULL,
                                    metadata_sheet    = 1,
                                    sample_col        = "SampleID",
                                    n_annotation_cols = NULL,
                                    verbose           = TRUE) {

  if (!file.exists(ms_file))
    stop("ms_file not found: ", ms_file)

  ext <- tolower(sub(".*\\.", "", basename(ms_file)))

  # ---------------------------------------------------------------------------
  # Read intensity file
  # ---------------------------------------------------------------------------

  if (verbose) cat("[ELEUTHIA] Reading mass spec file:", ms_file, "\n")

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Package 'readxl' is required to read .xlsx files.\n",
           "Install with: install.packages('readxl')")
    raw <- as.data.frame(readxl::read_excel(ms_file, sheet = sheet))
  } else if (ext %in% c("tsv", "txt")) {
    raw <- utils::read.delim(ms_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    raw <- utils::read.csv(ms_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported format: '.", ext, "'. Use .xlsx, .tsv, .txt, or .csv.")
  }

  if (verbose)
    cat("[ELEUTHIA]   Rows:", nrow(raw), " | Columns:", ncol(raw), "\n")

  # ---------------------------------------------------------------------------
  # Identify annotation columns vs sample columns
  # ---------------------------------------------------------------------------

  if (!is.null(n_annotation_cols)) {
    annot_idx  <- seq_len(n_annotation_cols)
    sample_idx <- seq(n_annotation_cols + 1L, ncol(raw))
  } else {
    # Auto-detect: sample columns have file-path separators or .d extension
    is_sample <- grepl("[/\\\\]|\\.d$", colnames(raw), ignore.case = TRUE)
    # Fallback: if no path-based columns detected, assume first 4 are annotation
    if (!any(is_sample)) {
      if (verbose)
        cat("[ELEUTHIA]   No file-path columns detected; assuming first 4 are annotations\n")
      is_sample <- c(rep(FALSE, min(4L, ncol(raw))),
                     rep(TRUE,  max(0L, ncol(raw) - 4L)))
    }
    annot_idx  <- which(!is_sample)
    sample_idx <- which(is_sample)
  }

  if (length(sample_idx) == 0L)
    stop("No sample columns detected. Check n_annotation_cols or file format.")

  if (verbose) {
    cat("[ELEUTHIA]   Annotation columns:", length(annot_idx),
        " | Sample columns:", length(sample_idx), "\n")
  }

  # ---------------------------------------------------------------------------
  # Parse sample IDs from column names (file paths → clean ID)
  # ---------------------------------------------------------------------------

  raw_col_names <- colnames(raw)[sample_idx]

  # Parse SampleID and PlateID from standardised column names.
  # Expected format: {SampleID}_Plate {N}
  # e.g. "1027_Plate 1", "SEP001-3_Plate 7", "pool1_Plate 12"
  has_plate  <- grepl("_Plate [0-9]+$", raw_col_names)
  plate_ids  <- ifelse(has_plate,
                       sub(".*_(Plate [0-9]+)$", "\\1", raw_col_names),
                       NA_character_)
  parsed_ids <- ifelse(has_plate,
                       sub("_Plate [0-9]+$", "", raw_col_names),
                       raw_col_names)

  if (verbose && any(!has_plate))
    cat("[ELEUTHIA]   Warning:", sum(!has_plate),
        "column(s) do not match the '{SampleID}_Plate N' convention — PlateID set to NA\n")

  if (anyDuplicated(parsed_ids)) {
    n_duped    <- sum(duplicated(parsed_ids))
    parsed_ids <- make.unique(parsed_ids, sep = "_")
    if (verbose)
      cat("[ELEUTHIA]   Deduplicated", n_duped,
          "repeated sample IDs with make.unique()\n")
  }

  # ---------------------------------------------------------------------------
  # Classify sample types
  # ---------------------------------------------------------------------------

  is_pool   <- grepl("^pool", parsed_ids, ignore.case = TRUE)
  is_sample <- grepl("^[0-9]+$", parsed_ids) |
               grepl("^SEP[0-9]+-", parsed_ids) |
               grepl("^SEP_CTR", parsed_ids)

  sample_type <- ifelse(is_pool,   "POOL",
                 ifelse(is_sample, "SAMPLE",
                                   "OTHER"))

  type_tbl <- table(sample_type)
  if (verbose) {
    cat("[ELEUTHIA]   Sample types: ",
        paste(names(type_tbl), type_tbl, sep = "=", collapse = ", "), "\n")
  }

  # ---------------------------------------------------------------------------
  # Drop rows with no primary protein identifier
  # DIA-NN appends summary-statistic rows (e.g. "Sum intensities",
  # "Number of valid values") with empty Protein.Group fields. These are
  # per-run QC statistics, not protein measurements, and must be excluded
  # before matrix construction. Gene ID is intentionally not used here:
  # valid proteins may lack a gene symbol (e.g. immunoglobulin constant
  # regions), but every real protein row has a Protein.Group accession.
  # ---------------------------------------------------------------------------

  protein_id_col <- colnames(raw)[annot_idx[1L]]
  valid_prot     <- !is.na(raw[[protein_id_col]]) &
                    nzchar(trimws(as.character(raw[[protein_id_col]])))
  if (any(!valid_prot)) {
    if (verbose)
      cat("[ELEUTHIA]   Dropped", sum(!valid_prot),
          "row(s) with no Protein.Group identifier (DIA-NN summary rows)\n")
    raw <- raw[valid_prot, , drop = FALSE]
  }

  # ---------------------------------------------------------------------------
  # Build wide matrix (proteins x samples), log2 transform
  # ---------------------------------------------------------------------------

  intensity_df <- raw[, sample_idx, drop = FALSE]
  wide_mat     <- suppressWarnings(
    matrix(as.numeric(as.matrix(intensity_df)),
           nrow     = nrow(raw),
           ncol     = length(sample_idx),
           dimnames = list(NULL, parsed_ids))
  )

  # Log2 transform: zeros and negatives → NA first
  wide_mat[!is.na(wide_mat) & wide_mat <= 0] <- NA
  wide_mat <- log2(wide_mat)

  # Rownames: use first annotation column as protein ID (Protein.Group)
  prot_ids <- as.character(raw[[protein_id_col]])
  if (anyDuplicated(prot_ids))
    prot_ids <- make.unique(prot_ids, sep = "_dup")
  rownames(wide_mat) <- prot_ids

  na_pct <- round(100 * mean(is.na(wide_mat)), 1)
  if (verbose)
    cat("[ELEUTHIA]   Matrix: ", nrow(wide_mat), " proteins x ",
        ncol(wide_mat), " samples | NA: ", na_pct, "%\n", sep = "")

  # ---------------------------------------------------------------------------
  # Build protein_meta
  # ---------------------------------------------------------------------------

  protein_meta              <- raw[, annot_idx, drop = FALSE]
  rownames(protein_meta)    <- prot_ids

  # ---------------------------------------------------------------------------
  # Build sample_meta
  # ---------------------------------------------------------------------------

  sample_meta <- data.frame(
    SampleID   = parsed_ids,
    SampleType = sample_type,
    PlateID    = plate_ids,
    stringsAsFactors = FALSE
  )
  colnames(sample_meta)[1L] <- sample_col

  # ---------------------------------------------------------------------------
  # Read and join metadata file
  # ---------------------------------------------------------------------------

  if (!is.null(metadata_file)) {
    if (!file.exists(metadata_file))
      stop("metadata_file not found: ", metadata_file)

    if (verbose) cat("[ELEUTHIA] Reading metadata:", metadata_file, "\n")

    meta_ext <- tolower(sub(".*\\.", "", basename(metadata_file)))

    if (meta_ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("Package 'readxl' is required for .xlsx metadata.")
      meta <- as.data.frame(readxl::read_excel(metadata_file,
                                                sheet = metadata_sheet))
    } else if (meta_ext == "csv") {
      meta <- utils::read.csv(metadata_file, stringsAsFactors = FALSE,
                               check.names = FALSE)
    } else {
      meta <- utils::read.delim(metadata_file, stringsAsFactors = FALSE,
                                 check.names = FALSE)
    }

    if (!sample_col %in% colnames(meta))
      stop("sample_col '", sample_col, "' not found in metadata_file.")

    meta[[sample_col]] <- as.character(meta[[sample_col]])
    meta <- meta[!is.na(meta[[sample_col]]) &
                   nzchar(trimws(meta[[sample_col]])), ]

    new_cols    <- setdiff(colnames(meta), colnames(sample_meta))
    meta_join   <- meta[, c(sample_col, new_cols), drop = FALSE]
    sample_meta <- merge(sample_meta, meta_join, by = sample_col,
                         all.x = TRUE, sort = FALSE)

    n_matched <- sum(sample_meta[[sample_col]] %in% meta[[sample_col]])
    if (verbose)
      cat("[ELEUTHIA]   Joined", length(new_cols), "metadata columns;",
          n_matched, "of", nrow(sample_meta), "samples matched\n")
  }

  rownames(sample_meta) <- sample_meta[[sample_col]]

  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------

  params <- list(
    ms_file       = ms_file,
    metadata_file = metadata_file,
    sample_col    = sample_col,
    sheet         = sheet
  )

  result <- list(
    wide         = wide_mat,
    sample_meta  = sample_meta,
    protein_meta = protein_meta,
    params       = params
  )
  class(result) <- c("massspec_data", "list")

  if (verbose)
    cat("[ELEUTHIA] massspec_data object ready\n")

  return(result)
}


#' @method print massspec_data
#' @export
print.massspec_data <- function(x, ...) {
  cat("Mass Spectrometry Data Object (DIA)\n")
  cat("-------------------------------------\n")
  cat("Proteins : ", nrow(x$wide), "\n", sep = "")
  cat("Samples  : ", ncol(x$wide), "\n", sep = "")

  if ("SampleType" %in% colnames(x$sample_meta)) {
    st_tbl <- table(x$sample_meta$SampleType)
    cat("SampleTypes: ",
        paste(names(st_tbl), st_tbl, sep = "=", collapse = ", "), "\n", sep = "")
  }

  na_pct <- round(100 * mean(is.na(x$wide)), 1)
  cat("NA (log2): ", na_pct, "%\n", sep = "")

  valid_prot <- sum(rowMeans(!is.na(x$wide)) > 0.5)
  cat("Proteins >50% valid: ", valid_prot, "\n", sep = "")

  rng <- range(x$wide, na.rm = TRUE)
  cat("log2 range: [", round(rng[1], 1), ", ", round(rng[2], 1), "]\n", sep = "")

  cat("\nSlots: $wide (matrix), $sample_meta, $protein_meta, $params\n")
  invisible(x)
}
