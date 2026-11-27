#' Load Mass Spectrometry Protein Quantities
#'
#' @description Reads a protein x sample quantification table and returns a
#' structured \code{massspec_data} object that the HADES, POSEIDON, ARTEMIS and
#' ELEUTHIA mass spectrometry functions work on. The loader knows nothing about
#' any particular acquisition or search tool: it expects a table that has been
#' prepared into a simple, explicit shape (see Details), and it stops with a
#' message naming the problem when the table does not fit that shape.
#'
#' The loader does only structural work: reading, choosing the sample columns,
#' naming samples, putting the values on the log2 scale, and joining sample
#' metadata. It does not filter, normalise, impute or map identifiers. Rows and
#' columns specific to one tool (for example MaxQuant's decoy and contaminant
#' flags, or a minimum peptide count) are removed by the user before loading;
#' because \code{ms_data} can be a data.frame, that takes a few lines of R and
#' no intermediate file. Worked examples for common tools are in the tutorial.
#'
#' @param ms_data Character path to a \code{.csv}, \code{.tsv}, \code{.txt} or
#'   \code{.xlsx} file, or a data.frame already read into R. One row per
#'   protein and one column per sample (proteins in rows, not columns).
#' @param scale Character, \strong{required}: \code{"log2"} if the quantities
#'   are already log2-transformed, or \code{"linear"} if they are raw
#'   intensities (zeros and negative values become \code{NA}, then
#'   \code{log2()} is applied). There is no default because a wrong choice
#'   silently corrupts every later step (a second log transform, or linear
#'   values given to imputation and limma as if they were log2).
#' @param id_col Name or position of the column holding the unique protein
#'   identifier. Default: the first column.
#' @param sample_regex Optional regular expression (Perl-compatible) that
#'   selects the sample columns by name. If it contains a capture group, the
#'   captured text becomes the sample ID; otherwise the whole column name is
#'   used. Every column that is not a sample column (other than the ID column)
#'   becomes protein annotation. If \code{NULL}, every column except the ID
#'   column and \code{annotation_cols} is a sample column. Examples:
#'   \code{"^LFQ intensity (.+)$"} (MaxQuant LFQ columns),
#'   \code{"^(.+)_Plate [0-9]+$"} (a sample ID followed by a plate suffix).
#' @param annotation_cols Optional character vector of columns to keep as
#'   protein annotation when \code{sample_regex} is \code{NULL}. Cannot be
#'   combined with \code{sample_regex}.
#' @param sheet Sheet to read when \code{ms_data} is an \code{.xlsx} file.
#'   Default: 1.
#' @param metadata Optional sample metadata: a path (\code{.csv}, \code{.tsv},
#'   \code{.txt}, \code{.xlsx}) or a data.frame. Joined onto \code{$sample_meta}
#'   by \code{sample_col}.
#' @param metadata_sheet Sheet to read when \code{metadata} is an \code{.xlsx}
#'   file. Default: 1.
#' @param sample_col Name of the sample ID column in \code{metadata}; it is also
#'   the name of the sample ID column in the returned \code{$sample_meta}.
#'   Default: \code{"SampleID"}.
#' @param pool_pattern Optional regular expression. Samples whose ID matches it
#'   are typed \code{"POOL"} (pooled QC injections, which usually have no
#'   metadata). Case-insensitive. Default: \code{NULL} (no pools).
#' @param verbose Logical. Print loading progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"massspec_data"} with:
#' \describe{
#'   \item{wide}{Numeric matrix, proteins (rows) x samples (columns), on the
#'     log2 scale. Rownames are the protein IDs, colnames the sample IDs.}
#'   \item{sample_meta}{Per-sample data.frame: the sample ID (column named
#'     \code{sample_col}), \code{SampleType}, \code{source_column} (the
#'     original column name in \code{ms_data}, kept so that information encoded
#'     in the name, such as a plate or run number, can be extracted
#'     afterwards), plus any columns joined from \code{metadata}.}
#'   \item{protein_meta}{Per-protein data.frame, the ID column first, then the
#'     annotation columns. Rownames match \code{rownames($wide)}.}
#'   \item{params}{List of the arguments used and counts of rows dropped (for
#'     provenance).}
#' }
#'
#' @details
#' \strong{What the table must look like.}
#' \itemize{
#'   \item Proteins in rows, samples in columns. Public matrices that have
#'     samples in rows must be transposed before loading.
#'   \item A unique protein identifier column. Rows with a blank identifier
#'     (for example summary rows some tools append) are dropped and counted.
#'     \strong{Duplicated identifiers are an error}: they usually mean the
#'     wrong column was chosen as the identifier (gene names, for example,
#'     repeat). Resolve them before loading.
#'   \item Sample columns that are numeric. \code{NA}, blank cells and the
#'     strings \code{"NaN"}, \code{"N/A"}, \code{"#N/A"} and \code{"Filtered"}
#'     are read as missing in the sample columns. Any other non-numeric text in
#'     a sample column is an error that names the column. Annotation and
#'     metadata text is left exactly as read (a blank cell stays blank).
#'   \item Unique sample IDs. Duplicates after \code{sample_regex} is applied
#'     are an error.
#' }
#'
#' \strong{Sample types.} In order: an ID matching \code{pool_pattern} is
#' \code{"POOL"}; otherwise, if \code{metadata} is supplied, a sample found in
#' it takes the value of its \code{SampleType} column if there is one
#' (missing values are an error) or \code{"SAMPLE"} if not; a sample absent from
#' the metadata is \code{"OTHER"}; with no metadata every other sample is
#' \code{"SAMPLE"}. \code{HADES_filter_massspec()} keeps SAMPLE and POOL by
#' default, so \code{"OTHER"} samples (no metadata) are dropped there.
#'
#' \strong{Batch information} (plate, run, date) is not parsed from the names.
#' Supply it through \code{metadata}, or derive it from
#' \code{$sample_meta$source_column} after loading.
#'
#' \strong{Scale check.} With \code{scale = "log2"} a warning is raised if the
#' median value is above 100, and with \code{scale = "linear"} if every value is
#' below 100, since either suggests the wrong scale was given. These
#' thresholds are heuristics with no published reference; verify the scale
#' yourself.
#'
#' \strong{The object.} \code{massspec_data} is a plain list of the four
#' components above with class \code{c("massspec_data", "list")}; the
#' downstream functions read \code{$wide}, \code{$sample_meta},
#' \code{$protein_meta} and \code{$params$sample_col}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # A prepared matrix (first column = protein ID, other columns = samples)
#' ms <- ELEUTHIA_load_massspec("proteins.csv", scale = "log2")
#'
#' # A MaxQuant proteinGroups.txt, filtered first, LFQ columns only
#' pg <- read.delim("proteinGroups.txt", check.names = FALSE)
#' pg <- pg[!(pg$`Potential contaminant` %in% "+") & !(pg$Decoy %in% "+"), ]
#' ms <- ELEUTHIA_load_massspec(pg, scale = "linear",
#'                              id_col = "Protein IDs",
#'                              sample_regex = "^LFQ intensity (.+)$",
#'                              metadata = "sample_meta.csv")
#'
#' # Sample columns named "<SampleID>_Plate <N>", with pooled QC samples
#' ms <- ELEUTHIA_load_massspec("prepped.csv", scale = "linear",
#'                              sample_regex = "^(.+)_Plate [0-9]+$",
#'                              pool_pattern = "^pool",
#'                              metadata = "meta.csv")
#' ms$sample_meta$PlateID <- sub(".*_(Plate [0-9]+)$", "\\1",
#'                               ms$sample_meta$source_column)
#' }
ELEUTHIA_load_massspec <- function(ms_data,
                                    scale,
                                    id_col          = 1,
                                    sample_regex    = NULL,
                                    annotation_cols = NULL,
                                    sheet           = 1,
                                    metadata        = NULL,
                                    metadata_sheet  = 1,
                                    sample_col      = "SampleID",
                                    pool_pattern    = NULL,
                                    verbose         = TRUE) {

  if (missing(scale))
    stop("'scale' is required: use scale = \"log2\" if the values are already ",
         "log2-transformed, or scale = \"linear\" for raw intensities.",
         call. = FALSE)
  if (!is.character(scale) || length(scale) != 1L || !scale %in% c("log2", "linear"))
    stop("'scale' must be \"log2\" or \"linear\".", call. = FALSE)

  if (!is.null(sample_regex) && !is.null(annotation_cols))
    stop("Use either 'sample_regex' or 'annotation_cols', not both: with a ",
         "sample_regex every non-sample column is kept as annotation.",
         call. = FALSE)

  if (verbose) cat("[ELEUTHIA] Loading mass spectrometry data\n")

  # ---------------------------------------------------------------------------
  # Read the quantification table
  # ---------------------------------------------------------------------------

  raw <- .massspec_read_table(ms_data, sheet = sheet, what = "ms_data")
  if (verbose)
    cat("    Input: ", if (is.data.frame(ms_data)) "data.frame" else ms_data,
        " (", nrow(raw), " rows x ", ncol(raw), " columns)\n", sep = "")

  # ---------------------------------------------------------------------------
  # Protein identifier column
  # ---------------------------------------------------------------------------

  id_idx <- .massspec_resolve_col(id_col, colnames(raw), "id_col")
  id_name <- colnames(raw)[id_idx]

  ids_chr <- trimws(as.character(raw[[id_idx]]))
  blank   <- is.na(ids_chr) | !nzchar(ids_chr)
  n_blank <- sum(blank)
  if (n_blank > 0L) {
    raw     <- raw[!blank, , drop = FALSE]
    ids_chr <- ids_chr[!blank]
    if (verbose)
      cat("    Dropped ", n_blank, " row(s) with a blank '", id_name, "'\n", sep = "")
  }
  if (nrow(raw) == 0L)
    stop("No rows left after dropping rows with a blank '", id_name, "'.", call. = FALSE)

  if (anyDuplicated(ids_chr)) {
    dups <- unique(ids_chr[duplicated(ids_chr)])
    stop("Duplicated protein identifier(s) in column '", id_name, "' (",
         length(dups), " distinct): ",
         paste(utils::head(dups, 10), collapse = ", "),
         if (length(dups) > 10L) ", ..." else "",
         ".\nChoose a column with unique values as 'id_col', or resolve the ",
         "duplicates before loading.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Sample columns (and the sample IDs they are named by)
  # ---------------------------------------------------------------------------

  other_idx <- setdiff(seq_len(ncol(raw)), id_idx)

  if (!is.null(sample_regex)) {
    if (!is.character(sample_regex) || length(sample_regex) != 1L)
      stop("'sample_regex' must be a single character string.", call. = FALSE)
    nms  <- colnames(raw)[other_idx]
    hits <- regmatches(nms, regexec(sample_regex, nms, perl = TRUE))
    matched <- lengths(hits) > 0L
    if (!any(matched))
      stop("'sample_regex' matched none of the column names. First columns: ",
           paste(utils::head(nms, 8), collapse = ", "), call. = FALSE)
    sample_idx <- other_idx[matched]
    # captured text if the regex has a capture group, otherwise the whole column name
    sample_ids <- vapply(seq_along(hits[matched]), function(i) {
      h <- hits[matched][[i]]
      if (length(h) >= 2L) h[2L] else nms[matched][i]
    }, character(1L))
    annot_idx  <- other_idx[!matched]
  } else {
    if (!is.null(annotation_cols)) {
      annot_idx <- vapply(annotation_cols, .massspec_resolve_col, integer(1L),
                          nms = colnames(raw), arg = "annotation_cols")
      annot_idx <- setdiff(unname(annot_idx), id_idx)
    } else {
      annot_idx <- integer(0L)
    }
    sample_idx <- setdiff(other_idx, annot_idx)
    sample_ids <- colnames(raw)[sample_idx]
  }

  if (length(sample_idx) == 0L)
    stop("No sample columns found.", call. = FALSE)

  bad_id <- is.na(sample_ids) | !nzchar(trimws(sample_ids))
  if (any(bad_id))
    stop("'sample_regex' produced an empty sample ID for column(s): ",
         paste(utils::head(colnames(raw)[sample_idx][bad_id], 5), collapse = ", "),
         call. = FALSE)
  if (anyDuplicated(sample_ids)) {
    dups <- unique(sample_ids[duplicated(sample_ids)])
    stop("Duplicated sample ID(s): ", paste(utils::head(dups, 10), collapse = ", "),
         if (length(dups) > 10L) ", ..." else "",
         ".\nSample IDs must be unique; adjust 'sample_regex' or the column ",
         "names.", call. = FALSE)
  }

  source_names <- colnames(raw)[sample_idx]

  # ---------------------------------------------------------------------------
  # Numeric matrix, then the log2 scale
  # ---------------------------------------------------------------------------

  num_cols <- lapply(sample_idx, function(j) .massspec_as_numeric(raw[[j]]))
  not_num  <- vapply(num_cols, is.null, logical(1L))
  if (any(not_num))
    stop("Sample column(s) that are not numeric: ",
         paste(utils::head(source_names[not_num], 10), collapse = ", "),
         if (sum(not_num) > 10L) ", ..." else "",
         ".\nIf they are annotation, list them in 'annotation_cols' (or select ",
         "the sample columns with 'sample_regex'); otherwise remove the text ",
         "values before loading.", call. = FALSE)

  wide <- matrix(unlist(num_cols, use.names = FALSE),
                 nrow = nrow(raw), ncol = length(sample_idx),
                 dimnames = list(ids_chr, sample_ids))

  n_nonpositive <- 0L
  if (scale == "linear") {
    nonpos        <- !is.na(wide) & wide <= 0
    n_nonpositive <- sum(nonpos)
    wide[nonpos]  <- NA
    if (all(is.na(wide)))
      stop("No positive values in the sample columns; is 'scale = \"linear\"' correct?",
           call. = FALSE)
    if (max(wide, na.rm = TRUE) < 100)
      warning("scale = \"linear\", but every value is below 100: the data may ",
              "already be log-transformed. Check the scale (heuristic check).",
              call. = FALSE)
    wide <- log2(wide)
  } else {
    if (any(is.finite(wide)) && stats::median(wide, na.rm = TRUE) > 100)
      warning("scale = \"log2\", but the median value is above 100: the data may ",
              "be on a linear scale. Check the scale (heuristic check).",
              call. = FALSE)
  }

  na_pct <- round(100 * mean(is.na(wide)), 1)
  if (verbose) {
    cat("    Matrix: ", nrow(wide), " proteins x ", ncol(wide), " samples | NA: ",
        na_pct, "%\n", sep = "")
    cat("    Scale: ", scale,
        if (scale == "linear")
          paste0(" -> log2 (", n_nonpositive, " zero/negative value(s) set to NA)"),
        "\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # protein_meta: the ID column first, then the annotation columns
  # ---------------------------------------------------------------------------

  protein_meta <- raw[, c(id_idx, annot_idx), drop = FALSE]
  protein_meta[[1L]] <- ids_chr
  rownames(protein_meta) <- ids_chr

  # ---------------------------------------------------------------------------
  # sample_meta: type, source column, and joined metadata
  # ---------------------------------------------------------------------------

  is_pool <- if (!is.null(pool_pattern)) {
    grepl(pool_pattern, sample_ids, ignore.case = TRUE, perl = TRUE)
  } else {
    rep(FALSE, length(sample_ids))
  }

  sample_type <- rep("SAMPLE", length(sample_ids))
  meta <- NULL

  if (!is.null(metadata)) {
    meta <- .massspec_read_table(metadata, sheet = metadata_sheet, what = "metadata")
    if (!sample_col %in% colnames(meta))
      stop("sample_col '", sample_col, "' not found in the metadata. Columns: ",
           paste(utils::head(colnames(meta), 10), collapse = ", "), call. = FALSE)

    meta[[sample_col]] <- trimws(as.character(meta[[sample_col]]))
    meta <- meta[!is.na(meta[[sample_col]]) & nzchar(meta[[sample_col]]), , drop = FALSE]
    if (anyDuplicated(meta[[sample_col]])) {
      dups <- unique(meta[[sample_col]][duplicated(meta[[sample_col]])])
      stop("Duplicated sample ID(s) in the metadata column '", sample_col, "': ",
           paste(utils::head(dups, 10), collapse = ", "), call. = FALSE)
    }
    if ("source_column" %in% colnames(meta))
      warning("The metadata has a 'source_column' column; it is not joined ",
              "because $sample_meta$source_column holds the original ",
              "ms_data column names.", call. = FALSE)

    m_idx   <- match(sample_ids, meta[[sample_col]])
    matched <- !is.na(m_idx)

    if ("SampleType" %in% colnames(meta)) {
      mtype <- trimws(as.character(meta$SampleType[m_idx[matched]]))
      bad   <- is.na(mtype) | !nzchar(mtype)
      if (any(bad))
        stop("Missing 'SampleType' in the metadata for sample(s): ",
             paste(utils::head(sample_ids[matched][bad], 10), collapse = ", "),
             call. = FALSE)
      sample_type[matched] <- mtype
    }
    sample_type[!matched] <- "OTHER"
  }
  sample_type[is_pool] <- "POOL"

  sample_meta <- data.frame(sample_ids, sample_type, source_names,
                            stringsAsFactors = FALSE)
  colnames(sample_meta) <- c(sample_col, "SampleType", "source_column")

  if (!is.null(meta)) {
    join_cols <- setdiff(colnames(meta), c(sample_col, "SampleType", "source_column"))
    if (length(join_cols) > 0L)
      sample_meta <- cbind(sample_meta,
                           meta[m_idx, join_cols, drop = FALSE],
                           stringsAsFactors = FALSE)
    if (verbose)
      cat("    Metadata: joined ", length(join_cols), " column(s); ",
          sum(matched), " of ", length(sample_ids), " samples matched\n", sep = "")
  }
  rownames(sample_meta) <- sample_ids

  if (verbose) {
    type_tbl <- table(sample_meta$SampleType)
    cat("    Sample types: ",
        paste(names(type_tbl), type_tbl, sep = "=", collapse = ", "), "\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # Assemble
  # ---------------------------------------------------------------------------

  params <- list(
    ms_data         = if (is.data.frame(ms_data)) "<data.frame>" else ms_data,
    metadata        = if (is.null(metadata)) NULL
                      else if (is.data.frame(metadata)) "<data.frame>" else metadata,
    id_col          = id_name,
    sample_regex    = sample_regex,
    annotation_cols = annotation_cols,
    scale           = scale,
    sample_col      = sample_col,
    pool_pattern    = pool_pattern,
    n_blank_id_rows_dropped      = n_blank,
    n_nonpositive_values_set_na  = n_nonpositive
  )

  result <- list(
    wide         = wide,
    sample_meta  = sample_meta,
    protein_meta = protein_meta,
    params       = params
  )
  class(result) <- c("massspec_data", "list")

  if (verbose) cat("[ELEUTHIA] massspec_data object ready\n")

  result
}


# --- Internal helpers ---------------------------------------------------------

# Strings treated as missing in the SAMPLE columns of a mass spectrometry table
# (besides NA). Text columns (annotation, metadata) are left exactly as read, so a
# blank annotation cell stays "" rather than becoming NA.
#' @noRd
.massspec_na_strings <- function() c("", "NA", "NaN", "N/A", "#N/A", "Filtered")

# A table from a path or a data.frame. `what` names the argument in messages.
#' @noRd
.massspec_read_table <- function(x, sheet = 1, what = "ms_data") {
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE,
                                             check.names = FALSE))
  if (!is.character(x) || length(x) != 1L)
    stop("'", what, "' must be a file path or a data.frame.", call. = FALSE)
  if (!file.exists(x)) stop(what, " file not found: ", x, call. = FALSE)

  ext <- tolower(sub(".*\\.", "", basename(x)))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Package 'readxl' is required to read .xlsx files.\n",
           "Install with: install.packages('readxl')", call. = FALSE)
    as.data.frame(readxl::read_excel(x, sheet = sheet, .name_repair = "minimal"))
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported format for '", what, "': '.", ext,
         "'. Use .csv, .tsv, .txt or .xlsx, or pass a data.frame.", call. = FALSE)
  }
}

# Position of a column given by name or number; errors name the argument.
#' @noRd
.massspec_resolve_col <- function(col, nms, arg) {
  if (length(col) != 1L)
    stop("'", arg, "' must be a single column name or position.", call. = FALSE)
  if (is.numeric(col)) {
    if (is.na(col) || col < 1 || col > length(nms) || col != round(col))
      stop("'", arg, "' = ", col, " is outside the ", length(nms), " columns.",
           call. = FALSE)
    return(as.integer(col))
  }
  idx <- match(col, nms)
  if (is.na(idx))
    stop("'", arg, "' column '", col, "' not found. Columns: ",
         paste(utils::head(nms, 10), collapse = ", "),
         if (length(nms) > 10L) ", ..." else "", call. = FALSE)
  idx
}

# A numeric vector from a column, or NULL if it holds text that is not a number
# or a recognised missing-value string.
#' @noRd
.massspec_as_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x) && all(is.na(x))) return(as.numeric(x))
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% .massspec_na_strings()] <- NA
  num <- suppressWarnings(as.numeric(x_chr))
  if (any(is.na(num) & !is.na(x_chr))) return(NULL)
  num
}


#' @method print massspec_data
#' @export
print.massspec_data <- function(x, ...) {
  cat("Mass Spectrometry Data Object\n")
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
