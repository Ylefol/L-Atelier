#' Eleuthia - Sample Sheet Functions
#'
#' @description Functions for validating and processing sample sheets
#' used in multi-omics analysis pipelines.


#' Validate Sample Sheet
#'
#' @description Validates a sample sheet for use in multi-omics analysis.
#' Checks for required columns, valid values, file existence, and generates
#' sample IDs programmatically.
#'
#' @param sample_sheet Either a data.frame containing sample information, or
#'   a character string path to a CSV file.
#' @param check_files Logical. If TRUE, verifies that data files exist at
#'   specified paths (default = TRUE).
#' @param check_beds Logical. If TRUE, verifies that BED files for quantification
#'   exist where required by format (peaks, bed). Only checked if check_files = TRUE
#'   (default = TRUE).
#' @param verbose Logical. If TRUE, prints validation progress and warnings
#'   (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{sample_sheet}{The validated and cleaned sample sheet with added
#'     sample_id column}
#'   \item{valid}{Logical. TRUE if validation passed with no errors}
#'   \item{errors}{Character vector of error messages (empty if valid)}
#'   \item{warnings}{Character vector of warning messages}
#' }
#'
#' @details
#' The function operates in two modes depending on whether a \code{sample_id}
#' column is present:
#'
#' \strong{Experimental mode} (no \code{sample_id} column): all of
#' \code{bio_rep}, \code{tech_rep}, and \code{batch} are required, and
#' \code{sample_id} is auto-generated as
#' \code{<omics>_<group>_<bio_rep><tech_rep>_b<batch>}.
#'
#' \strong{Cohort mode} (\code{sample_id} column present): \code{bio_rep},
#' \code{tech_rep}, and \code{batch} are optional. The provided
#' \code{sample_id} values are used directly and validated for uniqueness.
#' Suitable for public datasets (e.g., TCGA) where each sample is a unique
#' individual rather than a replicate within an experimental design.
#'
#' Required columns in both modes:
#' \itemize{
#'   \item file_loc: Directory containing the data file
#'   \item file_name: Name of the data file
#'   \item group: Experimental group (e.g., "WT", "SMUG1_KO", "male", "female")
#'   \item omics: Type of omics data (free-form label, e.g., "ATACseq", "CHIPseq",
#'     "RNAseq", "5hmu", "CUT&RUN", etc.)
#'   \item format: File format ("peaks", "bed", "counts")
#'   \item bed_loc: Path to fragment BED file for quantification (required for
#'     peaks and bed formats). Typically created with bedtools bamtobed.
#' }
#'
#' Additional required columns in experimental mode only:
#' \itemize{
#'   \item bio_rep: Biological replicate identifier
#'   \item tech_rep: Technical replicate identifier
#'   \item batch: Batch number
#' }
#'
#' The function performs the following validations:
#' \enumerate{
#'   \item Checks all required columns are present
#'   \item Removes empty rows
#'   \item Validates format values against allowed types
#'   \item Checks data files exist (if check_files = TRUE)
#'   \item Checks quantification BED files exist for peak/bed formats (if check_beds = TRUE)
#'   \item Generates unique sample_id for each row
#' }
#'
#' Note: The omics column accepts any string value. This allows flexibility for
#' custom omics types (e.g., "5hmu", "CUT&RUN", "MeDIP") that may be processed
#' similarly to standard types. The format column determines how data is loaded.
#'
#' In experimental mode, sample IDs are generated as:
#' \code{<omics>_<group>_<bio_rep><tech_rep>_b<batch>}
#'
#' If a \code{timepoint} column is present, it is appended for uniqueness:
#' \code{<omics>_<group>_<bio_rep><tech_rep>_b<batch>_t<timepoint>}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From CSV file
#' result <- ELEUTHIA_validate_sample_sheet("data/sample_sheet.csv")
#'
#' # From data.frame
#' result <- ELEUTHIA_validate_sample_sheet(my_sample_df)
#'
#' # Check validation status
#' if (result$valid) {
#'   sample_sheet <- result$sample_sheet
#' } else {
#'   stop(paste(result$errors, collapse = "\n"))
#' }
#'
#' # Skip file checks for quick validation
#' result <- ELEUTHIA_validate_sample_sheet(df, check_files = FALSE)
#'
#' }
ELEUTHIA_validate_sample_sheet <- function(sample_sheet,
                                            check_files = TRUE,
                                            check_beds = TRUE,
                                            verbose = TRUE) {

  errors <- character(0)
  warnings <- character(0)

  # ---------------------------------------------------------------------------
  # Load sample sheet if path provided
  # ---------------------------------------------------------------------------
  if (is.character(sample_sheet) && length(sample_sheet) == 1) {
    if (!file.exists(sample_sheet)) {
      return(list(
        sample_sheet = NULL,
        valid = FALSE,
        errors = paste("Sample sheet file not found:", sample_sheet),
        warnings = character(0)
      ))
    }
    if (verbose) cat("[ELEUTHIA] Loading sample sheet from:", sample_sheet, "\n")
    sample_sheet <- read.csv(sample_sheet, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(sample_sheet)) {
    return(list(
      sample_sheet = NULL,
      valid = FALSE,
      errors = "sample_sheet must be a data.frame or path to CSV file",
      warnings = character(0)
    ))
  }

  # ---------------------------------------------------------------------------
  # Detect cohort mode (user-supplied sample_id → bio_rep/tech_rep/batch optional)
  # ---------------------------------------------------------------------------
  cohort_mode <- "sample_id" %in% colnames(sample_sheet)
  if (verbose) {
    if (cohort_mode) {
      cat("[ELEUTHIA] Cohort mode detected (sample_id column present).",
          "bio_rep/tech_rep/batch are optional.\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Check required columns
  # ---------------------------------------------------------------------------
  base_required <- c("file_loc", "file_name", "group", "omics", "format", "bed_loc")
  extra_required <- if (cohort_mode) character(0) else c("bio_rep", "tech_rep", "batch")
  required_cols  <- c(base_required, extra_required)

  missing_cols <- setdiff(required_cols, colnames(sample_sheet))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required columns:",
                              paste(missing_cols, collapse = ", ")))
    return(list(
      sample_sheet = sample_sheet,
      valid = FALSE,
      errors = errors,
      warnings = warnings
    ))
  }

  if (verbose) cat("[ELEUTHIA] All required columns present.\n")

  # ---------------------------------------------------------------------------
  # Remove empty rows
  # ---------------------------------------------------------------------------
  # A row is empty if all required fields are NA or empty string
  # Use only columns that are actually present (guards against optional cols)
  check_cols <- intersect(required_cols, colnames(sample_sheet))
  empty_row_check <- apply(sample_sheet[, check_cols, drop = FALSE], 1, function(row) {
    all(is.na(row) | trimws(as.character(row)) == "")
  })

  n_empty <- sum(empty_row_check)
  if (n_empty > 0) {
    sample_sheet <- sample_sheet[!empty_row_check, , drop = FALSE]
    if (verbose) cat("[ELEUTHIA] Removed", n_empty, "empty row(s).\n")
  }

  if (nrow(sample_sheet) == 0) {
    errors <- c(errors, "Sample sheet contains no valid rows after removing empty rows")
    return(list(
      sample_sheet = sample_sheet,
      valid = FALSE,
      errors = errors,
      warnings = warnings
    ))
  }

  # Reset row names
  rownames(sample_sheet) <- NULL

  # ---------------------------------------------------------------------------
  # Validate format values
  # ---------------------------------------------------------------------------
  valid_formats <- c("peaks", "bed", "counts")
  invalid_formats <- unique(sample_sheet$format[!sample_sheet$format %in% valid_formats])

  if (length(invalid_formats) > 0) {
    errors <- c(errors, paste("Invalid format values:",
                              paste(invalid_formats, collapse = ", "),
                              "\nAllowed values:", paste(valid_formats, collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # Check for missing values in required fields
  # ---------------------------------------------------------------------------
  check_na_cols <- c("file_loc", "file_name", "group", "omics", "format")
  if (!cohort_mode) {
    check_na_cols <- c(check_na_cols, "bio_rep", "tech_rep", "batch")
  } else {
    # In cohort mode, sample_id itself must not be missing
    sid_na <- sum(is.na(sample_sheet$sample_id) | trimws(sample_sheet$sample_id) == "")
    if (sid_na > 0) {
      errors <- c(errors, paste("Column sample_id has", sid_na, "missing value(s)"))
    }
  }
  for (col in check_na_cols) {
    na_count <- sum(is.na(sample_sheet[[col]]) | trimws(as.character(sample_sheet[[col]])) == "")
    if (na_count > 0) {
      errors <- c(errors, paste("Column", col, "has", na_count, "missing value(s)"))
    }
  }

  # bed_loc can be NA for counts format
  needs_bed <- sample_sheet$format %in% c("peaks", "bed")
  bed_missing <- needs_bed & (is.na(sample_sheet$bed_loc) |
                               trimws(sample_sheet$bed_loc) == "" |
                               sample_sheet$bed_loc == "NA")
  if (any(bed_missing)) {
    errors <- c(errors, paste("bed_loc is required for peaks/bed formats but missing for",
                              sum(bed_missing), "row(s)"))
  }

  # ---------------------------------------------------------------------------
  # Check file existence
  # ---------------------------------------------------------------------------
  if (check_files && length(errors) == 0) {
    if (verbose) cat("[ELEUTHIA] Checking data file paths...\n")

    for (i in seq_len(nrow(sample_sheet))) {
      file_path <- file.path(sample_sheet$file_loc[i], sample_sheet$file_name[i])
      if (!file.exists(file_path)) {
        warnings <- c(warnings, paste("Data file not found:", file_path))
      }
    }

    if (length(warnings) > 0 && verbose) {
      cat("[ELEUTHIA] Warning:", length(warnings), "data file(s) not found.\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Check quantification BED file existence
  # ---------------------------------------------------------------------------
  if (check_files && check_beds && length(errors) == 0) {
    if (verbose) cat("[ELEUTHIA] Checking quantification BED file paths...\n")

    bed_rows <- which(sample_sheet$format %in% c("peaks", "bed"))
    bed_missing_files <- character(0)

    for (i in bed_rows) {
      bed_path <- sample_sheet$bed_loc[i]
      if (!is.na(bed_path) && bed_path != "NA" && !file.exists(bed_path)) {
        bed_missing_files <- c(bed_missing_files, bed_path)
      }
    }

    if (length(bed_missing_files) > 0) {
      # Unique BED paths that are missing
      unique_missing <- unique(bed_missing_files)
      warnings <- c(warnings, paste("Quantification BED file not found:",
                                    paste(unique_missing, collapse = "\n  ")))
      if (verbose) {
        cat("[ELEUTHIA] Warning:", length(unique_missing), "unique quantification BED file(s) not found.\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Generate or validate sample IDs
  # ---------------------------------------------------------------------------
  if (length(errors) == 0) {
    if (cohort_mode) {
      if (verbose) cat("[ELEUTHIA] Using provided sample_id column.\n")
      # Validate uniqueness of user-supplied IDs
      dup_ids <- sample_sheet$sample_id[duplicated(sample_sheet$sample_id)]
      if (length(dup_ids) > 0) {
        errors <- c(errors, paste("Duplicate sample_id values found:",
                                  paste(unique(dup_ids), collapse = ", ")))
      }
    } else {
      if (verbose) cat("[ELEUTHIA] Generating sample IDs...\n")

      # Base sample ID: {omics}_{group}_{bio_rep}{tech_rep}_b{batch}
      sample_sheet$sample_id <- paste0(
        sample_sheet$omics, "_",
        sample_sheet$group, "_",
        sample_sheet$bio_rep,
        sample_sheet$tech_rep, "_b",
        sample_sheet$batch
      )

      # If timepoint column exists, append _t{timepoint} for uniqueness
      if ("timepoint" %in% colnames(sample_sheet)) {
        if (verbose) cat("[ELEUTHIA] Timepoint column detected - including in sample IDs\n")
        sample_sheet$sample_id <- paste0(
          sample_sheet$sample_id, "_t",
          sample_sheet$timepoint
        )
      }

      # Check for duplicate sample IDs
      dup_ids <- sample_sheet$sample_id[duplicated(sample_sheet$sample_id)]
      if (length(dup_ids) > 0) {
        errors <- c(errors, paste("Duplicate sample IDs generated:",
                                  paste(unique(dup_ids), collapse = ", "),
                                  "\nThis indicates duplicate entries in the sample sheet."))
      }
    }
  }
  # Set rownames before returning
  rownames(sample_sheet) <- sample_sheet$sample_id
  
  # ---------------------------------------------------------------------------
  # Final validation status
  # ---------------------------------------------------------------------------
  valid <- length(errors) == 0

  if (verbose) {
    if (valid) {
      cat("\n[ELEUTHIA] Validation PASSED.\n")
      cat("    Total samples:", nrow(sample_sheet), "\n")
      cat("    Omics types:", paste(unique(sample_sheet$omics), collapse = ", "), "\n")
      cat("    Groups:", paste(unique(sample_sheet$group), collapse = ", "), "\n")
      if (length(warnings) > 0) {
        cat("    Warnings:", length(warnings), "\n")
      }
    } else {
      cat("\n[ELEUTHIA] Validation FAILED.\n")
      cat("    Errors:", length(errors), "\n")
      for (err in errors) {
        cat("    -", err, "\n")
      }
    }
  }
  
  return(list(
    sample_sheet = sample_sheet,
    valid = valid,
    errors = errors,
    warnings = warnings
  ))
}


#' Get Subset of Sample Sheet by Omics Type
#'
#' @description Extracts samples of a specific omics type from a validated
#' sample sheet.
#'
#' @param sample_sheet A validated sample sheet data.frame (with sample_id column).
#' @param omics Character string. The omics type to extract (e.g., "ATACseq",
#'   "CHIPseq", "RNAseq", "5hmu", or any custom omics label).
#'
#' @return A data.frame containing only rows matching the specified omics type.
#'   Returns empty data.frame if no matches found.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' atac_samples <- ELEUTHIA_get_omics_subset(result$sample_sheet, "ATACseq")
#' custom_samples <- ELEUTHIA_get_omics_subset(result$sample_sheet, "5hmu")
#'
#' }
ELEUTHIA_get_omics_subset <- function(sample_sheet, omics) {

  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]
  rownames(subset_df) <- NULL

  if (nrow(subset_df) == 0) {
    warning("No samples found for omics type: ", omics)
  }

  return(subset_df)
}


#' Summarize Sample Sheet Structure
#'
#' @description Prints a summary of the sample sheet structure showing
#' the distribution of samples across omics types, groups, batches, etc.
#'
#' @param sample_sheet A validated sample sheet data.frame.
#'
#' @return Invisibly returns a list with summary statistics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' ELEUTHIA_summarize_sample_sheet(result$sample_sheet)
#'
#' }
ELEUTHIA_summarize_sample_sheet <- function(sample_sheet) {

  cat("================================================================================\n")
  cat("[ELEUTHIA] SAMPLE SHEET SUMMARY\n")
  cat("================================================================================\n\n")

  cat("    Total samples:", nrow(sample_sheet), "\n\n")

  # By omics
  cat("    By Omics Type:\n")
  omics_table <- table(sample_sheet$omics)
  for (om in names(omics_table)) {
    cat(sprintf("    %-10s: %d samples\n", om, omics_table[om]))
  }
  cat("\n")

  # By group
  cat("    By Experimental Group:\n")
  group_table <- table(sample_sheet$group)
  for (grp in names(group_table)) {
    cat(sprintf("    %-15s: %d samples\n", grp, group_table[grp]))
  }
  cat("\n")

  # By batch
  cat("    By Batch:\n")
  batch_table <- table(sample_sheet$batch)
  for (b in names(batch_table)) {
    cat(sprintf("    Batch %s: %d samples\n", b, batch_table[b]))
  }
  cat("\n")

  # Cross-tabulation: omics x group
  cat("    Omics x Group:\n")
  cross_table <- table(sample_sheet$omics, sample_sheet$group)
  print(cross_table)
  cat("\n")

  # Biological replicates per group per omics
  cat("    Biological Replicates:\n")
  for (om in unique(sample_sheet$omics)) {
    om_data <- sample_sheet[sample_sheet$omics == om, ]
    cat(sprintf("    %s:\n", om))
    for (grp in unique(om_data$group)) {
      grp_data <- om_data[om_data$group == grp, ]
      n_bio <- length(unique(grp_data$bio_rep))
      cat(sprintf("    %s: %d biological rep(s)\n", grp, n_bio))
    }
  }

  cat("\n================================================================================\n")

  # Return summary stats invisibly
  invisible(list(
    total = nrow(sample_sheet),
    by_omics = as.list(omics_table),
    by_group = as.list(group_table),
    by_batch = as.list(batch_table)
  ))
}


#' Load Quantification BED Files from Sample Sheet
#'
#' @description Loads fragment BED files for quantification from the bed_loc
#' column of a validated sample sheet.
#'
#' @param sample_sheet A validated sample sheet data.frame with bed_loc column.
#' @param omics Character string. Omics type to load ("ATACseq" or "CHIPseq").
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A named list of BED data.frames, with names corresponding to sample_id.
#'
#' @details
#' This function reads fragment BED files (typically created with bedtools bamtobed)
#' for use with ELEUTHIA_quantify_bed(). It uses the bed_loc column from the
#' sample sheet rather than file_loc/file_name.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' atac_frags <- ELEUTHIA_load_quant_beds(result$sample_sheet, "ATACseq")
#' chip_frags <- ELEUTHIA_load_quant_beds(result$sample_sheet, "CHIPseq")
#'
#' }
ELEUTHIA_load_quant_beds <- function(sample_sheet,
                                      omics,
                                      verbose = TRUE) {

  # Check bed_loc column exists
  if (!"bed_loc" %in% colnames(sample_sheet)) {
    stop("sample_sheet must have a 'bed_loc' column")
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

  if (verbose) {
    cat("[ELEUTHIA] Loading", nrow(subset_df), omics, "quantification BED files...\n")
  }

  bed_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    bed_path <- subset_df$bed_loc[i]

    if (verbose) {
      cat("  Loading:", sample_id, "\n")
    }

    if (!file.exists(bed_path)) {
      stop("BED file not found: ", bed_path)
    }

    # Load BED file (simple 3-column minimum)
    bed <- read.table(
      bed_path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      comment.char = "#"
    )

    # Assign column names
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

    bed_list[[sample_id]] <- bed
  }

  if (verbose) {
    total_frags <- sum(sapply(bed_list, nrow))
    cat("[ELEUTHIA] Loaded", length(bed_list), "BED files.\n")
    cat("[ELEUTHIA] Total fragments:", format(total_frags, big.mark = ","), "\n")
  }

  return(bed_list)
}
