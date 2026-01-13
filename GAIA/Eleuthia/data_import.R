###############################################################################
########### Data Import Functions ###########
###############################################################################

#' Load Count Data from Directory
#'
#' @description Loads count files from a directory and merges them into a single
#' data frame. Expects files with gene_id in first column and counts in second column.
#' Handles file naming conventions with numeric prefixes.
#'
#' @param count_path Character string. Path to directory containing count files.
#'   Each file should be a tab-delimited text file with gene_id and count columns.
#' @param file_pattern Optional regex pattern to filter files (default = NULL, uses all files)
#' @param gene_col Name or index of gene ID column (default = 1)
#' @param count_col Name or index of count column (default = 2)
#' @param clean_sample_names Logical, should sample names be cleaned? (default = TRUE)
#'   Removes numeric prefixes like "123_456-" from filenames.
#'
#' @return Data frame with genes as rows and samples as columns.
#'   Row names are gene IDs.
#'
#' @details
#' This function:
#' - Reads all files in the specified directory (or files matching file_pattern)
#' - Extracts sample names from filenames (removes file extensions)
#' - Optionally cleans sample names by removing numeric prefixes
#' - Merges all count files by gene_id
#' - Returns a matrix with genes (rows) × samples (columns)
#'
#' @export
#'
#' @examples
#' # Load count data from directory
#' counts <- ELEUTHIA_load_counts_from_dir("path/to/counts/")
#'
#' # With file pattern filtering
#' counts <- ELEUTHIA_load_counts_from_dir(
#'   "path/to/counts/",
#'   file_pattern = "\\.txt$"
#' )
#'
ELEUTHIA_load_counts_from_dir <- function(count_path,
                                          file_pattern = NULL,
                                          gene_col = 1,
                                          count_col = 2,
                                          clean_sample_names = TRUE) {
  
  # Validate path
  if (!dir.exists(count_path)) {
    stop(sprintf("Directory does not exist: %s", count_path))
  }

  # Get file list
  all_files <- list.files(count_path)

  if (length(all_files) == 0) {
    stop(sprintf("No files found in directory: %s", count_path))
  }

  # Filter by pattern if provided
  if (!is.null(file_pattern)) {
    all_files <- all_files[grepl(file_pattern, all_files)]
    if (length(all_files) == 0) {
      stop(sprintf("No files matching pattern '%s' found in: %s",
                   file_pattern, count_path))
    }
  }

  cat(sprintf("Loading %d files from %s\n", length(all_files), count_path))

  # Initialize
  full_data <- data.frame(NULL)

  # Load and merge files
  for (i in all_files) {
    file_path <- paste0(count_path, i)

    # Read file
    temp_dta <- tryCatch({
      read.table(file_path, header = FALSE)
    }, error = function(e) {
      warning(sprintf("Failed to read file %s: %s", i, e$message))
      return(NULL)
    })

    if (is.null(temp_dta)) next

    # Extract sample name from filename (remove extension)
    name <- unlist(strsplit(i, '\\.'))[c(TRUE, FALSE)]

    # Clean sample name if requested
    if (clean_sample_names) {
      clean_name <- gsub("^\\d+_\\d+-", "", name)
    } else {
      clean_name <- name
    }

    # Set column names
    colnames(temp_dta) <- c('gene_id', clean_name)

    # Merge with existing data
    if (nrow(full_data) == 0) {
      full_data <- temp_dta
    } else {
      full_data <- merge(full_data, temp_dta, by = 'gene_id')
    }
  }

  # Check if any data was loaded
  if (nrow(full_data) == 0) {
    stop("No data was successfully loaded")
  }

  # Set row names and remove gene_id column
  row.names(full_data) <- full_data$gene_id
  full_data <- full_data[, -1, drop = FALSE]

  cat(sprintf("Loaded %d genes across %d samples\n",
              nrow(full_data), ncol(full_data)))

  return(full_data)
}

#' Load Multiple Count Datasets from Directories
#'
#' @description Convenience wrapper to load multiple count datasets from
#' different directories. Useful for loading multiple omics types.
#'
#' @param count_paths Named list of directory paths.
#'   Names will be used as dataset identifiers.
#' @param ... Additional arguments passed to ELEUTHIA_load_counts_from_dir
#'
#' @return Named list of data frames, one per directory
#'
#' @export
#'
#' @examples
#' # Load RNA-seq and ATAC-seq data
#' datasets <- ELEUTHIA_load_multi_counts(
#'   count_paths = list(
#'     rnaseq = "path/to/rnaseq/",
#'     atacseq = "path/to/atacseq/"
#'   )
#' )
#'
ELEUTHIA_load_multi_counts <- function(count_paths, ...) {

  if (!is.list(count_paths)) {
    stop("count_paths must be a named list")
  }

  if (is.null(names(count_paths)) || any(names(count_paths) == "")) {
    stop("count_paths must be a named list with non-empty names")
  }

  results <- list()

  for (dataset_name in names(count_paths)) {
    cat(sprintf("\n=== Loading %s ===\n", dataset_name))
    results[[dataset_name]] <- ELEUTHIA_load_counts_from_dir(
      count_paths[[dataset_name]],
      ...
    )
  }

  return(results)
}
