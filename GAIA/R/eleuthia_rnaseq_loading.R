#' Eleuthia - RNA-seq Loading Functions
#'
#' @description Functions for loading RNA-seq count data.


#' Load a Single Count File
#'
#' @description Reads a simple two-column count file (gene_id, count).
#'
#' @param file_path Character string. Path to the count file.
#' @param gene_col Integer. Column index for gene IDs (default = 1).
#' @param count_col Integer. Column index for counts (default = 2).
#' @param header Logical. Does the file have a header? (default = FALSE).
#'
#' @return A named numeric vector of counts, with gene IDs as names.
#'
#' @export
#'
ELEUTHIA_load_count_file <- function(file_path,
                                      gene_col = 1,
                                      count_col = 2,
                                      header = FALSE) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  counts <- read.table(
    file_path,
    header = header,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
  )

  # Extract gene IDs and counts
  result <- counts[[count_col]]
  names(result) <- counts[[gene_col]]

  return(result)
}


#' Load RNA-seq Counts from Sample Sheet
#'
#' @description Loads all RNA-seq count files from a validated sample sheet
#' and combines them into a count matrix.
#'
#' @param sample_sheet A validated sample sheet data.frame.
#' @param omics Character string. Omics type to load (default = "RNAseq").
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{counts}{Matrix of counts (genes x samples)}
#'   \item{targets}{Data.frame with sample metadata}
#' }
#'
#' @details
#' All count files must have the same genes (rows). The function will
#' use the intersection of genes if there are differences, with a warning.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
#' rna_data <- ELEUTHIA_load_rnaseq_from_sheet(result$sample_sheet)
#' count_matrix <- rna_data$counts
#'
#' }
ELEUTHIA_load_rnaseq_from_sheet <- function(sample_sheet,
                                             omics = "RNAseq",
                                             verbose = TRUE) {

  # Get subset for this omics type
  subset_df <- sample_sheet[sample_sheet$omics == omics, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No samples found for omics type: ", omics)
  }

  if (verbose) {
    cat("Loading", nrow(subset_df), omics, "count files...\n")
  }

  # Load all count files
  count_list <- list()

  for (i in seq_len(nrow(subset_df))) {
    sample_id <- subset_df$sample_id[i]
    file_path <- file.path(subset_df$file_loc[i], subset_df$file_name[i])

    if (verbose) {
      cat("  Loading:", sample_id, "\n")
    }

    count_list[[sample_id]] <- ELEUTHIA_load_count_file(file_path)
  }

  # Check gene consistency
  all_genes <- lapply(count_list, names)
  gene_sets_equal <- all(sapply(all_genes[-1], function(g) identical(g, all_genes[[1]])))

  if (!gene_sets_equal) {
    # Find common genes
    common_genes <- Reduce(intersect, all_genes)
    warning("Gene sets differ across samples. Using intersection of ",
            length(common_genes), " common genes.")

    # Subset to common genes
    count_list <- lapply(count_list, function(x) x[common_genes])
  }

  # Combine into matrix
  counts <- do.call(cbind, count_list)
  rownames(counts) <- names(count_list[[1]])
  colnames(counts) <- names(count_list)

  # Create targets data.frame
  targets <- subset_df[, c("sample_id", "group", "bio_rep", "tech_rep", "batch")]
  rownames(targets) <- targets$sample_id

  if (verbose) {
    cat("\nRNA-seq loading complete:\n")
    cat("  Genes:", nrow(counts), "\n")
    cat("  Samples:", ncol(counts), "\n")
    cat("  Total counts:", format(sum(counts), big.mark = ","), "\n")
  }

  return(list(
    counts = counts,
    targets = targets
  ))
}


#' Filter Low-Expressed Genes
#'
#' @description Filters genes with low expression across samples.
#' This is a convenience wrapper that works with RNA-seq data structure.
#'
#' @param rna_data A list from ELEUTHIA_load_rnaseq_from_sheet().
#' @param min_count Integer. Minimum count threshold (default = 10).
#' @param min_samples Integer. Minimum samples meeting threshold (default = 2).
#' @param verbose Logical. Print progress (default = TRUE).
#'
#' @return Filtered rna_data list.
#'
#' @details
#' Note: This is a simple count-based filter. For more sophisticated
#' filtering (e.g., CPM-based), use edgeR::filterByExpr() or similar.
#'
#' @export
#'
ELEUTHIA_filter_low_expression <- function(rna_data,
                                            min_count = 10,
                                            min_samples = 2,
                                            verbose = TRUE) {

  # This uses the same logic as POSEIDON_filter_low_counts
  # but is provided here for convenience with RNA-seq specific naming

  counts <- rna_data$counts

  n_above_threshold <- rowSums(counts >= min_count)
  keep <- n_above_threshold >= min_samples

  n_before <- nrow(counts)
  n_after <- sum(keep)

  if (verbose) {
    cat("Filtering low-expression genes:\n")
    cat("  Before:", n_before, "genes\n")
    cat("  After:", n_after, "genes\n")
    cat("  Removed:", n_before - n_after, "genes\n")
    cat("  (min_count =", min_count, ", min_samples =", min_samples, ")\n")
  }

  rna_data$counts <- counts[keep, , drop = FALSE]

  return(rna_data)
}
