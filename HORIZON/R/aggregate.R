#' Aggregate per-sample count files into a single count matrix
#'
#' Reads the per-sample count files produced by \code{\link{HORIZON_run_count}}
#' and combines them into a single gene x sample count matrix. Samples missing
#' a count file are skipped with a warning.
#'
#' Output is written to \code{<output_root>/aggregated/counts/}:
#' \itemize{
#'   \item \code{count_matrix.csv} — genes x samples count matrix
#'   \item \code{count_matrix.rds} — same matrix as an R object (if
#'     \code{save_rds = TRUE})
#'   \item \code{sample_metadata.csv} — sample sheet columns carried through
#'     as metadata (if \code{save_metadata = TRUE})
#' }
#'
#' @param sample_sheet Validated sample sheet data.frame from
#'   \code{\link{HORIZON_validate_sample_sheet}}.
#' @param output_root Character or NULL. Root output directory — the parent of
#'   the per-sample folders. If NULL, inferred from the first row's
#'   \code{output_dir}. Provide explicitly if samples have different
#'   \code{output_dir} values and you want the aggregated folder somewhere
#'   specific.
#' @param save_rds Logical. Save the count matrix as an RDS file in addition to
#'   CSV. Default TRUE.
#' @param save_metadata Logical. Save the sample sheet as
#'   \code{sample_metadata.csv} alongside the count matrix. Default TRUE.
#'
#' @return data.frame. The aggregated count matrix (genes x samples), invisibly.
#' @export
HORIZON_aggregate_counts <- function(sample_sheet,
                                     output_root    = NULL,
                                     save_rds       = TRUE,
                                     save_metadata  = TRUE) {
  if (is.null(output_root)) {
    output_root <- sample_sheet$output_dir[1]
    cat("output_root inferred from first sample: ", output_root)
  }

  agg_dir <- file.path(output_root, "aggregated", "counts")
  dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)

  count_list <- list()

  for (sid in sample_sheet$sample_id) {
    sample_out_dir <- sample_sheet$output_dir[sample_sheet$sample_id == sid]
    counts_file    <- file.path(sample_out_dir, sid, "counts",
                                paste0(sid, "_counts.txt"))

    if (!file.exists(counts_file)) {
      warning("Count file not found for sample '", sid, "' — skipped.")
      next
    }

    df              <- read.table(counts_file, header = TRUE, sep = "\t",
                                  stringsAsFactors = FALSE)
    count_list[[sid]] <- setNames(df$count, df$gene_id)
  }

  if (length(count_list) == 0) {
    stop("No count files found. Run HORIZON_run_count() for at least one sample first.")
  }

  # Align all samples on the same gene set (should be identical if same
  # annotation was used, but union handles any edge cases gracefully)
  all_genes    <- Reduce(union, lapply(count_list, names))
  count_matrix <- do.call(cbind, lapply(count_list, function(x) x[all_genes]))
  rownames(count_matrix) <- all_genes

  n_genes   <- nrow(count_matrix)
  n_samples <- ncol(count_matrix)
  cat("Aggregated count matrix: ", n_genes, " genes x ", n_samples, " samples")

  out_csv <- file.path(agg_dir, "count_matrix.csv")
  write.csv(count_matrix, out_csv)
  cat("Count matrix written to: ", out_csv)

  if (save_rds) {
    saveRDS(count_matrix, file.path(agg_dir, "count_matrix.rds"))
  }

  if (save_metadata) {
    # Only keep rows for samples that were successfully aggregated
    meta <- sample_sheet[sample_sheet$sample_id %in% colnames(count_matrix), ]
    write.csv(meta, file.path(agg_dir, "sample_metadata.csv"), row.names = FALSE)
    cat("Sample metadata written to: ", file.path(agg_dir, "sample_metadata.csv"))
  }

  invisible(count_matrix)
}
