#' @title HORIZON: Upstream Bioinformatics Processing
#'
#' @description
#' Upstream preprocessing pipeline for bulk RNA-seq data. Wraps Rfastp,
#' Rsubread, and Rsamtools into a consistent, RAM-aware interface designed
#' to feed downstream analysis in GAIA and CAULDRON.
#'
#' Typical workflow:
#' \enumerate{
#'   \item \code{\link{HORIZON_create_sample_sheet}} — generate a sample sheet template
#'   \item \code{\link{HORIZON_validate_sample_sheet}} — validate paths and columns
#'   \item \code{\link{HORIZON_build_index}} — build genome index (once per reference)
#'   \item \code{\link{HORIZON_run_qc_trim}} — QC + adapter trimming (per sample)
#'   \item \code{\link{HORIZON_run_align}} — splice-aware alignment (per sample)
#'   \item \code{\link{HORIZON_sort_index_bam}} — sort and index BAM (per sample)
#'   \item \code{\link{HORIZON_run_count}} — feature counting (per sample)
#'   \item \code{\link{HORIZON_aggregate_counts}} — combine into one count matrix
#' }
#'
#' @keywords internal
"_PACKAGE"
