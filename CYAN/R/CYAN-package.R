#' CYAN: Niche Genomic Analyses
#'
#' CYAN handles specialised genomic analyses — RNA editing detection via
#' REDItools2 and somatic variant calling via GATK MuTect2 — positioned
#' between upstream BAM processing (HORIZON) and expression-level analysis
#' (GAIA/CAULDRON).
#'
#' Named after CYAN from Horizon Zero Dawn: a focused, independent analytical
#' unit rather than a large orchestrating system.
#'
#' @section Package design:
#' \itemize{
#'   \item All public functions are prefixed \code{CYAN_}
#'   \item One sample per call; cohort looping is left to the user
#'   \item Large intermediates written to \code{work_dir}; final outputs to
#'     \code{output_dir}
#'   \item Python tools (REDItools2) run through a managed basilisk environment
#'   \item Java tools (GATK/MuTect2) invoked via \code{system2("gatk", ...)}
#' }
#'
#' @importFrom basilisk BasiliskEnvironment basiliskRun
#' @importFrom reticulate py_exe
#' @importFrom methods is
#'
#' @docType package
#' @name CYAN-package
"_PACKAGE"
