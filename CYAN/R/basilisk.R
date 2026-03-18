#' Basilisk Python environments for CYAN
#'
#' @description
#' Managed conda environments used internally by CYAN functions.
#' Not exported; accessed only through \code{basiliskRun()} calls.
#'
#' @name cyan-basilisk-envs
#' @keywords internal
NULL

#' REDItools2 basilisk environment
#'
#' Python dependencies for \code{CYAN_run_reditools()}.
#' Packages pinned to validated, compatible versions.
#'
#' @keywords internal
.cyan_reditools_env <- basilisk::BasiliskEnvironment(
  envname  = "cyan_reditools_env",
  pkgname  = "CYAN",
  packages = c(
    "python==3.10.14",
    "pysam==0.22.1",
    "sortedcontainers==2.4.0",
    "psutil==6.0.0",
    "netifaces==0.11.0"
  ),
  channels = c("bioconda", "conda-forge", "defaults")
)
