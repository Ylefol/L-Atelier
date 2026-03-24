# HORIZON environment setup
#
# All external tools are accessed through a single user-managed conda
# environment.  The user creates this environment once and registers its
# path via HORIZON_set_conda_env() at the start of each session.
#
# Required tools in the conda environment:
#   samtools  >= 1.21
#   bedtools  >= 2.31
#   macs3     >= 3.0
#   deeptools >= 3.5
#
# Recommended setup (run once in a terminal):
#   conda create -n horizon_cli -c bioconda -c conda-forge \
#     samtools=1.23.1 bedtools=2.31.1 macs3=3.0.2 deeptools=3.5.5
#
# Then in R at the start of each session:
#   HORIZON_set_conda_env("~/miniconda3/envs/horizon_cli")


# ---------------------------------------------------------------------------
# HORIZON_set_conda_env
# ---------------------------------------------------------------------------

#' Set the conda environment used by HORIZON for all external tools
#'
#' Registers the path to a conda environment containing \code{samtools},
#' \code{bedtools}, \code{macs3}, and \code{deeptools}.  Must be called
#' once per session before any pipeline function that invokes these tools.
#'
#' Recommended setup (run once in a terminal):
#' \preformatted{
#' conda create -n horizon_cli -c bioconda -c conda-forge \
#'   samtools=1.23.1 bedtools=2.31.1 macs3=3.0.2 deeptools=3.5.5
#' }
#'
#' @param path Character. Path to the conda environment directory
#'   (the directory that contains a \code{bin/} subdirectory).
#'
#' @return The resolved path, invisibly.
#' @export
HORIZON_set_conda_env <- function(path) {
  path <- path.expand(path)
  if (!dir.exists(path))
    stop("Conda environment directory not found: ", path, call. = FALSE)
  bin_dir <- file.path(path, "bin")
  if (!dir.exists(bin_dir))
    stop("No bin/ directory found in: ", path,
         "\n  Ensure 'path' points to a conda environment root.", call. = FALSE)
  options(horizon.conda_env = path)
  message("[HORIZON] Conda env set: ", path)
  invisible(path)
}


# ---------------------------------------------------------------------------
# .horizon_run_cli (internal)
#
# Runs any tool from the user-registered conda environment via system2().
# Resolves the binary as <conda_env>/bin/<tool>.
#
# @param tool   Character. Executable name (e.g. "samtools", "macs3").
# @param args   Character vector of arguments passed to system2().
# @param stdout Passed to system2() stdout. "" = print; TRUE = capture.
# @param stderr Passed to system2() stderr. "" = print; TRUE = capture.
#
# @return Integer exit code, or character vector when stdout = TRUE.
# ---------------------------------------------------------------------------
.horizon_run_cli <- function(tool, args, stdout = "", stderr = "") {
  env_path <- getOption("horizon.conda_env", default = NULL)
  if (is.null(env_path))
    stop("HORIZON conda environment is not set.\n",
         "  Call HORIZON_set_conda_env('/path/to/conda/envs/horizon_cli') ",
         "before running pipeline functions.", call. = FALSE)
  bin <- file.path(env_path, "bin", tool)
  if (!file.exists(bin))
    stop("Tool '", tool, "' not found in conda env at: ",
         file.path(env_path, "bin"),
         "\n  Ensure the tool is installed: conda install -n horizon_cli -c bioconda ",
         tool, call. = FALSE)
  system2(bin, args = args, stdout = stdout, stderr = stderr, wait = TRUE)
}
