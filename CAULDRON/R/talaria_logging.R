# ==============================================================================
# TALARIA - Logging
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# Console output logging via sink() with split=TRUE — duplicates stdout to a
# log file while keeping it visible in the console.  Warnings and errors are
# not captured (they remain console-only for immediate visibility).
#
# Typical usage:
#   log_file <- TALARIA_generate_log_filename("results/logs", prefix = "first_steps")
#   TALARIA_start_log(log_file)
#   # ... analysis code ...
#   TALARIA_stop_log()
#
# All functions prefixed: TALARIA_
# ==============================================================================


# Package-level environment storing logging state (one active log at a time).
.talaria_log_env <- new.env(parent = emptyenv())
.talaria_log_env$active     <- FALSE
.talaria_log_env$log_file   <- NULL
.talaria_log_env$output_con <- NULL
.talaria_log_env$start_time <- NULL


#' Generate a timestamp-based log filename
#'
#' Creates a log filename using the current date and time, suitable for
#' automatic naming so log files sort chronologically on disk.
#'
#' @param log_dir Character or \code{NULL}.  Directory in which to place the
#'   file.  \code{NULL} (default) returns just the filename with no path.
#' @param prefix Character.  Filename prefix.  Default \code{"log"}.
#' @param extension Character.  File extension.  Default \code{".log"}.
#' @param format Character.  Timestamp style:
#'   \describe{
#'     \item{\code{"compact"}}{Default. \code{YYYY_MM_DD__HH_MM} — filesystem-safe,
#'       no seconds (e.g. \code{log_2026_03_05__14_30.log}).}
#'     \item{\code{"readable"}}{\code{YYYY-MM-DD_HH-MM-SS} — includes seconds
#'       (e.g. \code{log_2026-03-05_14-30-52.log}).}
#'   }
#'
#' @return Character string — full file path if \code{log_dir} is provided,
#'   otherwise just the filename.
#' @export
TALARIA_generate_log_filename <- function(log_dir   = NULL,
                                           prefix    = "log",
                                           extension = ".log",
                                           format    = c("compact", "readable")) {

  format <- match.arg(format)
  ts     <- Sys.time()

  time_str <- if (format == "compact") {
    paste0(base::format(ts, "%Y_%m_%d"), "__", base::format(ts, "%H_%M"))
  } else {
    base::format(ts, "%Y-%m-%d_%H-%M-%S")
  }

  filename <- paste0(prefix, "_", time_str, extension)

  if (!is.null(log_dir)) file.path(log_dir, filename) else filename
}


#' Start console logging
#'
#' Begins duplicating standard output (\code{stdout}) to a log file using
#' \code{sink(split = TRUE)}, so output appears in both the console and the
#' file simultaneously.  Warnings and errors are intentionally not captured —
#' they remain console-only for immediate visibility.
#'
#' The log file opens with a header recording the start timestamp, R version,
#' platform, and working directory.
#'
#' Must be paired with \code{\link{TALARIA_stop_log}} to close the connection
#' cleanly and append session information.
#'
#' @param log_file Character.  Path to the log file to create (or append to).
#'   Parent directories are created automatically if they do not exist.
#' @param append Logical.  Append to an existing file (\code{TRUE}) or
#'   overwrite (\code{FALSE}, default).
#'
#' @return Invisibly returns \code{log_file}.
#' @export
TALARIA_start_log <- function(log_file, append = FALSE) {

  if (.talaria_log_env$active) {
    warning("Logging is already active. ",
            "Call TALARIA_stop_log() before starting a new log.",
            call. = FALSE)
    return(invisible(NULL))
  }

  if (!is.character(log_file) || length(log_file) != 1L)
    stop("log_file must be a single character string.", call. = FALSE)

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir))
    dir.create(log_dir, recursive = TRUE)

  .talaria_log_env$start_time <- Sys.time()
  .talaria_log_env$log_file   <- log_file
  .talaria_log_env$output_con <- file(log_file, open = if (append) "a" else "w")

  con <- .talaria_log_env$output_con
  cat("================================================================================\n", file = con)
  cat("LOG START:", base::format(.talaria_log_env$start_time, "%Y-%m-%d %H:%M:%S"), "\n", file = con)
  cat("================================================================================\n", file = con)
  cat("R version       :", R.version.string, "\n", file = con)
  cat("Platform        :", R.version$platform, "\n", file = con)
  cat("Working directory:", getwd(), "\n", file = con)
  cat("================================================================================\n\n", file = con)

  sink(file = con, append = TRUE, split = TRUE)

  .talaria_log_env$active <- TRUE

  cat(sprintf("Logging started: %s\n", log_file))

  invisible(log_file)
}


#' Stop console logging
#'
#' Stops duplicating output to the log file, appends a session information
#' block (end timestamp, duration, \code{sessionInfo()} output), and closes
#' the file connection cleanly.
#'
#' @return Invisibly returns the path to the closed log file, or \code{NULL}
#'   if no log was active.
#' @export
TALARIA_stop_log <- function() {

  if (!.talaria_log_env$active) {
    warning("No active log to stop.", call. = FALSE)
    return(invisible(NULL))
  }

  cat("\nLogging stopped.\n")
  sink()

  end_time <- Sys.time()
  duration <- difftime(end_time, .talaria_log_env$start_time, units = "auto")
  con      <- .talaria_log_env$output_con

  cat("\n================================================================================\n", file = con, append = TRUE)
  cat("SESSION INFORMATION\n",                                                              file = con, append = TRUE)
  cat("================================================================================\n", file = con, append = TRUE)
  cat("Log end  :", base::format(end_time, "%Y-%m-%d %H:%M:%S"), "\n",                    file = con, append = TRUE)
  cat("Duration :", base::format(duration), "\n\n",                                        file = con, append = TRUE)
  cat("R Session Info:\n",                                                                 file = con, append = TRUE)
  cat("--------------------------------------------------------------------------------\n", file = con, append = TRUE)
  writeLines(capture.output(sessionInfo()), con = con)
  cat("\n================================================================================\n", file = con, append = TRUE)
  cat("LOG END\n",                                                                          file = con, append = TRUE)
  cat("================================================================================\n", file = con, append = TRUE)

  close(con)

  log_file_path <- .talaria_log_env$log_file

  .talaria_log_env$active     <- FALSE
  .talaria_log_env$log_file   <- NULL
  .talaria_log_env$output_con <- NULL
  .talaria_log_env$start_time <- NULL

  invisible(log_file_path)
}


#' Check whether logging is active
#'
#' @return Logical — \code{TRUE} if a log is currently running.
#' @export
TALARIA_is_logging <- function() {
  .talaria_log_env$active
}


#' Get the current log file path
#'
#' @return Character path to the active log file, or \code{NULL} if no log
#'   is running.
#' @export
TALARIA_get_log_file <- function() {
  if (.talaria_log_env$active) .talaria_log_env$log_file else NULL
}
