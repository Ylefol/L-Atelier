#' HORIZON - Logging Functions
#'
#' @description Functions for managing console output logging to files.
#' Duplicates stdout to both screen and log file (does not redirect).
#' Modelled on GAIA's ELEUTHIA logging.

# Package environment to store log state
.log_env <- new.env(parent = emptyenv())
.log_env$active     <- FALSE
.log_env$log_file   <- NULL
.log_env$output_con <- NULL
.log_env$start_time <- NULL


#' Generate Timestamp-Based Log Filename
#'
#' @description Creates a log filename based on current date and time.
#'
#' @param log_dir Character. Directory for the log file. NULL returns filename only.
#' @param prefix Character. Optional filename prefix (default = "log").
#' @param extension Character. File extension (default = ".log").
#' @param format Character. "compact" = YYYY_MM_DD__HH_MM (default),
#'   "readable" = YYYY-MM-DD_HH-MM-SS.
#'
#' @return Character string with the full log file path or filename.
#'
#' @export
HORIZON_generate_log_filename <- function(log_dir   = NULL,
                                          prefix    = "log",
                                          extension = ".log",
                                          format    = "compact") {

  timestamp <- Sys.time()

  if (format == "compact") {
    date_str <- format(timestamp, "%Y_%m_%d")
    time_str <- format(timestamp, "%H_%M")
    time_str <- paste0(date_str, "__", time_str)
  } else if (format == "readable") {
    time_str <- format(timestamp, "%Y-%m-%d_%H-%M-%S")
  } else {
    stop("format must be 'compact' or 'readable'")
  }

  filename <- paste0(prefix, "_", time_str, extension)

  if (!is.null(log_dir)) {
    filename <- file.path(log_dir, filename)
  }

  filename
}


#' Start Console Logging
#'
#' @description Begins duplicating stdout to a log file while still displaying
#' output in the console. Standard output (cat, print) appears in both
#' locations; messages and warnings appear in console only.
#'
#' @param log_file Character. Path to the log file to create/append to.
#' @param append Logical. Append to existing file (default = FALSE).
#'
#' @return Invisibly returns the log file path.
#'
#' @details
#' Records: start timestamp, R version, working directory, and all subsequent
#' standard output. Must be paired with HORIZON_stop_log().
#'
#' Note: \code{message()} output (stderr) goes to the console only and is not
#' captured in the log file.  Explicit \code{cat()} step banners are captured.
#'
#' @export
HORIZON_start_log <- function(log_file, append = FALSE) {

  if (.log_env$active) {
    warning("Logging already active. Stop current log before starting a new one.")
    return(invisible(NULL))
  }

  if (!is.character(log_file) || length(log_file) != 1) {
    stop("log_file must be a single character string")
  }

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  .log_env$start_time <- Sys.time()
  .log_env$log_file   <- log_file
  .log_env$output_con <- file(log_file, open = if (append) "a" else "w")

  cat("================================================================================\n",
      file = .log_env$output_con)
  cat("LOG START:", format(.log_env$start_time, "%Y-%m-%d %H:%M:%S"), "\n",
      file = .log_env$output_con)
  cat("================================================================================\n",
      file = .log_env$output_con)
  cat("R version:", R.version.string, "\n",
      file = .log_env$output_con)
  cat("Platform:", R.version$platform, "\n",
      file = .log_env$output_con)
  cat("Working directory:", getwd(), "\n",
      file = .log_env$output_con)
  cat("================================================================================\n\n",
      file = .log_env$output_con)

  sink(file = .log_env$output_con, append = TRUE, split = TRUE)

  .log_env$active <- TRUE
  cat(sprintf("Logging started: %s\n", log_file))

  invisible(log_file)
}


#' Stop Console Logging
#'
#' @description Stops stdout duplication and appends session information.
#'
#' @return Invisibly returns the log file path, or NULL if no log was active.
#'
#' @export
HORIZON_stop_log <- function() {

  if (!.log_env$active) {
    warning("No active log to stop.")
    return(invisible(NULL))
  }

  cat("\nLogging stopped.\n")
  sink()

  end_time <- Sys.time()
  duration <- difftime(end_time, .log_env$start_time, units = "auto")

  cat("\n================================================================================\n",
      file = .log_env$output_con, append = TRUE)
  cat("SESSION INFORMATION\n",
      file = .log_env$output_con, append = TRUE)
  cat("================================================================================\n",
      file = .log_env$output_con, append = TRUE)
  cat("Log end:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n",
      file = .log_env$output_con, append = TRUE)
  cat("Duration:", format(duration), "\n\n",
      file = .log_env$output_con, append = TRUE)

  cat("R Session Info:\n",
      file = .log_env$output_con, append = TRUE)
  cat("--------------------------------------------------------------------------------\n",
      file = .log_env$output_con, append = TRUE)
  writeLines(capture.output(sessionInfo()), con = .log_env$output_con)

  cat("\n================================================================================\n",
      file = .log_env$output_con, append = TRUE)
  cat("LOG END\n",
      file = .log_env$output_con, append = TRUE)
  cat("================================================================================\n",
      file = .log_env$output_con, append = TRUE)

  close(.log_env$output_con)

  log_file_path <- .log_env$log_file

  .log_env$active     <- FALSE
  .log_env$log_file   <- NULL
  .log_env$output_con <- NULL
  .log_env$start_time <- NULL

  invisible(log_file_path)
}


#' Check if Logging is Active
#'
#' @return Logical.
#' @export
HORIZON_is_logging <- function() {
  .log_env$active
}


#' Get Current Log File Path
#'
#' @return Character log file path, or NULL.
#' @export
HORIZON_get_log_file <- function() {
  if (.log_env$active) .log_env$log_file else NULL
}
