#' Eleuthia - Logging Functions
#'
#' @description Functions for managing console output logging to files.
#' Duplicates console output to both screen and log file (does not redirect).


# Package environment to store log state
.log_env <- new.env(parent = emptyenv())
.log_env$active <- FALSE
.log_env$log_file <- NULL
.log_env$output_con <- NULL
.log_env$message_con <- NULL
.log_env$start_time <- NULL


#' Generate Timestamp-Based Log Filename
#'
#' @description Creates a log filename based on current date and time for
#' automatic log file naming and tracking.
#'
#' @param log_dir Character string. Directory path where log will be saved.
#'   If NULL, returns just the filename (default = NULL).
#' @param prefix Character string. Optional prefix for the log filename
#'   (default = "log").
#' @param extension Character string. File extension (default = ".log").
#' @param format Character string. Timestamp format style:
#'   "compact" = YYYY_MM_DD__HH_MM (default, e.g., log_2026_01_14__14_30.log),
#'   "readable" = YYYY-MM-DD_HH-MM-SS (e.g., log_2026-01-14_14-30-52.log)
#'
#' @return Character string with full log file path or just filename.
#'
#' @details
#' Generates filenames that:
#' \itemize{
#'   \item Sort chronologically when listed
#'   \item Are filesystem-safe (no problematic characters)
#'   \item Include full timestamp for uniqueness
#'   \item Indicate when the log was created
#' }
#'
#' @export
#'
#' @examples
#' # Just the filename
#' ELEUTHIA_generate_log_filename()
#' # "log_2026_01_14__14_30.log"
#'
#' # With directory path
#' ELEUTHIA_generate_log_filename(log_dir = "results/logs")
#' # "results/logs/log_2026_01_14__14_30.log"
#'
#' # With custom prefix
#' ELEUTHIA_generate_log_filename(prefix = "scca_cv", log_dir = "output")
#' # "output/scca_cv_2026_01_14__14_30.log"
#'
#' # Readable format (with seconds)
#' ELEUTHIA_generate_log_filename(format = "readable")
#' # "log_2026-01-14_14-30-52.log"
#'
#' # Use with start_log
#' log_file <- ELEUTHIA_generate_log_filename(log_dir = "results")
#' ELEUTHIA_start_log(log_file)
#'
ELEUTHIA_generate_log_filename <- function(log_dir = NULL,
                                            prefix = "log",
                                            extension = ".log",
                                            format = "compact") {

  # Get current timestamp
  timestamp <- Sys.time()

  # Format timestamp based on style
  if (format == "compact") {
    # YYYY_MM_DD__HH_MM format (underscores, no seconds)
    date_str <- format(timestamp, "%Y_%m_%d")
    time_str <- format(timestamp, "%H_%M")
    time_str <- paste0(date_str, "__", time_str)
  } else if (format == "readable") {
    # YYYY-MM-DD_HH-MM-SS format (more readable, with seconds)
    time_str <- format(timestamp, "%Y-%m-%d_%H-%M-%S")
  } else {
    stop("format must be either 'compact' or 'readable'")
  }

  # Construct filename
  filename <- paste0(prefix, "_", time_str, extension)

  # Add directory if provided
  if (!is.null(log_dir)) {
    filename <- file.path(log_dir, filename)
  }

  return(filename)
}


#' Start Console Logging
#'
#' @description Begins duplicating console output (stdout) to a log file while
#' still displaying in the console. Standard output appears in both locations,
#' while warnings and errors appear only in the console for immediate visibility.
#'
#' @param log_file Character string. Path to the log file to create/append to.
#' @param append Logical. If TRUE, append to existing log file. If FALSE,
#'   overwrite existing file (default = FALSE).
#'
#' @return Invisibly returns the log file path.
#'
#' @details
#' This function uses sink() with split=TRUE to duplicate stdout rather than
#' redirect it. Console output continues to appear normally while also being
#' written to the log file.
#'
#' The function records:
#' \itemize{
#'   \item Start timestamp
#'   \item R version and platform
#'   \item Working directory
#'   \item All subsequent standard output (print, cat, etc.)
#' }
#'
#' Note: Warnings and errors appear only in the console (not in log file) for
#' immediate visibility during analysis.
#'
#' Must be paired with ELEUTHIA_stop_log() to properly close the log and
#' append session information.
#'
#' @export
#'
#' @examples
#' # Start logging
#' ELEUTHIA_start_log("analysis.log")
#'
#' # Your analysis code here
#' print("This appears in both console and log file")
#'
#' # Stop logging and append session info
#' ELEUTHIA_stop_log()
#'
ELEUTHIA_start_log <- function(log_file, append = FALSE) {

  # Check if logging is already active
  if (.log_env$active) {
    warning("Logging is already active. Stop current log before starting a new one.")
    return(invisible(NULL))
  }

  # Validate log_file path
  if (!is.character(log_file) || length(log_file) != 1) {
    stop("log_file must be a single character string")
  }

  # Create directory if it doesn't exist
  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  # Record start time
  .log_env$start_time <- Sys.time()
  .log_env$log_file <- log_file

  # Open log file connection
  .log_env$output_con <- file(log_file, open = if (append) "a" else "w")

  # Write header to log file
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

  # Start duplicating stdout to log file
  # Warnings and errors remain in console only for immediate visibility
  sink(file = .log_env$output_con, append = TRUE, split = TRUE)

  # Mark logging as active
  .log_env$active <- TRUE

  # Print confirmation to console (which will also go to log)
  cat(sprintf("Logging started: %s\n", log_file))

  invisible(log_file)
}


#' Stop Console Logging
#'
#' @description Stops duplicating console output to the log file and appends
#' session information before closing.
#'
#' @return Invisibly returns the log file path, or NULL if no log was active.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Stops stdout duplication (closes sink)
#'   \item Appends session information (packages, versions, etc.)
#'   \item Records end timestamp and total duration
#'   \item Closes file connections cleanly
#' }
#'
#' Session information includes:
#' \itemize{
#'   \item End timestamp and total log duration
#'   \item Loaded packages and versions
#'   \item Locale settings
#'   \item Full sessionInfo() output
#' }
#'
#' Note: Only standard output is logged. Warnings and errors appear only
#' in the console for immediate visibility.
#'
#' @export
#'
#' @examples
#' # Start logging
#' ELEUTHIA_start_log("analysis.log")
#'
#' # Your analysis code here
#' print("Analysis output")
#'
#' # Stop logging (appends session info)
#' ELEUTHIA_stop_log()
#'
ELEUTHIA_stop_log <- function() {

  # Check if logging is active
  if (!.log_env$active) {
    warning("No active log to stop.")
    return(invisible(NULL))
  }

  # Print completion message (will appear in both console and log)
  cat("\n")
  cat("Logging stopped.\n")

  # Close stdout sink
  sink()

  # Record end time
  end_time <- Sys.time()
  duration <- difftime(end_time, .log_env$start_time, units = "auto")

  # Append session information directly to log file
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

  # Capture session info to log
  cat("R Session Info:\n",
      file = .log_env$output_con, append = TRUE)
  cat("--------------------------------------------------------------------------------\n",
      file = .log_env$output_con, append = TRUE)

  # Write sessionInfo() output to log
  session_output <- capture.output(sessionInfo())
  writeLines(session_output, con = .log_env$output_con)

  cat("\n================================================================================\n",
      file = .log_env$output_con, append = TRUE)
  cat("LOG END\n",
      file = .log_env$output_con, append = TRUE)
  cat("================================================================================\n",
      file = .log_env$output_con, append = TRUE)

  # Close file connection
  close(.log_env$output_con)

  # Store log file path before resetting
  log_file_path <- .log_env$log_file

  # Reset log environment
  .log_env$active <- FALSE
  .log_env$log_file <- NULL
  .log_env$output_con <- NULL
  .log_env$message_con <- NULL
  .log_env$start_time <- NULL

  invisible(log_file_path)
}


#' Check if Logging is Active
#'
#' @description Helper function to check if console logging is currently active.
#'
#' @return Logical. TRUE if logging is active, FALSE otherwise.
#'
#' @export
#'
#' @examples
#' ELEUTHIA_is_logging()
#' # FALSE
#'
#' ELEUTHIA_start_log("test.log")
#' ELEUTHIA_is_logging()
#' # TRUE
#'
#' ELEUTHIA_stop_log()
#' ELEUTHIA_is_logging()
#' # FALSE
#'
ELEUTHIA_is_logging <- function() {
  return(.log_env$active)
}


#' Get Current Log File Path
#'
#' @description Returns the path to the currently active log file, or NULL
#' if no log is active.
#'
#' @return Character string with log file path, or NULL if logging is inactive.
#'
#' @export
#'
#' @examples
#' ELEUTHIA_get_log_file()
#' # NULL
#'
#' ELEUTHIA_start_log("my_analysis.log")
#' ELEUTHIA_get_log_file()
#' # "my_analysis.log"
#'
#' ELEUTHIA_stop_log()
#'
ELEUTHIA_get_log_file <- function() {
  if (.log_env$active) {
    return(.log_env$log_file)
  } else {
    return(NULL)
  }
}
