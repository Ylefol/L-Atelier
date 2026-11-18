# Tests for the console logging helpers. Logging uses sink(), so every test
# makes sure the log is stopped again, even if an expectation fails.

stop_if_logging <- function() {
  if (HORIZON_is_logging()) suppressWarnings(HORIZON_stop_log())
}

test_that("generate_log_filename builds the expected names", {
  f <- HORIZON_generate_log_filename(prefix = "run", extension = ".txt")
  expect_match(f, "^run_\\d{4}_\\d{2}_\\d{2}__\\d{2}_\\d{2}\\.txt$")

  f2 <- HORIZON_generate_log_filename(log_dir = "logs", format = "readable")
  expect_match(f2, "^logs/log_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}\\.log$")

  expect_error(HORIZON_generate_log_filename(format = "nonsense"),
               "compact.*readable")
})

test_that("logging is off by default and reports no log file", {
  stop_if_logging()
  expect_false(HORIZON_is_logging())
  expect_null(HORIZON_get_log_file())
})

test_that("start_log / stop_log capture output and record session info", {
  stop_if_logging()
  on.exit(stop_if_logging(), add = TRUE)

  log_file <- file.path(tempfile("logdir_"), "test.log")
  on.exit(unlink(dirname(log_file), recursive = TRUE), add = TRUE)

  HORIZON_start_log(log_file)
  expect_true(HORIZON_is_logging())
  expect_equal(HORIZON_get_log_file(), log_file)
  cat("hello from the test\n")
  HORIZON_stop_log()

  expect_false(HORIZON_is_logging())
  expect_true(file.exists(log_file))
  txt <- readLines(log_file)
  expect_true(any(grepl("LOG START", txt)))
  expect_true(any(grepl("hello from the test", txt)))
  # Session info is what supports reproducibility, so check it is recorded
  expect_true(any(grepl("R Session Info", txt)))
  expect_true(any(grepl("LOG END", txt)))
})

test_that("start_log validates its input and warns when already active", {
  stop_if_logging()
  on.exit(stop_if_logging(), add = TRUE)

  expect_error(HORIZON_start_log(c("a.log", "b.log")), "single character string")

  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)
  HORIZON_start_log(log_file)
  expect_warning(HORIZON_start_log(log_file), "already active")
  HORIZON_stop_log()
})

test_that("stop_log warns when there is no active log", {
  stop_if_logging()
  expect_warning(HORIZON_stop_log(), "No active log")
})
