# Generate Timestamp-Based Log Filename

Creates a log filename based on current date and time.

## Usage

``` r
HORIZON_generate_log_filename(
  log_dir = NULL,
  prefix = "log",
  extension = ".log",
  format = "compact"
)
```

## Arguments

- log_dir:

  Character. Directory for the log file. NULL returns filename only.

- prefix:

  Character. Optional filename prefix (default = "log").

- extension:

  Character. File extension (default = ".log").

- format:

  Character. "compact" = YYYY_MM_DD\_\_HH_MM (default), "readable" =
  YYYY-MM-DD_HH-MM-SS.

## Value

Character string with the full log file path or filename.
