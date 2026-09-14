# Generate Timestamp-Based Log Filename

Creates a log filename based on current date and time for automatic log
file naming and tracking.

## Usage

``` r
ELEUTHIA_generate_log_filename(
  log_dir = NULL,
  prefix = "log",
  extension = ".log",
  format = "compact"
)
```

## Arguments

- log_dir:

  Character string. Directory path where log will be saved. If NULL,
  returns just the filename (default = NULL).

- prefix:

  Character string. Optional prefix for the log filename (default =
  "log").

- extension:

  Character string. File extension (default = ".log").

- format:

  Character string. Timestamp format style: "compact" =
  YYYY_MM_DD\_\_HH_MM (default, e.g., log_2026_01_14\_\_14_30.log),
  "readable" = YYYY-MM-DD_HH-MM-SS (e.g., log_2026-01-14_14-30-52.log)

## Value

Character string with full log file path or just filename.

## Details

Generates filenames that:

- Sort chronologically when listed

- Are filesystem-safe (no problematic characters)

- Include full timestamp for uniqueness

- Indicate when the log was created

## Examples

``` r
if (FALSE) { # \dontrun{
# Just the filename
ELEUTHIA_generate_log_filename()
# "log_2026_01_14__14_30.log"

# With directory path
ELEUTHIA_generate_log_filename(log_dir = "results/logs")
# "results/logs/log_2026_01_14__14_30.log"

# With custom prefix
ELEUTHIA_generate_log_filename(prefix = "scca_cv", log_dir = "output")
# "output/scca_cv_2026_01_14__14_30.log"

# Readable format (with seconds)
ELEUTHIA_generate_log_filename(format = "readable")
# "log_2026-01-14_14-30-52.log"

# Use with start_log
log_file <- ELEUTHIA_generate_log_filename(log_dir = "results")
ELEUTHIA_start_log(log_file)

} # }
```
