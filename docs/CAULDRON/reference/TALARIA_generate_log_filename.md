# Generate a timestamp-based log filename

Creates a log filename using the current date and time, suitable for
automatic naming so log files sort chronologically on disk.

## Usage

``` r
TALARIA_generate_log_filename(
  log_dir = NULL,
  prefix = "log",
  extension = ".log",
  format = c("compact", "readable")
)
```

## Arguments

- log_dir:

  Character or `NULL`. Directory in which to place the file. `NULL`
  (default) returns just the filename with no path.

- prefix:

  Character. Filename prefix. Default `"log"`.

- extension:

  Character. File extension. Default `".log"`.

- format:

  Character. Timestamp style:

  `"compact"`

  :   Default. `YYYY_MM_DD__HH_MM` — filesystem-safe, no seconds (e.g.
      `log_2026_03_05__14_30.log`).

  `"readable"`

  :   `YYYY-MM-DD_HH-MM-SS` — includes seconds (e.g.
      `log_2026-03-05_14-30-52.log`).

## Value

Character string — full file path if `log_dir` is provided, otherwise
just the filename.
