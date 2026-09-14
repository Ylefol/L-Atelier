# Start console logging

Begins duplicating standard output (`stdout`) to a log file using
`sink(split = TRUE)`, so output appears in both the console and the file
simultaneously. Warnings and errors are intentionally not captured —
they remain console-only for immediate visibility.

## Usage

``` r
TALARIA_start_log(log_file, append = FALSE)
```

## Arguments

- log_file:

  Character. Path to the log file to create (or append to). Parent
  directories are created automatically if they do not exist.

- append:

  Logical. Append to an existing file (`TRUE`) or overwrite (`FALSE`,
  default).

## Value

Invisibly returns `log_file`.

## Details

The log file opens with a header recording the start timestamp, R
version, platform, and working directory.

Must be paired with
[`TALARIA_stop_log`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_stop_log.md)
to close the connection cleanly and append session information.
