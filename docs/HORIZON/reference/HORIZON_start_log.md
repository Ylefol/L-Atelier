# Start Console Logging

Begins duplicating stdout to a log file while still displaying output in
the console. Standard output (cat, print) appears in both locations;
messages and warnings appear in console only.

## Usage

``` r
HORIZON_start_log(log_file, append = FALSE)
```

## Arguments

- log_file:

  Character. Path to the log file to create/append to.

- append:

  Logical. Append to existing file (default = FALSE).

## Value

Invisibly returns the log file path.

## Details

Records: start timestamp, R version, working directory, and all
subsequent standard output. Must be paired with HORIZON_stop_log().

Note: [`message()`](https://rdrr.io/r/base/message.html) output (stderr)
goes to the console only and is not captured in the log file. Explicit
[`cat()`](https://rdrr.io/r/base/cat.html) step banners are captured.
