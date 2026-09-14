# Stop console logging

Stops duplicating output to the log file, appends a session information
block (end timestamp, duration,
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) output), and
closes the file connection cleanly.

## Usage

``` r
TALARIA_stop_log()
```

## Value

Invisibly returns the path to the closed log file, or `NULL` if no log
was active.
