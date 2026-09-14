# Start Console Logging

Begins duplicating console output (stdout) to a log file while still
displaying in the console. Standard output appears in both locations,
while warnings and errors appear only in the console for immediate
visibility.

## Usage

``` r
ELEUTHIA_start_log(log_file, append = FALSE)
```

## Arguments

- log_file:

  Character string. Path to the log file to create/append to.

- append:

  Logical. If TRUE, append to existing log file. If FALSE, overwrite
  existing file (default = FALSE).

## Value

Invisibly returns the log file path.

## Details

This function uses sink() with split=TRUE to duplicate stdout rather
than redirect it. Console output continues to appear normally while also
being written to the log file.

The function records:

- Start timestamp

- R version and platform

- Working directory

- All subsequent standard output (print, cat, etc.)

Note: Warnings and errors appear only in the console (not in log file)
for immediate visibility during analysis.

Must be paired with ELEUTHIA_stop_log() to properly close the log and
append session information.

## Examples

``` r
if (FALSE) { # \dontrun{
# Start logging
ELEUTHIA_start_log("analysis.log")

# Your analysis code here
print("This appears in both console and log file")

# Stop logging and append session info
ELEUTHIA_stop_log()

} # }
```
