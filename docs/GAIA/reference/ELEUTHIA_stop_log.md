# Stop Console Logging

Stops duplicating console output to the log file and appends session
information before closing.

## Usage

``` r
ELEUTHIA_stop_log()
```

## Value

Invisibly returns the log file path, or NULL if no log was active.

## Details

This function:

- Stops stdout duplication (closes sink)

- Appends session information (packages, versions, etc.)

- Records end timestamp and total duration

- Closes file connections cleanly

Session information includes:

- End timestamp and total log duration

- Loaded packages and versions

- Locale settings

- Full sessionInfo() output

Note: Only standard output is logged. Warnings and errors appear only in
the console for immediate visibility.

## Examples

``` r
if (FALSE) { # \dontrun{
# Start logging
ELEUTHIA_start_log("analysis.log")

# Your analysis code here
print("Analysis output")

# Stop logging (appends session info)
ELEUTHIA_stop_log()

} # }
```
