# Get Current Log File Path

Returns the path to the currently active log file, or NULL if no log is
active.

## Usage

``` r
ELEUTHIA_get_log_file()
```

## Value

Character string with log file path, or NULL if logging is inactive.

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_get_log_file()
# NULL

ELEUTHIA_start_log("my_analysis.log")
ELEUTHIA_get_log_file()
# "my_analysis.log"

ELEUTHIA_stop_log()

} # }
```
