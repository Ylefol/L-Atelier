# Check if Logging is Active

Helper function to check if console logging is currently active.

## Usage

``` r
ELEUTHIA_is_logging()
```

## Value

Logical. TRUE if logging is active, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_is_logging()
# FALSE

ELEUTHIA_start_log("test.log")
ELEUTHIA_is_logging()
# TRUE

ELEUTHIA_stop_log()
ELEUTHIA_is_logging()
# FALSE

} # }
```
