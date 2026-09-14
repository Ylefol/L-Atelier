# Execute REDItools2 via basilisk + system2

Activates the managed Python environment, retrieves the Python binary,
and calls the REDItools2 script as a subprocess. Python stdout/stderr
are suppressed unless `debug = TRUE`.

## Usage

``` r
.cyan_reditools_run(script_path, cli_args, verbose, debug)
```
