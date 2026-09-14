# Resolve the GATK executable path

Resolves the GATK executable in the following order:

1.  Explicit `gatk_path` argument

2.  `GATK_PATH` environment variable

3.  `gatk` on the system `PATH` (via `Sys.which`)

Stops with an informative error if GATK cannot be found.

## Usage

``` r
CYAN_resolve_gatk(gatk_path = NULL)
```

## Arguments

- gatk_path:

  Character or `NULL`. Explicit path to the GATK executable (e.g.
  `"/opt/gatk-4.4/gatk"`). If `NULL` the function falls back to the
  environment variable and then `PATH`.

## Value

Character. Resolved path to the GATK executable.
