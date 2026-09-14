# Validate a PARSE Biosciences sample sheet

Reads a PARSE sample sheet CSV (or accepts a data.frame) and validates
it before use. Checks required columns, unique run IDs, valid chemistry
and kit values, and optionally whether all paths exist on disk.

## Usage

``` r
HORIZON_validate_parse_sheet(path, check_files = TRUE)
```

## Arguments

- path:

  Character or data.frame. Path to a CSV, or an already-loaded
  data.frame.

- check_files:

  Logical. Whether to verify that FASTQ files and genome directories
  exist on disk. Default TRUE. Set FALSE if paths are on a drive not
  currently mounted.

## Value

The validated sample sheet as a data.frame.
