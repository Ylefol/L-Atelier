# Validate and load a sample sheet

Reads a HORIZON sample sheet CSV (or accepts a data.frame) and validates
it before use. Checks for required columns, unique sample IDs, valid
strandedness values, and optionally whether all FASTQ paths exist on
disk.

## Usage

``` r
HORIZON_validate_sample_sheet(path, check_files = TRUE)
```

## Arguments

- path:

  Character or data.frame. Path to a sample sheet CSV, or an
  already-loaded data.frame.

- check_files:

  Logical. Whether to verify that FASTQ file paths exist on disk.
  Default TRUE. Set to FALSE if files are on a drive not currently
  mounted.

## Value

The validated sample sheet as a data.frame.
