# Load a BED File

Reads a BED file into a data.frame with standardized column names.

## Usage

``` r
.eleuthia_load_seacr_peaks(file_path)
```

## Arguments

- file_path:

  Character string. Path to the BED file.

## Value

A data.frame with at minimum columns: chr, start, end. Additional
columns depend on BED format (BED3, BED6, etc.).
