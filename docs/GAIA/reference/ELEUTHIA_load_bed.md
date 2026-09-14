# Load a BED File

Reads a BED format file (BED3, BED6, narrowPeak, broadPeak, etc.) into a
data.frame with standardized column names. The first three columns are
always named `chr`, `start`, `end`; additional columns are named `name`,
`score`, `strand` where present.

## Usage

``` r
ELEUTHIA_load_bed(file_path)
```

## Arguments

- file_path:

  Character string. Path to the BED file.

## Value

A data.frame with at minimum columns: `chr`, `start`, `end`. Coordinates
are returned as-read (BED 0-based half-open intervals).
