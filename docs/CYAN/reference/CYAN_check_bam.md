# Validate a BAM file

Checks that a BAM file exists and has an accompanying index (`.bai`).
Stops with an informative message if either is missing.

## Usage

``` r
CYAN_check_bam(bam)
```

## Arguments

- bam:

  Character. Path to the BAM file.

## Value

Invisibly returns `bam` (normalised absolute path).
