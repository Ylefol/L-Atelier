# Validate reference files

Checks that a reference FASTA exists and has an accompanying FASTA index
(`.fai`). Optionally validates additional VCF or BED paths.

## Usage

``` r
CYAN_check_references(reference, ...)
```

## Arguments

- reference:

  Character. Path to the reference FASTA.

- ...:

  Additional file paths to check for existence (e.g. VCFs, BED files).
  Checked for existence only, not for indexing.

## Value

Invisibly returns `reference` (normalised absolute path).
