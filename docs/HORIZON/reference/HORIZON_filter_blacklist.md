# Remove reads overlapping a genomic blacklist

Filters a processed BAM file using `bedtools intersect -v` to exclude
reads that overlap a blacklist BED file (e.g. ENCODE blacklist regions).
For assemblies without an ENCODE blacklist (such as T2T-CHM13), pass
`blacklist_bed = NULL` to skip this step.

## Usage

``` r
HORIZON_filter_blacklist(
  sample_sheet,
  sample_id,
  blacklist_bed,
  threads = 4L,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame.

- sample_id:

  Character. Sample ID to process.

- blacklist_bed:

  Character. Path to blacklist BED file. Pass `NULL` to skip filtering
  (a warning is issued).

- threads:

  Integer. Threads for re-indexing the filtered BAM. Default 4.

- force:

  Logical. If `FALSE` (default), skip filtering when
  `_blacklist_filtered.bam` already exists. Set `TRUE` to re-filter and
  overwrite.

## Value

Character. Path to the blacklist-filtered BAM file (or the unfiltered
processed BAM if `blacklist_bed = NULL`), invisibly.
