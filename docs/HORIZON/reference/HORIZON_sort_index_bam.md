# Sort and index a BAM file

Coordinate-sorts the BAM produced by
[`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md)
and creates a BAI index file. Sorting is required before feature
counting.

## Usage

``` r
HORIZON_sort_index_bam(
  sample_sheet,
  sample_id,
  threads = 4,
  memory_per_thread = 2000,
  remove_unsorted = TRUE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- sample_id:

  Character. Sample ID to process.

- threads:

  Integer. Number of sort threads. Default 4.

- memory_per_thread:

  Integer. Memory per sort thread in MB. Total RAM used is approximately
  `memory_per_thread * threads` MB. Default 2000.

- remove_unsorted:

  Logical. Delete the unsorted BAM after sorting to free disk space.
  Default TRUE.

## Value

Character. Path to the sorted, indexed BAM file, invisibly.

## Details

By default, the unsorted BAM is deleted after sorting to avoid doubling
disk usage (set `remove_unsorted = FALSE` to keep it).

RAM note: total RAM used during sorting is approximately
`memory_per_thread * threads` MB. Reduce `threads` or
`memory_per_thread` on memory-constrained machines.
