# Compute insert size distribution with bamPEFragmentSize

Wraps `deeptools bamPEFragmentSize` to produce an insert-size histogram.
For ATAC-seq data, the nucleosomal periodicity (peaks near ~200 bp, ~400
bp, ~600 bp) is visible in this plot.

## Usage

``` r
HORIZON_run_fragment_size(
  sample_sheet,
  sample_id,
  threads = 4L,
  max_length = 1000L,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame.

- sample_id:

  Character. Sample ID to process.

- threads:

  Integer. Number of processors (`--numberOfProcessors`). Default 4.

- max_length:

  Integer. Maximum fragment length shown on the x-axis
  (`--maxFragmentLength`). Default 1000.

- force:

  Logical. If `FALSE` (default), skip when the output PNG already
  exists. Set `TRUE` to rerun.

## Value

Character. Path to the output PNG, invisibly.
