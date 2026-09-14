# Compute Coverage for a Genomic Region Directly from a BED File

Extracts fragments overlapping a specific genomic region from a
fragment-level BED file and computes binned coverage — without
generating a full genome-wide BigWig file. Intended as the efficient
backend for
[`AETHER_plot_coverage_tracks()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_coverage_tracks.md)
when BigWig files are not pre-computed.

## Usage

``` r
ELEUTHIA_bed_region_coverage(
  bed_file,
  region,
  chrom_sizes,
  normalize = "CPM",
  bin_size = 10L,
  verbose = FALSE
)
```

## Arguments

- bed_file:

  Character. Path to the fragment-level BED file (0-based half-open
  coordinates: chr, start, end).

- region:

  Genomic region. One of:

  - Character string: `"chr1:1000000-2000000"`

  - Named character vector:
    `c(chr="chr1", start="1000000", end="2000000")`

  - A `GRanges` object (single range)

- chrom_sizes:

  Exact chromosome lengths. Accepts either:

  - A named numeric vector: `c(chr1 = 248956422L, chr2 = 242193529L)`

  - A data.frame with columns `chr` and `size`, as returned by
    [`APOLLO_get_chromosome_sizes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_chromosome_sizes.md).

- normalize:

  Character. `"CPM"` (default) scales by counts per million total
  fragments in the file; `"raw"` returns absolute counts.

- bin_size:

  Integer. Bin width in base pairs. Default `10`.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A data.frame with columns `start`, `end`, `score`, `midpoint` covering
the requested region. Same format as
[`rtracklayer::import.bw()`](https://rdrr.io/pkg/rtracklayer/man/BigWigFile.html)
output, making it interchangeable in downstream plotting code.

## Details

The entire BED file is read once (using `data.table`) to obtain both the
total fragment count (for CPM) and the region-filtered subset. Only bins
within the requested region are computed — no genome-wide tiling.

Requires: `data.table`, `GenomicRanges`, `IRanges`, `GenomeInfoDb`.
