# Convert Fragment-Level BED to BigWig

Converts a fragment-level BED file (one row per fragment) into a BigWig
coverage file suitable for genome browser visualization or use with
[`AETHER_plot_coverage_tracks()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_coverage_tracks.md).
Coverage is computed at a configurable bin resolution and optionally
normalized to CPM.

## Usage

``` r
ELEUTHIA_bed_to_bigwig(
  bed_file,
  output_path,
  chrom_sizes,
  normalize = "CPM",
  bin_size = 10L,
  chromosomes = NULL,
  verbose = TRUE
)
```

## Arguments

- bed_file:

  Character. Path to the fragment-level BED file (0-based half-open
  coordinates: chr, start, end).

- output_path:

  Character. Output file path. Extension should be `.bw` or `.bigwig`;
  format is inferred by `rtracklayer` from the extension.

- chrom_sizes:

  Exact chromosome lengths. Accepts either:

  - A named numeric vector: `c(chr1 = 248956422L, chr2 = 242193529L)`

  - A data.frame with columns `chr` and `size`, as returned by
    [`APOLLO_get_chromosome_sizes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_chromosome_sizes.md).
    Names/values must match chromosome names in the BED file exactly. No
    inferring from data.

- normalize:

  Character. Normalization method: `"CPM"` (default) scales signal by
  counts per million total fragments; `"raw"` leaves coverage as
  absolute fragment counts.

- bin_size:

  Integer. Bin width in base pairs for coverage aggregation. Larger
  values produce smoother signal and smaller files. Default `10`.

- chromosomes:

  Character vector. Restrict processing to these chromosomes. Default
  `NULL` processes all chromosomes present in both the BED file and
  `chrom_sizes`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

Invisibly returns `output_path`.

## Details

BED files use 0-based half-open coordinates. These are converted to
1-based closed intervals for GRanges internally.

Requires Bioconductor packages: `GenomicRanges`, `IRanges`,
`GenomeInfoDb`, `rtracklayer`. Install via
[`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).
