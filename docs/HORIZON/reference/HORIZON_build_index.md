# Build a Rsubread genome index

Wrapper around
[`buildindex`](https://rdrr.io/pkg/Rsubread/man/buildindex.html) for
constructing the genome index required for alignment. This is a one-time
operation per reference genome; the index can be stored on any
accessible drive (including external storage).

## Usage

``` r
HORIZON_build_index(
  basename,
  reference,
  memory = 8000,
  index_split = TRUE,
  gapped = TRUE,
  force = FALSE,
  ...
)
```

## Arguments

- basename:

  Character. Path prefix for the index files (e.g.,
  `"/data/refs/hg38/hg38_subread"`). The directory must already exist.

- reference:

  Character. Path to the reference genome FASTA file.

- memory:

  Integer. RAM to use for index building in MB. Default 8000 (~8 GB).
  Reduce on memory-constrained machines.

- index_split:

  Logical. Split the index into chunks to reduce RAM usage during
  alignment. Recommended for large genomes. Default TRUE.

- gapped:

  Logical. Build a gapped (splice-aware) index suitable for RNA-seq
  alignment. Default TRUE.

- force:

  Logical. If TRUE, rebuild the index even if one already exists at
  `basename`. Default FALSE.

- ...:

  Additional arguments passed to
  [`buildindex`](https://rdrr.io/pkg/Rsubread/man/buildindex.html).

## Value

Character. The index basename path, invisibly.

## Details

RAM note: index building is the most RAM-intensive step. The `memory`
parameter directly controls how much RAM Rsubread uses during
construction. With `index_split = TRUE` (default), the resulting index
is stored in chunks, which significantly reduces per-alignment RAM usage
at the cost of slightly slower alignment.
