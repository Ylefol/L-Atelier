# Get chromosome sizes

Returns a two-column data.frame of chromosome names and sizes. Accepts
three input types so callers can use whichever file is most convenient:

## Usage

``` r
HORIZON_get_chrom_sizes(
  genome,
  chromosomes = NULL,
  cache_dir = NULL,
  cache_name = NULL,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- genome:

  Character. Path to a FASTA file, `.fai` index, or pre-made two-column
  chromosome sizes file.

- chromosomes:

  Character vector or `NULL`. If supplied, only rows whose `chr` value
  is in this vector are returned.

- cache_dir:

  Character or `NULL`. Directory for the cached RDS file. Defaults to a
  `data/` subdirectory relative to the current working directory.

- cache_name:

  Character or `NULL`. Stem of the cache filename (without extension).
  Defaults to the input filename stem with `_chrom_sizes` appended.

- force:

  Logical. If `TRUE`, regenerate the cache even when it already exists.
  Default `FALSE`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `data.frame` with columns `chr` (character) and `size` (integer), one
row per chromosome.

## Details

- FASTA (`.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`):

  Runs `samtools faidx` via the registered conda environment to generate
  a `.fai` index next to the FASTA, then reads the first two columns.
  Requires
  [`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md)
  to have been called first.

- FAI (`.fai`):

  Reads the samtools fasta index directly. Columns 1 and 2 are the
  chromosome name and length.

- Two-column text file:

  Any tab- or whitespace-separated file whose first column is the
  chromosome name and second column is the size. Typical extensions:
  `.sizes`, `.chrom.sizes`, `.tsv`, `.txt`.

The return value uses the same `data.frame(chr, size)` format as
`APOLLO_get_chromosome_sizes()` in the GAIA package so that the two
outputs are interchangeable.
