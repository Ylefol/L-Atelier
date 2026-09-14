# Build a Bowtie2 genome index

Wraps
[`Rbowtie2::bowtie2_build()`](https://rdrr.io/pkg/Rbowtie2/man/bowtie2-build.html)
to create a Bowtie2 index from one or more reference FASTA files. The
index only needs to be built once per reference genome; the resulting
`index_name` path is passed directly to
[`HORIZON_run_bowtie2`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md).

## Usage

``` r
HORIZON_build_bowtie2_index(
  reference,
  index_dir,
  index_name = NULL,
  threads = 4L,
  force = FALSE,
  ...
)
```

## Arguments

- reference:

  Character vector. Path(s) to reference FASTA file(s). Multiple FASFTAs
  are concatenated by bowtie2-build.

- index_dir:

  Character. Directory in which to write the index files. Created
  automatically if it does not exist.

- index_name:

  Character or `NULL`. Basename for the index files (without extension).
  Defaults to the filename of `reference[1]` without its extension.

- threads:

  Integer. Number of threads passed to `--threads`. Default 4.

- force:

  Logical. If `TRUE`, rebuild the index even when it already exists.
  Default `FALSE`.

- ...:

  Additional arguments passed to
  [`Rbowtie2::bowtie2_build`](https://rdrr.io/pkg/Rbowtie2/man/bowtie2-build.html).

## Value

Character. The full index basename path (i.e.
`<index_dir>/<index_name>`), invisibly. Pass this to `index` in
[`HORIZON_run_bowtie2`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md).
