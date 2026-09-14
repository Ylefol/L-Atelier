# Align reads to the genome using Rsubread

Wraps [`align`](https://rdrr.io/pkg/Rsubread/man/align.html) for
splice-aware RNA-seq alignment against a pre-built Rsubread index. By
default, uses the trimmed FASTQ files produced by
[`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md).
Raw FASTQs from the sample sheet can be used instead by setting
`use_trimmed = FALSE`.

## Usage

``` r
HORIZON_run_align(
  sample_sheet,
  sample_id,
  index,
  use_trimmed = TRUE,
  threads = 4,
  max_mismatches = 3,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- sample_id:

  Character. Sample ID to process.

- index:

  Character. Path to the Rsubread index basename as passed to
  [`HORIZON_build_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_index.md).

- use_trimmed:

  Logical. Use trimmed FASTQs from the `qc/` subfolder (output of
  [`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md)).
  Default TRUE. Set to FALSE to align the raw FASTQs listed in the
  sample sheet.

- threads:

  Integer. Number of alignment threads. Default 4.

- max_mismatches:

  Integer. Maximum mismatches allowed per read. Default 3.

- ...:

  Additional arguments passed to
  [`align`](https://rdrr.io/pkg/Rsubread/man/align.html) (e.g.,
  `minFragLength`, `maxFragLength`, `unique`, `nBestLocations`).

## Value

Character. Path to the output BAM file, invisibly.

## Details

Output BAM is written to `<output_dir>/<sample_id>/aligned/`. The BAM is
unsorted at this stage; pass it to
[`HORIZON_sort_index_bam`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_sort_index_bam.md)
before counting.

RAM note: Rsubread's RAM usage during alignment is largely determined by
the index size. Using a split index (`index_split = TRUE` in
[`HORIZON_build_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_index.md))
substantially reduces peak RAM. Thread count has a modest additional
effect.
