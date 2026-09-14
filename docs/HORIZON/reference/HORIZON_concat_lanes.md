# Concatenate multi-lane FASTQ files into one file per read direction

Some aligners/tools accept multiple FASTQ files per sample natively
(e.g.
[`HORIZON_run_salmon`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_salmon.md)'s
`fastq_r1`/`fastq_r2` vectors), but most of the pipeline
([`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md),
[`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md))
expects exactly one FASTQ per sample per read direction. When a sample
was sequenced across multiple flowcells/lanes, this function merges
those lane-level FASTQs into a single R1 (and R2, for paired-end data)
file per sample, so the rest of the pipeline can treat it like any
single-lane sample.

## Usage

``` r
HORIZON_concat_lanes(
  fastq_r1,
  fastq_r2 = NULL,
  output_dir,
  sample_id,
  force = FALSE
)
```

## Arguments

- fastq_r1:

  Character vector. Paths to R1 FASTQ files across lanes, for one
  sample.

- fastq_r2:

  Character vector or `NULL`. Paths to R2 FASTQ files across lanes,
  positionally matched to `fastq_r1` (same length, same lane order —
  `fastq_r1[i]` and `fastq_r2[i]` must be the same lane). `NULL` for
  single-end data.

- output_dir:

  Character. Directory to write the concatenated FASTQ(s) into. Created
  if it doesn't exist.

- sample_id:

  Character. Sample ID; output files are named `<sample_id>_R1.fastq.gz`
  / `<sample_id>_R2.fastq.gz`.

- force:

  Logical. If `FALSE` (default), skip concatenation when the output
  file(s) already exist and return their paths directly. Set to `TRUE`
  to re-run and overwrite.

## Value

Named list with elements `fastq_r1` and `fastq_r2` (the latter `NULL`
for single-end data) giving the paths to the concatenated FASTQ(s),
invisibly.

## Details

Concatenation is done at the byte level (raw `.fastq.gz` bytes appended
in sequence), never decompressing/recompressing. A gzip file is a series
of independent "members" (RFC 1952); concatenating gzip files this way
produces another valid gzip file that decompresses to the full,
correctly-ordered read set — this is the standard way multi-lane FASTQs
are merged (equivalent to
`cat lane1.fastq.gz lane2.fastq.gz > merged.fastq.gz`), and every
downstream reader here (Rfastp, Rsubread, zcat) handles multi-member
gzip transparently.
