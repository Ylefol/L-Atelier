# Run REDItools2 on multiple samples

Convenience wrapper around
[`CYAN_run_reditools`](https://ylefol.github.io/L-Atelier/CYAN/reference/CYAN_run_reditools.md)
that iterates over a vector of BAM files. If the vector is named, each
name is forwarded as the `sample_name` argument, controlling output file
prefixes. If unnamed, sample names are derived from the BAM filenames as
usual.

## Usage

``` r
CYAN_run_reditools_batch(bam, reference, output_dir, ..., verbose = TRUE)
```

## Arguments

- bam:

  Named or unnamed character vector of BAM file paths. Names, if
  present, are used as `sample_name` for each sample.

- reference:

  Character. Path to the reference FASTA (shared across all samples).

- output_dir:

  Character. Directory for all result files (created if absent). Each
  sample writes its own prefixed files here.

- ...:

  Additional arguments passed to `CYAN_run_reditools` (e.g.
  `strand_mode`, `min_coverage`, `keep_tmp`).

- verbose:

  Logical. Print `[CYAN]` progress messages including a per-sample
  summary on completion (default `TRUE`).

## Value

An object of class `cyan_batch_result`: a named list with one element
per BAM file. Each element is either a `cyan_editing_result` (success)
or `NULL` (failure). Failed samples are identified in the printed
summary.

## Details

Samples are processed sequentially. A failed sample issues a warning and
is recorded as `NULL` in the returned list; remaining samples continue.
Because `CYAN_run_reditools` checkpoints its raw output, re-running
after a partial failure skips already-completed samples automatically.

## See also

[`CYAN_run_reditools`](https://ylefol.github.io/L-Atelier/CYAN/reference/CYAN_run_reditools.md)
