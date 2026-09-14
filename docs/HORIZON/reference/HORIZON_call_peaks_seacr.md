# Call peaks with SEACR (pure-R reimplementation)

A native R reimplementation of the SEACR (Sparse Enrichment Analysis for
CUT\\RUN) peak-calling method published by Meers, Tenenbaum and Bhatt
(2019,
[doi:10.1186/s13072-019-0287-4](https://doi.org/10.1186/s13072-019-0287-4)
). **This function reimplements their published algorithm; it is not an
independent method.** Results should be validated against the original
SEACR tool (<https://github.com/FredHutch/SEACR>) before use in
publications.

## Usage

``` r
HORIZON_call_peaks_seacr(
  treatment_bed,
  control_bed = NULL,
  output_dir,
  sample_name,
  chrom_sizes,
  threshold = 0.05,
  normalize = c("non", "norm"),
  stringency = c("stringent", "relaxed"),
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- treatment_bed:

  Character. Path to the treatment fragment BED file (from
  [`HORIZON_bam_to_bed`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bed.md)).

- control_bed:

  Character or `NULL`. Path to an IgG / input fragment BED file. When
  `NULL` (default), threshold-based calling is used (`threshold`).

- output_dir:

  Character. Directory for the output peak file.

- sample_name:

  Character. Prefix applied to the output filename.

- chrom_sizes:

  A `data.frame` with columns `chr` and `size` (from
  [`HORIZON_get_chrom_sizes`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_get_chrom_sizes.md)).

- threshold:

  Numeric between 0 and 1. Fraction of blocks to retain by total AUC
  when no control is provided. Default `0.05` (top 5\\ Ignored when
  `control_bed` is supplied.

- normalize:

  Character. `"non"` (default): use raw fragment counts. `"norm"`: scale
  coverage by total fragments (per million) before computing AUC,
  matching SEACR's `norm` mode.

- stringency:

  Character. `"stringent"` (default): a block must exceed the threshold
  on both total AUC and max signal. `"relaxed"`: a block passes if
  either criterion is met.

- force:

  Logical. Re-run even if the output file already exists. Default
  `FALSE`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

Character. Path to the output BED file
(`<sample_name>_seacr_peaks.bed`), invisibly. Columns: `peak_id`, `chr`,
`start`, `end`, `total_signal`, `max_signal`.

## Details

**Algorithm summary (Meers et al. 2019):**

1.  Build per-base fragment coverage from the treatment BED file.

2.  Identify contiguous blocks of non-zero signal.

3.  Compute the total signal (AUC) and maximum signal height per block.

4.  If a control BED is provided: subtract the per-block control AUC
    from the treatment AUC; retain blocks with positive net signal.

5.  If no control is provided: rank blocks by AUC and retain the top
    `threshold` fraction.

6.  Apply stringency: `"stringent"` requires a block to exceed the
    threshold on both total AUC and max signal; `"relaxed"` requires
    either criterion.

The interface mirrors
[`HORIZON_call_peaks`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_call_peaks.md)
(MACS3) so the two functions can be called with the same core arguments
in a comparison workflow.
