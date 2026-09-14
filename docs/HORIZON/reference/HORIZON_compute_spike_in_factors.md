# Compute per-sample scale factors from spike-in read counts

Compute per-sample scale factors from spike-in read counts

## Usage

``` r
HORIZON_compute_spike_in_factors(
  spikein_bams,
  output_file = NULL,
  threads = 4L
)
```

## Arguments

- spikein_bams:

  Named character vector or list. Names are sample IDs; values are paths
  to spike-in BAMs from
  [`HORIZON_separate_spike_in`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_separate_spike_in.md).

- output_file:

  Character or `NULL`. Path to write a CSV of all factors (e.g.
  `output/spikein_factors.csv`). The file can be read back later to
  supply factors to DESeq2 without re-running this function. Default
  `NULL` (no file written).

- threads:

  Integer. Threads for `samtools view`. Default 4.

## Value

A data frame with columns `sample_id`, `spikein_reads`,
`scale_factor_bw`, `size_factor_deseq2`, invisibly. Samples with zero
spike-in reads receive `NA` factors with a warning.

## Details

Counts aligned read pairs in each spike-in BAM (produced by
[`HORIZON_separate_spike_in`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_separate_spike_in.md))
and returns a data frame containing two normalization factors per sample
derived from the same underlying spike-in read counts, for use at
different stages of the analysis:

- `scale_factor_bw`:

  For
  [`HORIZON_bam_to_bigwig`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bigwig.md).
  Formula: \\\min(N) / N_i\\. Values \\\leq 1\\; the least-sequenced
  sample gets 1.0 and all others are scaled down. Passed to
  `bamCoverage --scaleFactor`.

- `size_factor_deseq2`:

  For `DESeq2::sizeFactors()`. Formula: \\N_i / \min(N)\\. Values \\\geq
  1\\; DESeq2 divides raw counts by this value, so the least-sequenced
  sample (factor = 1.0) is unchanged and more deeply sequenced samples
  are scaled down. Assign via
  `sizeFactors(dds) <- factors$size_factor_deseq2` to bypass DESeq2's
  internal normalization and use spike-in calibration instead. This is
  conceptually equivalent to RNA-seq spike-in normalization but computed
  from whole-genome read counts rather than spike-in gene rows in a
  count matrix.
