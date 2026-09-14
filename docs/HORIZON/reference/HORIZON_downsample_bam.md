# Downsample a BAM file by a spike-in derived scale factor

Randomly subsamples the most downstream BAM for a sample using
`samtools view -s`, producing a `<sample_id>_downsampled.bam` that is
picked up automatically by all downstream pipeline functions via
`.latest_bam()`.

## Usage

``` r
HORIZON_downsample_bam(
  sample_sheet,
  sample_id,
  scale_factor,
  seed = 42L,
  threads = 4L,
  remove_input = FALSE,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame.

- sample_id:

  Character. Sample ID to process.

- scale_factor:

  Numeric in `(0, 1]`. Downsampling fraction. Supply the
  `scale_factor_bw` value from
  [`HORIZON_compute_spike_in_factors`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md).

- seed:

  Integer. Random seed for `samtools view -s`. Using the same seed
  across samples ensures reproducibility. Default `42L`.

- threads:

  Integer. Threads for samtools. Default `4L`.

- remove_input:

  Logical. Remove the input BAM (and its index) after downsampling.
  Default `FALSE`.

- force:

  Logical. If `FALSE` (default), skip downsampling when
  `_downsampled.bam` already exists. Set `TRUE` to redo.

## Value

Character. Path to the downsampled BAM file, invisibly.

## Details

The `scale_factor` should come from the `scale_factor_bw` column of
[`HORIZON_compute_spike_in_factors`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md):
it is the ratio \\\min(N) / N_i\\ (\\\leq 1\\), where the sample with
the fewest spike-in reads receives a factor of 1 (no downsampling) and
all others are scaled down proportionally.

When `scale_factor = 1`, no reads are discarded. A copy of the input BAM
is written under the `_downsampled.bam` name so that all samples in a
batch have a consistent output file and `.latest_bam()` behaves
identically for every sample.
