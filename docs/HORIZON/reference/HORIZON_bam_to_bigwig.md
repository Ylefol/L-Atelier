# Convert a processed BAM to a BigWig using bamCoverage

Wraps the deeptools `bamCoverage` command-line tool. Supports two
normalization modes:

- **RPKM** (default) — suitable for ATAC-seq and ChIP-seq without a
  spike-in control.

- **Spike-in** — pass a numeric `scale_factor` from
  [`HORIZON_compute_spike_in_factors`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md).
  Uses `--normalizeUsing None --scaleFactor <value>`, which preserves
  absolute differences in chromatin occupancy across samples.

## Usage

``` r
HORIZON_bam_to_bigwig(
  sample_sheet = NULL,
  sample_id = NULL,
  bam_file = NULL,
  output_dir = NULL,
  effective_genome_size,
  scale_factor = NULL,
  threads = 4L,
  bin_size = 10L,
  ignore_for_normalization = c("chrX", "chrM"),
  ignore_duplicates = TRUE,
  extend_reads = TRUE,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).
  Required in sample-sheet mode; must be `NULL` when `bam_file` is
  supplied.

- sample_id:

  Character. Sample ID to process (sample-sheet mode). Optional in
  direct BAM mode — inferred from the filename if omitted.

- bam_file:

  Character or `NULL`. Path to a coordinate-sorted, indexed BAM file.
  When supplied, `sample_sheet` must be `NULL`. Default `NULL`.

- output_dir:

  Character or `NULL`. Output directory for direct BAM mode. Ignored in
  sample-sheet mode. Defaults to the directory containing `bam_file`.

- effective_genome_size:

  Integer, numeric, or character. The effective (mappable) genome size
  passed to `--effectiveGenomeSize`. Pass a numeric value directly, or
  one of the following assembly names:

  - `"GRCh38"` — 2913022398

  - `"GRCm38"` — 2652783500

  - `"WBcel235"` — 100272607

  - `"T2TCHM13"` — 3117292070

- scale_factor:

  Numeric or `NULL`. When `NULL` (default), RPKM normalisation is used.
  When a numeric value is supplied (e.g. from
  [`HORIZON_compute_spike_in_factors`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md)),
  spike-in calibration normalisation is used instead
  (`--normalizeUsing None --scaleFactor <value>`).

- threads:

  Integer. Number of parallel threads (`-p`). Default 4.

- bin_size:

  Integer. Bin width in base pairs (`--binSize`). Default 10.

- ignore_for_normalization:

  Character vector. Chromosome names excluded from total read count
  during RPKM normalisation (`--ignoreForNormalization`). Ignored when
  `scale_factor` is set. Default `c("chrX", "chrM")`.

- ignore_duplicates:

  Logical. Skip duplicate reads (`--ignoreDuplicates`). Default `TRUE`.

- extend_reads:

  Logical. Extend reads to fragment size (`--extendReads`). Default
  `TRUE`.

- force:

  Logical. If `FALSE` (default), skip BigWig generation when the output
  file already exists. Set `TRUE` to regenerate.

## Value

Character. Path to the BigWig file, invisibly.

## Details

**Two usage modes:**

1.  **Sample-sheet mode** (default): provide `sample_sheet` and
    `sample_id`. The most downstream available BAM is used automatically
    via `.latest_bam()` (priority: blacklist-filtered \> processed \>
    host \> mapq-filtered). Output is written to the standard HORIZON
    directory layout.

2.  **Direct BAM mode**: provide `bam_file` (path to a
    coordinate-sorted, indexed BAM). The BAM index (`.bai`) must exist.
    `sample_id` is inferred from the filename; output goes to
    `output_dir` (defaults to the directory containing the BAM).

`bamCoverage` is provided via the user-managed conda environment
registered with
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).
