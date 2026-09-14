# Compute TSS Enrichment Score from BigWig files

Computes per-sample TSS enrichment scores from BigWig files. The TSS
enrichment score is a standard ATAC-seq QC metric (ENCODE guidelines)
that quantifies how specifically Tn5 inserted at open chromatin near
transcription start sites relative to flanking background signal.

The score is computed by:

1.  Aggregating normalised BigWig signal across all TSS windows into a
    single mean profile.

2.  Normalising by the mean signal in the outermost `flank_bins` bins on
    each side of the window (background estimate).

3.  Reporting the maximum of the normalised profile as the enrichment
    score (typically the centre bin over the TSS).

ENCODE recommends a score \>= 6 for human ATAC-seq data as a passing
threshold. This value was derived empirically and may not be universally
appropriate — verify for your cell type and coverage depth (Corces et
al., 2018, Nature Methods, doi:10.1038/s41592-018-0105-7).

## Usage

``` r
HADES_tss_enrichment(
  sample_sheet = NULL,
  bw_file = NULL,
  txdb,
  window = 3000L,
  bin_size = 10L,
  flank_bins = 10L,
  threshold = 6,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  data.frame. Must contain at minimum columns `sample_id` and `bw_loc`
  (path to the BigWig file for each sample). An optional `condition`
  column is carried through to the result and used for colouring in
  [`AETHER_plot_tss_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_tss_enrichment.md).
  Exactly one of `sample_sheet` or `bw_file` must be provided.

- bw_file:

  Character. Path to a single BigWig file. Use when processing one
  sample without a sample sheet. Exactly one of `sample_sheet` or
  `bw_file` must be provided.

- txdb:

  TxDb object. Used to extract TSS positions. Create with
  [`APOLLO_make_txdb()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_make_txdb.md)
  or load a pre-built TxDb annotation package.

- window:

  Integer. Base pairs upstream and downstream of each TSS to include in
  the profile window. Default `3000L`.

- bin_size:

  Integer. Approximate width of each bin in bp. The window is divided
  into `2 * floor(window / bin_size) + 1` equal bins. Default `10L`.

- flank_bins:

  Integer. Number of bins at each edge of the window used to estimate
  background signal. Background is the mean of the outermost
  `flank_bins` on each side. At the default settings (bin_size = 10,
  flank_bins = 10) this represents 100 bp of background on each flank.
  Default `10L`.

- threshold:

  Numeric. Score threshold used to flag samples as passing QC. The
  ENCODE recommended value of 6 is a guideline for human data; no
  universally established threshold exists — treat this as a reference
  point, not an absolute standard. Default `6`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `hades_tss_enrichment` S3 object containing:

- `$scores`:

  data.frame with columns `sample_id`, `tss_enrichment_score`, `pass`
  (score \>= threshold), and `condition` if present in the sample sheet.

- `$profiles`:

  Named list of normalised mean signal profiles (one numeric vector per
  sample, length = `n_bins`).

- `$params`:

  List of run parameters: `window`, `bin_size`, `flank_bins`, `n_bins`,
  `n_tss`, `threshold`.

## See also

[`AETHER_plot_tss_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_tss_enrichment.md)
