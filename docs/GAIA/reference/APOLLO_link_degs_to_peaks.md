# Link Differentially Expressed Genes to Nearby ATAC Peaks

Identifies ATAC-seq peaks located within a defined distance of the TSS
of differentially expressed genes (DEGs). Uses a window-based approach
that captures all peak-gene pairs within range (many-to-many), not just
the single nearest gene per peak.

Significance filtering should be applied to `deg_result` by the caller
before passing it to this function.

## Usage

``` r
APOLLO_link_degs_to_peaks(
  deg_result,
  peaks,
  gtf_file = NULL,
  txdb = NULL,
  distance = 50000,
  gene_col = "gene_name",
  chr_mapping = NULL,
  verbose = TRUE
)
```

## Arguments

- deg_result:

  DEA result. Either an `artemis_dea` S3 object from
  `ARTEMIS_perform_dea()`, or a plain data.frame. The data.frame must
  contain at minimum a gene name column (see `gene_col`).

- peaks:

  Peak data. Either a data.frame with at least three columns (chr,
  start, end — BED-style 0-based starts) or a `GRanges` object.

- gtf_file:

  Character. Path to GTF/GFF annotation file. Used to extract TSS
  positions. Either `gtf_file` or `txdb` must be provided. When both are
  given, `gtf_file` takes precedence.

- txdb:

  A TxDb object (from
  [`APOLLO_make_txdb()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_make_txdb.md)).
  Alternative to `gtf_file`. Gene IDs in the TxDb must match the values
  in `gene_col` of `deg_result` (e.g. both use gene symbols).

- distance:

  Numeric. Maximum distance in bp from a peak (any edge) to a TSS for
  the pair to be included. Default `50000` (50 kb).

- gene_col:

  Character. Column in `deg_result` holding gene symbols. Default
  `"gene_name"`.

- chr_mapping:

  Character or named vector. Chromosome renaming to apply to gene
  coordinates from the GTF so they match peak coordinates. Pass `"T2T"`
  for the T2T-CHM13 NCBI→UCSC mapping. Ignored when using `txdb`
  (mapping is applied at TxDb creation time). Default `NULL` (no
  renaming).

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

An S3 object of class `"apollo_deg_peaks"` (a named list):

- `$linked`:

  data.frame. One row per peak × DEG pair within `distance`. Contains
  all columns from `deg_result`, plus peak coordinates, TSS position,
  and `distance_to_tss` (signed: negative = peak is upstream of TSS
  relative to transcription direction, positive = downstream).

- `$deg_summary`:

  data.frame. One row per input DEG. Includes all `deg_result` columns
  plus `has_peak` (logical), `n_peaks` (integer), and
  `nearest_peak_distance` (bp).

- `$peak_summary`:

  data.frame. One row per input peak. Includes peak coordinates plus
  `near_deg` (logical), `n_degs`, `nearest_deg` (gene name), and
  `nearest_deg_distance`.

- `$stats`:

  Named list of summary counts.

- `$params`:

  Named list of parameters used.

## Details

**Signed distance convention** (strand-aware):

- \+ strand: `distance_to_tss = peak_midpoint - TSS_position`

- \- strand: `distance_to_tss = TSS_position - peak_midpoint`

Negative values = peak is upstream (regulatory region); positive =
downstream (gene body or beyond).

**TSS extraction:** When using `gtf_file`, only gene-level records are
read, making import fast even for large annotation files. The most 5'
TSS per gene symbol is used (in transcription direction). Gene name is
read from the `gene_name` or `gene` GTF attribute (tried in that order).
