# Extract summit-centered BED files from narrowPeak files

Reads MACS2/MACS3 narrowPeak files and writes fixed-width BED files
centered on each peak summit (column 10: 0-based offset from peak
start). The output is a named vector of BED paths that feeds directly
into
[`APOLLO_homer_motif_enrichment_batch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_homer_motif_enrichment_batch.md).

## Usage

``` r
APOLLO_extract_summits(
  peak_files,
  output_dir,
  half_window = 100L,
  chr_mapping = NULL,
  skip_missing = TRUE,
  verbose = TRUE
)
```

## Arguments

- peak_files:

  Named character vector of narrowPeak file paths. Names are used as set
  identifiers and as the base of output BED filenames.

- output_dir:

  Directory to save summit BED files. Created if it doesn't exist.

- half_window:

  Integer. Half-width (bp) around each summit. Total window =
  `2 * half_window`. Default: 100 (200 bp total).

- chr_mapping:

  Optional chromosome name conversion applied to the chr column before
  writing the BED file. Useful when peak files use UCSC-style names
  (chr1, chr2) but the HOMER genome was indexed with NCBI accessions
  (e.g., NC_060925.1 for T2T). Accepts:

  - A genome name string (e.g., `"T2T"`) — uses the built-in mapping
    from
    [`APOLLO_get_chr_mapping()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_chr_mapping.md),
    automatically reversed to UCSC → NCBI direction.

  - A named character vector where names are the source chr names and
    values are the target chr names.

  - `NULL` (default) — no conversion.

  Peaks whose chromosome is not found in the mapping are dropped with a
  warning.

- skip_missing:

  Logical. If TRUE, missing peak files are skipped with a warning rather
  than causing an error. Default: TRUE.

- verbose:

  Logical. Print per-set progress. Default: TRUE.

## Value

Named character vector of output BED file paths (one per input set).
Silently drops entries for files that were skipped or contained no valid
summits.

## Details

NarrowPeak column 10 is the 0-based distance from the peak start to the
summit. Rows where this value is `-1` (as in broadPeak-style output with
no called summit) are skipped. Coordinates that would extend below 0 are
clamped to 0; no upper-boundary clamping is applied (chromosome sizes
are not required).

The output BED is 6-column (chr, start, end, name, score, strand). Name
and score are taken from the narrowPeak where available.

## Examples

``` r
if (FALSE) { # \dontrun{
peak_files <- c(
  KO1a = "data/ATACseq/KO1a_sorted_peaks.narrowPeak",
  WT1a = "data/ATACseq/WT1a_sorted_peaks.narrowPeak"
)

# Standard hg38
summit_beds <- APOLLO_extract_summits(peak_files, "results/summits/")

# T2T: convert chr1/chr2/... to NC_060925.1/NC_060926.1/...
summit_beds <- APOLLO_extract_summits(peak_files, "results/summits/", chr_mapping = "T2T")

# Feed directly into HOMER batch
batch <- APOLLO_homer_motif_enrichment_batch(summit_beds, "T2T", "results/motifs/")
} # }
```
