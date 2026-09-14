# Quantify BED Fragments Against Regions

Counts the number of BED file fragments overlapping each region.
Memory-efficient: loads one BED file at a time from the sample sheet.
Works for any omics type (ATACseq, CHIPseq, etc.).

## Usage

``` r
ELEUTHIA_quantify_bed(
  sample_sheet,
  regions,
  omics,
  min_overlap = 1,
  stranded = FALSE,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame with bed_loc column.

- regions:

  A regions data.frame with chr, start, end, and peak_id columns (from
  ELEUTHIA_call_regions_from_fragments or
  ELEUTHIA_create_consensus_peaks).

- omics:

  Character string. Omics type to quantify ("ATACseq", "CHIPseq", etc.).

- min_overlap:

  Integer. Minimum bp overlap required to count a fragment (default = 1,
  any overlap counts).

- stranded:

  Logical. If TRUE, only count a fragment toward a region when the
  fragment's strand (BED column 6) matches the region's strand (requires
  a `strand` column in `regions`, values `"+"`/ `"-"`). Use this for
  strand-specific peak sets (e.g. DRIPc-seq R-loop peaks, which are
  called separately per strand) where a region at the same coordinates
  on the opposite strand is a distinct feature. Default `FALSE`
  (strand-blind, matches prior behavior).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A list containing:

- counts:

  Matrix of fragment counts (regions x samples)

- annotation:

  Data.frame with region annotations

- targets:

  Data.frame with sample metadata

## Details

This function quantifies fragment overlaps with genomic regions using
BED files specified in the sample sheet's bed_loc column. It is
memory-efficient because it loads and processes one BED file at a time,
discarding each before loading the next.

A fragment is counted if it overlaps the region by at least
`min_overlap` base pairs. Each fragment is counted at most once per
region.

This function works for any omics type - it simply filters the sample
sheet by the specified omics value and loads BED files from bed_loc.

## Examples

``` r
if (FALSE) { # \dontrun{
# ATAC-seq: quantify against consensus peaks
atac_peaks <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
consensus <- ELEUTHIA_create_consensus_peaks(atac_peaks, min_overlap = 2)
atac_counts <- ELEUTHIA_quantify_bed(sample_sheet, consensus, "ATACseq")

# ChIP-seq: quantify against called regions
chip_beds <- ELEUTHIA_load_bed_from_sheet(sample_sheet, "CHIPseq")
chip_regions <- ELEUTHIA_call_regions_from_fragments(chip_beds, min_fragments = 10)
chip_counts <- ELEUTHIA_quantify_bed(sample_sheet, chip_regions, "CHIPseq")

} # }
```
