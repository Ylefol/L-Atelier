# Artemis - Utility Functions

Miscellaneous utility functions for genomic analysis that don't fit into
more specific categories. Calculate Distance to Nearest Feature

For each query region, finds the nearest subject region and reports the
distance. Wraps GenomicRanges::distanceToNearest() with a convenient
data.frame interface.

## Usage

``` r
ARTEMIS_distance_to_nearest(
  query_regions,
  subject_regions,
  query_id_col = NULL,
  subject_id_col = NULL,
  verbose = TRUE
)
```

## Arguments

- query_regions:

  Data.frame with at minimum: chr, start, end. These are the regions you
  want to annotate (e.g., ChIP anchors).

- subject_regions:

  Data.frame with at minimum: chr, start, end. These are the reference
  features to find distances to (e.g., ATAC peaks).

- query_id_col:

  Character or NULL. Column name to use as query IDs. If NULL,
  auto-detects from common names (peak_id, name, region_id) or generates
  sequential IDs.

- subject_id_col:

  Character or NULL. Column name to use as subject IDs. Same
  auto-detection logic as query_id_col.

- verbose:

  Logical. Print summary statistics (default = TRUE).

## Value

A data.frame with columns:

- query_id:

  ID of the query region

- query_chr, query_start, query_end:

  Query region coordinates

- nearest_subject_id:

  ID of the nearest subject region

- subject_chr, subject_start, subject_end:

  Nearest subject coordinates

- distance:

  Distance in bp (0 if overlapping)

- overlaps:

  Logical. TRUE if distance == 0

## Details

Distance is calculated as the minimum number of base pairs separating
the two regions. Overlapping regions have distance = 0.

This function requires the GenomicRanges package from Bioconductor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find distance from ChIP anchors to nearest ATAC peak
atac_peaks <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
consensus_peaks <- ELEUTHIA_create_consensus_peaks(atac_peaks)

distances <- ARTEMIS_distance_to_nearest(
  query_regions = chip_anchors,
  subject_regions = consensus_peaks
)

# Check overlap rate
mean(distances$overlaps)  # Fraction of anchors at ATAC peaks

# Examine distance distribution for non-overlapping
hist(distances$distance[!distances$overlaps])

} # }
```
