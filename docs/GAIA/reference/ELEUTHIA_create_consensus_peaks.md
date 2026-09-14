# Create Consensus Peak Set

Merges overlapping peaks across multiple samples to create a consensus
peak set for quantification.

## Usage

``` r
ELEUTHIA_create_consensus_peaks(
  peak_list,
  min_overlap = 1,
  merge_distance = 0,
  sample_groups = NULL,
  verbose = TRUE
)
```

## Arguments

- peak_list:

  Named list of peak data.frames (output from
  ELEUTHIA_load_peaks_from_sheet).

- min_overlap:

  Integer. Minimum number of samples a peak must appear in to be
  included in consensus (default = 1, include all).

- merge_distance:

  Integer. Maximum distance (bp) between peaks to merge them into a
  single region (default = 0, only merge overlapping peaks).

- sample_groups:

  Named character vector or NULL. Maps sample IDs (names) to group
  labels (values). When provided, a `groups` column is added to the
  output containing sorted, comma-separated unique group labels that
  contributed to each consensus peak (default = NULL, no group column
  added).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A data.frame with columns: chr, start, end, peak_id, n_samples, and
optionally groups (when `sample_groups` is supplied). Each row
represents a consensus peak region.

## Details

The function:

1.  Combines all peaks from all samples

2.  Sorts by chromosome and position

3.  Merges overlapping/adjacent peaks (within merge_distance)

4.  Counts how many samples contributed to each merged region

5.  Filters by min_overlap threshold

6.  Assigns unique peak_id to each consensus region

## Examples

``` r
if (FALSE) { # \dontrun{
peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
consensus <- ELEUTHIA_create_consensus_peaks(peak_list, min_overlap = 2)

} # }
```
