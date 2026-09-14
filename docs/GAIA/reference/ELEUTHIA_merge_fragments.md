# Merge Fragments into Candidate Regions

Pools and merges fragments from BED files into candidate regions with a
baseline noise filter. This is the first step of a two-step workflow
where you can then examine the distribution and choose a final
threshold.

## Usage

``` r
ELEUTHIA_merge_fragments(
  bed_files,
  merge_distance = 150,
  min_fragments = 10,
  extend = 75,
  chrom_sizes = NULL,
  pool_samples = TRUE,
  verbose = TRUE
)
```

## Arguments

- bed_files:

  Character vector of paths to BED files, or a named list of loaded BED
  data.frames.

- merge_distance:

  Integer. Maximum distance (bp) between fragments to merge into a
  single region (default = 150).

- min_fragments:

  Integer. Minimum fragments to retain a candidate region (default =
  10). This is a noise filter - regions below this are discarded during
  merging. Set to 1 to keep all regions.

- extend:

  Integer. Extend each fragment by this many bp on each side before
  merging (default = 75). Helps connect nearby fragments.

- chrom_sizes:

  Optional. Chromosome sizes for clamping extended coordinates. Prevents
  regions from extending beyond chromosome boundaries. Can be:

  - A named numeric vector where names are chromosome names and values
    are chromosome lengths

  - A data.frame with columns 'chr' and 'size'

  - NULL (default) - no clamping, extended coordinates may exceed
    chromosome length

- pool_samples:

  Logical. If TRUE, pool all samples before merging (recommended for
  sparse data). If FALSE, not yet implemented.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A data.frame with columns: chr, start, end, n_fragments, width. Contains
candidate regions passing the min_fragments noise filter. Use
ELEUTHIA_select_regions() for further quantile-based filtering.

## Details

This function performs the merging step with a baseline noise filter:

1.  Pools fragments from all samples (if pool_samples = TRUE)

2.  Extends each fragment by `extend` bp on each side

3.  Merges overlapping/nearby extended fragments (within
    `merge_distance`)

4.  Counts fragments in each merged region

5.  Discards regions with fewer than `min_fragments` (noise filter)

The min_fragments parameter removes obvious noise (regions with very few
fragments) to make the quantile distribution more meaningful. This is
different from the final selection threshold in
ELEUTHIA_select_regions().

## See also

[`ELEUTHIA_select_regions`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_select_regions.md)
for quantile-based filtering,
[`ELEUTHIA_plot_region_distribution`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_plot_region_distribution.md)
for visualization,
[`ELEUTHIA_call_regions_from_fragments`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_call_regions_from_fragments.md)
for one-step workflow

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-step workflow:
# 1. Merge fragments into candidate regions (with noise filter)
candidates <- ELEUTHIA_merge_fragments(bed_files, min_fragments = 10)

# 2. Examine distribution of candidates
ELEUTHIA_plot_region_distribution(candidates)

# 3. Select top regions based on distribution
regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)

} # }
```
