# Call Regions from Fragment BED Files (Sparse Data)

Convenience wrapper that combines ELEUTHIA_merge_fragments() and
ELEUTHIA_select_regions() into a single call. For more control, use the
two-step workflow with those functions directly.

## Usage

``` r
ELEUTHIA_call_regions_from_fragments(
  bed_files,
  merge_distance = 150,
  min_fragments = 5,
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

  Integer. Minimum number of fragments required in a region to call it
  as a "peak" (default = 5).

- extend:

  Integer. Extend each fragment by this many bp on each side before
  merging (default = 75). Helps connect nearby fragments.

- chrom_sizes:

  Optional. Chromosome sizes for clamping extended coordinates.

- pool_samples:

  Logical. If TRUE, pool all samples before calling regions (recommended
  for sparse data). (default = TRUE).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A data.frame with columns: chr, start, end, peak_id, n_fragments, width.
Each row represents a called region.

## Details

This is a convenience function that wraps the two-step workflow:

1.  ELEUTHIA_merge_fragments() - merge fragments into candidate regions

2.  ELEUTHIA_select_regions() - filter by min_fragments threshold

For more control (e.g., quantile-based selection, examining
distribution), use the two-step workflow directly.

## See also

[`ELEUTHIA_merge_fragments`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_merge_fragments.md),
[`ELEUTHIA_select_regions`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_select_regions.md),
[`ELEUTHIA_plot_region_distribution`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_plot_region_distribution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# One-step workflow (quick)
regions <- ELEUTHIA_call_regions_from_fragments(bed_files, min_fragments = 100)

# Two-step workflow (more control)
candidates <- ELEUTHIA_merge_fragments(bed_files)
ELEUTHIA_plot_region_distribution(candidates)  # examine distribution
regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)

} # }
```
