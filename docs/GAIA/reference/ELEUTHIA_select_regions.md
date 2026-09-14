# Select Regions Based on Fragment Count or Density Threshold

Filters candidate regions from ELEUTHIA_merge_fragments() based on
fragment count or fragment density, with optional maximum width and
minimum sample support filters.

## Usage

``` r
ELEUTHIA_select_regions(
  candidate_regions,
  min_fragments = NULL,
  quantile_threshold = NULL,
  max_width = NULL,
  min_samples = NULL,
  by_density = FALSE,
  min_density = NULL,
  verbose = TRUE
)
```

## Arguments

- candidate_regions:

  Data.frame from ELEUTHIA_merge_fragments() with columns: chr, start,
  end, n_fragments, n_samples.

- min_fragments:

  Integer or NULL. Absolute minimum fragment count to keep a region.
  Ignored when `by_density = TRUE`; use `min_density` instead. If NULL,
  only quantile_threshold is used.

- quantile_threshold:

  Numeric between 0 and 1, or NULL. Keep regions above this quantile.
  When `by_density = FALSE` (default) the quantile is computed on raw
  fragment count; when `by_density = TRUE` it is computed on fragment
  density (fragments per bp). E.g., 0.90 keeps the top 10% of regions.
  If NULL, only the absolute threshold is used.

- max_width:

  Integer or NULL. Maximum region width in bp. Regions wider than this
  are removed before any other filter. Useful for discarding broad,
  diffuse accumulations that inflate raw fragment counts without genuine
  local enrichment. Default NULL (no cap).

- min_samples:

  Integer or NULL. Minimum number of samples that must have contributed
  at least one fragment to a region for it to be retained. Requires the
  `n_samples` column produced by
  [`ELEUTHIA_merge_fragments()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_merge_fragments.md).
  This filter is applied after `max_width` and before the fragment
  count/density threshold, so it can substantially reduce the candidate
  pool before quantile computation. Default NULL (no sample support
  requirement).

- by_density:

  Logical. When TRUE, threshold filtering is performed on fragment
  density (n_fragments / width in bp) rather than raw fragment count.
  This normalises for region width so that narrow, concentrated peaks
  are preferred over broad, diffuse ones with equivalent total counts.
  Default FALSE (backward-compatible behaviour).

- min_density:

  Numeric or NULL. Absolute minimum fragment density (fragments per bp)
  to keep a region. Only used when `by_density = TRUE`. If NULL, only
  quantile_threshold is used.

- verbose:

  Logical. Print selection summary (default = TRUE).

## Value

A data.frame with columns: chr, start, end, peak_id, n_fragments,
n_samples, fragment_density, width. Only regions passing all filters are
included. `fragment_density` is always returned for inspection.

## Details

Filters are applied in the following order:

1.  `max_width`: remove regions wider than the cap

2.  `min_samples`: require fragment contribution from at least N samples

3.  Fragment count or density threshold (whichever of the absolute and
    quantile thresholds is more stringent)

Applying `min_samples` before the quantile threshold is intentional: the
quantile is then computed only on the biologically supported subset,
preventing low-support regions from diluting the distribution.

**Choosing between count and density mode:** Raw fragment count
(`by_density = FALSE`) reflects total signal accumulation and is
appropriate when region widths are comparable. Fragment density
(`by_density = TRUE`) normalises for region width and better matches
RPKM-normalised BigWig tracks, making it preferable when the merging
step produces regions of highly variable size.

## See also

[`ELEUTHIA_merge_fragments`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_merge_fragments.md)
for the merging step,
[`ELEUTHIA_plot_region_distribution`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_plot_region_distribution.md)
to visualize before selecting

## Examples

``` r
if (FALSE) { # \dontrun{
# Merge first
candidates <- ELEUTHIA_merge_fragments(bed_files)

# Select by absolute fragment count (original behaviour)
regions <- ELEUTHIA_select_regions(candidates, min_fragments = 100)

# Require signal in at least 3 of N samples + top 10% by density
regions <- ELEUTHIA_select_regions(candidates,
                                    quantile_threshold = 0.90,
                                    min_samples        = 3,
                                    by_density         = TRUE)

# Full combination: width cap + sample support + density
regions <- ELEUTHIA_select_regions(candidates,
                                    quantile_threshold = 0.90,
                                    max_width          = 2000,
                                    min_samples        = 3,
                                    by_density         = TRUE)
} # }
```
