# Detect Sample Outliers in Olink Data

Flags multivariate outliers in an `olink_data` object using PCA-based
Euclidean distance from the sample centroid. Outliers are annotated in
`$sample_meta` via two new columns (`outlier`, `outlier_distance`) and
optionally removed.

Mean imputation and per-protein scaling are applied internally for the
outlier computation only — the NPX values in `$data` and `$wide` are
never modified.

## Usage

``` r
HADES_detect_outliers_olink(
  olink_data,
  n_pcs = 5L,
  threshold = 3,
  method = c("sd", "mad"),
  filter = FALSE,
  verbose = TRUE
)
```

## Arguments

- olink_data:

  An `olink_data` object from
  [`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md)
  or
  [`HADES_filter_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_olink.md).

- n_pcs:

  Integer. Number of principal components to use for the distance
  computation. Default: 5.

- threshold:

  Numeric. Multiplier applied to the spread measure (SD or MAD) above
  the mean/median distance. Samples whose distance exceeds
  `centre + threshold * spread` are flagged. Default: 3.

- method:

  Character. Spread measure: `"sd"` uses `mean + threshold * sd`;
  `"mad"` uses `median + threshold * mad`. MAD is more robust when
  multiple outliers are present. Default: `"sd"`.

- filter:

  Logical. If `TRUE`, flagged samples are removed from `$data`, `$wide`,
  and `$sample_meta`. If `FALSE` (default), samples are annotated but
  retained.

- verbose:

  Logical. Print detection summary. Default: TRUE.

## Value

An `olink_data` object. `$sample_meta` gains two columns:

- outlier_distance:

  Euclidean distance from the centroid in PC space.

- outlier:

  Logical. `TRUE` for samples exceeding the threshold.

When `filter = TRUE`, flagged samples are also removed from all data
slots.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- HADES_filter_olink(ol_raw)
ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)

# Detect and annotate only (inspect before deciding)
ol <- HADES_detect_outliers_olink(ol)
ol$sample_meta[ol$sample_meta$outlier, c("SampleID", "outlier_distance")]

# Remove detected outliers
ol <- HADES_detect_outliers_olink(ol, filter = TRUE)
} # }
```
