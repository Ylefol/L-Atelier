# Detect Sample Outliers in Mass Spectrometry Data

Flags multivariate outliers in a `massspec_data` object using PCA-based
Euclidean distance from the sample centroid. Mean imputation and
per-protein scaling are applied internally for the outlier computation
only — log2 values in `$wide` are never modified.

## Usage

``` r
HADES_detect_outliers_massspec(
  massspec_data,
  n_pcs = 5L,
  threshold = 3,
  method = c("sd", "mad"),
  filter = FALSE,
  verbose = TRUE
)
```

## Arguments

- massspec_data:

  A `massspec_data` object after filtering and normalization.

- n_pcs:

  Integer. Number of principal components for the distance computation.
  Default: 5.

- threshold:

  Numeric. Multiplier applied to the spread measure above the centre.
  Default: 3.

- method:

  Character. `"sd"` (mean + threshold × sd) or `"mad"` (median +
  threshold × mad). Default: `"sd"`.

- filter:

  Logical. If TRUE, flagged samples are removed from `$wide` and
  `$sample_meta`. Default: FALSE (annotate only).

- verbose:

  Logical. Print detection summary. Default: TRUE.

## Value

A `massspec_data` object. `$sample_meta` gains `outlier_distance` and
`outlier` columns.

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- HADES_detect_outliers_massspec(ms, filter = FALSE)
ms$sample_meta[ms$sample_meta$outlier, c("SampleID", "outlier_distance")]
} # }
```
