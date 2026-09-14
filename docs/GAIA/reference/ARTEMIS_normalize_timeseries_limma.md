# Normalize Time Series Data for limma

Validates and stores a log2-scale matrix with sample metadata as the
first step in the limma-based time series workflow. Unlike the DESeq2
workflow (`ARTEMIS_normalize_timeseries`), no transformation is applied
here: the matrix is assumed to be already on a log2 scale (Olink NPX,
log2 mass spec intensities, microarray log2 values). This step aligns
samples, validates required columns, and packages the inputs into an
`artemis_ts_norm_limma` object consumed by
[`ARTEMIS_timeseries_conditional_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_conditional_limma.md)
and
[`ARTEMIS_timeseries_temporal_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_temporal_limma.md).

## Usage

``` r
ARTEMIS_normalize_timeseries_limma(
  matrix,
  targets,
  group_col = "group",
  time_col = "timepoint",
  batch_col = NULL,
  verbose = TRUE
)
```

## Arguments

- matrix:

  Numeric matrix, features x samples (log2-scale). Rownames required.
  Columns must match `rownames(targets)`.

- targets:

  Data.frame with sample metadata. Rownames must match
  `colnames(matrix)`. Must include columns for group and timepoint.

- group_col:

  Character. Column in targets for group labels. Default: "group".

- time_col:

  Character. Column in targets for timepoint labels. Default:
  "timepoint".

- batch_col:

  Character or NULL. Batch column in targets. When provided, it is
  automatically included as a covariate in all downstream comparisons
  (prepended to any `covariates` passed to conditional/temporal
  functions). Default: NULL.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_norm_limma"` containing:

- matrix:

  The aligned log2-scale matrix (features x samples)

- targets:

  The aligned targets data.frame

- parameters:

  List: group_col, time_col, batch_col, n_samples, n_features

## Details

**Olink NPX**: values are already log2-normalized by the Olink platform.
Pass `ol$wide` directly — no prior transformation needed.

**Mass spectrometry**:
[`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md)
stores log2-transformed intensities in `massspec_data$wide`. Verify that
sample-level normalization (e.g., median centering) has been applied
upstream if required by your experiment design before passing the matrix
here — unlike Olink, DIA-NN output does not guarantee cross-sample
comparability without an explicit normalization step.

## Examples

``` r
if (FALSE) { # \dontrun{
# Olink longitudinal
ts_norm <- ARTEMIS_normalize_timeseries_limma(
  ol$wide, ol$sample_meta,
  group_col = "Group", time_col = "Timepoint"
)

# Mass spec with batch
ts_norm <- ARTEMIS_normalize_timeseries_limma(
  ms$wide, ms$sample_meta,
  group_col = "Group", time_col = "Visit", batch_col = "Plate"
)
} # }
```
