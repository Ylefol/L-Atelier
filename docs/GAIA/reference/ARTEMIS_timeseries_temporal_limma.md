# Time Series Temporal Differential Expression (limma)

Compares timepoints across ALL samples, pooling groups together, to
extract the pure temporal effect — features that change over time
regardless of condition. For example, TP3 vs TP1 pools all TP3 samples
against all TP1 samples; the condition variation averages out as noise.
When subjects are measured at multiple timepoints, use `block_col` (e.g.
SubjectID) to account for repeated measures via
[`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html).

## Usage

``` r
ARTEMIS_timeseries_temporal_limma(
  ts_norm,
  comparisons = "consecutive",
  block_col = NULL,
  covariates = NULL,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- ts_norm:

  An `artemis_ts_norm_limma` object from
  [`ARTEMIS_normalize_timeseries_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries_limma.md).

- comparisons:

  Character. "consecutive" (TP2 vs TP1, TP3 vs TP2) or "all" (all
  pairwise timepoint combinations). Default: "consecutive".

- block_col:

  Character or NULL. Column in targets for the blocking variable for
  repeated measures (e.g. SubjectID when the same individuals are
  measured at multiple timepoints). When provided,
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  estimates the within-subject correlation per comparison. Default:
  NULL.

- covariates:

  Character vector or NULL. Additional covariate columns from targets.
  If `batch_col` was set in
  [`ARTEMIS_normalize_timeseries_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries_limma.md),
  it is automatically prepended. Default: NULL.

- alpha:

  Numeric. FDR threshold for counting significant features. Default:
  0.05.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_de"` containing:

- results:

  Named list of per-comparison result lists, each with `$results`
  (data.frame: feature_id, log2FoldChange, AveExpr, t, pvalue, padj, B,
  sig), `$fit`, `$summary`

- summary:

  Data.frame with comparison details and significance counts

- type:

  "temporal"

- comparison:

  Character description

- norm_counts:

  The log2 matrix

- parameters:

  List of parameters used

## Examples

``` r
if (FALSE) { # \dontrun{
ts_norm <- ARTEMIS_normalize_timeseries_limma(ol$wide, ol$sample_meta,
  group_col = "Group", time_col = "Timepoint")

# Consecutive comparisons, blocking on subject (paired design)
temp_de <- ARTEMIS_timeseries_temporal_limma(
  ts_norm, block_col = "SubjectID"
)

# All pairwise timepoint comparisons
temp_de <- ARTEMIS_timeseries_temporal_limma(
  ts_norm, comparisons = "all", block_col = "SubjectID"
)
} # }
```
