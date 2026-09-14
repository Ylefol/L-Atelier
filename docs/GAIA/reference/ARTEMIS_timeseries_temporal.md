# Time Series Temporal Differential Expression

Compares timepoints across ALL samples, pooling groups together. This
extracts the pure temporal effect - genes that change over time
regardless of condition. For example, TP3 vs TP1 compares all TP3
samples (IgM + LPS) against all TP1 samples (IgM + LPS). The
group/condition variation becomes biological noise that averages out,
leaving only the time-dependent changes.

## Usage

``` r
ARTEMIS_timeseries_temporal(
  ts_norm,
  comparisons = "consecutive",
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- ts_norm:

  An `artemis_ts_norm` object from
  [`ARTEMIS_normalize_timeseries()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries.md).

- comparisons:

  Character. "consecutive" (TP2 vs TP1, TP3 vs TP2) or "all" (all
  pairwise timepoint combinations). Default: "consecutive".

- alpha:

  Numeric. FDR threshold for counting significant genes. Default: 0.05.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_de"` containing:

- results:

  Named list of per-comparison DESeq2 results

- summary:

  Data.frame with comparison details and significance counts

- type:

  "temporal"

- norm_counts:

  Normalized count matrix from ts_norm

- parameters:

  List of parameters used

## Details

Uses size factors from the pre-normalized dataset to ensure consistent
normalization across all comparisons.

Unlike conditional DEA (which compares groups at each timepoint),
temporal DEA pools all groups together and uses timepoint as the
"condition" for DESeq2. This design captures genes with consistent
temporal behavior across all experimental conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
# First normalize
ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)

# Consecutive timepoint comparisons (TP2 vs TP1, TP3 vs TP2, etc.)
temp_de <- ARTEMIS_timeseries_temporal(ts_norm)

# All pairwise timepoint comparisons
temp_de <- ARTEMIS_timeseries_temporal(ts_norm, comparisons = "all")

} # }
```
