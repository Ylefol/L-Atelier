# Normalize Time Series Data with DESeq2

Creates a DESeq2 object with size factors estimated on the FULL dataset.
This is a REQUIRED first step before running time series DEA.
Normalizing all samples together ensures consistent size factors across
all subsequent comparisons.

## Usage

``` r
ARTEMIS_normalize_timeseries(
  counts,
  targets,
  group_col = "group",
  time_col = "timepoint",
  batch_col = NULL,
  verbose = TRUE
)
```

## Arguments

- counts:

  Integer matrix of raw counts (genes x samples). Rownames required.

- targets:

  Data.frame with sample metadata. Rownames must match colnames of
  counts. Must include columns for group and timepoint.

- group_col:

  Character. Column in targets for group labels. Default: "group".

- time_col:

  Character. Column in targets for timepoint labels. Default:
  "timepoint".

- batch_col:

  Character or NULL. Optional batch column to include in design.
  Default: NULL.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_norm"` containing:

- dds:

  DESeqDataSet with size factors estimated on full dataset

- norm_counts:

  Normalized count matrix

- targets:

  The targets data.frame

- size_factors:

  Named vector of size factors per sample

- parameters:

  List of parameters used

## Details

This function:

1.  Creates a DESeqDataSet from the full count matrix

2.  Estimates size factors using all samples together

3.  Stores the normalized counts and size factors

The resulting object is then passed to
[`ARTEMIS_timeseries_conditional()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_conditional.md)
or
[`ARTEMIS_timeseries_temporal()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_temporal.md)
which will subset and run DESeq on subsets while preserving the global
size factors.

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1: Normalize the full dataset
ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)

# Step 2: Run conditional DEA (uses pre-computed size factors)
cond_de <- ARTEMIS_timeseries_conditional(
  ts_norm, reference = "Control", experiment = "Treatment"
)

} # }
```
