# Time Series Conditional Differential Expression

Compares two groups at each timepoint independently. For example,
Treatment vs Control at TP1, TP2, TP3. Uses size factors from the
pre-normalized dataset to ensure consistent normalization.

## Usage

``` r
ARTEMIS_timeseries_conditional(
  ts_norm,
  reference,
  experiment,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- ts_norm:

  An `artemis_ts_norm` object from
  [`ARTEMIS_normalize_timeseries()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries.md).
  This ensures all comparisons use consistent size factors estimated on
  the full dataset.

- reference:

  Character. Reference/baseline group (denominator in FC).

- experiment:

  Character. Experimental group (numerator in FC).

- alpha:

  Numeric. FDR threshold. Default: 0.05.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_de"` containing:

- results:

  Named list of per-timepoint DESeq2 result data.frames

- summary:

  Data.frame: timepoint, n_tested, n_sig_up, n_sig_down

- type:

  "conditional"

- comparison:

  Character description

- norm_counts:

  Normalized count matrix from ts_norm

- parameters:

  List of parameters used

## Examples

``` r
if (FALSE) { # \dontrun{
# First normalize
ts_norm <- ARTEMIS_normalize_timeseries(counts, targets)

# Then run conditional DEA
cond_de <- ARTEMIS_timeseries_conditional(
  ts_norm, reference = "Control", experiment = "Treatment"
)

} # }
```
