# Normalize a Log2 Mass Spectrometry Intensity Matrix

Normalizes a log2-scale DIA mass spectrometry intensity matrix using
either Variance Stabilizing Normalization (VSN, default) or per-sample
median centering.

## Usage

``` r
POSEIDON_normalize_massspec(
  wide,
  method = c("vsn", "median_centering"),
  verbose = TRUE
)
```

## Arguments

- wide:

  Numeric matrix, proteins (rows) × samples (cols). Values must be
  log2-transformed (as produced by
  [`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md)).
  Both methods expect log2 input; VSN back-transforms to raw intensities
  internally before fitting.

- method:

  Character. Normalization method: `"vsn"` (default) or
  `"median_centering"`.

- verbose:

  Logical. Print normalization summary. Default: TRUE.

## Value

A numeric matrix of the same dimensions, normalized. NA values are
preserved in their original positions.

## Details

**VSN** (Huber et al. 2002, Bioinformatics): fits a generalized log
(glog/arsinh) calibration per sample via robust iterative regression on
the raw intensity scale. Simultaneously normalizes cross-sample shifts
AND stabilizes variance across the full intensity range
(heteroscedasticity correction) — proteins at low and high abundance are
treated equally. Output is on a log2-like scale (the "h" scale). Because
the pipeline stores log2 data in `ms$wide`, VSN mode back-transforms
(`2^wide`) before fitting and re-applies the original NA mask after,
since [`vsn::justvsn()`](https://rdrr.io/pkg/vsn/man/justvsn.html)
fitting is driven by observed values only. Requires the `vsn`
Bioconductor package: `BiocManager::install("vsn")`.

**Median centering**: subtracts the per-sample median (non-NA values)
from each column. Removes run-to-run location shifts on the log2 scale
but does not correct heteroscedasticity. Retained as a lightweight
fallback or for comparison.

For the analysis pipeline:

    ms$wide <- POSEIDON_normalize_massspec(ms$wide)

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
ms$wide <- POSEIDON_normalize_massspec(ms$wide, method = "vsn")
} # }
```
