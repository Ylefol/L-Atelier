# QRILC Imputation for Mass Spectrometry Data

Imputes missing values in a `massspec_data` object using Quantile
Regression Imputation of Left-Censored data (QRILC, Lazar et al. 2016).
Designed for DIA proteomics data where missingness is predominantly MNAR
(missing-not-at-random) due to values falling below the instrument
detection limit.

## Usage

``` r
POSEIDON_impute_massspec(massspec_data, verbose = TRUE)
```

## Arguments

- massspec_data:

  A `massspec_data` S3 object as produced by
  [`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md)
  and processed through
  [`POSEIDON_normalize_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_normalize_massspec.md)
  and
  [`POSEIDON_correct_batch_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_correct_batch_massspec.md).
  The `$wide` slot must be a log-scale matrix (log2 or VSN glog) with
  NAs marking missing proteins.

- verbose:

  Logical. Print imputation summary. Default: TRUE.

## Value

The input `massspec_data` object with three modifications:

- `$wide`:

  Fully imputed matrix — no NAs.

- `$na_mask`:

  Logical matrix (same dimensions as `$wide`), TRUE at positions that
  were originally NA. Retained so downstream code can distinguish real
  measurements from imputed values if needed.

- `$params$imputation`:

  List recording method, total imputed count, and percentage.

## Details

QRILC assumes that missing values exist because the true protein
abundance fell below the detection limit of the instrument
(left-censoring). For each sample, it fits a quantile regression model
using the observed intensity distribution to estimate the shape of the
unobserved left tail, then draws imputed values from that estimated
tail. The censoring threshold is inferred per sample from its own
missing-value proportion — no explicit threshold parameter is required.

This function should be called **after** normalization and batch
correction, so that imputed values land on the final processed scale.
The original NA positions are preserved in `$na_mask` so that
differential analysis can be made aware of which values are synthetic.

Requires the `imputeLCMD` package: `install.packages("imputeLCMD")`. If
not on CRAN, install from GitHub:
`remotes::install_github("cran/imputeLCMD")`.

## References

Lazar, C. et al. (2016). Accounting for the Multiple Natures of Missing
Values in Label-Free Quantitative Proteomics Data Sets to Compare
Imputation Strategies. *Journal of Proteome Research*, 15(4), 1116–1125.
