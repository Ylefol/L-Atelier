# Filter Mass Spectrometry Proteins by Missingness

Removes proteins (rows) from a `massspec_data` object whose missingness
fraction across samples exceeds a threshold. DIA data typically has high
per-protein missingness, so a higher threshold than used for Olink data
is appropriate.

Also removes rows with no primary protein identifier (the first column
of `$protein_meta` is `NA`) before applying the missingness threshold.
[`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md)
already drops DIA-NN summary rows (e.g. "Sum intensities", "Number of
valid values") at load time, so this check acts as a secondary safety
net for `massspec_data` objects constructed outside of
[`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md).

## Usage

``` r
HADES_filter_massspec_proteins(
  massspec_data,
  max_na_fraction = 0.5,
  verbose = TRUE
)
```

## Arguments

- massspec_data:

  A `massspec_data` object, usually after
  [`HADES_filter_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec.md).

- max_na_fraction:

  Numeric (0–1). Maximum allowable proportion of NA values across
  samples. Proteins above this threshold are removed. Default: `0.5`
  (retain proteins detected in at least 50% of samples).

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

A filtered `massspec_data` object. `$sample_meta` is unchanged.

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- HADES_filter_massspec(ms_raw)
ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
} # }
```
