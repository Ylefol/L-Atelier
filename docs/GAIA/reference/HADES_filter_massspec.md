# Filter Mass Spectrometry Samples

Filters a `massspec_data` object by sample type and/or per-sample
missingness. Updates `$wide` and `$sample_meta` to reflect the retained
samples.

## Usage

``` r
HADES_filter_massspec(
  massspec_data,
  keep_types = c("SAMPLE", "POOL"),
  max_sample_na_fraction = 0.9,
  verbose = TRUE
)
```

## Arguments

- massspec_data:

  A `massspec_data` object from
  [`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md).

- keep_types:

  Character vector. SampleType values to retain. Default keeps
  `"SAMPLE"` (all biological samples — numeric IDs, SEP###-#
  longitudinal, and SEP_CTR\_\* controls; their distinction is carried
  by the `Group` metadata column) and `"POOL"` (pooled QC injections
  used for batch drift monitoring). `"OTHER"` is excluded by default:
  these are columns whose SampleID did not match any known pattern in
  [`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md),
  most commonly samples not renamed via `prep_MS_data.R` or absent from
  the metadata file. They have no metadata and cannot be used in any
  downstream analysis.

- max_sample_na_fraction:

  Numeric (0–1). Samples whose log2 NA fraction across all proteins
  exceeds this threshold are removed. Default: `0.9` (remove only
  extreme outliers; DIA data normally has ~70% overall missingness).

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

A filtered `massspec_data` object. `$protein_meta` is unchanged.

## Examples

``` r
if (FALSE) { # \dontrun{
ms_raw <- ELEUTHIA_load_massspec("prepped_ms_data.csv",
                                  metadata_file = "meta.csv")
ms <- HADES_filter_massspec(ms_raw)
} # }
```
