# Filter Olink Proteins by AssayQC Warn Fraction

Removes proteins (rows) from an `olink_data` object whose AssayQC warn
fraction across samples exceeds a threshold. A high warn fraction
indicates systematic assay instability for that protein.

`warn_fraction` is computed at load time by
[`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md)
and stored in `$assay_meta`. It reflects the proportion of samples where
that protein received `AssayQC == "WARN"`.

## Usage

``` r
HADES_filter_olink_assay_warn(
  olink_data,
  max_warn_fraction = 0.2,
  verbose = TRUE
)
```

## Arguments

- olink_data:

  An `olink_data` object from
  [`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md).

- max_warn_fraction:

  Numeric (0–1). Maximum allowable proportion of samples with
  `AssayQC == "WARN"` for a protein to be retained. Proteins exceeding
  this threshold are removed. Default: `0.2` (20%). No established
  universal standard exists; verify this threshold against your panel
  and study design.

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

A filtered `olink_data` object. `$wide`, `$data`, and `$assay_meta`
contain only retained proteins. `$sample_meta` is unchanged.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- HADES_filter_olink(ol_raw)
ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
ol <- HADES_filter_olink_lod(ol, max_lod_fraction = 0.5)
ol <- HADES_filter_olink_assay_warn(ol, max_warn_fraction = 0.2)
} # }
```
