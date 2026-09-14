# Filter Olink Proteins by Missingness

Removes proteins (rows) from an `olink_data` object whose NPX
missingness fraction across samples exceeds a threshold. High-NA
proteins are dropped because missingness cannot be reliably attributed
to a true zero signal versus technical absence.

## Usage

``` r
HADES_filter_olink_proteins(olink_data, max_na_fraction = 0.2, verbose = TRUE)
```

## Arguments

- olink_data:

  An `olink_data` object from
  [`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md)
  or
  [`HADES_filter_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_olink.md).

- max_na_fraction:

  Numeric (0-1). Maximum allowable proportion of NA NPX values across
  samples for a protein to be retained. Proteins with NA fraction above
  this threshold are removed. Default: 0.2 (20%).

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

A filtered `olink_data` object. `$wide`, `$data`, and `$assay_meta`
contain only the retained proteins.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- HADES_filter_olink(ol_raw)
ol <- HADES_filter_olink_proteins(ol, max_na_fraction = 0.2)
} # }
```
