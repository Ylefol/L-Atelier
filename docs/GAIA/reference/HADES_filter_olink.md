# Filter Olink Samples

Applies QC-based filtering to an `olink_data` object: optionally removes
non-biological control samples and/or samples that failed Olink's
SampleQC. Updates `$data`, `$wide`, and `$sample_meta` to reflect the
filtered set.

## Usage

``` r
HADES_filter_olink(
  olink_data,
  keep_controls = FALSE,
  filter_qc = TRUE,
  verbose = TRUE
)
```

## Arguments

- olink_data:

  An `olink_data` object from
  [`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md).

- keep_controls:

  Logical. Retain non-SAMPLE SampleTypes (e.g. NEGATIVE_CONTROL,
  PLATE_CONTROL) in the output. Default: FALSE (removes them, keeping
  only SampleType == "SAMPLE").

- filter_qc:

  Logical. Remove samples where SampleQC == "FAIL". Default: TRUE.

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

A filtered `olink_data` object of the same class. `$data`, `$wide`, and
`$sample_meta` contain only the retained samples. `$assay_meta` is
unchanged (warn_fraction is computed from the original load and reflects
all samples).

## Examples

``` r
if (FALSE) { # \dontrun{
ol_raw <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
ol <- HADES_filter_olink(ol_raw, keep_controls = FALSE, filter_qc = TRUE)
print(ol)
} # }
```
