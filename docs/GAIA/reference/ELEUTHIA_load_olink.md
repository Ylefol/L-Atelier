# Load Olink NPX Data

Reads an Olink Explore NPX file (.parquet or .csv), optionally joins
plate layout metadata, and returns a structured `olink_data` object
containing all samples and proteins as-loaded.

This function performs no QC filtering. To remove failed samples or
non-biological controls, pass the result to
[`HADES_filter_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_olink.md).

## Usage

``` r
ELEUTHIA_load_olink(
  npx_file,
  metadata_file = NULL,
  metadata_sheet = 1,
  sample_col = "SampleID",
  npx_col = "NPX",
  restrict_to_metadata = FALSE,
  verbose = TRUE
)
```

## Arguments

- npx_file:

  Character. Path to the Olink NPX file (.parquet or .csv).

- metadata_file:

  Character or NULL. Path to plate layout metadata (.xlsx, .csv, or
  .tsv). Joined onto data by `sample_col`. Default: NULL.

- metadata_sheet:

  Integer or character. Sheet to read from .xlsx files. Default: 1.

- sample_col:

  Character. Column name present in both the NPX data and metadata used
  as the join key. Default: "SampleID".

- npx_col:

  Character. Which NPX representation to use as the primary value. One
  of `"NPX"`, `"ExtNPX"`, or `"PCNormalizedNPX"`. Default: `"NPX"`. NPX
  is log2-scale and approximately normally distributed — do NOT
  log-transform again before analysis.

- restrict_to_metadata:

  Logical. If `TRUE` and a `metadata_file` is provided, removes any
  samples from the NPX data that have no corresponding entry in the
  metadata. Useful when the NPX file contains additional samples (e.g.
  validation runs, reference samples) that are intentionally absent from
  the study metadata. Default: `FALSE`.

- verbose:

  Logical. Print loading progress. Default: TRUE.

## Value

An S3 object of class `"olink_data"` containing:

- data:

  Long-format data.frame: SampleID, OlinkID, Assay, UniProt, Panel, NPX
  (the selected NPX column, renamed to "NPX"), AssayQC. Metadata columns
  are not replicated here — join from `$sample_meta` when needed.

- wide:

  Numeric matrix: proteins (rows) x samples (cols). Rownames = OlinkID.
  Colnames = SampleID. NPX NAs propagated as-is.

- sample_meta:

  Per-sample data.frame for all samples. Contains SampleID, SampleType,
  PlateID, WellID, SampleQC, plus any joined metadata columns.

- assay_meta:

  Per-protein data.frame: OlinkID, Assay, UniProt, Panel, warn_fraction
  (proportion of samples with AssayQC == "WARN").

- params:

  List of loader parameters used (for provenance).

## Details

Olink process control probes (AssayType != "assay") are excluded during
parsing — these are internal instrument rows (ext_ctrl, inc_ctrl,
amp_ctrl), not proteins, and are not meaningful for analysis.

All biological samples and QC control samples are retained in the
returned object. Use
[`HADES_filter_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_olink.md)
to remove SampleQC failures or restrict to biological samples only.

AssayQC == "WARN" proteins are flagged but retained. `$assay_meta`
carries a `warn_fraction` column (proportion of samples with WARN).

Requires `arrow` for .parquet input; requires `readxl` for .xlsx
metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load only
ol_raw <- ELEUTHIA_load_olink("data.parquet")

# Load with plate layout metadata, then filter
ol_raw <- ELEUTHIA_load_olink(
  npx_file      = "data.parquet",
  metadata_file = "plate_layout.xlsx"
)
ol <- HADES_filter_olink(ol_raw, keep_controls = FALSE, filter_qc = TRUE)
} # }
```
