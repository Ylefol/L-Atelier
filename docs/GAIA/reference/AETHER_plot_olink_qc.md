# Olink QC Plots

Generates QC outputs for an `olink_data` object: (1) per-plate NPX
boxplot, (2) SampleQC pass/fail summary (table), and optionally (3) a
histogram of AssayQC warn fractions (plot).

## Usage

``` r
AETHER_plot_olink_qc(
  olink_data,
  plate_col = "PlateID",
  plot_warn_proteins = FALSE,
  warn_threshold = 0.1,
  font_size = 9,
  title = NULL
)
```

## Arguments

- olink_data:

  An `olink_data` object from
  [`ELEUTHIA_load_olink()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_olink.md).

- plate_col:

  Character. Column in `$sample_meta` (or `$data`) to group the
  distribution plot by (typically PlateID). Falls back to PlateID →
  SampleID if not found. Default: `"PlateID"`.

- plot_warn_proteins:

  Logical. Include a histogram of AssayQC warn fractions across all
  proteins, with a vertical line at `warn_threshold`. Default: FALSE.

- warn_threshold:

  Numeric (0–1). Reference line drawn on the warn fraction histogram.
  Default: 0.1.

- font_size:

  Numeric. Base font size for all plots. Default: 9.

- title:

  Character or NULL. Optional prefix prepended to each plot title.
  Default: NULL.

## Value

A named list with the following elements:

- npx_distributions:

  ggplot. Boxplot of NPX values per plate (or group).

- sample_qc_table:

  data.frame. PASS/FAIL counts per SampleType (or overall if SampleType
  is absent). More legible than a bar chart when most samples pass.

- warn_proteins:

  (Only if `plot_warn_proteins = TRUE`.) ggplot. Histogram of
  warn_fraction across all proteins, with a dashed reference line at
  `warn_threshold`.

## Details

SampleQC data is drawn from `$sample_meta`, which retains all samples
(biological + controls) regardless of the `keep_controls` setting used
at load time. This gives a complete picture of plate-level QC.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")

# Basic QC plots
qc <- AETHER_plot_olink_qc(ol)
print(qc$npx_distributions)
print(qc$sample_qc_table)

# Include warn protein barplot
qc <- AETHER_plot_olink_qc(ol, plot_warn_proteins = TRUE, warn_threshold = 0.05)
print(qc$warn_proteins)
} # }
```
