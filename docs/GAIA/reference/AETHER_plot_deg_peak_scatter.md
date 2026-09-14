# Scatter Plot of DEG LFC vs Distance to Nearest ATAC Peak

Visualises the output of
[`APOLLO_link_degs_to_peaks()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_link_degs_to_peaks.md)
as a scatter plot of log2 fold change against signed distance to the TSS
(negative = upstream, positive = downstream). Each point is one peak-DEG
pair. Point size encodes statistical significance (-log10 padj); colour
encodes expression direction.

## Usage

``` r
AETHER_plot_deg_peak_scatter(
  deg_peaks,
  color_up = "#D6604D",
  color_down = "#2166AC",
  alpha = 0.7,
  size_range = c(1, 5),
  label_top = 10,
  label_col = NULL,
  show_tss_line = TRUE,
  show_lfc_line = TRUE,
  title = NULL,
  x_limits = NULL
)
```

## Arguments

- deg_peaks:

  An `apollo_deg_peaks` object from
  [`APOLLO_link_degs_to_peaks()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_link_degs_to_peaks.md).

- color_up:

  Character. Colour for upregulated DEGs. Default `"#D6604D"` (red).

- color_down:

  Character. Colour for downregulated DEGs. Default `"#2166AC"` (blue).

- alpha:

  Numeric. Point transparency. Default `0.7`.

- size_range:

  Numeric vector of length 2. Min and max point sizes (mapped to -log10
  padj). Default `c(1, 5)`.

- label_top:

  Integer. Number of top genes (by padj) to label with `ggrepel`.
  Default `10`. Set to `0` to suppress labels.

- label_col:

  Character. Column in `$linked` used for labels. Default `NULL` (uses
  the `gene_col` from the linking step).

- show_tss_line:

  Logical. Draw a vertical dashed line at distance = 0 (the TSS).
  Default `TRUE`.

- show_lfc_line:

  Logical. Draw a horizontal dashed line at LFC = 0. Default `TRUE`.

- title:

  Character. Plot title. Default `NULL` (auto-generated).

- x_limits:

  Numeric vector of length 2 or `NULL`. X-axis limits in bp. Default
  `NULL` (auto from data).

## Value

A `ggplot` object.
