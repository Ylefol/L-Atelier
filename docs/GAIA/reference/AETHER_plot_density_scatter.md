# 2D Density Scatter Plot

Plots a 2D histogram where color encodes the number of data points in
each bin. Designed to accept output from
[`ARTEMIS_prepare_bw_comparison`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_prepare_bw_comparison.md)
but works with any data frame containing two numeric columns.

## Usage

``` r
AETHER_plot_density_scatter(
  data,
  x_col = "x",
  y_col = "y",
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  nbins = 100L,
  log_axes = FALSE,
  log_density = FALSE,
  fill_cap_quantile = 0.99,
  palette = "PiYG",
  palette_direction = -1L,
  color_low = NULL,
  color_high = NULL,
  show_zero_lines = FALSE,
  show_diagonal = FALSE,
  quadrant_labels = NULL
)
```

## Arguments

- data:

  data.frame with at least two numeric columns.

- x_col:

  Character. Column name for x-axis values. Default `"x"`.

- y_col:

  Character. Column name for y-axis values. Default `"y"`.

- x_label:

  Character or `NULL`. X-axis label. Falls back to
  `attr(data, "x_label")` then `x_col` if not supplied.

- y_label:

  Character or `NULL`. Y-axis label. Falls back to
  `attr(data, "y_label")` then `y_col` if not supplied.

- title:

  Character or `NULL`. Plot title. Default `NULL`.

- nbins:

  Integer. Number of bins in each dimension. Default 100.

- log_axes:

  Logical. If `TRUE`, apply `log1p` to x and y values before plotting.
  Useful for direct signal comparisons where a few high-signal outliers
  compress the distribution into a corner. Axis labels are automatically
  wrapped in `log1p(...)`. Default `FALSE`.

- log_density:

  Logical. If `TRUE`, apply `log1p` to the bin count colour scale.
  Useful when a few bins dominate the colour range. Default `FALSE`.

- fill_cap_quantile:

  Numeric (0–1) or `NULL`. Caps the fill colour scale at this quantile
  of non-empty bin counts, preventing a single dominant bin (e.g. the
  zero-signal spike at the origin) from compressing the rest of the
  palette. Bins above the cap are painted the maximum colour
  ([`scales::squish`](https://scales.r-lib.org/reference/oob.html)).
  Default `0.99`. Set to `NULL` to disable capping.

- palette:

  Character or `NULL`. Named RColorBrewer palette for the fill scale.
  Default `"PiYG"`. Set to `NULL` to use `color_low`/`color_high`
  instead.

- palette_direction:

  Integer. `-1` reverses the palette order (default); `1` uses the
  palette as-is. For `"PiYG"` the RColorBrewer order is pink→green, so
  `-1` gives green (low) → pink (high). For sequential palettes like
  `"Blues"` (light→dark), `1` is typically preferred.

- color_low:

  Character or `NULL`. Low-density fill colour. Only used when
  `palette = NULL`. Default `NULL` (falls back to `"white"`).

- color_high:

  Character or `NULL`. High-density fill colour. Only used when
  `palette = NULL`. Default `NULL` (falls back to `"darkblue"`).

- show_zero_lines:

  Logical. If `TRUE`, draw dashed reference lines at x = 0 and y = 0.
  Useful for log2 ratio plots. Default `FALSE`.

- show_diagonal:

  Logical. If `TRUE`, draw a dotted y = x reference line. Useful for
  direct signal comparisons where deviation from the diagonal indicates
  change. Default `FALSE`.

- quadrant_labels:

  `NULL`, `"auto"`, or a 4-element character vector. When `"auto"`,
  labels the four quadrants created by the zero reference lines using
  subject names parsed from the axis labels, e.g. `"K27ac(+) ATAC(+)"`
  (Q1 top-right), `"K27ac(-) ATAC(+)"` (Q2 top-left),
  `"K27ac(-) ATAC(-)"` (Q3 bottom-left), `"K27ac(+) ATAC(-)"` (Q4
  bottom-right). Subject names are extracted by stripping `" log2(...)"`
  from the axis labels. Supply a 4-element vector to use fully custom
  labels in the same corner order. Intended for use with
  `show_zero_lines = TRUE`. Default `NULL` (no labels).

## Value

A ggplot2 object.
