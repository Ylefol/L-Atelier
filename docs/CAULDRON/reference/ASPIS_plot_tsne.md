# tSNE embedding plot

Plots the tSNE embedding stored in `reducedDims(sce)[["tSNE"]]`,
coloured by a `colData` column or gene expression value.

## Usage

``` r
ASPIS_plot_tsne(
  sce,
  colour_by = "cluster",
  point_size = 0.8,
  point_alpha = 0.6,
  palette = NULL,
  title = NULL,
  label_clusters = FALSE,
  label_size = 4,
  assay_name = "logcounts",
  ncol = NULL,
  style = c("points", "contour", "both"),
  n_levels = 8L,
  contour_alpha = 0.7,
  contour_adjust = 1,
  show_legend = TRUE,
  boundary_pad = 0.1,
  legend_point_size = 4
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `"tSNE"` in `reducedDims` (run
  [`TALOS_run_tsne`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_tsne.md)
  first).

- colour_by:

  Character scalar or vector. One or more `colData` column names or gene
  names present in `rownames(sce)`. A single value returns a `ggplot`; a
  vector produces one panel per element arranged in a grid (requires
  `gridExtra`). Default `"cluster"`.

- point_size:

  Numeric. Point size. Default `0.8`.

- point_alpha:

  Numeric. Point transparency (0–1). Default `0.6`.

- palette:

  Character vector or `NULL`. For discrete variables: a vector of
  colours (recycled as needed). For continuous variables or gene
  expression: a 2-element vector `c(low, high)`. `NULL` (default) uses
  built-in palettes.

- title:

  Character or `NULL`. Plot title. `NULL` auto-generates
  `"tSNE — <colour_by>"`. Ignored when `colour_by` is a vector.

- label_clusters:

  Logical. Overlay cluster centroid labels. Only applied when
  `colour_by` resolves to a discrete variable. Default `FALSE`.

- label_size:

  Numeric. Size of centroid labels. Default `4`.

- assay_name:

  Character. Assay used when `colour_by` is a gene. Default
  `"logcounts"`.

- ncol:

  Integer or `NULL`. Number of columns in the panel grid when
  `colour_by` is a vector. `NULL` (default) uses
  `min(length(colour_by), 3)`.

- style:

  Character. Visual style: `"points"` (default, classic scatter),
  `"contour"` (density ring lines only, no points), or `"both"` (points
  with density ring lines overlaid). Topographic styles require a
  discrete `colour_by`.

- n_levels:

  Integer or `"dynamic"`. Number of contour levels when using a fixed
  spacing. Pass `"dynamic"` to compute per-cluster normalised KDE
  contours. Default `8`.

- contour_alpha:

  Numeric. Opacity of the contour lines. Default `0.7`.

- contour_adjust:

  Numeric. Bandwidth multiplier for the KDE used to draw contours.
  Values above `1` produce smoother, more spread-out silhouettes; values
  below `1` tighten the contours. Default `1`.

- show_legend:

  Logical. Whether to display the colour legend. Default `TRUE`.

- boundary_pad:

  Numeric. Fractional expansion of the plot boundaries beyond the data
  range (e.g. `0.10` = 10\\ lines are clipped at the edges; decrease if
  there is too much whitespace. Default `0.10`.

## Value

A `ggplot` (single `colour_by`) or a `gtable` grid (multiple
`colour_by`).
