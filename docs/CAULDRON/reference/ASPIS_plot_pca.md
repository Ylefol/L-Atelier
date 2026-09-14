# PCA scatter plot

Plots two principal components from a `reducedDims` slot, coloured by a
`colData` column or gene expression value. Axis labels include the
percentage of variance explained when a `percentVar` attribute is
present on the reduced-dimension matrix (set automatically by
[`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md),
or by the caller for custom PCA objects).

## Usage

``` r
ASPIS_plot_pca(
  sce,
  colour_by = "cluster",
  dimred = "PCA",
  point_size = 0.8,
  point_alpha = 0.6,
  palette = NULL,
  title = NULL,
  label_clusters = FALSE,
  label_size = 4,
  assay_name = "logcounts",
  ncol = NULL,
  show_legend = TRUE,
  legend_point_size = 4
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with the target reduced dimension populated.

- colour_by:

  Character scalar or vector. One or more `colData` column names or gene
  names present in `rownames(sce)`. A single value returns a `ggplot`; a
  vector produces one panel per element arranged in a grid (requires
  `gridExtra`). Default `"cluster"`.

- dimred:

  Character. Name of the `reducedDims` slot to plot. Default `"PCA"`.
  Set to a custom name (e.g. `"PCA_singler"`) when plotting a
  non-standard PCA — attach a `percentVar` numeric vector as an
  attribute on that matrix to show variance on the axes.

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
  `"PCA — <colour_by>"`.

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
  `colour_by` is a vector. `NULL` uses `min(length(colour_by), 3)`.

- show_legend:

  Logical. Whether to display the colour legend. Default `TRUE`.

- legend_point_size:

  Numeric. Size of the coloured point in the discrete legend. Default
  `4`.

## Value

A `ggplot` (single `colour_by`) or a `gtable` grid (multiple
`colour_by`).
