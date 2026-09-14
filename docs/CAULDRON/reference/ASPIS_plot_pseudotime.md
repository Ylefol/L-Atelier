# Plot pseudotime as a colour gradient on a 2D embedding

Colours cells by their Monocle3 pseudotime value on a UMAP or tSNE
embedding. Cells with `NA` pseudotime (those in partitions without a
root) are drawn in `na_colour`.

## Usage

``` r
ASPIS_plot_pseudotime(
  sce,
  dimred = NULL,
  pseudotime_col = "monocle_pseudotime",
  na_colour = "grey80",
  palette = NULL,
  point_size = 0.5,
  point_alpha = 0.7,
  title = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with pseudotime results from
  [`TRIPODES_run_monocle`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_monocle.md).

- dimred:

  Character. Embedding to plot on. `NULL` (default) uses the embedding
  recorded in `metadata(sce)$monocle_dimred`; otherwise specify any name
  present in `reducedDimNames(sce)`.

- pseudotime_col:

  Character. `colData` column containing pseudotime values. Default
  `"monocle_pseudotime"`.

- na_colour:

  Character. Colour for cells with `NA` pseudotime. Default `"grey80"`.

- palette:

  Character vector of colours for the gradient, or `NULL` for viridis
  (default).

- point_size:

  Numeric. Cell point size. Default `0.5`.

- point_alpha:

  Numeric. Cell point opacity. Default `0.7`.

- title:

  Character or `NULL`. Plot title.

## Value

A `ggplot` object.
