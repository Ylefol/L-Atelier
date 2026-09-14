# Plot Monocle3 principal graph trajectory on a 2D embedding

Overlays the Monocle3 principal graph edges on a cell embedding, with
cells optionally coloured by pseudotime, cluster, or any other feature.
The principal graph node coordinates stored by
[`TRIPODES_run_monocle`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_monocle.md)
are used directly, so the graph is always aligned to the correct
embedding.

## Usage

``` r
ASPIS_plot_trajectory(
  sce,
  dimred = NULL,
  colour_by = "monocle_pseudotime",
  na_colour = "grey80",
  edge_colour = "black",
  edge_size = 0.8,
  edge_alpha = 0.8,
  point_size = 0.5,
  point_alpha = 0.7,
  palette = NULL,
  assay_name = "logcounts",
  title = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with trajectory results from
  [`TRIPODES_run_monocle`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_monocle.md).

- dimred:

  Character. Embedding to plot on. `NULL` (default) uses the embedding
  recorded in `metadata(sce)$monocle_dimred`. If a different embedding
  is specified, a warning is issued because the principal graph node
  coordinates were computed for the original embedding.

- colour_by:

  Character. `colData` column or gene name for cell colour. Default
  `"monocle_pseudotime"`.

- na_colour:

  Character. Colour for cells with `NA` values (e.g. unordered cells
  when `colour_by = "monocle_pseudotime"`). Default `"grey80"`.

- edge_colour:

  Character. Colour of principal graph edges. Default `"black"`.

- edge_size:

  Numeric. Line width of graph edges. Default `0.8`.

- edge_alpha:

  Numeric. Opacity of graph edges. Default `0.8`.

- point_size:

  Numeric. Cell point size. Default `0.5`.

- point_alpha:

  Numeric. Cell point opacity. Default `0.7`.

- palette:

  Character vector or `NULL` for defaults.

- assay_name:

  Character. Assay for gene expression when `colour_by` is a gene.
  Default `"logcounts"`.

- title:

  Character or `NULL`. Plot title.

## Value

A `ggplot` object.
