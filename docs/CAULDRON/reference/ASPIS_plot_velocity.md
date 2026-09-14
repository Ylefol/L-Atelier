# Plot RNA velocity with coloured per-cell arrows and origin circles

Overlays per-cell velocity arrows on a UMAP or tSNE embedding. Each
cell's arrow and its semi-transparent origin circle share the cell's
cluster/feature colour, replicating the `scv.pl.velocity_embedding`
style. For smooth streamlines see
[`ASPIS_plot_velocity_stream`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_velocity_stream.md).

## Usage

``` r
ASPIS_plot_velocity(
  sce,
  dimred = "UMAP",
  colour_by = "cluster",
  arrow_scale = 1,
  arrow_density = 0.3,
  arrow_alpha = 0.9,
  arrow_size = 0.4,
  arrow_length = 0.08,
  circle_size = 3,
  circle_alpha = 0.2,
  point_size = 0.5,
  point_alpha = 0.8,
  palette = NULL,
  assay_name = "logcounts",
  title = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with velocity results from
  [`TRIPODES_run_velocity`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_velocity.md).

- dimred:

  Character. Embedding to use. Default `"UMAP"`.

- colour_by:

  Character. `colData` column or gene name for cell colour. Default
  `"cluster"`.

- arrow_scale:

  Numeric. Multiplier on arrow length. Default `1`.

- arrow_density:

  Numeric in (0,1\]. Fraction of cells that get arrows. Reduce for large
  datasets. Default `0.3`.

- arrow_alpha:

  Numeric. Arrow opacity. Default `0.9`.

- arrow_size:

  Numeric. Arrow line width. Default `0.4`.

- arrow_length:

  Numeric. Arrowhead length in cm. Default `0.08`.

- circle_size:

  Numeric. Size of the semi-transparent origin circle. Default `3`.

- circle_alpha:

  Numeric. Opacity of the origin circle. Default `0.2`.

- point_size:

  Numeric. Cell point size. Default `0.5`.

- point_alpha:

  Numeric. Cell point opacity. Default `0.8`.

- palette:

  Character vector or `NULL` for defaults.

- assay_name:

  Character. Assay for gene expression. Default `"logcounts"`.

- title:

  Character or `NULL`. Plot title.

## Value

A `ggplot` object.
