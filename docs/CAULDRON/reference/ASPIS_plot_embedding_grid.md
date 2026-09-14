# Grid of embedding plots from a parameter sweep

Takes the output of
[`TALOS_tune_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_umap.md)
or
[`TALOS_tune_tsne`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_tsne.md)
and plots every parameter combination as a small embedding, arranged in
a grid. The best combination (by composite score) is marked with a `★`
in its panel title.

## Usage

``` r
ASPIS_plot_embedding_grid(
  sweep,
  sce,
  colour_by = "cluster",
  point_size = 0.3,
  point_alpha = 0.5,
  palette = NULL,
  assay_name = "logcounts"
)
```

## Arguments

- sweep:

  A `talos_embedding_sweep` object produced by
  [`TALOS_tune_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_umap.md)
  or
  [`TALOS_tune_tsne`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_tsne.md).

- sce:

  A `SingleCellExperiment` used only for `colData` and assay access when
  resolving `colour_by`.

- colour_by:

  Character. A `colData` column name (e.g. `"cluster"`, `"region"`,
  `"condition"`) or a gene name. Default `"cluster"`.

- point_size:

  Numeric. Point size. Smaller values work better in a dense grid.
  Default `0.3`.

- point_alpha:

  Numeric. Point transparency. Default `0.5`.

- palette:

  Character vector or `NULL`. Passed to the colour scale — see
  [`ASPIS_plot_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_umap.md)
  for details. Default `NULL` uses built-in palettes.

- assay_name:

  Character. Assay used when `colour_by` is a gene. Default
  `"logcounts"`.

## Value

A `ggplot` object.

## Details

Embeddings are read directly from the sweep object (pre-computed during
the tuning run) so no re-computation is needed.
