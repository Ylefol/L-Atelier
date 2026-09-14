# Harmony batch integration

Corrects batch effects in the PCA embedding using Harmony, while
preserving genuine biological variation. The harmonised embedding is
stored in `reducedDims(sce)[["Harmony"]]` and is designed to be passed
directly to
[`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
via `use_rep = "Harmony"`.

## Usage

``` r
TALOS_run_harmony(
  sce,
  batch_col = "sample_id",
  n_pcs = 30L,
  theta = 2,
  max_iter = 10L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `reducedDims(sce)[["PCA"]]` populated by
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md).

- batch_col:

  Character. Column in `colData(sce)` identifying the batch variable.
  Default `"sample_id"`.

- n_pcs:

  Integer. Number of PCA dimensions to pass to Harmony. Default `30L`.

- theta:

  Numeric. Diversity penalty — higher values enforce stronger mixing
  across batches. Default `2`.

- max_iter:

  Integer. Maximum Harmony iterations. Default `10L`.

- seed:

  Integer. Random seed. Default `42L`.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

SCE with `reducedDims(sce)[["Harmony"]]` populated (cells x `n_pcs`
matrix). Downstream call: `TALOS_build_graph(sce, use_rep = "Harmony")`.

## Details

Harmony iteratively finds soft cell-type clusters and removes the
component of each cluster's centroid that varies across batches, leaving
genuine between-condition differences intact.
