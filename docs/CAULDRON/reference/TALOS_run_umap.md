# UMAP embedding

Computes a 2-dimensional UMAP embedding from a reduced-dimension
representation. Defaults to PCA; pass `use_rep = "scVI"` to embed from
the scVI latent space produced by
[`TALOS_run_scvi`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md).
Results are stored in `reducedDims(sce)[["UMAP"]]`.

## Usage

``` r
TALOS_run_umap(
  sce,
  use_rep = "PCA",
  n_pcs = 30L,
  n_neighbors = 15L,
  min_dist = 0.1,
  use_graph = TRUE,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with the chosen reduced dimension already
  populated.

- use_rep:

  Character. Name of the `reducedDims` slot to use as input. Default
  `"PCA"`. Use `"scVI"` for the batch-corrected latent embedding from
  [`TALOS_run_scvi`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md).

- n_pcs:

  Integer. Number of PCA components to use when `use_rep = "PCA"`.
  Ignored for other representations (all dims are used). Default `30`.

- n_neighbors:

  Integer. Number of nearest neighbours for the UMAP graph. Higher
  values produce a more global view; lower values emphasise local
  structure. Default `15`.

- min_dist:

  Numeric. Minimum distance between points in the embedding. Smaller
  values pack points more tightly. Default `0.1`.

- seed:

  Integer. Random seed for reproducibility. Default `42L`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `reducedDims(sce)[["UMAP"]]` populated.
