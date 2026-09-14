# tSNE embedding

Computes a 2-dimensional tSNE embedding from the PCA reduced dimensions.
Requires
[`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)
to have been run first. Results are stored in
`reducedDims(sce)[["tSNE"]]`.

## Usage

``` r
TALOS_run_tsne(
  sce,
  n_pcs = 30L,
  perplexity = 30,
  max_iter = 1000L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"PCA"` reduced dimension.

- n_pcs:

  Integer. Number of PCA components to use. Default `30`.

- perplexity:

  Numeric. tSNE perplexity parameter — roughly the number of effective
  nearest neighbours. Typical range 5–50. Default `30`.

- max_iter:

  Integer. Maximum number of tSNE iterations. Default `1000`.

- seed:

  Integer. Random seed. Default `42L`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `reducedDims(sce)[["tSNE"]]` populated.
