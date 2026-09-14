# Build a shared nearest-neighbour graph

Constructs a SNN graph from a reduced-dimension representation using
[`scran::buildSNNGraph`](https://rdrr.io/pkg/scran/man/buildSNNGraph.html).
Defaults to PCA; pass `use_rep = "scVI"` to build the graph in the scVI
latent space produced by
[`TALOS_run_scvi`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md).
The graph is stored in `metadata(sce)$snn_graph` for use by
[`TALOS_cluster`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_cluster.md).

## Usage

``` r
TALOS_build_graph(
  sce,
  use_rep = "PCA",
  n_pcs = 30L,
  k = 20L,
  type = c("rank", "number", "jaccard"),
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with the chosen reduced dimension already
  populated.

- use_rep:

  Character. Name of the `reducedDims` slot to use. Default `"PCA"`. Use
  `"scVI"` for the batch-corrected latent embedding from
  [`TALOS_run_scvi`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md).

- n_pcs:

  Integer. Number of PCA components to use when `use_rep = "PCA"`.
  Ignored for other representations (all dims are used). Default `30`.

- k:

  Integer. Number of nearest neighbours. Higher values produce broader,
  more robust clusters. Default `20`.

- type:

  Character. Edge weight scheme: `"rank"` (default), `"number"`, or
  `"jaccard"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `metadata(sce)$snn_graph` populated.
