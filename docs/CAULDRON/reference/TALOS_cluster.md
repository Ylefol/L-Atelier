# Graph-based clustering

Applies Leiden or Louvain community detection to the SNN graph built by
[`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md).
Cluster labels are stored in `colData(sce)[[cluster_col]]`.

## Usage

``` r
TALOS_cluster(
  sce,
  method = c("leiden", "louvain"),
  resolution = 1,
  cluster_col = "cluster",
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `metadata(sce)$snn_graph` populated
  (i.e.
  [`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
  must have been run first).

- method:

  Character. `"leiden"` (default) or `"louvain"`. Leiden generally
  produces better-connected, more reproducible partitions.

- resolution:

  Numeric. Resolution parameter controlling cluster granularity. Higher
  values yield more, smaller clusters. Default `1.0`. Applied to both
  methods.

- cluster_col:

  Character. `colData` column name for storing labels. Default
  `"cluster"`.

- seed:

  Integer. Random seed. Default `42L`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with cluster labels in `colData(sce)[[cluster_col]]`.
