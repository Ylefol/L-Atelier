# Sweep UMAP parameters

UMAP has two parameters that most strongly shape the embedding.
`n_neighbors` governs how many local neighbours each cell considers:
small values emphasise fine local structure; large values preserve more
global topology at the cost of local detail. `min_dist` controls how
tightly points are packed: small values produce compact, well-separated
clusters; large values spread points more evenly, which can reveal
continuous structure. This function tests all combinations of both
parameters and scores each embedding on two complementary metrics:

## Usage

``` r
TALOS_tune_umap(
  sce,
  n_neighbors_range = c(10L, 15L, 20L, 30L, 50L),
  min_dist_range = c(0.05, 0.1, 0.3, 0.5),
  use_rep = "PCA",
  n_pcs = 30L,
  knn_k = 15L,
  cluster_col = "cluster",
  use_graph = TRUE,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `"PCA"` in `reducedDims` (run
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)
  first).

- n_neighbors_range:

  Integer vector of n_neighbors values to test. Default
  `c(10, 15, 20, 30, 50)`.

- min_dist_range:

  Numeric vector of min_dist values to test. Default
  `c(0.05, 0.1, 0.3, 0.5)`.

- n_pcs:

  Integer. Number of PCA components used as UMAP input and for KNN
  reference computation. Default `30`.

- knn_k:

  Integer. Number of nearest neighbours used for the KNN overlap metric.
  Default `15`.

- cluster_col:

  Character. `colData` column containing cluster labels used for
  silhouette computation. Default `"cluster"`.

- seed:

  Integer. Random seed for UMAP reproducibility. Default `42L`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `talos_embedding_sweep` list with:

- `results`:

  Data frame with one row per n_neighbors × min_dist combination:
  `n_neighbors`, `min_dist`, `knn_overlap`, `mean_sil` (`NA` if cluster
  labels absent), `composite`.

- `embeddings`:

  Named list of cells × 2 matrices — one UMAP embedding per parameter
  combination, in the same order as `results`. Passed directly to
  [`ASPIS_plot_embedding_grid()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_embedding_grid.md)
  for visual comparison.

- `plot`:

  Tile heatmap of composite score (or KNN overlap if no clusters) across
  the n_neighbors × min_dist grid, with the best combination outlined in
  red.

- `best_params`:

  Named list (`n_neighbors`, `min_dist`) for the combination with the
  highest composite score.

- `params`:

  List of sweep parameters for reproducibility.

## Details

- KNN overlap:

  Mean Jaccard similarity between each cell's `k` nearest neighbours in
  PCA space and in the UMAP embedding. Measures how faithfully local
  neighbourhood structure from PCA is preserved in the 2D layout.
  Computed once from a shared KNN reference in PCA space
  ([`BiocNeighbors::findKNN`](https://rdrr.io/pkg/BiocNeighbors/man/findKNN-methods.html)),
  which is reused across all embedding combinations for efficiency.
  Ranges from 0 (no overlap) to 1 (identical neighbourhoods).

- Mean silhouette width:

  Mean silhouette width computed in 2D UMAP space using existing cluster
  labels from `cluster_col`. Measures how visually separated the
  clusters are in the final embedding. Ranges from \\-1\\ (misassigned)
  to \\+1\\ (well-separated). Skipped when `cluster_col` is not present
  in `colData` — only KNN overlap is used for selection in that case.

Both metrics are normalised to \\\[0, 1\]\\ and averaged into a
composite score. The parameter combination with the highest composite
score is suggested as the optimum.
