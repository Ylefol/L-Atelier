# Sweep k values for SNN graph construction

The `k` parameter in
[`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
controls how many nearest neighbours each cell shares edges with. Small
`k` values produce fine-grained, loosely connected graphs that tend to
over-fragment; large values produce coarser graphs that may merge
distinct populations. This function sweeps a range of `k` values, builds
an SNN graph and applies community detection at each, and scores the
resulting partitions on two complementary metrics:

## Usage

``` r
TALOS_tune_k(
  sce,
  k_range = seq(5L, 50L, by = 5L),
  use_rep = "PCA",
  n_pcs = 30L,
  type = c("rank", "number", "jaccard"),
  method = c("leiden", "louvain"),
  resolution = 1,
  compute_silhouette = TRUE,
  max_cells_sil = 10000L,
  n_runs = 5L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with the chosen reduced dimension already
  populated.

- k_range:

  Integer vector of k values to test. Default `seq(5, 50, by = 5)`.

- use_rep:

  Character. Name of the `reducedDims` slot to use for graph
  construction and silhouette distances. Default `"PCA"`. Use `"scVI"`
  for the batch-corrected latent embedding from
  [`TALOS_run_scvi`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md).

- n_pcs:

  Integer. Number of PCA components to use when `use_rep = "PCA"`.
  Ignored for other representations (all dims are used). Default `30`.

- type:

  Character. Edge weight scheme for
  [`scran::buildSNNGraph`](https://rdrr.io/pkg/scran/man/buildSNNGraph.html):
  `"rank"` (default), `"number"`, or `"jaccard"`.

- method:

  Character. Clustering algorithm: `"leiden"` (default) or `"louvain"`.

- resolution:

  Numeric. Clustering resolution (held fixed during sweep). Default
  `1.0`.

- compute_silhouette:

  Logical. Compute mean silhouette width. Automatically disabled when
  `ncol(sce) > max_cells_sil`. Default `TRUE`.

- max_cells_sil:

  Integer. Cell count above which silhouette is skipped. Default
  `10000`.

- n_runs:

  Integer. Number of community-detection runs per k value. Results are
  averaged over runs (modularity mean; representative membership closest
  to the mean used for silhouette). Increasing `n_runs` dampens
  Leiden/Louvain stochasticity at the cost of proportionally longer
  runtime. Default `5L`.

- seed:

  Integer. Random seed for the first run; subsequent runs use
  `seed + 1`, `seed + 2`, … Default `42L`.

- verbose:

  Logical. Print per-k progress and a final summary. Default `TRUE`.

## Value

A `talos_k_sweep` list with:

- `results`:

  Data frame with one row per k: `k`, `n_clusters`, `modularity` (mean
  over `n_runs`), `mean_sil` (`NA` if silhouette was skipped).

- `plot`:

  Multi-panel line plot (one panel per metric) with a red dashed
  vertical line marking the suggested `k`.

- `best_k`:

  Suggested k — the value with the highest mean silhouette (or highest
  modularity when silhouette is skipped).

- `params`:

  List of sweep parameters for reproducibility.

## Details

- Modularity:

  Graph modularity of the partition — measures how cleanly communities
  separate relative to a random null graph. Higher values indicate more
  distinct cluster boundaries. Always computed.

- Mean silhouette width:

  Average silhouette width in the representation space — for each cell,
  how much more similar it is to its own cluster than to the nearest
  other cluster. Ranges from \\-1\\ (misassigned) to \\+1\\
  (well-separated). The `k` with the highest mean silhouette is
  suggested. Automatically skipped when `ncol(sce) > max_cells_sil`
  (pairwise distances are \\O(N^2)\\); modularity is used instead.

Resolution and clustering method are held fixed during the sweep — they
are not targets here. Use
[`TALOS_tune_resolution`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_resolution.md)
afterwards to optimise those parameters with the chosen `k`.
