# Sweep number of PCA components

Choosing the number of PCA components (`n_pcs`) affects every downstream
step: too few discard biological signal; too many introduce technical
noise into the graph and embedding. This function sweeps a range of
`n_pcs` values and scores each choice on two complementary metrics:

## Usage

``` r
TALOS_tune_pc(
  sce,
  n_pcs_range = seq(5L, 50L, by = 5L),
  k = 20L,
  type = c("rank", "number", "jaccard"),
  method = c("leiden", "louvain"),
  resolution = 1,
  n_runs = 5L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `"PCA"` in `reducedDims` (run
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)
  first). The stored `percentVar` attribute is used for cumulative
  variance — the PCA result is not recomputed.

- n_pcs_range:

  Integer vector of component counts to test. Values exceeding the
  number of computed PCs are silently dropped. Default
  `seq(5, 50, by = 5)`.

- k:

  Integer. Number of nearest neighbours for SNN graph construction (held
  fixed). Default `20`.

- type:

  Character. Edge weight scheme for
  [`scran::buildSNNGraph`](https://rdrr.io/pkg/scran/man/buildSNNGraph.html):
  `"rank"` (default), `"number"`, or `"jaccard"`.

- method:

  Character. Clustering algorithm: `"leiden"` (default) or `"louvain"`.

- resolution:

  Numeric. Clustering resolution (held fixed). Default `1.0`.

- n_runs:

  Integer. Number of community-detection runs per `n_pcs` value.
  Modularity is averaged over all runs to dampen Leiden/Louvain
  stochasticity; the representative run (closest to mean modularity) is
  used for `n_clusters`. Default `5L`.

- seed:

  Integer. Random seed for the first run; subsequent runs use
  `seed + 1`, `seed + 2`, … Default `42L`.

- verbose:

  Logical. Print per-step progress and a final summary. Default `TRUE`.

## Value

A `talos_pc_sweep` list with:

- `results`:

  Data frame with one row per `n_pcs` value: `n_pcs`, `cum_var`
  (cumulative \\ `modularity` (mean over `n_runs`), `n_clusters`.

- `plot`:

  Two-panel line plot (cumulative variance explained and graph
  modularity vs `n_pcs`), with a red dashed vertical line marking the
  suggested value.

- `best_n_pcs`:

  Suggested number of components — the elbow of the modularity curve.

- `params`:

  List of sweep parameters for reproducibility.

## Details

- Cumulative variance explained:

  The percentage of total gene-expression variance captured by the first
  `n_pcs` components. Read directly from the stored PCA result — no
  recomputation required. Provides context for how much information each
  cutoff retains, but does not on its own indicate when additional
  components stop helping downstream tasks.

- Graph modularity:

  The SNN graph is built from the first `n_pcs` PCA components (at fixed
  `k`, `type`, `method`, and `resolution`) and community detection is
  applied. Graph modularity measures how cleanly the resulting partition
  separates relative to a random null. This is the primary optimisation
  target: the `n_pcs` at the elbow of the modularity curve is the point
  of diminishing returns — adding more components no longer meaningfully
  improves cluster separation. K, method, and resolution are held fixed;
  use
  [`TALOS_tune_k`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_k.md)
  and
  [`TALOS_tune_resolution`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_resolution.md)
  afterwards to optimise those.

The elbow in the modularity curve is detected via the kneedle method
(point of maximum perpendicular distance from the line connecting the
first and last values).
