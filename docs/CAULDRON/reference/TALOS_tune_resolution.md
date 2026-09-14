# Sweep clustering resolution

Resolution controls the granularity of Leiden and Louvain community
detection: higher values push the algorithm to accept smaller, more
numerous clusters; lower values favour larger, coarser partitions. This
function sweeps a range of resolution values on the pre-built SNN graph
and scores each partition on four complementary metrics:

## Usage

``` r
TALOS_tune_resolution(
  sce,
  resolution_range = seq(0.1, 2, by = 0.1),
  method = c("both", "leiden", "louvain"),
  use_rep = "PCA",
  n_pcs = 30L,
  compute_silhouette = TRUE,
  max_cells_sil = 10000L,
  n_runs = 5L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `metadata(sce)$snn_graph` populated (run
  [`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
  first). PCA must also be present for silhouette computation.

- resolution_range:

  Numeric vector of resolution values to test. Default
  `seq(0.1, 2.0, by = 0.1)`.

- method:

  Character. `"leiden"`, `"louvain"`, or `"both"` (default).

- n_pcs:

  Integer. PCA components used for the silhouette distance matrix.
  Default `30`.

- compute_silhouette:

  Logical. Compute mean silhouette. Automatically disabled when
  `ncol(sce) > max_cells_sil`. Default `TRUE`.

- max_cells_sil:

  Integer. Cell count above which silhouette is skipped. Default
  `10000`.

- n_runs:

  Integer. Number of community-detection runs per method × resolution
  combination. Modularity is averaged over all runs; the representative
  membership (closest to mean modularity) is used for silhouette and
  cluster-size metrics. Default `5L`.

- seed:

  Integer. Random seed for the first run; subsequent runs use
  `seed + 1`, `seed + 2`, … Default `42L`.

- verbose:

  Logical. Print per-combination progress. Default `TRUE`.

## Value

A `talos_resolution_sweep` list with:

- `results`:

  Data frame with one row per method × resolution combination: `method`,
  `resolution`, `n_clusters`, `modularity` (mean over `n_runs`),
  `mean_sil` (`NA` if skipped), `min_cl_size`.

- `plot`:

  Four-panel line plot (n_clusters, modularity, mean_sil, min_cl_size vs
  resolution), one coloured line per method, with dashed vertical lines
  at the best resolution per method.

- `best_params`:

  Named list (`method`, `resolution`) for the single combination with
  the highest score across all methods.

- `best_per_method`:

  Data frame with the best-scoring row per clustering algorithm — useful
  when `method = "both"` to compare Leiden and Louvain optima side by
  side.

- `params`:

  List of sweep parameters for reproducibility.

## Details

- n_clusters:

  Number of distinct clusters produced. Rises monotonically with
  resolution; useful for anchoring biological expectations (e.g. known
  number of major cell types).

- modularity:

  Graph modularity of the partition — how cleanly community boundaries
  separate relative to a random null. Higher is better, but very high
  resolution can inflate modularity by fragmenting real populations.

- mean_sil:

  Mean silhouette width in PCA space — how confidently each cell is
  assigned to its cluster vs the nearest alternative. Ranges from \\-1\\
  (misassigned) to \\+1\\ (well-separated). The resolution per method
  with the highest mean silhouette is suggested. Skipped when
  `ncol(sce) > max_cells_sil`; modularity is used instead.

- min_cl_size:

  Size of the smallest cluster. Near-zero values flag pathological
  resolutions that fragment cells into singletons or micro-clusters,
  typically a sign the resolution is too high.

When `method = "both"`, Leiden and Louvain are run at every resolution
for direct comparison. The PCA distance matrix for silhouette is
computed once before the loop, so runtime scales with the number of
resolution values, not the number of cells per iteration.
