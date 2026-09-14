# Sweep tSNE perplexity

Perplexity is tSNE's primary parameter: it loosely controls how many
neighbours each cell considers when constructing its local probability
distribution. Low perplexity values emphasise very local structure and
can fragment continuous populations into disconnected islands; high
values integrate broader context but may merge distinct clusters.
Typical effective ranges are 5–50 for most single-cell datasets. This
function tests a range of perplexity values and scores each embedding on
two complementary metrics:

## Usage

``` r
TALOS_tune_tsne(
  sce,
  perplexity_range = c(10, 20, 30, 50, 100),
  use_rep = "PCA",
  n_pcs = 30L,
  max_iter = 1000L,
  knn_k = 15L,
  cluster_col = "cluster",
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with `"PCA"` in `reducedDims` (run
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)
  first).

- perplexity_range:

  Numeric vector of perplexity values to test. Values exceeding \\N/3\\
  are removed automatically. Default `c(10, 20, 30, 50, 100)`.

- n_pcs:

  Integer. PCA components used for tSNE input and KNN reference. Default
  `30`.

- max_iter:

  Integer. tSNE iterations per run (fixed; not a tuning target). Default
  `1000`.

- knn_k:

  Integer. Nearest neighbours for the KNN overlap metric. Default `15`.

- cluster_col:

  Character. `colData` column for silhouette labels. Default
  `"cluster"`.

- seed:

  Integer. Random seed. Default `42L`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `talos_embedding_sweep` list with:

- `results`:

  Data frame with one row per perplexity value: `perplexity`,
  `knn_overlap`, `mean_sil` (`NA` if cluster labels absent),
  `composite`.

- `embeddings`:

  List of cells × 2 matrices — one tSNE embedding per perplexity value,
  in the same order as `results`. Passed directly to
  [`ASPIS_plot_embedding_grid()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_embedding_grid.md)
  for visual comparison.

- `plot`:

  Multi-panel line plot (KNN overlap, mean silhouette if available,
  composite vs perplexity), with a red dashed vertical line at the
  suggested value.

- `best_params`:

  Named list (`perplexity`) for the value with the highest composite
  score.

- `params`:

  List of sweep parameters for reproducibility.

## Details

- KNN overlap:

  Mean Jaccard similarity between each cell's `k` nearest neighbours in
  PCA space and in the tSNE embedding. Measures how faithfully local
  neighbourhood structure is preserved in the 2D layout. A shared KNN
  reference in PCA space is computed once and reused across all
  perplexity values. Ranges from 0 (no overlap) to 1 (identical
  neighbourhoods).

- Mean silhouette width:

  Mean silhouette width in 2D tSNE space using existing cluster labels
  from `cluster_col`. Measures how visually well-separated the clusters
  are in the final embedding. Ranges from \\-1\\ (misassigned) to \\+1\\
  (well-separated). Skipped when `cluster_col` is absent from `colData`.

Both metrics are normalised to \\\[0, 1\]\\ and averaged into a
composite score. The perplexity with the highest composite score is
suggested. Perplexity values exceeding \\N / 3\\ (where \\N\\ is the
number of cells) are invalid for tSNE and are silently removed before
sweeping.

This function returns the same `talos_embedding_sweep` class as
[`TALOS_tune_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_umap.md),
so `print` and downstream helpers work identically for both.
