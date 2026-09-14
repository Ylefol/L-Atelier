# Sweep HVG count to find a data-driven optimum

Tests a range of `n_hvgs` values and, for each, scores the resulting
gene set on two complementary metrics:

## Usage

``` r
PYRI_tune_hvg(
  sce,
  n_hvgs_range = c(500L, 1000L, 2000L, 3000L, 5000L, 7500L),
  n_pcs = 20L,
  scale = FALSE,
  block_col = NULL,
  assay_name = "logcounts",
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"logcounts"` assay
  ([`PYRI_normalize`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_normalize.md)
  must have been run first).

- n_hvgs_range:

  Integer vector of candidate HVG counts to test. Values exceeding
  `nrow(sce)` are silently dropped. Default
  `c(500, 1000, 2000, 3000, 5000, 7500)`.

- n_pcs:

  Integer. Number of PCA components used to measure variance explained
  at each step. Smaller values are faster; `20` is usually sufficient
  for a sweep. Default `20`.

- scale:

  Logical. Whether to scale genes to unit variance before PCA, matching
  [`TALOS_run_pca`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)'s
  `scale` argument. Default `FALSE`.

- block_col:

  Character. Column in `colData` for per-sample blocking in the
  mean-variance model (matches
  [`PYRI_select_hvg`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_select_hvg.md)).
  `NULL` (default) fits one trend across all cells.

- assay_name:

  Character. Assay to model. Default `"logcounts"`.

- seed:

  Integer. Random seed for irlba reproducibility. Default `42L`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `pyri_hvg_sweep` list with:

- `results`:

  Data frame: `n_hvgs`, `cum_var`, `mean_bio_var`.

- `plot`:

  Two-panel ggplot (cum var + mean bio var vs n_hvgs).

- `best_n_hvgs`:

  Suggested HVG count (elbow in cum_var).

- `params`:

  Sweep parameters.

## Details

1.  **Cumulative variance explained** — the percentage of total
    gene-expression variance captured by the top `n_pcs` principal
    components computed from the selected HVGs. Increasing `n_hvgs`
    initially raises this sharply, then flattens as noise genes dilute
    the signal. The elbow of this curve (detected via the kneedle
    method) is used as the suggested optimum.

2.  **Mean biological variance** — the mean of the scran biological
    variance component across selected genes. This falls monotonically
    as lower-variance genes enter the set and serves as a sanity check:
    the suggested `n_hvgs` should lie in the region where mean bio
    variance is still meaningfully above zero.

The mean-variance model is fitted **once** across all genes; the sweep
only changes the cutoff on the pre-ranked list, so runtime is dominated
by the repeated lightweight PCA steps, not by modelling. PCA is
performed with
[`irlba::irlba`](https://rdrr.io/pkg/irlba/man/irlba.html) directly
(BPCells-aware — no full matrix materialisation required).
