# Test for differences in cell-type proportions across samples (propeller)

Wraps
[`speckle::propeller()`](https://rdrr.io/pkg/speckle/man/propeller.html)
to test whether cell-type (cluster) proportions differ between two or
more experimental groups, using sample-level replication. Unlike
differential expression, this only reads `colData(sce)` – no assay
matrix is touched, so the function runs in a fraction of a second
regardless of dataset size.

## Usage

``` r
KERAUNOS_propeller_proportions(
  sce,
  cluster_col = "cluster",
  sample_col,
  group_col,
  cluster_colors = NULL,
  group_colors = NULL,
  transform = c("logit", "asin"),
  trend = FALSE,
  robust = TRUE,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- cluster_col:

  Character. `colData(sce)` column naming the cluster/cell-type
  assignment. Default `"cluster"`.

- sample_col:

  Character. `colData(sce)` column naming the biological
  replicate/sample each cell belongs to. Required – this is the unit
  proportions are computed over.

- group_col:

  Character. `colData(sce)` column naming the experimental
  group/condition to compare (e.g. genotype). Must be constant within
  each `sample_col` level.

- cluster_colors:

  Named character vector, cluster label -\> hex colour, one entry for
  every level of `cluster_col` present in `sce`. Optional – if supplied,
  carried through in `$params` so
  [`ASPIS_plot_composition`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_composition.md)
  uses these colours instead of assigning its own, keeping cluster
  colours consistent with however they're coloured elsewhere (e.g. on a
  UMAP). Default `NULL` (the plot functions fall back to their own
  default palette).

- group_colors:

  Named character vector, group label -\> hex colour, one entry for
  every level of `group_col` present in `sce`. Same purpose as
  `cluster_colors`, consumed by
  [`ASPIS_plot_propeller`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_propeller.md)
  and
  [`ASPIS_plot_composition`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_composition.md)'s
  facet strips. Default `NULL`.

- transform:

  Character. Variance-stabilising transform applied to proportions
  before testing: `"logit"` (default, matches
  [`speckle::propeller()`](https://rdrr.io/pkg/speckle/man/propeller.html)'s
  own default) or `"asin"` (arcsine-square-root).

- trend:

  Logical. Fit a mean-variance trend on the transformed proportions
  (passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)).
  Default `FALSE`, matching the `propeller` default.

- robust:

  Logical. Use robust empirical Bayes shrinkage of variances. Default
  `TRUE`, matching the `propeller` default. Note `propeller` itself
  forces this to `FALSE` internally when fewer than 3 clusters are
  tested.

- fdr_threshold:

  Numeric. FDR cutoff used to populate `$significant`. Default `0.05`.

- verbose:

  Logical. Print `[KERAUNOS]` progress output, including a low-replicate
  warning when any group has fewer than 3 samples. Default `TRUE`.

## Value

A classed list `keraunos_propeller`:

- results:

  Tidy data.frame, one row per cluster: `cluster`, `baseline_prop`,
  per-group mean proportion (`PropMean.<group>`),
  `PropRatio`/`Tstatistic` (two groups) or `Fstatistic` (more than two
  groups), `P.Value`, `FDR`. Ordered by `P.Value`. `baseline_prop` is
  normalised here to a single consistent name – speckle itself calls it
  `BaselineProp.Freq` in the two-group case and `BaselineProp` in the
  \>2-group case.

- significant:

  Subset of `results` with `FDR < fdr_threshold`.

- proportions:

  Tidy data.frame, one row per sample x cluster: `sample`, `cluster`,
  `proportion`, `group`. Feeds
  [`ASPIS_plot_propeller`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_propeller.md)
  so the plot never needs `sce` directly.

- summary:

  One-row data.frame: `n_clusters`, `n_sig`.

- params:

  List of parameters used, including per-group sample counts
  (`n_per_group`) and any supplied `cluster_colors`/ `group_colors`.

## Details

propeller compares proportions at the level of biological samples (one
proportion vector per `sample_col` level), not individual cells –
treating cells as the unit of replication would pseudoreplicate. It
variance-stabilises sample-level proportions (arcsine-square-root or
logit transform), then fits an empirical-Bayes-moderated t-test (two
groups) or F-test/ANOVA (more than two groups) via `limma`, borrowing
variance information across clusters the same way limma borrows
information across genes.

As with any test on few biological replicates, statistical power is
limited when there are only 2-3 samples per group – propeller's variance
moderation makes this more reliable than an unmoderated per-cluster
test, but does not manufacture power the replicate count doesn't
support; treat marginal p-values from very small designs cautiously.
