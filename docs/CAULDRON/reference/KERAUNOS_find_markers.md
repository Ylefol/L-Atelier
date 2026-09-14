# Identify cluster marker genes

Finds differentially expressed marker genes for each cluster using
pairwise statistical tests via
[`scran::findMarkers()`](https://rdrr.io/pkg/scran/man/findMarkers.html).
By default runs pairwise Wilcoxon rank-sum tests and returns genes that
are significantly upregulated in at least one pairwise comparison
(`pval_type = "any"`).

## Usage

``` r
KERAUNOS_find_markers(
  sce,
  cluster_col = "cluster",
  assay_name = "logcounts",
  test_type = c("wilcox", "t", "binom"),
  pval_type = c("any", "some", "all"),
  direction = c("up", "down", "any"),
  min_prop = 0.5,
  block_col = NULL,
  top_n = 10L,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `logcounts` (or other) assay.

- cluster_col:

  Character. `colData` column with cluster labels. Default `"cluster"`.

- assay_name:

  Character. Assay to test. Default `"logcounts"`.

- test_type:

  Character. Test to use: `"wilcox"` (default), `"t"`, or `"binom"`.

- pval_type:

  Character. How to combine pairwise p-values: `"any"` (default — gene
  significant in any comparison), `"some"` (at least `min_prop`
  comparisons), or `"all"` (all comparisons).

- direction:

  Character. Direction of effect: `"up"` (default), `"down"`, or
  `"any"`.

- min_prop:

  Numeric. Minimum proportion of comparisons required when
  `pval_type = "some"`. Default `0.5`.

- block_col:

  Character or `NULL`. `colData` column to use as a blocking factor
  (e.g. sample ID). When provided, each pairwise comparison is
  stratified by block, controlling for batch or sample effects. Default
  `NULL` (no blocking).

- top_n:

  Integer. Maximum number of top markers per cluster in the tidy summary
  table. Default `10L`.

- fdr_threshold:

  Numeric. FDR cut-off for the tidy summary table. Default `0.05`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `keraunos_markers` object (list) with:

- `$markers` — named list of DataFrames from
  [`scran::findMarkers()`](https://rdrr.io/pkg/scran/man/findMarkers.html),
  one per cluster.

- `$top` — tidy data.frame of top markers (columns: cluster, gene, Top,
  p.value, FDR, mean_effect).

- `$params` — list of parameters used.

## Details

The result object contains both the full ranked gene list per cluster
(`$markers`) and a tidy top-N summary table (`$top`).
