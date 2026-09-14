# Cell-level contrast between groups (quick path)

Tests for differentially expressed genes between two (or more) groups at
the single-cell level using
[`scran::findMarkers()`](https://rdrr.io/pkg/scran/man/findMarkers.html).
Can run globally or within each cluster independently.

## Usage

``` r
KERAUNOS_contrast_groups(
  sce,
  group_col,
  cluster_col = NULL,
  assay_name = "logcounts",
  test_type = c("wilcox", "t", "binom"),
  direction = c("any", "up", "down"),
  fdr_threshold = 0.05,
  top_n = 10L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `logcounts` (or other) assay.

- group_col:

  Character. `colData` column defining the groups to contrast (must have
  \\\geq 2\\ levels).

- cluster_col:

  Character or `NULL`. If `NULL` (default), a single global contrast is
  performed. If a `colData` column name is provided, the contrast is run
  independently within each cluster.

- assay_name:

  Character. Assay to test. Default `"logcounts"`.

- test_type:

  Character. Test: `"wilcox"` (default), `"t"`, or `"binom"`.

- direction:

  Character. Direction of effect: `"any"` (default), `"up"`, or
  `"down"`.

- fdr_threshold:

  Numeric. FDR cut-off for the tidy summary. Default `0.05`.

- top_n:

  Integer. Maximum top hits per group/cluster in the tidy table. Default
  `10L`.

- verbose:

  Logical. Print a summary and pseudoreplication warning. Default
  `TRUE`.

## Value

A `keraunos_contrast` object (list) with:

- `$results` — named list of findMarkers DataFrames (one per group if
  global; named list `[[cluster]][[group]]` if clustered).

- `$top` — tidy data.frame of top hits.

- `$params` — list of parameters used.

## Details

**Important:** This is a cell-level test. With multiple biological
samples per group, individual cells from the same sample are not
independent — this is the pseudoreplication problem. For a statistically
rigorous multi-sample comparison, use
[`KERAUNOS_de_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_de_pseudobulk.md)
instead. Use this function for exploratory analysis or single-sample
datasets.
