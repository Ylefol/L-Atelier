# Pseudobulk differential expression via DESeq2

Aggregates raw counts per sample per cluster (pseudobulk), then runs
`DESeq2` on each cluster's pseudobulk matrix to identify genes
differentially expressed between two conditions.

## Usage

``` r
KERAUNOS_de_pseudobulk(
  sce,
  condition_col,
  sample_col,
  cluster_col = "cluster",
  assay_name = "counts",
  reference_level = NULL,
  min_cells = 10L,
  min_samples = 2L,
  fdr_threshold = 0.05,
  pval_threshold = NULL,
  lfc_threshold = 0,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a raw `counts` assay.

- condition_col:

  Character. `colData` column with condition labels (must have exactly 2
  levels).

- sample_col:

  Character. `colData` column with sample/replicate IDs used for
  pseudobulk aggregation (one pseudobulk column per sample).

- cluster_col:

  Character. `colData` column with cluster labels. Default `"cluster"`.

- assay_name:

  Character. Assay to aggregate. Must be raw integer counts. Default
  `"counts"`.

- reference_level:

  Character or `NULL`. Which condition level is the reference
  (denominator) in the DESeq2 contrast. When `NULL` (default), the
  alphabetically first level is used.

- min_cells:

  Integer. Minimum cells a sample must contribute to a cluster to be
  included as a pseudobulk replicate. Samples below this threshold are
  dropped. Default `10L`.

- min_samples:

  Integer. Minimum number of replicates required per condition per
  cluster after `min_cells` filtering. Clusters not meeting this
  requirement are skipped. Default `2L`.

- fdr_threshold:

  Numeric or `NULL`. FDR (adjusted p-value) threshold used to define
  significant genes in the summary table and for DESeq2 independent
  filtering. Default `0.05`. Set to `NULL` when using `pval_threshold`
  instead.

- pval_threshold:

  Numeric or `NULL`. Raw p-value threshold. Only used when
  `fdr_threshold = NULL`. Intended for niche cases where FDR adjustment
  is not appropriate (e.g. very few tests). Providing both
  `fdr_threshold` and `pval_threshold` raises an error. Default `NULL`.

- lfc_threshold:

  Numeric. Absolute log2 fold-change threshold passed to
  [`DESeq2::results()`](https://rdrr.io/pkg/DESeq2/man/results.html).
  Default `0` (no LFC filter).

- verbose:

  Logical. Print per-cluster progress and a summary. Default `TRUE`.

## Value

A `keraunos_pseudobulk` object (list) with:

- `$results` — named list of data.frames, one per cluster (columns:
  gene, cluster, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj).

- `$dds` — named list of fitted `DESeqDataSet` objects (for downstream
  use, e.g.
  [`DESeq2::lfcShrink()`](https://rdrr.io/pkg/DESeq2/man/lfcShrink.html)).

- `$summary` — data.frame: cluster \| n_tested \| n_sig \| n_up \|
  n_down.

- `$skipped` — character vector of cluster names that were skipped.

- `$params` — list of parameters used.

## Details

This is the statistically rigorous multi-sample approach: each
biological replicate contributes one pseudobulk column, eliminating
pseudoreplication. Requires at least `min_samples` replicates per
condition per cluster after cell-count filtering.
