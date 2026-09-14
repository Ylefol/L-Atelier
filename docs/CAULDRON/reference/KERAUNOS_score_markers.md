# Annotation-focused cluster marker scoring

Scores cluster markers using pairwise effect sizes via
[`scran::scoreMarkers()`](https://rdrr.io/pkg/scran/man/scoreMarkers.html).
Unlike
[`KERAUNOS_find_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_find_markers.md),
this function does not compute p-values; instead it returns AUC and
Cohen's d summaries across all pairwise cluster comparisons. This is the
recommended approach when the goal is cell type annotation rather than
hypothesis testing, as p-values from single-cell tests are heavily
inflated due to non-independence of cells.

## Usage

``` r
KERAUNOS_score_markers(
  sce,
  cluster_col = "cluster",
  assay_name = "logcounts",
  block_col = NULL,
  restrict_to = NULL,
  min_auc = 0.5,
  top_n = 10L,
  full_stats = FALSE,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `logcounts` (or other) assay.

- cluster_col:

  Character. `colData` column with cluster labels. Default `"cluster"`.

- assay_name:

  Character. Assay to score. Default `"logcounts"`.

- block_col:

  Character or `NULL`. `colData` column to use as a blocking factor
  (e.g. sample ID). Stratifies each pairwise comparison by block,
  controlling for batch or sample effects. Default `NULL`.

- restrict_to:

  Character vector or `NULL`. When provided, only genes present in this
  set are considered when building the tidy `$top` table. The full
  `$markers` list is unaffected (all genes are scored). Use
  [`KERAUNOS_fetch_marker_genes`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_fetch_marker_genes.md)
  to obtain a validated whitelist from MSigDB C8 or PanglaoDB, which
  filters out uninformative features such as lncRNAs, ribosomal genes,
  and ubiquitous housekeeping genes. Default `NULL` (no restriction).

- min_auc:

  Numeric. Minimum `median.AUC` for a gene to appear in the tidy summary
  table. AUC of 0.5 = no better than random; 1.0 = perfect classifier.
  Default `0.5`.

- top_n:

  Integer. Maximum number of top markers per cluster in the tidy summary
  table, ranked by `rank.AUC`. Default `10L`.

- full_stats:

  Logical. Passed through to
  [`scran::scoreMarkers()`](https://rdrr.io/pkg/scran/man/scoreMarkers.html)'s
  `full.stats` argument. When `TRUE`, each per-cluster DataFrame in
  `$markers` additionally carries a `full.AUC` (and `full.logFC.cohen`,
  `full.logFC.detected`) nested column giving the per-gene effect size
  against every *individual* other cluster, not just the
  `median`/`mean`/`min`/`max` summary across all of them. Required by
  [`KERAUNOS_rank_pairwise_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_rank_pairwise_markers.md)
  to build a ranked gene list for one specific pair of clusters (e.g.
  for GSEA between two clusters without needing pseudobulk replicate
  samples). Increases result object size roughly in proportion to the
  number of clusters. Default `FALSE` (unchanged from previous
  behaviour).

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `keraunos_markers` object (list) with:

- `$markers` — named list of DataFrames from
  [`scran::scoreMarkers()`](https://rdrr.io/pkg/scran/man/scoreMarkers.html),
  one per cluster. Each DataFrame contains per-gene effect size
  summaries (`median.AUC`, `rank.AUC`, `median.logFC.cohen`, etc.), plus
  `full.*` pairwise columns when `full_stats = TRUE`.

- `$top` — tidy data.frame of top markers (columns: cluster, gene, Top,
  median_auc, mean_effect, FDR). `Top = rank.AUC`; `FDR = NA` (no
  p-values).

- `$params` — list of parameters used.

## Details

Genes are ranked by `rank.AUC` (the minimum rank across all pairwise AUC
comparisons — lower is better) and filtered by `min_auc`
(`median.AUC >= min_auc` across pairwise comparisons, where 0.5 = random
and 1.0 = perfect classifier).

The returned object is a `keraunos_markers` and is fully compatible with
[`ASPIS_plot_marker_dotplot()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_marker_dotplot.md)
and
[`ASPIS_plot_marker_heatmap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_marker_heatmap.md).
