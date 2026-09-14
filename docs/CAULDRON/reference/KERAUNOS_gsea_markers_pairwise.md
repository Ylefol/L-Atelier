# Run GSEA between pairs of clusters using marker effect sizes

A marker-score-based alternative to
[`KERAUNOS_gsea_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea_pseudobulk.md)
for comparing specific clusters to each other rather than testing a
condition within each cluster. Ranking statistics come from
[`scran::scoreMarkers()`](https://rdrr.io/pkg/scran/man/scoreMarkers.html)
pairwise effect sizes (via
[`KERAUNOS_rank_pairwise_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_rank_pairwise_markers.md)),
computed at single-cell resolution — unlike DESeq2 pseudobulk DE, this
does not require \\\ge2\\ biological replicate samples per cluster,
making it suitable for comparing two specific clusters directly (e.g.
two developmentally-related or visually-adjacent clusters) even when
sample composition is skewed or a cluster is dominated by one
sample/clone.

## Usage

``` r
KERAUNOS_gsea_markers_pairwise(
  markers_result,
  cluster_pairs = NULL,
  stat = c("AUC", "logFC.cohen", "logFC.detected"),
  gene_sets = NULL,
  species = "Homo sapiens",
  db_species = NULL,
  collection = c("H", "C2", "C5"),
  min_ranked_genes = 10L,
  min_size = 15L,
  max_size = 500L,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- markers_result:

  A `keraunos_markers` object from
  `KERAUNOS_score_markers(sce, ..., full_stats = TRUE)`.

- cluster_pairs:

  `NULL` (default) or a list of length-2 character vectors, e.g.
  `list(c("1","2"), c("8","13"))`, restricting the run to specific
  pairs. `NULL` runs every pairwise combination of clusters present in
  `markers_result`.

- stat:

  Character. Which pairwise effect size to rank on — see
  [`KERAUNOS_rank_pairwise_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_rank_pairwise_markers.md).
  Default `"AUC"`.

- gene_sets:

  Named list of character vectors, or `NULL` to fetch from MSigDB via
  msigdbr (fetched once and reused across all pairs).

- species:

  Character. Species for MSigDB gene set retrieval. Default
  `"Homo sapiens"`.

- db_species:

  Character or `NULL`. MSigDB database to query. `NULL` (default): human
  gene sets; `"MM"`: native mouse (use `M`-prefixed collection codes,
  e.g. `c("MH", "M2", "M5")`). See
  [`KERAUNOS_list_gene_sets`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_list_gene_sets.md)
  to browse available codes.

- collection:

  Character vector. MSigDB collection codes. Default
  `c("H", "C2", "C5")`. Subcollection can be appended with a colon, e.g.
  `"C2:CP:REACTOME"`.

- min_ranked_genes:

  Integer. Minimum number of non-`NA` ranked genes required to run GSEA
  for a pair. Pairs below this threshold are skipped. Default `10L`.

- min_size:

  Integer. Minimum gene set size after overlap with ranked genes.
  Default `15L`.

- max_size:

  Integer. Maximum gene set size. Default `500L`.

- fdr_threshold:

  Numeric. FDR threshold for `$significant`. Default `0.05`.

- verbose:

  Logical. Print per-pair progress and a summary table. Default `TRUE`.

## Value

A `keraunos_gsea_pairwise` object (also inherits `keraunos_gsea_multi`,
so it is directly compatible with
[`TALARIA_export_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gsea.md))
with:

- `$results` — named list of `keraunos_gsea` objects, one per pair run
  (names are `"<cluster1>_vs_<cluster2>"`).

- `$significant` — named list of significant-only data.frames, one per
  pair.

- `$summary` — data.frame: comparison / cluster1 / cluster2 / n_ranked /
  n_sets / n_sig / n_enriched / n_depleted.

- `$skipped` — character vector of pair labels skipped.

- `$params` — parameters used.
