# Run GSEA across all clusters from pseudobulk DE results

A convenience wrapper around
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)
that accepts a `keraunos_pseudobulk` object and runs GSEA for every
cluster automatically. Gene sets are fetched **once** from MSigDB (not
per cluster), making multi-cluster runs efficient.

## Usage

``` r
KERAUNOS_gsea_pseudobulk(
  de_result,
  gene_sets = NULL,
  species = "Homo sapiens",
  db_species = NULL,
  collection = c("H", "C2", "C5"),
  rank_by = c("stat", "lfc_pval"),
  min_ranked_genes = 10L,
  min_size = 15L,
  max_size = 500L,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- de_result:

  A `keraunos_pseudobulk` object (output of
  [`KERAUNOS_de_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_de_pseudobulk.md)).

- gene_sets:

  Named list of character vectors, or `NULL` (default) to fetch from
  MSigDB via `msigdbr`. Fetched once and reused across all clusters.

- species:

  Character. Species for MSigDB gene set retrieval. Default
  `"Homo sapiens"`.

- db_species:

  Character or `NULL`. MSigDB database to query. Use `"MM"` for native
  mouse gene sets (requires M-prefixed collection codes, e.g.
  `c("MH", "M2", "M5")`). See
  [`KERAUNOS_list_gene_sets`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_list_gene_sets.md)
  to browse available codes.

- collection:

  Character vector. MSigDB collection codes. Default
  `c("H", "C2", "C5")`. Subcategory can be appended with a colon, e.g.
  `"C2:CP:REACTOME"`.

- rank_by:

  Character. Ranking statistic to derive from DESeq2 results. `"stat"`
  (default) or `"lfc_pval"`.

- min_ranked_genes:

  Integer. Minimum number of non-`NA` ranked genes required to run GSEA
  for a cluster. Clusters below this threshold are skipped. Default
  `10L`.

- min_size:

  Integer. Minimum gene set size after overlap with ranked genes.
  Default `15L`.

- max_size:

  Integer. Maximum gene set size. Default `500L`.

- fdr_threshold:

  Numeric. FDR threshold for `$significant`. Default `0.05`.

- verbose:

  Logical. Print per-cluster progress and a summary table. Default
  `TRUE`.

## Value

A `keraunos_gsea_multi` object (list) with:

- `$results` — named list of `keraunos_gsea` objects, one per cluster
  that was run.

- `$significant` — named list of significant-only data.frames, one per
  cluster.

- `$summary` — data.frame: cluster / n_ranked / n_sets / n_sig /
  n_enriched / n_depleted.

- `$skipped` — character vector of cluster names skipped.

- `$params` — parameters used.

## Details

The ranking statistic is derived from the DESeq2 result columns:

- `"stat"`:

  DESeq2 Wald test statistic (default). Most informative ranking;
  captures both magnitude and significance.

- `"lfc_pval"`:

  \\\mathrm{sign}(\log_2 FC) \times -\log\_{10}(p\text{-value})\\.
  Useful when `stat` has many `NA`s (e.g. genes filtered by DESeq2
  independent filtering).

Genes with `NA` ranking values are removed before the GSEA call.
Clusters with fewer than `min_ranked_genes` valid ranking values are
skipped with a message.
