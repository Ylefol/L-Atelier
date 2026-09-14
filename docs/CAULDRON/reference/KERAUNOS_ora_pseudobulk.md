# Run ORA across all clusters from pseudobulk DE results

A convenience wrapper around
[`KERAUNOS_ora`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_ora.md)
that accepts a `keraunos_pseudobulk` object and runs ORA separately on
upregulated and downregulated significant genes for every cluster.

## Usage

``` r
KERAUNOS_ora_pseudobulk(
  de_result,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  lfc_threshold = 0,
  background = c("tested", "genome"),
  fdr_threshold = 0.05,
  min_sig_genes = 5L,
  verbose = TRUE
)
```

## Arguments

- de_result:

  A `keraunos_pseudobulk` object (output of
  [`KERAUNOS_de_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_de_pseudobulk.md)).

- organism:

  Character. gprofiler2 organism code. Default `"hsapiens"`. Use
  `"mmusculus"` for mouse.

- sources:

  Character vector. Databases to query. Default
  `c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")`.

- lfc_threshold:

  Numeric. Minimum absolute \\\log_2\\ fold-change for a gene to be
  considered directionally significant. Genes with \\\|LFC\| \leq
  lfc\\threshold\\ are excluded from both up and down sets. Default `0`
  (all significant genes included).

- background:

  Character. Background gene set to pass to gprofiler2. `"tested"`
  (default) — all genes tested in the cluster (statistically correct for
  pseudobulk); `"genome"` — gprofiler2 reference genome (no custom
  background).

- fdr_threshold:

  Numeric. FDR threshold for `$significant` ORA terms. Default `0.05`.

- min_sig_genes:

  Integer. Minimum number of significant DE genes required to run ORA
  for a given cluster/direction. Default `5L`.

- verbose:

  Logical. Print per-cluster progress and a summary table. Default
  `TRUE`.

## Value

A `keraunos_ora_multi` object (list) with:

- `$results` — named list (cluster → list with `$up` and `$down`, each a
  `keraunos_ora` object or `NULL` if skipped).

- `$significant` — named list (cluster → list with `$up` and `$down`
  significant-only data.frames).

- `$summary` — data.frame: cluster / n_sig_genes / n_up_genes /
  n_down_genes / n_sig_terms_up / n_sig_terms_down.

- `$skipped` — character vector of cluster names where neither direction
  had enough genes.

- `$params` — parameters used.

## Details

Significance is determined using whichever threshold was applied when
[`KERAUNOS_de_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_de_pseudobulk.md)
was run (FDR or raw p-value, stored in `de_result$params`). Direction is
determined by `log2FoldChange`: genes with \\LFC \> lfc\\threshold\\ are
"up"; genes with \\LFC \< -lfc\\threshold\\ are "down".

Clusters, or directions within a cluster, with fewer than
`min_sig_genes` significant genes are skipped with a message.
