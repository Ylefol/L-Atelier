# Gene set enrichment analysis via fgsea

Runs
[`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html)
on a pre-ranked gene list against a collection of gene sets. Gene sets
can be supplied as a named list or fetched automatically from MSigDB via
`msigdbr`.

## Usage

``` r
KERAUNOS_gsea(
  ranked_genes,
  gene_sets = NULL,
  species = "Homo sapiens",
  db_species = NULL,
  collection = c("H", "C2", "C5"),
  min_size = 15L,
  max_size = 500L,
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- ranked_genes:

  Named numeric vector. Names are gene symbols (or IDs matching the gene
  sets); values are the ranking statistic (e.g. signed
  \\-\log\_{10}(p)\\ \\\times\\ sign(LFC), or the DESeq2 `stat` column).
  Must have no missing names.

- gene_sets:

  Named list of character vectors (gene set name → member genes), or
  `NULL` to fetch from MSigDB via `msigdbr`.

- species:

  Character. Species name passed to
  [`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  when `gene_sets = NULL`. Default `"Homo sapiens"`.

- db_species:

  Character or `NULL`. MSigDB species database to query. `NULL`
  (default) lets msigdbr use its default (human-centric, with orthologs
  for non-human species). Use `"MM"` to query the native mouse database
  (gene sets defined in mouse studies with native Mus musculus gene
  symbols, not human orthologs). Only relevant when `gene_sets = NULL`.

- collection:

  Character vector. MSigDB collection codes to fetch when
  `gene_sets = NULL`. Default `c("H", "C2", "C5")`. Subcategory can be
  appended with a colon (e.g. `"C2:CP:REACTOME"`).

- min_size:

  Integer. Minimum gene set size (after overlap with ranked genes).
  Default `15L`.

- max_size:

  Integer. Maximum gene set size. Default `500L`.

- fdr_threshold:

  Numeric. FDR threshold for `$significant`. Default `0.05`.

- verbose:

  Logical. Print progress and summary. Default `TRUE`.

## Value

A `keraunos_gsea` object (list) with:

- `$results` — full fgsea result data.frame (ordered by padj).

- `$significant` — filtered to `fdr_threshold`.

- `$params` — list of parameters used.
