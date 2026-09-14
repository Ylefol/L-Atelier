# Gene Set Enrichment Analysis (fgsea)

Runs
[`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html)
on a pre-ranked gene list against a collection of gene sets. Gene sets
can be supplied as a named list or fetched automatically from MSigDB via
`msigdbr`. Use
[`APOLLO_rank_from_de()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_rank_from_de.md)
to build `ranked_genes` from an `artemis_limma` or `artemis_ts_de`
result.

## Usage

``` r
APOLLO_gsea(
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
  sets); values are the ranking statistic (e.g. limma's `t`, DESeq2's
  `stat`, or signed -log10(p)). Must have names.

- gene_sets:

  Named list of character vectors (gene set name -\> member genes), or
  `NULL` to fetch from MSigDB via `msigdbr`.

- species:

  Character. Species name passed to
  [`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  when `gene_sets = NULL`. Default: `"Homo sapiens"`.

- db_species:

  Character or `NULL`. MSigDB species database to query. `NULL`
  (default) lets msigdbr use its default (human-centric, with orthologs
  for non-human species). Use `"MM"` to query the native mouse database
  (gene sets defined in mouse studies with native Mus musculus gene
  symbols, not human orthologs). Only relevant when `gene_sets = NULL`.

- collection:

  Character vector. MSigDB collection codes to fetch when
  `gene_sets = NULL`. Default: `c("H", "C2", "C5")`. A subcollection can
  be appended with a colon (e.g. `"C2:CP:REACTOME"`, `"C5:GO:BP"`).
  Mouse-native collections (`db_species = "MM"`) use `M`-prefixed codes
  (e.g. `"MH"`, `"M2"`, `"M5"`).

- min_size:

  Integer. Minimum gene set size (after overlap with ranked genes).
  Default: 15.

- max_size:

  Integer. Maximum gene set size. Default: 500.

- fdr_threshold:

  Numeric. FDR threshold for `$significant`. Default: 0.05.

- verbose:

  Logical. Print progress and summary. Default: TRUE.

## Value

An S3 object of class `"apollo_gsea"` containing:

- results:

  Full fgsea result data.frame (ordered by padj): pathway, pval, padj,
  ES, NES, size, leadingEdge (collapsed string), collection (if gene
  sets came from msigdbr).

- significant:

  `results` filtered to `fdr_threshold`.

- ranked_genes:

  The input vector, NA-removed and sorted descending.

- gene_sets:

  Named list of only the gene sets actually tested (passed min/max
  size) - used by
  [`AETHER_plot_gsea_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_gsea_enrichment.md).

- params:

  List of parameters used.

## Details

Requires the `fgsea` package (Suggests). When `gene_sets = NULL`, also
requires `msigdbr` (Suggests).

## Examples

``` r
if (FALSE) { # \dontrun{
ranked <- APOLLO_rank_from_de(de, rank_by = "t")
gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))

# Custom gene sets
my_sets <- list(my_pathway = c("TP53", "BRCA1", "EGFR"))
gsea_result <- APOLLO_gsea(ranked, gene_sets = my_sets)
} # }
```
