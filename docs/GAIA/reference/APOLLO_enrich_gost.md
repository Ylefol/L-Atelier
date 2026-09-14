# Run enrichment analysis using gprofiler2

Performs functional enrichment analysis on gene lists using
gprofiler2::gost(). Supports GO, KEGG, Reactome, WikiPathways, and more
in a single call. Accepts multiple ID types directly (ENSEMBL, SYMBOL,
etc.) - no conversion needed.

## Usage

``` r
APOLLO_enrich_gost(
  gene_lists,
  organism = "hsapiens",
  sources = c("GO:BP", "KEGG", "REAC"),
  user_threshold = 0.05,
  correction_method = "g_SCS",
  domain_scope = "annotated",
  custom_bg = NULL,
  min_term_size = 1,
  max_term_size = 1e+06,
  max_query_size = 10000,
  significant = TRUE,
  exclude_iea = FALSE,
  verbose = TRUE
)
```

## Arguments

- gene_lists:

  Named list of gene vectors. Can be ENSEMBL IDs (with or without
  version numbers), gene symbols, or other identifiers.

- organism:

  gprofiler2 organism code. Common values: "hsapiens" (human),
  "mmusculus" (mouse), "rnorvegicus" (rat). Default: "hsapiens".

- sources:

  Character vector of data sources to query. Options include: "GO:BP",
  "GO:MF", "GO:CC" (Gene Ontology), "KEGG" (KEGG pathways), "REAC"
  (Reactome), "WP" (WikiPathways), "TF" (TRANSFAC transcription
  factors), "MIRNA" (miRTarBase miRNA targets), "CORUM" (protein
  complexes), "HP" (Human Phenotype Ontology). Default: c("GO:BP",
  "KEGG", "REAC").

- user_threshold:

  Significance threshold for term filtering. Default: 0.05.

- correction_method:

  Multiple testing correction method: "g_SCS" (gprofiler's default),
  "fdr", "bonferroni". Default: "g_SCS".

- domain_scope:

  Background for statistical test: "annotated" (genes with any
  annotation), "known" (all known genes), "custom" (provide custom_bg).
  Default: "annotated".

- custom_bg:

  Custom background gene set (when domain_scope = "custom").

- min_term_size:

  Minimum term/pathway size. Default: 1 (no filtering, matches gprofiler
  online defaults).

- max_term_size:

  Maximum term/pathway size. Default: 1e6 (no effective upper limit,
  matches gprofiler online defaults).

- max_query_size:

  Maximum number of genes per query. Gene lists exceeding this are
  skipped. Prevents slow, uninformative enrichment on very large gene
  sets (e.g., WGCNA grey module). Default: 10000. Set to NULL to
  disable.

- significant:

  Logical. Only include significant results in `combined`/`summary`. All
  terms are always fetched from gprofiler2 in a single API call
  regardless of this setting (filtering is applied client-side), so
  toggling this does not incur extra requests. Default: TRUE.

- exclude_iea:

  Logical. Exclude GO terms inferred from electronic annotation.
  Default: FALSE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

A list with class "gost_enrichment" containing:

- results:

  Named list of FULL (unfiltered) gost result objects, one per gene list
  — includes every evaluated term regardless of significance, for use
  with
  [`AETHER_plot_gost_full()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_gost_full.md)

- combined:

  Combined data.frame of `significant`-filtered results with 'query'
  column indicating source

- summary:

  Summary data.frame with counts per gene list and source

- metadata:

  List of analysis parameters

## Details

gprofiler2 automatically detects the input ID type in most cases. For
best results with Ensembl IDs that have version numbers (e.g.,
ENSG00000141510.16), use APOLLO_extract_cluster_genes() with
strip_version = TRUE before calling this function.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with gene symbols
gene_lists <- list(
  module1 = c("TP53", "BRCA1", "EGFR", "MYC"),
  module2 = c("IL6", "TNF", "IL1B", "CXCL8")
)
results <- APOLLO_enrich_gost(gene_lists)

# From WGCNA modules (Ensembl IDs)
module_genes <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE)
results <- APOLLO_enrich_gost(module_genes, sources = c("GO:BP", "GO:MF", "KEGG", "REAC"))

# Mouse data
results <- APOLLO_enrich_gost(gene_lists, organism = "mmusculus")

} # }
```
