# KEGG Pathway Enrichment Analysis

Performs KEGG pathway enrichment analysis on genes associated with
annotated peaks using clusterProfiler.

## Usage

``` r
APOLLO_enrich_kegg(
  annotated_peaks,
  org_db = "org.Hs.eg.db",
  organism = "hsa",
  gene_id_type = "ENTREZID",
  pval_cutoff = 0.05,
  qval_cutoff = 0.1,
  verbose = TRUE
)
```

## Arguments

- annotated_peaks:

  Data.frame from APOLLO_annotate_peaks() with gene_id column.
  Alternatively, a character vector of gene IDs.

- org_db:

  Character string or OrgDb object. The organism annotation database for
  ID conversion. Default is "org.Hs.eg.db" for human.

- organism:

  Character. KEGG organism code (default = "hsa" for human).

- gene_id_type:

  Character. Type of gene IDs provided. One of "ENTREZID", "ENSEMBL",
  "SYMBOL", "REFSEQ" (default = "ENTREZID").

- pval_cutoff:

  Numeric. P-value cutoff for enrichment (default = 0.05).

- qval_cutoff:

  Numeric. Adjusted p-value cutoff (default = 0.1).

- verbose:

  Logical. Print summary (default = TRUE).

## Value

An enrichResult object from clusterProfiler.

## Examples

``` r
if (FALSE) { # \dontrun{
annotated <- APOLLO_annotate_peaks(sig_peaks, txdb)
kegg_results <- APOLLO_enrich_kegg(annotated)
dotplot(kegg_results)

} # }
```
