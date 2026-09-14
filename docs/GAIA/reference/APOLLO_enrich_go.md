# GO Enrichment Analysis for Annotated Peaks

Performs Gene Ontology enrichment analysis on genes associated with
annotated peaks using clusterProfiler.

## Usage

``` r
APOLLO_enrich_go(
  annotated_peaks,
  org_db = "org.Hs.eg.db",
  ont = "BP",
  gene_id_type = "ENTREZID",
  pval_cutoff = 0.05,
  qval_cutoff = 0.1,
  min_gs_size = 10,
  max_gs_size = 500,
  verbose = TRUE
)
```

## Arguments

- annotated_peaks:

  Data.frame from APOLLO_annotate_peaks() with gene_id column.
  Alternatively, a character vector of gene IDs.

- org_db:

  Character string or OrgDb object. The organism annotation database.
  Default is "org.Hs.eg.db" for human.

- ont:

  Character. GO ontology: "BP" (Biological Process), "MF" (Molecular
  Function), "CC" (Cellular Component), or "ALL" (default = "BP").

- gene_id_type:

  Character. Type of gene IDs provided. One of "ENTREZID", "ENSEMBL",
  "SYMBOL", "REFSEQ" (default = "ENTREZID").

- pval_cutoff:

  Numeric. P-value cutoff for enrichment (default = 0.05).

- qval_cutoff:

  Numeric. Adjusted p-value cutoff (default = 0.1).

- min_gs_size:

  Integer. Minimum gene set size (default = 10).

- max_gs_size:

  Integer. Maximum gene set size (default = 500).

- verbose:

  Logical. Print summary (default = TRUE).

## Value

An enrichResult object from clusterProfiler. Key methods:

- as.data.frame():

  Convert to data.frame of results

- dotplot():

  Create dot plot of top terms

- barplot():

  Create bar plot of top terms

- cnetplot():

  Gene-concept network plot

## Details

If gene IDs are not ENTREZID, the function attempts to convert them
using the org.Db. This is necessary because GO enrichment requires
ENTREZ IDs.

For T2T-CHM13 or other genomes, gene symbols should still work with
org.Hs.eg.db since the gene names are the same.

## Examples

``` r
if (FALSE) { # \dontrun{
# From annotated peaks
annotated <- APOLLO_annotate_peaks(sig_peaks, txdb)
go_results <- APOLLO_enrich_go(annotated)

# View results
head(as.data.frame(go_results))

# Plot
dotplot(go_results, showCategory = 20)

# From gene vector
genes <- c("TP53", "BRCA1", "EGFR")
go_results <- APOLLO_enrich_go(genes, gene_id_type = "SYMBOL")

} # }
```
