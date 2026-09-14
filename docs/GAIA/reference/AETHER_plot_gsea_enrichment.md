# Running-Score (Mountain) Plot for a Single GSEA Pathway

Wraps
[`fgsea::plotEnrichment()`](https://rdrr.io/pkg/fgsea/man/plotEnrichment.html)
to draw the running enrichment score curve for one gene set from a GSEA
result
([`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md)).

## Usage

``` r
AETHER_plot_gsea_enrichment(gsea_result, pathway, title = NULL)
```

## Arguments

- gsea_result:

  An `apollo_gsea` object from
  [`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md).
  Must retain `$gene_sets` and `$ranked_genes` (both included by
  default).

- pathway:

  Character. Name of the gene set to plot (must be present in
  `gsea_result$gene_sets`).

- title:

  Character or `NULL`. Plot title. `NULL` (default) uses the pathway
  name.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
gsea_result <- APOLLO_gsea(ranked, collection = "H")
p <- AETHER_plot_gsea_enrichment(gsea_result, "HALLMARK_INFLAMMATORY_RESPONSE")
} # }
```
