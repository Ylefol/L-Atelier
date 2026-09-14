# Identify hub genes in modules

Identifies hub genes - genes that are highly connected within their
module and/or highly correlated with traits of interest.

## Usage

``` r
ARTEMIS_wgcna_hub_genes(
  modules,
  gene_sig = NULL,
  trait_name = NULL,
  n_top = 10,
  mm_threshold = 0.8,
  gs_threshold = 0.2,
  verbose = TRUE
)
```

## Arguments

- modules:

  A wgcna_modules object from ARTEMIS_wgcna_detect_modules().

- gene_sig:

  Optional. A wgcna_gene_sig object. If provided, uses GS to identify
  trait-relevant hubs.

- trait_name:

  Trait name to use for trait-relevant hub identification. Required if
  gene_sig is provided and has multiple traits.

- n_top:

  Number of top hub genes to return per module. Default: 10. Set to NULL
  to return all genes above mm_threshold (useful for enrichment analysis
  where a kME threshold is more appropriate than a fixed count).

- mm_threshold:

  Module membership threshold for hub genes. Default: 0.8.

- gs_threshold:

  Gene significance threshold (if using trait). Default: 0.2.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_hubs" containing:

- hub_genes:

  Data.frame of all hub genes across modules

- by_module:

  List of hub genes per module

- hub_summary:

  Summary of hub genes per module

- criteria:

  Criteria used for hub identification

## Details

Hub genes are identified based on:

- Module Membership (MM): How well the gene's expression pattern matches
  the module eigengene. High MM = gene is central to module.

- Gene Significance (GS, optional): Correlation with trait of interest.
  High GS + High MM = biologically relevant hub.

## Examples

``` r
if (FALSE) { # \dontrun{
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
hubs <- ARTEMIS_wgcna_hub_genes(modules)

# With trait information
gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)
hubs <- ARTEMIS_wgcna_hub_genes(modules, gene_sig, trait_name = "age")

} # }
```
