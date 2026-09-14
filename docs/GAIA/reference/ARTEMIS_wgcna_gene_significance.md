# Calculate gene significance for traits

Computes gene significance (GS) and module membership (MM) for each
gene. Gene significance measures the correlation of each gene with a
trait. Module membership measures how well each gene fits its assigned
module.

## Usage

``` r
ARTEMIS_wgcna_gene_significance(
  modules,
  traits,
  trait_name = NULL,
  cor_method = "pearson",
  verbose = TRUE
)
```

## Arguments

- modules:

  A wgcna_modules object from ARTEMIS_wgcna_detect_modules().

- traits:

  Data.frame of traits, or wgcna_data object.

- trait_name:

  Name of specific trait to analyze. If NULL, analyzes all traits.

- cor_method:

  Correlation method: "pearson" or "spearman". Default: "pearson".

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_gene_sig" containing:

- gene_info:

  Data.frame with gene significance and module membership

- by_module:

  List of data.frames, one per module

- top_genes_per_module:

  Top significant genes for each module

- trait_names:

  Names of traits analyzed

## Details

For each gene, calculates:

- Gene Significance (GS): correlation between gene expression and trait

- Module Membership (MM): correlation between gene expression and module
  eigengene

High GS + High MM = genes that are both important for the trait AND
central to their module (good biomarker candidates).

## Examples

``` r
if (FALSE) { # \dontrun{
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data)

# For specific trait
gene_sig <- ARTEMIS_wgcna_gene_significance(modules, wgcna_data,
                                             trait_name = "age")

} # }
```
