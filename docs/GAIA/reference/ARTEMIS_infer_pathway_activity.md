# Infer pathway activity

Convenience wrapper for pathway activity inference using PROGENy
signatures.

## Usage

``` r
ARTEMIS_infer_pathway_activity(
  mat,
  method,
  database = "progeny",
  organism = "human",
  top = 500,
  minsize = 5,
  verbose = TRUE,
  ...
)
```

## Arguments

- mat:

  Numeric matrix of gene expression data. Genes as rows, samples as
  columns.

- method:

  Character. Statistical method to use. REQUIRED - no default.

- database:

  Character. Pathway database. Currently only "progeny" supported.
  Default: "progeny".

- organism:

  Character. Organism: "human" or "mouse". Default: "human".

- top:

  Integer. Number of top responsive genes per pathway. Default: 500.

- minsize:

  Integer. Minimum genes per pathway. Default: 5.

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Additional arguments passed to ARTEMIS_run_decoupler().

## Value

A decoupler_result object with pathway activity scores.

## Details

Note: This provides pathway activity via a different approach than
enrichment analysis (gprofiler2). PROGENy uses pathway-responsive genes
derived from perturbation experiments, while enrichment uses pathway
membership. Both approaches are complementary.

## Examples

``` r
if (FALSE) { # \dontrun{
pathway_activity <- ARTEMIS_infer_pathway_activity(expr_matrix, method = "mlm")

} # }
```
