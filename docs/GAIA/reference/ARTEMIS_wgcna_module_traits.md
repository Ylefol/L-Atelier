# Correlate modules with traits

Calculates correlations between module eigengenes and sample traits.
Identifies which modules are associated with which clinical/phenotypic
variables.

## Usage

``` r
ARTEMIS_wgcna_module_traits(
  modules,
  traits,
  cor_method = "pearson",
  p_adjust = "BH",
  verbose = TRUE
)
```

## Arguments

- modules:

  A wgcna_modules object from ARTEMIS_wgcna_detect_modules().

- traits:

  Data.frame of traits, or wgcna_data object (uses datTraits). Must have
  samples as rows matching the module data.

- cor_method:

  Correlation method: "pearson" or "spearman". Default: "pearson".

- p_adjust:

  Method for p-value adjustment: "BH", "bonferroni", "none", etc.
  Default: "BH" (Benjamini-Hochberg).

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_trait_cor" containing:

- cor_matrix:

  Matrix of correlations (modules x traits)

- pvalue_matrix:

  Matrix of p-values

- padj_matrix:

  Matrix of adjusted p-values

- significant_associations:

  Data.frame of significant module-trait pairs

- n_significant:

  Number of significant associations

## Examples

``` r
if (FALSE) { # \dontrun{
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = 6)
trait_cor <- ARTEMIS_wgcna_module_traits(modules, wgcna_data)

} # }
```
