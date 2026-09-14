# List species available in MSigDB

Wraps
[`msigdbr::msigdbr_species()`](https://igordot.github.io/msigdbr/reference/msigdbr_species.html)
and prints a formatted table of available species names and their
two-letter `db_species` codes. Use this to find the correct `species`
and `db_species` values for
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)
and
[`KERAUNOS_list_gene_sets`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_list_gene_sets.md).

## Usage

``` r
KERAUNOS_list_species(verbose = TRUE)
```

## Arguments

- verbose:

  Logical. Print the formatted table. Default `TRUE`.

## Value

Invisibly returns the species `data.frame` from
[`msigdbr::msigdbr_species()`](https://igordot.github.io/msigdbr/reference/msigdbr_species.html).
