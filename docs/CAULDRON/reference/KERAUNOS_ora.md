# Over-representation analysis via gprofiler2

Tests whether a set of genes is enriched for annotated biological terms
using
[`gprofiler2::gost()`](https://rdrr.io/pkg/gprofiler2/man/gost.html).
Supports GO terms, KEGG, Reactome, and other databases available through
g:Profiler.

## Usage

``` r
KERAUNOS_ora(
  genes,
  background = NULL,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  fdr_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- genes:

  Character vector of gene symbols (or Ensembl IDs, Entrez IDs, etc. —
  whatever gprofiler2 accepts for the chosen organism).

- background:

  Character vector of background genes, or `NULL` (default) to use the
  gprofiler2 reference genome as background.

- organism:

  Character. gprofiler2 organism code. Default `"hsapiens"`. Use
  `"mmusculus"` for mouse, etc. See
  [`gprofiler2::gost()`](https://rdrr.io/pkg/gprofiler2/man/gost.html)
  for full list.

- sources:

  Character vector. Databases to query. Default
  `c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")`.

- fdr_threshold:

  Numeric. FDR threshold for `$significant`. Default `0.05`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `keraunos_ora` object (list) with:

- `$results` — full gprofiler2 result data.frame (ordered by p-value).

- `$significant` — filtered to `fdr_threshold`.

- `$params` — list of parameters used.
