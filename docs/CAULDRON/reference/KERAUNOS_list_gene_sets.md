# List available MSigDB gene set collections

Wraps
[`msigdbr::msigdbr_collections()`](https://igordot.github.io/msigdbr/reference/msigdbr_collections.html)
and prints a formatted table grouped by database (human vs mouse). Use
this to discover which collection codes to pass to
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md).

## Usage

``` r
KERAUNOS_list_gene_sets(db_species = NULL, verbose = TRUE)
```

## Arguments

- db_species:

  Character or `NULL`. Filter to a specific database: `"HS"` (human),
  `"MM"` (mouse), or `NULL` (default — show both).

- verbose:

  Logical. Print the formatted table. Default `TRUE`.

## Value

Invisibly returns the full collections `data.frame` from
[`msigdbr::msigdbr_collections()`](https://igordot.github.io/msigdbr/reference/msigdbr_collections.html),
with an added `database` column (`"HS"` or `"MM"`).

## Details

Human collections use plain codes (`H`, `C2`, `C5`, …). Mouse
collections (`db_species = "MM"`) use `M`-prefixed codes (`MH`, `M2`,
`M5`, …).
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)
adjusts these automatically when `db_species = "MM"` is set, so you can
still pass `collection = c("H", "C2")` and the correct codes will be
used.
