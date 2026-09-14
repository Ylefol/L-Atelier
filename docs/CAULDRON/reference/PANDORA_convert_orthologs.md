# Convert gene symbols to orthologs

Converts a character vector of gene symbols from one species to another
using orthogene. By default the returned vector is the same length as
the input, with `NA` for genes that could not be mapped. Set
`drop_na = TRUE` to silently remove unmapped entries (used internally by
[`PANDORA_annotate_sctype`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_sctype.md)).

## Usage

``` r
PANDORA_convert_orthologs(
  genes,
  from = "human",
  to,
  method = "homologene",
  drop_na = FALSE,
  verbose = TRUE
)
```

## Arguments

- genes:

  Character vector of gene symbols to convert.

- from:

  Character. Source species. Default `"human"`.

- to:

  Character. Target species (e.g. `"mouse"`, `"Mus musculus"`,
  `"zebrafish"`). Any species string accepted by
  [`orthogene::convert_orthologs()`](https://rdrr.io/pkg/orthogene/man/convert_orthologs.html)
  is valid.

- method:

  Character. Ortholog database passed to
  [`orthogene::convert_orthologs()`](https://rdrr.io/pkg/orthogene/man/convert_orthologs.html).
  Default `"homologene"`. Use `"gprofiler"` for broader species
  coverage.

- drop_na:

  Logical. If `FALSE` (default), unmapped genes are returned as `NA`
  (output length equals input length). If `TRUE`, unmapped genes are
  dropped from the result.

- verbose:

  Logical. Print mapping statistics. Default `TRUE`.

## Value

A named character vector. Names are the input gene symbols; values are
the corresponding ortholog symbols (or `NA` when unmapped and
`drop_na = FALSE`).
