# List available annotation references

Prints a formatted summary of supported `celldex` reference shorthands
for
[`PANDORA_annotate_singler`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_singler.md)
and available tissue types from the bundled ScTypeDB for
[`PANDORA_annotate_sctype`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_sctype.md).
Use this to quickly confirm whether either approach covers your dataset
before committing to an annotation strategy.

## Usage

``` r
PANDORA_list_references(source = c("all", "singler", "sctype"), db = NULL)
```

## Arguments

- source:

  Character. Which sources to list: `"all"` (default), `"singler"`, or
  `"sctype"`.

- db:

  `NULL` (use bundled `ScTypeDB_full.xlsx`), a file path to a custom
  Excel or CSV, or a `data.frame` already loaded in R. Only relevant
  when `source` includes `"sctype"`.

## Value

Invisibly returns a named list with `$singler` (reference info
`data.frame`) and `$sctype` (character vector of tissue types).
