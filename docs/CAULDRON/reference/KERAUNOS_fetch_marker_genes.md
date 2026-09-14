# Fetch known cell type marker genes from a curated database

Retrieves the union of all gene symbols annotated as cell type markers
across a curated database. The resulting character vector is intended to
be passed to the `restrict_to` parameter of
[`KERAUNOS_score_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_score_markers.md),
limiting reported markers to biologically validated genes and filtering
out uninformative features such as long non-coding RNAs, ribosomal
genes, and ubiquitously expressed housekeeping genes.

## Usage

``` r
KERAUNOS_fetch_marker_genes(
  source = c("C8", "PanglaoDB"),
  species = "Homo sapiens",
  tissue = NULL,
  path = NULL,
  verbose = TRUE
)
```

## Arguments

- source:

  Character. Database to query: `"C8"` (default) or `"PanglaoDB"`.

- species:

  Character. Species for gene symbol resolution. For `source = "C8"`:
  full species name passed to `msigdbr` (e.g. `"Homo sapiens"`,
  `"Mus musculus"`). For `source = "PanglaoDB"`: `"Hs"` (human,
  default), `"Mm"` (mouse), or `"both"` (no species filter). Default
  `"Homo sapiens"`.

- tissue:

  Character or `NULL`. For `source = "PanglaoDB"` only: restrict to
  markers annotated for a specific tissue/organ (e.g. `"Brain"`).
  Case-insensitive partial match against the `organ` column. `NULL`
  (default) returns all tissues.

- path:

  Character or `NULL`. For `source = "PanglaoDB"` only: path to the
  local PanglaoDB markers TSV file. When `NULL` (default), looks for
  `data/PanglaoDB_markers.tsv` relative to the current working
  directory.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A sorted character vector of unique gene symbols known to be cell type
markers in the chosen database.

## Sources

- `"C8"` (default):

  MSigDB collection C8 — cell type signature gene sets aggregated from
  PanglaoDB, CellMarker, DICE, Blueprint, Monaco, HPCA, and others
  (~15,000 unique genes in human). Requires the `msigdbr` package, which
  is already used by
  [`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md).

- `"PanglaoDB"`:

  PanglaoDB marker database (~8,000 marker associations across \>1,100
  cell types in human and mouse). Reads from a local TSV file (download
  from <https://panglaodb.se/markers.html>). Default path:
  `data/PanglaoDB_markers.tsv` relative to the current working
  directory. Override with the `path` argument. Supports tissue
  filtering via the `tissue` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
known_markers <- KERAUNOS_fetch_marker_genes(source = "C8",
                                              species = "Homo sapiens")
markers <- KERAUNOS_score_markers(sce, restrict_to = known_markers)

# PanglaoDB — brain markers only
brain_markers <- KERAUNOS_fetch_marker_genes(source = "PanglaoDB",
                                              tissue = "Brain")
} # }
```
