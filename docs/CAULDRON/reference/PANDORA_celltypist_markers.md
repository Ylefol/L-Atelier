# CellTypist's own classifier-driving genes for a cell type

CellTypist classifies cells with a linear model trained on a reference
atlas. This extracts the genes with the largest coefficients for one or
more of that model's cell types — i.e. the genes the model itself relied
on to define the population in its training atlas.

This answers a different question from
[`KERAUNOS_find_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_find_markers.md):
that function tests for genes differentially expressed between clusters
in the CURRENT dataset (a statistical marker test), whereas this reports
what CellTypist's classifier was actually trained on (a property of the
reference atlas). The two are complementary and worth cross-checking
against each other, but neither substitutes for the other — a gene the
classifier relies on heavily need not be significantly differential in
every dataset it's applied to, and vice versa.

## Usage

``` r
PANDORA_celltypist_markers(
  cell_types,
  model,
  top_n = 10L,
  only_positive = TRUE,
  force_update = FALSE,
  verbose = TRUE
)
```

## Arguments

- cell_types:

  Character vector. Cell type label(s), must match an entry in the
  model's own `cell_types` exactly — an error lists the offending
  label(s) if not.

- model:

  Character. Model filename, e.g. `"Mouse_Whole_Brain.pkl"` — same
  models listed by
  [`PANDORA_list_celltypist_models`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_celltypist_models.md).

- top_n:

  Integer. Number of driving genes per cell type. Default `10`.

- only_positive:

  Logical. Only genes with a positive coefficient (higher expression
  pushes toward this class). Set `FALSE` to also include negative
  markers (genes whose LOW expression is diagnostic). Default `TRUE`.

- force_update:

  Logical. Re-download the model even if already cached. Default
  `FALSE`.

- verbose:

  Logical. Print the result. Default `TRUE`.

## Value

Named list, cell type -\> character vector of driving genes, ordered by
decreasing coefficient magnitude.

## Examples

``` r
if (FALSE) { # \dontrun{
PANDORA_celltypist_markers(
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  model = "Mouse_Whole_Brain.pkl", top_n = 20
)
} # }
```
