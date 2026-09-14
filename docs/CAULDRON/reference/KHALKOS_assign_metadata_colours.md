# Assign a stable colour to every level of qualifying colData columns

CAULDRON's plotting functions otherwise assign discrete colours
positionally, fresh on every call (see
[`KHALKOS_default_palette`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_default_palette.md)
and how e.g. `ASPIS_plot_umap` uses it) – so the same category (e.g.
cluster `"4"`, genotype `"KO"`) can end up a different colour in two
different plots if the set/order of levels present differs between
calls. This builds one colour map per qualifying `colData(sce)` column
instead, meant to be computed once and threaded through to every
plotting call in a report (e.g. `ASPIS_plot_umap(palette = ...)`,
`KERAUNOS_propeller_proportions(cluster_colors = ..., group_colors = ...)`)
so the same category renders identically everywhere.

## Usage

``` r
KHALKOS_assign_metadata_colours(
  sce,
  palette = NULL,
  max_levels = 30,
  overrides = list(),
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- palette:

  Character vector of hex colours to draw from, recycled if a column has
  more levels than colours supplied. Default `NULL` uses
  [`KHALKOS_default_palette`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_default_palette.md).

- max_levels:

  Integer. Only `colData` columns with at most this many unique values
  are assigned colours. Default `30`.

- overrides:

  Named list, one entry per `colData` column to override (need not cover
  every column). Each entry is itself a named character vector (category
  label -\> hex colour); only the categories named are overwritten – the
  rest of that column keeps its palette-assigned colour. Default: empty
  list (no overrides).

- verbose:

  Logical. Print which columns were included/skipped. Default `TRUE`.

## Value

Named list, one entry per qualifying `colData` column, each a named
character vector (category label -\> hex colour). A column absent from
the result (excluded by type or `max_levels`) simply isn't a name in
this list – indexing it (e.g. `colours$gene_count`) returns `NULL`,
which every consumer here already treats as "use the default palette
instead".

## Details

Only character/factor columns with at most `max_levels` unique values
qualify – this excludes continuous metadata (e.g. a UMI count column) by
type, and high-cardinality identifier columns (e.g. cell barcodes) by
level count, without needing either named explicitly.
