# Select highly variable genes

Models the mean-variance relationship across genes using
[`scran::modelGeneVar`](https://rdrr.io/pkg/scran/man/modelGeneVar.html)
and selects the top `n_hvgs` genes by biological variance component.
Results are stored in two places so downstream functions can find them
easily:

## Usage

``` r
PYRI_select_hvg(
  sce,
  n_hvgs = 2000L,
  block_col = NULL,
  assay_name = "logcounts",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"logcounts"` assay (i.e.
  [`PYRI_normalize`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_normalize.md)
  must have been run first).

- n_hvgs:

  Integer. Number of top HVGs to select. Default `2000`.

- block_col:

  Character. Column in `colData` to use as a blocking factor (typically
  sample or animal ID). `NULL` (default) fits one trend across all
  cells.

- assay_name:

  Character. Assay to model variance from. Default `"logcounts"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

The input SCE with `rowData(sce)$is_hvg` and `metadata(sce)$hvg`
populated.

## Details

- `rowData(sce)$is_hvg` — logical vector flagging selected genes

- `metadata(sce)$hvg` — character vector of selected gene names

For multi-sample datasets it is strongly recommended to set `block_col`
to a `colData` column identifying sample of origin (e.g. `"animal_id"`).
This models the mean-variance trend per sample and combines the results,
preventing any single sample with an unusual depth from dominating the
selection.
