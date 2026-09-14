# Safe assay assignment

A thin wrapper around `assay(sce, assay_name) <- value` that returns the
modified SCE, making it pipe-friendly.

## Usage

``` r
KHALKOS_set_assay(sce, assay_name, value)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- assay_name:

  Character scalar. Name of the assay to set (created if absent,
  replaced if present).

- value:

  Matrix to store.

## Value

SCE with the assay updated.
