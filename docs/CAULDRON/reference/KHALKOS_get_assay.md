# Safe assay retrieval

Returns `assay(sce, assay_name)` with a clear error message listing
available assays when the requested name is not found.

## Usage

``` r
KHALKOS_get_assay(sce, assay_name, error_if_missing = TRUE)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- assay_name:

  Character scalar. Name of the assay to retrieve.

- error_if_missing:

  Logical. If `TRUE` (default), throw an error when the assay is absent.
  If `FALSE`, return `NULL` silently.

## Value

The assay matrix, or `NULL` if missing and `error_if_missing = FALSE`.
