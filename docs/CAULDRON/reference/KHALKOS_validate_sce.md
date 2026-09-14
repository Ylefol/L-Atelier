# Validate a SingleCellExperiment object

Checks that an object is a `SingleCellExperiment` and optionally that
required assays, reduced dimensions, `colData` columns, and `metadata`
keys are present. All failures are collected and reported together in a
single informative error.

## Usage

``` r
KHALKOS_validate_sce(
  sce,
  require_assays = NULL,
  require_reduceddims = NULL,
  require_coldata = NULL,
  require_metadata = NULL,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` to validate.

- require_assays:

  Character vector of assay names that must be present.

- require_reduceddims:

  Character vector of `reducedDim` names that must be present.

- require_coldata:

  Character vector of `colData` column names that must be present.

- require_metadata:

  Character vector of `metadata` keys that must be present.

- verbose:

  Logical. Print a short confirmation on success. Default `TRUE`.

## Value

Invisibly returns `TRUE` if all checks pass; throws an informative error
otherwise.
