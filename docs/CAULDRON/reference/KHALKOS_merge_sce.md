# Merge multiple SingleCellExperiment objects

Column-binds a list of `SingleCellExperiment` objects into a single
combined SCE. Optionally adds a `colData` column identifying which input
object each cell came from.

## Usage

``` r
KHALKOS_merge_sce(
  sce_list,
  sample_col = "sample_id",
  sample_names = NULL,
  verbose = TRUE
)
```

## Arguments

- sce_list:

  A named or unnamed list of `SingleCellExperiment` objects. Must
  contain at least 2 elements.

- sample_col:

  Character scalar or `NULL`. If non-`NULL`, a new `colData` column with
  this name is added to each SCE before merging, recording which sample
  it came from. Default `"sample_id"`.

- sample_names:

  Character vector of length `length(sce_list)`. Labels for each SCE.
  Defaults to `names(sce_list)` or `"sample_1"`, `"sample_2"`, … if
  unnamed.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A merged `SingleCellExperiment`.

## Details

All SCEs must share the same genes (rows) in the same order. If assay
names differ across objects, `cbind` will fill missing assays with `NA`
— ensure all objects have been processed to the same point before
merging.
