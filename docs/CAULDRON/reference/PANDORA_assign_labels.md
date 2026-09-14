# Manually assign cell type labels to clusters

Maps cluster IDs to user-defined cell type labels using a named
character vector, storing the result as a new `colData` column. Clusters
without a matching entry in `labels` receive `NA`.

## Usage

``` r
PANDORA_assign_labels(
  sce,
  labels,
  cluster_col = "cluster",
  label_col,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- labels:

  Named character vector. Names are cluster IDs (must match levels in
  `cluster_col`); values are cell type labels. Example:
  `c("0" = "T cell", "1" = "B cell", "2" = "Myeloid")`.

- cluster_col:

  Character. `colData` column containing cluster IDs. Default
  `"cluster"`.

- label_col:

  Character. Name of the new `colData` column to create. Must be
  supplied explicitly; no default to prevent accidental overwrites.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `colData(sce)[[label_col]]` populated.
