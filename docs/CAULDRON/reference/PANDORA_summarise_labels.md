# Summarise cell type label composition

Returns overall label frequencies and optionally a breakdown by cluster
and/or sample. Prints a formatted composition table to the console.

## Usage

``` r
PANDORA_summarise_labels(
  sce,
  label_col,
  cluster_col = NULL,
  sample_col = NULL,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- label_col:

  Character. `colData` column with cell type labels (output of any
  PANDORA annotation function or
  [`PANDORA_assign_labels`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_assign_labels.md)).

- cluster_col:

  Character or `NULL`. If provided, also tabulates label composition
  within each cluster. Default `NULL`.

- sample_col:

  Character or `NULL`. If provided, also tabulates label composition
  within each sample. Default `NULL`.

- verbose:

  Logical. Print the composition table. Default `TRUE`.

## Value

Invisibly returns a named list:

- `$overall` — data.frame of overall label frequencies.

- `$by_cluster` — data.frame or `NULL`.

- `$by_sample` — data.frame or `NULL`.
