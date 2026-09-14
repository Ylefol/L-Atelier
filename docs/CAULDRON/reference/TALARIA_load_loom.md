# Load a Loom file into a SingleCellExperiment

Reads a `.loom` file using LoomExperiment and returns a
`SingleCellExperiment`. Loom is a legacy format; prefer `.h5ad` for new
projects.

## Usage

``` r
TALARIA_load_loom(path, verbose = TRUE)
```

## Arguments

- path:

  Character. Path to the `.loom` file.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `SingleCellExperiment`.
