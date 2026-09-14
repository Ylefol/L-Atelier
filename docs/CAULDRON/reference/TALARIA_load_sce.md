# Load a SingleCellExperiment object from disk

Reads an RDS file previously saved with
[`TALARIA_save_sce`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_save_sce.md)
and validates that it contains a `SingleCellExperiment`.

## Usage

``` r
TALARIA_load_sce(path, verbose = TRUE)
```

## Arguments

- path:

  Character. Path to the `.rds` file.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `SingleCellExperiment` object.
