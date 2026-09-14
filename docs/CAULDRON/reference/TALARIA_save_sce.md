# Save a SingleCellExperiment object to disk

Serialises the SCE to an RDS file. Use
[`TALARIA_load_sce`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_sce.md)
to reload.

## Usage

``` r
TALARIA_save_sce(sce, output_path, verbose = TRUE)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- output_path:

  Character. Full path for the `.rds` file.

- verbose:

  Logical. Print confirmation. Default `TRUE`.

## Value

`output_path`, invisibly.
