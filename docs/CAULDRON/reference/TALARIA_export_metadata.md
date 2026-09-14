# Export cell metadata to CSV

Writes `colData(sce)` to a CSV file, with cell barcodes as the first
column.

## Usage

``` r
TALARIA_export_metadata(sce, output_path, verbose = TRUE)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- output_path:

  Character. Full path for the output CSV file.

- verbose:

  Logical. Print confirmation. Default `TRUE`.

## Value

`output_path`, invisibly.
