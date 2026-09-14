# Export a count matrix to disk

Writes the specified assay to disk in CSV or 10x MEX format.

## Usage

``` r
TALARIA_export_matrix(
  sce,
  output_dir,
  format = c("csv", "mex"),
  assay_name = "counts",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- output_dir:

  Character. Directory to write output into. Created if it does not
  exist.

- format:

  Character. `"csv"` (default) or `"mex"`.

- assay_name:

  Character. Assay to export. Default `"counts"`.

- verbose:

  Logical. Print confirmation. Default `TRUE`.

## Value

Output path(s), invisibly.

## Details

**Note:** CSV export calls
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) on the assay, which
loads the full matrix into RAM. For large datasets prefer `"mex"`
(requires DropletUtils).
