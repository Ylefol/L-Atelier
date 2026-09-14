# Convert a SingleCellExperiment to a Seurat object

Wraps
[`Seurat::as.Seurat()`](https://satijalab.github.io/seurat-object/reference/as.Seurat.html)
to convert a `SingleCellExperiment` to a `Seurat` object, transferring
counts, normalised data, dimensionality reductions, and cell metadata.

## Usage

``` r
TALARIA_to_seurat(sce, counts = "counts", data = NULL)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- counts:

  Character. Assay to use as the raw counts slot in Seurat. Default
  `"counts"`.

- data:

  Character. Assay to use as the normalised data slot. Default
  `"logcounts"` if present, otherwise `NULL`.

## Value

A `Seurat` object.
