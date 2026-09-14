# Convert a Seurat object to a SingleCellExperiment

Wraps
[`Seurat::as.SingleCellExperiment()`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html)
to convert a `Seurat` object to a `SingleCellExperiment`, preserving
assays, embeddings, and cell metadata.

## Usage

``` r
TALARIA_from_seurat(sobj, assay = "RNA", verbose = TRUE)
```

## Arguments

- sobj:

  A `Seurat` object.

- assay:

  Character. Which Seurat assay to convert. Default `"RNA"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `SingleCellExperiment`.
