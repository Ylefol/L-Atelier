# Load a 10x Genomics HDF5 file into a SingleCellExperiment

Reads a `.h5` file produced by Cell Ranger using DropletUtils.

## Usage

``` r
TALARIA_load_h5(
  path,
  sample_names = NULL,
  backend = c("memory", "bpcells"),
  bpcells_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- path:

  Character. Path to the `.h5` file.

- sample_names:

  Character. Sample name assigned to the `colData` `Sample` column.
  Defaults to the file name without extension.

- backend:

  Character. `"memory"` (default) or `"bpcells"`. See
  [`TALARIA_load_h5ad`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_h5ad.md)
  for details.

- bpcells_dir:

  Character. Path for the BPCells on-disk store directory. Required when
  `backend = "bpcells"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A `SingleCellExperiment` with a `"counts"` assay.
