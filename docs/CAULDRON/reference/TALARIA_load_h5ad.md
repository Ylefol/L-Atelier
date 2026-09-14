# Load an AnnData H5AD file into a SingleCellExperiment

Reads a `.h5ad` file produced by Python's AnnData / Scanpy and returns a
`SingleCellExperiment` object. Conversion is handled by zellkonverter,
which maps:

## Usage

``` r
TALARIA_load_h5ad(
  path,
  use_hdf5 = FALSE,
  reader = c("R", "python"),
  backend = c("memory", "bpcells"),
  bpcells_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- path:

  Character. Path to the `.h5ad` file.

- use_hdf5:

  Logical. If `TRUE`, count matrices are kept on-disk as `HDF5Array`
  objects rather than loaded into RAM. Useful for very large datasets.
  Default `FALSE`.

- reader:

  Character. `"R"` (default, pure R via rhdf5) or `"python"` (requires a
  working Python / basilisk environment).

- backend:

  Character. `"memory"` (default) loads the counts assay into RAM as a
  standard sparse matrix. `"bpcells"` writes the counts assay to an
  on-disk BPCells store and replaces the in-memory matrix with a lazy
  `BPCells::IterableMatrix` that streams from disk. Requires BPCells.

- bpcells_dir:

  Character. Path for the BPCells on-disk store directory. Required (and
  created) when `backend = "bpcells"`; ignored otherwise. Choose a
  permanent location — the path is embedded in the returned SCE and must
  remain accessible after saving/reloading the object.

- verbose:

  Logical. Print a summary of the loaded object. Default `TRUE`.

## Value

A `SingleCellExperiment` object.

## Details

- `X` (and `raw.X` if present) → assay(s) of the SCE

- `layers` → additional named assays

- `obs` columns → `colData`

- `var` columns → `rowData`

- `obsm` entries (e.g. `X_pca`, `X_umap`) → `reducedDims`
