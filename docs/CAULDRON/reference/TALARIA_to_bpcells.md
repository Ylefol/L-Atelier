# Convert a SingleCellExperiment counts assay to an on-disk BPCells matrix

Writes the specified assay to a BPCells bitpacked directory store and
replaces the in-memory matrix with a lazy `IterableMatrix` that streams
from disk. Most downstream CAULDRON operations (QC metrics,
log-normalisation, HVG selection, PCA via irlba) are natively compatible
with BPCells and will not load the full matrix into RAM.

## Usage

``` r
TALARIA_to_bpcells(sce, bpcells_dir, assay_name = "counts", verbose = TRUE)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- bpcells_dir:

  Character. Directory path for the BPCells on-disk store. Created
  (recursively) if it does not exist. Must be a permanent path —
  temporary directories will break the object after session restart.

- assay_name:

  Character. Assay to convert. Default `"counts"`.

- verbose:

  Logical. Print a summary message. Default `TRUE`.

## Value

The input `SingleCellExperiment` with the specified assay replaced by a
BPCells `IterableMatrix`. `metadata(sce)$bpcells_dir` is set to
`bpcells_dir`.

## Details

**Note on `AEGIS_detect_doublets`:** scDblFinder may internally coerce
the BPCells matrix to a standard sparse matrix, loading the full counts
into memory for doublet detection. This is a limitation of the
scDblFinder implementation. For very large datasets consider running
doublet detection on a per-sample subset.

**Persistence:** The `bpcells_dir` path is stored in
`metadata(sce)$bpcells_dir`. When saving the SCE with
[`TALARIA_save_sce`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_save_sce.md)
and reloading with
[`TALARIA_load_sce`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_sce.md),
the lazy matrix is automatically re-connected — provided `bpcells_dir`
still exists at the same path. Use a permanent, session-independent
path.
