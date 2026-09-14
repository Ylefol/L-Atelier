# Principal component analysis

Runs PCA on the normalised `logcounts` assay, restricted to HVGs
selected by
[`PYRI_select_hvg`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_select_hvg.md)
when available. Results are stored in `reducedDims(sce)[["PCA"]]`.

## Usage

``` r
TALOS_run_pca(
  sce,
  n_pcs = 50L,
  use_hvg = TRUE,
  assay_name = "logcounts",
  scale = FALSE,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"logcounts"` assay.

- n_pcs:

  Integer. Number of principal components to compute. Default `50`.

- use_hvg:

  Logical. If `TRUE` (default), restricts to genes in
  `metadata(sce)$hvg`. Falls back to all genes with a warning if HVGs
  are not found.

- assay_name:

  Character. Assay to use. Default `"logcounts"`.

- scale:

  Logical. Scale genes to unit variance before PCA. Default `FALSE`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `reducedDims(sce)[["PCA"]]` populated.

## Details

PCA is computed via
[`irlba::irlba`](https://rdrr.io/pkg/irlba/man/irlba.html) directly for
all matrix backends (in-memory `dgCMatrix` and BPCells `IterableMatrix`
alike). This ensures identical numerical behaviour regardless of the
backend chosen at load time, making results comparable between the two
modes.
[`BiocGenerics::t()`](https://rdrr.io/pkg/BiocGenerics/man/t.html) is
used for the transpose so that S4 dispatch works correctly for all
matrix classes.
