# Detect doublets with scDblFinder

Identifies likely doublets (droplets containing two cells) using
[`scDblFinder`](https://plger.github.io/scDblFinder/reference/scDblFinder.html).
Doublet scores and classifications are appended to `colData` as
`scDblFinder.score` and `scDblFinder.class`.

## Usage

``` r
AEGIS_detect_doublets(
  sce,
  sample_col = NULL,
  clusters_col = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  seed = 42L,
  force_bpcells = FALSE,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- sample_col:

  Character. Column in `colData` identifying sample of origin (e.g.
  `"mouse_id"`). Doublet detection is then run per sample. `NULL`
  (default) treats the entire dataset as one sample.

- clusters_col:

  Character. Column in `colData` with cluster labels. Providing
  pre-computed clusters improves sensitivity. `NULL` (default) lets
  scDblFinder cluster internally.

- BPPARAM:

  A `BiocParallelParam` object controlling parallelisation. Default
  `SerialParam()` (single-threaded).

- seed:

  Integer. Random seed for reproducibility. Default `42L`.

- force_bpcells:

  Logical. When the counts assay is a BPCells-backed `IterableMatrix`,
  `scDblFinder` will internally coerce the full matrix into RAM, which
  may exhaust memory on large datasets. By default (`FALSE`), doublet
  detection is skipped and a message is printed. Set `TRUE` to run
  regardless — only do this if you are confident the dataset fits in
  available RAM.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

The input `SingleCellExperiment` with doublet annotations added to
`colData`, or the unchanged SCE if the step was skipped due to a BPCells
backend.

## Details

When the dataset contains multiple samples, pass the `sample_col`
argument so that doublet simulation is performed per sample, which
substantially improves accuracy.
