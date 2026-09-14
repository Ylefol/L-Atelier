# Normalise a SingleCellExperiment

Computes normalised log-counts and stores them as the `"logcounts"`
assay. Two methods are available:

## Usage

``` r
PYRI_normalize(
  sce,
  method = c("lognorm", "scran"),
  assay_name = "counts",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"counts"` assay.

- method:

  Character. `"lognorm"` (default) or `"scran"`.

- assay_name:

  Character. Input count assay. Default `"counts"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

The input SCE with a `"logcounts"` assay added and, for
`method = "scran"`, size factors stored in `sizeFactors(sce)`.

## Details

- `"lognorm"`:

  Library-size normalisation followed by a log1p transform
  ([`scater::logNormCounts`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html)).
  Fast and appropriate for most datasets. Each cell is scaled to a
  common size factor derived from its total count.

- `"scran"`:

  Pooling-based size factor estimation
  ([`scran::quickCluster`](https://rdrr.io/pkg/scran/man/quickCluster.html) +
  [`scran::computeSumFactors`](https://rdrr.io/pkg/scran/man/computeSumFactors.html))
  followed by
  [`scater::logNormCounts`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html).
  More robust for datasets with large differences in cell type
  composition or sequencing depth, at the cost of additional compute
  time.
