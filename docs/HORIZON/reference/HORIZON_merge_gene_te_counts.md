# Merge gene-level and TE-level count matrices

Combines a canonical gene-level count matrix (e.g. from
[`HORIZON_aggregate_counts`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_counts.md))
with a TE-subfamily count matrix (from
[`HORIZON_aggregate_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_tecount.md))
into a single feature x sample matrix suitable for
`ARTEMIS_normalize_counts()`. There is no prior HORIZON/GAIA convention
for combining two distinct feature sets (existing aggregators only ever
combine samples for a single feature set), so this fills that gap.

## Usage

``` r
HORIZON_merge_gene_te_counts(gene_counts, te_counts, verbose = TRUE)
```

## Arguments

- gene_counts:

  Matrix. Gene x sample count matrix, with gene IDs as rownames and
  sample IDs as colnames.

- te_counts:

  Matrix. TE subfamily x sample count matrix, with composite TE IDs as
  rownames and sample IDs as colnames.

- verbose:

  Logical. Print merge diagnostics (dimensions, any dropped samples).
  Default TRUE.

## Value

The combined matrix (genes + TEs) x sample, restricted to samples
present in both inputs.
