# Print method for ARTEMIS differential analysis results

Reprints the comparison and significance summary (and top 5 hits) from
an `artemis_dea` object without rerunning DESeq2 – the same summary
[`ARTEMIS_differential_counts`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
prints live when `verbose = TRUE`, but callable on the already-computed
result.

## Usage

``` r
# S3 method for class 'artemis_dea'
print(x, ...)
```

## Arguments

- x:

  An `artemis_dea` object from
  [`ARTEMIS_differential_counts`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md).

- ...:

  Additional arguments (ignored)
