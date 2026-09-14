# Print method for artemis_ts_de

Mirrors the per-comparison format of
[`print.artemis_dea`](https://ylefol.github.io/L-Atelier/GAIA/reference/print.artemis_dea.md)
(features tested, significance breakdown, top hits), looped once per
comparison bundled in this object – since `ARTEMIS_differential_counts`
and the time series DEA functions in this file (DESeq2- and limma-backed
alike) all answer the same question (differential expression between two
conditions), their console summaries read the same way, one comparison
at a time, regardless of which engine produced them.

## Usage

``` r
# S3 method for class 'artemis_ts_de'
print(x, ...)
```

## Arguments

- x:

  An artemis_ts_de object

- ...:

  Additional arguments (ignored)
