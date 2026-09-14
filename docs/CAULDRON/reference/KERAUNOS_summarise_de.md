# Summarise differential expression results

Produces a concise table (n_tested, n_sig, n_up, n_down per
cluster/group) from any KERAUNOS DE result object.

## Usage

``` r
KERAUNOS_summarise_de(de_result, fdr_threshold = 0.05, lfc_threshold = 0)
```

## Arguments

- de_result:

  A `keraunos_markers`, `keraunos_contrast`, or `keraunos_pseudobulk`
  object.

- fdr_threshold:

  Numeric. FDR cut-off applied to the results. Default `0.05`.

- lfc_threshold:

  Numeric. Minimum absolute log2 fold-change required (only applied to
  pseudobulk results, which have an LFC column). Default `0`.

## Value

A data.frame with one row per cluster/group.
