# Build a ranked gene vector for GSEA from a differential expression result

Extracts a named numeric vector (gene -\> ranking statistic) from an
`artemis_limma` or `artemis_ts_de` result, suitable for
[`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md).

## Usage

``` r
APOLLO_rank_from_de(
  de_result,
  rank_by = c("t", "log2FoldChange", "signed_neg_log10p"),
  gene_col = "feature_id",
  comparison = NULL,
  verbose = TRUE
)
```

## Arguments

- de_result:

  An `artemis_limma` object (from
  [`ARTEMIS_limma_de()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_limma_de.md)),
  an `artemis_ts_de` object (from
  [`ARTEMIS_timeseries_conditional()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_conditional.md)/`_temporal()`
  or their `_limma` counterparts), or a plain data.frame with the
  relevant columns.

- rank_by:

  Character. Ranking statistic:

  `"t"`

  :   The test statistic - limma's `t` column or DESeq2's `stat` column
      (whichever is present is used automatically). Accounts for both
      effect size and precision; the recommended default.

  `"log2FoldChange"`

  :   Fold change only - ignores variance.

  `"signed_neg_log10p"`

  :   `sign(log2FoldChange) * -log10(pvalue)`. Useful when no test
      statistic column is available.

  Default: `"t"`.

- gene_col:

  Character. Column in the results data.frame to use as gene names.
  Default: `"feature_id"`. Use `"original_id"` for Olink data when the
  gene sets are expected to match the matrix's original IDs (e.g.
  OlinkIDs) rather than the substituted `feature_id_col` symbols.

- comparison:

  Character, integer, or `NULL`. For `artemis_ts_de` objects (which
  contain multiple comparisons): the comparison name or index to
  extract. `NULL` uses the first comparison and prints a note. Ignored
  for `artemis_limma` and data.frame inputs.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Named numeric vector, sorted descending, NA values removed, duplicate
gene names removed (first occurrence kept).

## Examples

``` r
if (FALSE) { # \dontrun{
de <- ARTEMIS_limma_de(ol$wide, ol$sample_meta, group_col = "Group",
                        reference = "Control", experiment = "Sepsis")
ranked <- APOLLO_rank_from_de(de, rank_by = "t")

# From a timeseries result, specific comparison
ranked_tp2 <- APOLLO_rank_from_de(temp_de, rank_by = "t",
                                   comparison = "2_vs_1")
} # }
```
