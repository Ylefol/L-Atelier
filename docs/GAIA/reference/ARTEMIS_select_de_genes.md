# Select DE Genes Across Differential Expression Results

Collects the union (or intersection) of significant genes across one or
more differential expression comparisons. Works with time series DEA
results, standard DEA results from
[`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md),
or raw data.frames containing DE statistics.

## Usage

``` r
ARTEMIS_select_de_genes(
  de_results,
  l2fc_thresh = 1,
  p_col = "padj",
  p_thresh = 0.05,
  gene_col = NULL,
  union = TRUE,
  verbose = TRUE
)
```

## Arguments

- de_results:

  DE results in any of the following formats:

  - `artemis_ts_de` object (from time series analysis)

  - DEA result from
    [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)

  - A data.frame with gene IDs, log2FoldChange, and p-values

  - A list of any combination of the above

- l2fc_thresh:

  Numeric. Absolute log2 fold-change threshold. Only genes with
  \|log2FC\| \>= this value are selected. Default: 1.0.

- p_col:

  Character. P-value column to use: "padj" or "pvalue". Default: "padj".

- p_thresh:

  Numeric. P-value threshold. Default: 0.05.

- gene_col:

  Character or NULL. Column name containing gene IDs. If NULL,
  auto-detects from common names (gene_id, feature_id, gene, etc.).
  Default: NULL.

- union:

  Logical. If TRUE, take union of genes across comparisons. If FALSE,
  take intersection. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

A character vector of selected gene IDs with attributes:

- gene_summary:

  Data.frame showing which comparisons each gene was significant in

- per_comparison:

  Named list of gene sets per comparison

- parameters:

  Selection parameters used

## Examples

``` r
if (FALSE) { # \dontrun{
# From time series DEA
genes <- ARTEMIS_select_de_genes(list(cond_de, temp_de), l2fc_thresh = 2)

# From standard DEA (ARTEMIS_differential_counts)
genes <- ARTEMIS_select_de_genes(dea_result, l2fc_thresh = 1)

# From multiple DEA results
genes <- ARTEMIS_select_de_genes(
  list(dea_result1, dea_result2, dea_result3),
  l2fc_thresh = 1,
  union = TRUE
)

} # }
```
