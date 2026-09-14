# Volcano Plot for Differential Expression/Accessibility Results

Creates a four-category volcano plot (-log10 p-value vs log2 fold
change) with embedded counts in the legend and optional gene labeling.
Categories: up-regulated, down-regulated, low-regulation (significant
but below l2FC threshold), and non-significant.

## Usage

``` r
AETHER_plot_volcano(
  dea_result,
  label_col = NULL,
  genes_of_interest = NULL,
  show_non_sig_interest = TRUE,
  label_top_n = 0,
  filter_choice = "padj",
  l2fc_thresh = 1,
  p_thresh = 0.05,
  title = "Volcano Plot",
  colors = NULL,
  point_size = 0.8,
  point_alpha = 0.7
)
```

## Arguments

- dea_result:

  A data.frame of DEA results, or a list with a `$results` element (as
  returned by `ARTEMIS_perform_dea()`). Required columns:
  `log2FoldChange`, `pvalue`, and the column named by `filter_choice`.

- label_col:

  Character. Column to use for point labels. If NULL, auto-selects
  `gene_name` if present, otherwise `feature_id`.

- genes_of_interest:

  Character vector. Values in `label_col` to always label. Default:
  NULL.

- show_non_sig_interest:

  Logical. If FALSE, genes of interest that do not meet both thresholds
  are not labeled. Default: TRUE.

- label_top_n:

  Integer. Label the top N features by `filter_choice` regardless of
  direction. Default: 0 (disabled).

- filter_choice:

  Column for significance filtering: `"padj"` or `"pvalue"`. Default:
  `"padj"`.

- l2fc_thresh:

  Numeric. log2 fold change threshold. Default: 1.

- p_thresh:

  Numeric. Significance threshold. Default: 0.05.

- title:

  Character. Plot title. Default: "Volcano Plot".

- colors:

  Named character vector with colors for `"up"`, `"down"`, `"low_reg"`,
  `"non_sig"`. Default uses red/blue/green/gray.

- point_size:

  Numeric. Point size. Default: 0.8.

- point_alpha:

  Numeric. Point transparency. Default: 0.7.

## Value

A ggplot object.

## Details

The four categories are:

- **up-reg**: `log2FoldChange > l2fc_thresh` AND
  `filter_choice < p_thresh`

- **down-reg**: `log2FoldChange < -l2fc_thresh` AND
  `filter_choice < p_thresh`

- **low-regulation**: `|log2FoldChange| <= l2fc_thresh` AND
  `filter_choice < p_thresh`

- **non-significant**: `filter_choice >= p_thresh` or NA

A horizontal dashed line marks the pvalue of the least-significant gene
that still passes the filter threshold, showing the actual significance
boundary in the plotted space. Vertical dashed lines mark
`±l2fc_thresh`.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_volcano(dea_result, title = "KO vs WT")

p <- AETHER_plot_volcano(dea_result, label_top_n = 10,
                         genes_of_interest = c("SMUG1", "OGG1"))

p <- AETHER_plot_volcano(dea_result, filter_choice = "pvalue", p_thresh = 0.01)
} # }
```
