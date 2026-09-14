# MA Plot for Differential Expression/Accessibility Results

Creates a four-category MA plot (log2 fold change vs log2 mean
expression) with embedded counts in the legend and optional gene
labeling.

## Usage

``` r
AETHER_plot_ma(
  dea_result,
  label_col = NULL,
  genes_of_interest = NULL,
  label_top_n = 0,
  filter_choice = "padj",
  l2fc_thresh = 1,
  p_thresh = 0.05,
  title = "MA Plot",
  colors = NULL,
  point_size = 0.8,
  point_alpha = 0.7
)
```

## Arguments

- dea_result:

  A data.frame of DEA results, or a list with a `$results` element (as
  returned by `ARTEMIS_perform_dea()`). Required columns: `baseMean`,
  `log2FoldChange`, and the column named by `filter_choice`.

- label_col:

  Character. Column to use for point labels. If NULL, auto-selects
  `gene_name` if present, otherwise `feature_id`.

- genes_of_interest:

  Character vector. Values in `label_col` to label. Default: NULL.

- label_top_n:

  Integer. Label the top N features by `filter_choice`. Default: 0
  (disabled).

- filter_choice:

  Column for significance filtering: `"padj"` or `"pvalue"`. Default:
  `"padj"`.

- l2fc_thresh:

  Numeric. log2 fold change threshold. Default: 1.

- p_thresh:

  Numeric. Significance threshold. Default: 0.05.

- title:

  Character. Plot title. Default: "MA Plot".

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

X-axis is `log2(baseMean + 1)`. Y-axis is `log2FoldChange`. A solid red
line marks y = 0. Dashed black lines mark `±l2fc_thresh`.
Non-significant points are drawn first (background) with significant
points on top.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_ma(dea_result, title = "KO vs WT")

p <- AETHER_plot_ma(dea_result, label_top_n = 10,
                    genes_of_interest = c("SMUG1", "OGG1"))
} # }
```
