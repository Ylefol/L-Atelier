# Bar Plot for Enrichment Results

Creates a horizontal bar plot for GO/KEGG enrichment results. Simpler
alternative to dotplot, showing either gene count or odds ratio.

## Usage

``` r
AETHER_plot_enrichment_bar(
  enrich_result,
  show_category = 15,
  x_var = "Count",
  fill_by = "p.adjust",
  order_by = NULL,
  title = "Enrichment Analysis"
)
```

## Arguments

- enrich_result:

  An enrichResult object or data.frame.

- show_category:

  Integer. Number of top categories to display (default = 15).

- x_var:

  Character. Variable for bar length: "Count", "odds_ratio", or
  "GeneRatio" (default = "Count").

- fill_by:

  Character. Variable for bar fill: "p.adjust" (default), "pvalue", or a
  fixed color string.

- order_by:

  Character. Variable to order terms (default = same as x_var).

- title:

  Character. Plot title.

## Value

A ggplot object.
