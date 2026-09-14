# Dotplot for Enrichment Results

Creates a dotplot visualization for GO/KEGG enrichment results. Supports
coloring by p-value, adjusted p-value, or odds ratio.

## Usage

``` r
AETHER_plot_enrichment_dotplot(
  enrich_result,
  show_category = 20,
  color_by = "p.adjust",
  size_by = "Count",
  order_by = "p.adjust",
  title = "Enrichment Analysis",
  x_label = NULL,
  font_size = 10,
  color_low = "#2166ac",
  color_high = "#b2182b",
  color_midpoint = NULL
)
```

## Arguments

- enrich_result:

  An enrichResult object from clusterProfiler (output of
  APOLLO_enrich_go or APOLLO_enrich_kegg), or a data.frame with the
  required columns.

- show_category:

  Integer. Number of top categories to display (default = 20).

- color_by:

  Character. Variable for color scale: "p.adjust" (default), "pvalue",
  or "odds_ratio".

- size_by:

  Character. Variable for dot size: "Count" (default) or "GeneRatio".

- order_by:

  Character. Variable to order terms: "p.adjust" (default), "pvalue",
  "odds_ratio", or "Count".

- title:

  Character. Plot title (default = "Enrichment Analysis").

- x_label:

  Character. X-axis label. If NULL, auto-generated based on size_by
  parameter.

- font_size:

  Numeric. Base font size for term labels (default = 10).

- color_low:

  Character. Color for low values (default = "#2166ac" blue).

- color_high:

  Character. Color for high values (default = "#b2182b" red).

- color_midpoint:

  Numeric or NULL. Midpoint for diverging color scale. If NULL, uses
  sequential scale.

## Value

A ggplot object.

## Details

The dotplot shows:

- Y-axis: enriched terms (GO terms or KEGG pathways)

- X-axis: gene ratio or count

- Dot size: number of genes (or gene ratio)

- Dot color: significance (p-value/q-value) or effect size (odds ratio)

**Why use Odds Ratio for color?** P-values conflate effect size with
sample size. A pathway with OR=10 and p=0.01 may be more biologically
relevant than OR=1.5 with p=1e-10. Odds ratio directly measures the
strength of association:

- OR = 1: no association

- OR \> 1: enrichment (higher = stronger)

- OR \< 1: depletion

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard dotplot (color by adjusted p-value)
p <- AETHER_plot_enrichment_dotplot(go_results)

# Color by odds ratio for effect size
p <- AETHER_plot_enrichment_dotplot(go_results, color_by = "odds_ratio")

# Show more categories, order by odds ratio
p <- AETHER_plot_enrichment_dotplot(
  go_results,
  show_category = 30,
  color_by = "odds_ratio",
  order_by = "odds_ratio"
)

} # }
```
