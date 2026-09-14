# Aether - Statistical Comparison Plots

Visualisation functions for statistical comparison results, including
group-level comparisons of distributions and proportions. Plot Peak
Annotation Distribution Comparison Between Groups

Visualises the output of
[`ARTEMIS_compare_annotation_distribution()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_compare_annotation_distribution.md).
For each genomic feature, shows per-sample proportions as dots (coloured
by group) with a horizontal mean line, and overlays significance
brackets where the group comparison passes the chosen threshold.

## Usage

``` r
AETHER_plot_annotation_comparison(
  annotation_test,
  colors = NULL,
  title = "Peak Annotation Distribution by Group",
  p_threshold = 0.05,
  show_ns = TRUE,
  nrow = 2L
)
```

## Arguments

- annotation_test:

  List returned by
  [`ARTEMIS_compare_annotation_distribution()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_compare_annotation_distribution.md).

- colors:

  Named character vector mapping group labels to colours. If `NULL`,
  uses a default blue/red palette.

- title:

  Character. Plot title.

- p_threshold:

  Numeric. p.adj threshold for significance brackets. Default: `0.05`.

- show_ns:

  Logical. If `TRUE` (default), label non-significant features with
  "ns". If `FALSE`, omit brackets for non-significant features entirely.

- nrow:

  Integer. Number of rows in the facet grid. Default: `2`.

## Value

A ggplot object.

## Details

Significance is encoded as: \*\*\* p.adj \< 0.001, \*\* \< 0.01, \* \<
0.05, ns \>= 0.05.

The y-axis is free per facet so that features with very different
overall proportions (e.g., Intergenic ~40% vs 5' UTR ~1%) remain
legible.

## See also

[`ARTEMIS_compare_annotation_distribution`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_compare_annotation_distribution.md),
[`AETHER_plot_annotation_bar`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_annotation_bar.md).

## Examples

``` r
if (FALSE) { # \dontrun{
test_result <- ARTEMIS_compare_annotation_distribution(
  annotated_list = annotated_list,
  metadata       = meta,
  group_col      = "group"
)
p <- AETHER_plot_annotation_comparison(test_result)
ggsave("annotation_comparison.pdf", p, width = 14, height = 8)
} # }
```
