# Plot Full (Unfiltered) g:GOSt Manhattan Plot

Wraps
[`gprofiler2::gostplot()`](https://rdrr.io/pkg/gprofiler2/man/gostplot.html)
to render the traditional, interactive g:GOSt Manhattan-style plot for a
single gene list/module, showing every evaluated term — not just those
passing the significance threshold. Complements the
significance-filtered results from
[`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)
so a cluster's biology can still be inspected visually when no term
clears multiple-testing correction.

## Usage

``` r
AETHER_plot_gost_full(gost_result, capped = TRUE, save_html = NULL)
```

## Arguments

- gost_result:

  A single raw gost result object — one element of
  `enrich_result$results` (e.g. `enrich_result$results[["C1"]]`) as
  returned by
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md).
  This already holds the full, unfiltered term set.

- capped:

  Logical. Cap the -log10(p-value) y-axis at 16 for readability
  (gprofiler2 default behavior). Default: TRUE.

- save_html:

  Character or NULL. Path to save as a standalone interactive HTML file.
  Default: NULL.

## Value

A plotly htmlwidget object.

## Details

Interactive only — gprofiler2's Manhattan plot is plotly-native and has
no meaningful static (ggplot2) equivalent; hovering each point reveals
the term name, source, and p-value, which is what makes the
non-significant terms interpretable.

## Examples

``` r
if (FALSE) { # \dontrun{
enrich <- APOLLO_enrich_gost(cluster_genes)
AETHER_plot_gost_full(enrich$results[["C1"]], save_html = "C1_gostplot_full.html")

} # }
```
