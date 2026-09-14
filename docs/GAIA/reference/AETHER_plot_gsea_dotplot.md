# NES Dotplot for GSEA Results

Creates a dotplot showing the top enriched and depleted gene sets from a
GSEA result
([`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md)),
ranked by normalized enrichment score (NES).

## Usage

``` r
AETHER_plot_gsea_dotplot(
  gsea_result,
  top_n = 20,
  title = "GSEA",
  font_size = 8,
  max_label_length = 55
)
```

## Arguments

- gsea_result:

  An `apollo_gsea` object from
  [`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md),
  or a data.frame with columns `pathway`, `NES`, `padj`, `size`
  (typically `$results` or `$significant`).

- top_n:

  Integer. Total number of gene sets to show, split between top enriched
  (NES \> 0) and top depleted (NES \< 0) by padj. Default: 20.

- title:

  Character. Plot title. Default: "GSEA".

- font_size:

  Numeric. Base font size for pathway labels. Default: 8.

- max_label_length:

  Integer. Pathway names longer than this are word-wrapped onto multiple
  lines (never truncated/"..."-ed). Default: 55.

## Value

A ggplot object, or `NULL` (with a warning) if there are no results to
plot.

## Details

Dot size = number of genes in the leading edge (`size` column); dot
color = -log10(padj). A dashed vertical line marks NES = 0.

Long pathway names are wrapped (via
[`strwrap()`](https://rdrr.io/r/base/strwrap.html)), not truncated –
MSigDB-style names (e.g. `"GOBP_CHROMOSOME_SEGREGATION"`) have no spaces
to wrap on, so underscores are treated as break points the same way
spaces would be, with the underscore at the chosen break consumed by the
line break (matching how a normal word-wrap drops the space it breaks
on).

## Examples

``` r
if (FALSE) { # \dontrun{
gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))
p <- AETHER_plot_gsea_dotplot(gsea_result, top_n = 30)
} # }
```
