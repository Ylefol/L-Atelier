# Export an interactive gene-expression-by-group widget for R Markdown

Builds a self-contained, dependency-free HTML/JS widget that lets a
report viewer search any gene and see its mean expression across levels
of one or more grouping variables (e.g. cluster, genotype, cell label),
rendered as an SVG bar chart. Meant to be called directly inside an R
Markdown chunk: the return value knits as raw HTML automatically (see
[`asis_output`](https://rdrr.io/pkg/knitr/man/asis_output.html)) — no
`results='asis'` chunk option or manual
[`cat()`](https://rdrr.io/r/base/cat.html) needed.

## Usage

``` r
TALARIA_export_gene_expr_widget(
  sce,
  groupings = NULL,
  assay_name = "logcounts",
  widget_id = "gene-expr",
  digits = 4
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- groupings:

  Character vector of `colData(sce)` column names to expose as grouping
  options (e.g. `c("cluster", "genotype", "cell_label")`).

- assay_name:

  Character. Assay averaged per group. Default `"logcounts"`.

- widget_id:

  Character. DOM id prefix — must be unique if more than one widget is
  embedded in the same document. Default `"gene-expr"`.

- digits:

  Integer. Rounding applied to exported mean values. Default `4`.

## Value

An object of class `knit_asis` (see
[`asis_output`](https://rdrr.io/pkg/knitr/man/asis_output.html)).
Printing or auto-printing it inside an R Markdown chunk renders the
widget; outside of knitr it just prints as plain HTML text.

## Details

Expression is precomputed as a gene x group-level **mean** matrix, one
per grouping variable, not exported per cell. For a full transcriptome
this keeps the embedded payload on the order of a few MB rather than the
multi-GB a per-cell export would require, at the cost of showing each
group's mean rather than per-cell resolution — plotting a per-cell UMAP
colored by a per-cluster average would imply resolution the data doesn't
have, so a bar chart is used instead.
