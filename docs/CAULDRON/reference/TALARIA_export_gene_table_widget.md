# Export a searchable rowData(sce) table widget for R Markdown

Builds a self-contained, dependency-free HTML/JS widget that lets a
report viewer search genes by name (matched against `rownames(sce)`) and
see their associated `rowData(sce)` columns in a live-filtered table.
Mirrors
[`TALARIA_export_gene_expr_widget`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gene_expr_widget.md)'s
embedded-JSON + vanilla-JS approach so the report stays free of external
JS dependencies (a package like DT or crosstalk would each pull in their
own bundled JS/CSS) and stays portable as a single, self-contained HTML
file. Meant to be called directly inside an R Markdown chunk: the return
value knits as raw HTML automatically (see
[`asis_output`](https://rdrr.io/pkg/knitr/man/asis_output.html)) – no
`results='asis'` chunk option or manual
[`cat()`](https://rdrr.io/r/base/cat.html) needed.

## Usage

``` r
TALARIA_export_gene_table_widget(
  sce,
  columns = NULL,
  widget_id = "gene-table",
  max_rows = 200
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- columns:

  Character vector of `rowData(sce)` column names to expose. `NULL`
  (default) exposes every column.

- widget_id:

  Character. DOM id prefix – must be unique if more than one widget is
  embedded in the same document. Default `"gene-table"`.

- max_rows:

  Integer. Maximum number of matching rows rendered into the table at
  once. Default `200`.

## Value

An object of class `knit_asis` (see
[`asis_output`](https://rdrr.io/pkg/knitr/man/asis_output.html)).
Printing or auto-printing it inside an R Markdown chunk renders the
widget; outside of knitr it just prints as plain HTML text.

## Details

Only the currently matching rows are ever rendered into the DOM (capped
at `max_rows`), rather than every gene at once – with a full
transcriptome (tens of thousands of genes), building a table row for
every gene up front would bloat the page and make typing sluggish.
