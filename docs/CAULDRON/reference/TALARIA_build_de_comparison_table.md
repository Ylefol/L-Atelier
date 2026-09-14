# Build a browsable table of DE results across several comparisons

Returns a wide data.frame – one row per gene significant in at least one
comparison, with a log2FoldChange and padj column pair per comparison –
meant to be rendered with
[`DT::datatable()`](https://rdrr.io/pkg/DT/man/datatable.html) directly
above
[`TALARIA_export_de_comparison_widget`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_de_comparison_widget.md)
in the same report, so a viewer can browse (sort/search) the full list
of searchable genes before picking one to plot in the widget. Uses the
same significance/union logic as that function, so the two stay in sync
when called with matching `sig_col`/`sig_thresh`/`l2fc_thresh`.

## Usage

``` r
TALARIA_build_de_comparison_table(
  de_list,
  sig_col = "padj",
  sig_thresh = 0.05,
  l2fc_thresh = 0.5,
  digits = 3
)
```

## Arguments

- de_list:

  Named list of data.frames, one per comparison (names become the
  comparison shown in each column pair's name, in the order given). Each
  data.frame must have `gene`, `log2FoldChange`, and `padj` columns
  (plus `sig_col` itself, if different from `"padj"`) – the same shape
  written by
  [`TALARIA_export_de()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_de.md)'s
  per-cluster pseudobulk DE CSVs.

- sig_col:

  Character. Column used to decide which genes are included (i.e.
  significant in at least one comparison). Default `"padj"`. Set to
  `"pvalue"` for pipelines that call significance on the raw p-value
  instead of FDR.

- sig_thresh:

  Numeric. Cutoff applied to `sig_col`. Default `0.05`.

- l2fc_thresh:

  Numeric. Absolute log2 fold change cutoff used to decide which genes
  are included. Default `0.5`.

- digits:

  Integer. Rounding applied to log2FoldChange columns (padj columns are
  kept to 3 significant figures regardless). Default `3`.

## Value

A data.frame with one `Gene` column plus a `"<comparison> log2FC"` /
`"<comparison> padj"` column pair per comparison in `de_list`, in the
order given.

## Examples

``` r
if (FALSE) { # \dontrun{
tbl <- TALARIA_build_de_comparison_table(de_list, sig_col = "pvalue")
DT::datatable(tbl, rownames = FALSE)
} # }
```
