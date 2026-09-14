# Export an interactive DE cross-comparison bar-chart widget for R Markdown

Builds a self-contained, dependency-free HTML/JS widget that lets a
report viewer search any gene that was called significant in at least
one of several differential expression comparisons, and see a bar chart
of that gene's log2 fold change across *every* comparison supplied –
including comparisons where the gene did not reach significance – so all
comparisons are visible side by side for one gene. Each bar is labelled
with its log2 fold change and its adjusted p-value. Mirrors
[`TALARIA_export_gene_expr_widget`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gene_expr_widget.md)'s
embedded-JSON + vanilla-JS approach so the report stays free of external
JS dependencies (a package like plotly or crosstalk would each pull in
their own bundled JS/CSS) and stays portable as a single, self-contained
HTML file.

## Usage

``` r
TALARIA_export_de_comparison_widget(
  de_list,
  sig_col = "padj",
  sig_thresh = 0.05,
  l2fc_thresh = 0.5,
  widget_id = "de-comparison",
  digits = 4
)
```

## Arguments

- de_list:

  Named list of data.frames, one per comparison (names become the
  comparison labels shown on the x-axis, in the order given). Each
  data.frame must have `gene`, `log2FoldChange`, and `padj` columns
  (plus `sig_col` itself, if different from `"padj"`) – the same shape
  written by
  [`TALARIA_export_de()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_de.md)'s
  per-cluster pseudobulk DE CSVs.

- sig_col:

  Character. Column used to decide which genes are searchable (i.e.
  significant in at least one comparison). Default `"padj"`. Set to
  `"pvalue"` for pipelines that call significance on the raw p-value
  instead of FDR (e.g.
  `KERAUNOS_de_pseudobulk(..., fdr_threshold = NULL)`) – the bar label
  always shows `padj` regardless of this setting.

- sig_thresh:

  Numeric. Cutoff applied to `sig_col`. Default `0.05`.

- l2fc_thresh:

  Numeric. Absolute log2 fold change cutoff used only to decide which
  genes are searchable. Default `0.5`.

- widget_id:

  Character. DOM id prefix – must be unique if more than one widget is
  embedded in the same document. Default `"de-comparison"`.

- digits:

  Integer. Rounding applied to exported log2FoldChange values (padj is
  kept to 3 significant figures regardless). Default `4`.

## Value

An object of class `knit_asis` (see
[`asis_output`](https://rdrr.io/pkg/knitr/man/asis_output.html)).
Printing or auto-printing it inside an R Markdown chunk renders the
widget; outside of knitr it just prints as plain HTML text.

## Details

A gene not tested in a given comparison (e.g. filtered out upstream by
`min_cells`/`min_samples`) is shown as a "not tested" marker at the zero
line for that comparison rather than a bar, so its absence isn't
mistaken for a log2FoldChange of exactly zero.

## Examples

``` r
if (FALSE) { # \dontrun{
TALARIA_export_de_comparison_widget(
  de_list = list(
    "SNc DA Neuron: AST23 vs CL21" = de_results[["SNc_DA_Neuron_CL21_vs_AST23"]]$results[["SNc_DA_Neuron"]],
    "SNc DA Neuron: AST23 vs CL18" = de_results[["SNc_DA_Neuron_CL18_vs_AST23"]]$results[["SNc_DA_Neuron"]],
    "DA Neuron: AST23 vs CL21"     = de_results[["DA_Neuron_CL21_vs_AST23"]]$results[["DA_Neuron"]],
    "DA Neuron: AST23 vs CL18"     = de_results[["DA_Neuron_CL18_vs_AST23"]]$results[["DA_Neuron"]]
  ),
  sig_col = "pvalue", sig_thresh = 0.05, l2fc_thresh = 0.5
)
} # }
```
