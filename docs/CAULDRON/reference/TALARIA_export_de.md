# Export DE results to files

Writes CSV tables, an RDS object, and optional plots for a
`keraunos_markers`, `keraunos_contrast`, or `keraunos_pseudobulk`
object.

## Usage

``` r
TALARIA_export_de(
  de_result,
  output_dir,
  prefix = "de",
  plots = TRUE,
  top_n = 10L,
  l2fc_thresh = 1,
  plot_width = 8,
  plot_height = 6,
  verbose = TRUE
)
```

## Arguments

- de_result:

  A `keraunos_markers`, `keraunos_contrast`, or `keraunos_pseudobulk`
  object.

- output_dir:

  Character. Directory to write into (created if absent).

- prefix:

  Character. Filename prefix. Default `"de"`.

- plots:

  Logical. Generate volcano and MA plots (pseudobulk only). Default
  `TRUE`.

- top_n:

  Integer. Number of genes to label on each plot by significance.
  Requires `ggrepel` (Suggests). Default `10L`.

- l2fc_thresh:

  Numeric. \\\log_2\\ fold-change threshold for up/down categories on
  plots. Default `1`.

- plot_width:

  Numeric. PNG width in inches. Default `8`.

- plot_height:

  Numeric. PNG height in inches. Default `6`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

Invisibly, a character vector of file paths written.

## Details

Volcano and MA plots are generated for `keraunos_pseudobulk` only
(requires `log2FoldChange`, `pvalue`, and `baseMean` columns from
DESeq2). Marker and contrast results receive only CSV + RDS.
