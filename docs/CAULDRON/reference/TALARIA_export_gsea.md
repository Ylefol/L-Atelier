# Export GSEA results to files

Writes CSV tables, an RDS object, NES dotplots (one per gene set
collection), and optional per-pathway enrichment score plots for a
`keraunos_gsea` or `keraunos_gsea_multi` object. Dotplots are controlled
by `make_dotplots`; enrichment score plots are controlled by
`do_gsea_plots`. CSV/RDS export always happens regardless of either
setting.

## Usage

``` r
TALARIA_export_gsea(
  result,
  output_dir,
  prefix = "gsea",
  make_dotplots = TRUE,
  do_gsea_plots = TRUE,
  top_n = 20L,
  plot_width = 8,
  plot_height = 6,
  enrich_plot_height = 4,
  verbose = TRUE
)
```

## Arguments

- result:

  A `keraunos_gsea` or `keraunos_gsea_multi` object.

- output_dir:

  Character. Directory to write into (created if absent).

- prefix:

  Character. Filename prefix. Default `"gsea"`.

- make_dotplots:

  Logical. Generate NES dotplots (one per gene set collection). Set
  `FALSE` when exporting many comparisons at once (e.g. all pairwise
  cluster combinations) to avoid writing one dotplot per comparison —
  CSV/RDS results are still written either way, so dotplots can be
  produced later for a chosen subset by calling this function again with
  `make_dotplots = TRUE` on just those results. Default `TRUE`.

- do_gsea_plots:

  Logical. Generate per-pathway enrichment score plots (running-score
  curves via
  [`fgsea::plotEnrichment`](https://rdrr.io/pkg/fgsea/man/plotEnrichment.html)).
  Default `TRUE`.

- top_n:

  Integer. Pathways shown per dotplot (split equally between enriched
  and depleted) and number of enrichment score plots per cluster.
  Default `20L`.

- plot_width:

  Numeric. PNG width in inches. Default `8`.

- plot_height:

  Numeric. PNG height for dotplots in inches. Default `6`.

- enrich_plot_height:

  Numeric. PNG height for enrichment score plots. Default `4`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

Invisibly, a character vector of file paths written.

## Details

For `keraunos_gsea_multi` each cluster gets its own subdirectory.
Enrichment score plots require `$ranked_genes` and `$gene_sets` stored
in the object (present automatically when produced by
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)
or
[`KERAUNOS_gsea_pseudobulk`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea_pseudobulk.md)).
