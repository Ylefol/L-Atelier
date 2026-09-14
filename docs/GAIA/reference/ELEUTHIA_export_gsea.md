# Export GSEA results to files

Exports fgsea results (from
[`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md))
to CSV files and optionally generates an NES dotplot and per-pathway
running-score enrichment plots for significant pathways.

## Usage

``` r
ELEUTHIA_export_gsea(
  gsea_result,
  output_dir,
  prefix = "gsea",
  plot_top_n = 20,
  save_plots = TRUE,
  do_enrichment_plots = TRUE,
  max_enrichment_plots = 20,
  plot_format = "png",
  plot_width = 10,
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- gsea_result:

  An `apollo_gsea` object from
  [`APOLLO_gsea()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_gsea.md).

- output_dir:

  Output directory path.

- prefix:

  Prefix for output filenames. Default: "gsea".

- plot_top_n:

  Integer. Number of top pathways (by padj, split enriched/depleted)
  shown in the dotplot. Default: 20.

- save_plots:

  Logical. Generate dotplot and enrichment plots. Default: TRUE.

- do_enrichment_plots:

  Logical. Generate per-pathway running-score plots for significant
  pathways. Default: TRUE.

- max_enrichment_plots:

  Integer. Cap on the number of per-pathway plots generated (ordered by
  padj), to avoid producing hundreds of files when many pathways are
  significant. Default: 20.

- plot_format:

  Character. "png", "pdf", or "both". Default: "png".

- plot_width:

  Numeric. Dotplot width in inches. Default: 10.

- save_rds:

  Logical. Save full apollo_gsea object as RDS. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Creates the following files:

- \_results.csv: full fgsea results (all tested gene sets)

- \_significant.csv: results filtered to the FDR threshold used in
  APOLLO_gsea()

- \_plots/dotplot.: NES dotplot (top plot_top_n pathways)

- *plots/enrichment*.: running-score plot per significant pathway (up to
  max_enrichment_plots), if do_enrichment_plots = TRUE

- .rds: full apollo_gsea object (if save_rds = TRUE)

## Examples

``` r
if (FALSE) { # \dontrun{
ranked <- APOLLO_rank_from_de(de, rank_by = "t")
gsea_result <- APOLLO_gsea(ranked, collection = c("H", "C2:CP:REACTOME"))
ELEUTHIA_export_gsea(gsea_result, "results/gsea")
} # }
```
