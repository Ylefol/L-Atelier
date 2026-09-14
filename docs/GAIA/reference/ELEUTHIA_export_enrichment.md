# Export enrichment results to files

Exports gprofiler2 enrichment results (from APOLLO_enrich_gost) to CSV
files and optionally generates dotplots for each source.

## Usage

``` r
ELEUTHIA_export_enrichment(
  enrichment,
  output_dir,
  prefix = "enrichment",
  by_source = FALSE,
  by_module = TRUE,
  save_plots = TRUE,
  save_full_gostplot = TRUE,
  plot_top_n = 10,
  plot_format = "png",
  plot_width = 10,
  plot_height_per_module = 3,
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- enrichment:

  A gost_enrichment object from APOLLO_enrich_gost().

- output_dir:

  Output directory path.

- prefix:

  Prefix for output filenames. Default: "enrichment".

- by_source:

  Logical. Create separate files per source (GO:BP, KEGG, etc.).
  Default: FALSE.

- by_module:

  Logical. Create separate files per module/cluster (all sources
  combined per query). Unlike `combined`/`by_source`, these hold EVERY
  evaluated term (not just significance-filtered ones) with a
  `significant` column to filter on, since `enrichment$results` already
  stores the full unfiltered per-module gost output. Default: TRUE.

- save_plots:

  Logical. Generate dotplots per source. Default: TRUE.

- save_full_gostplot:

  Logical. Also generate one traditional, interactive g:GOSt Manhattan
  plot per module (saved as HTML), showing ALL evaluated terms
  regardless of significance. Unlike the significance-filtered dotplots,
  this is generated even for modules with zero significant terms.
  Default: TRUE.

- plot_top_n:

  Integer. Number of top terms per module in dotplots. Default: 10.

- plot_format:

  Character. Plot file format: "png", "pdf", or "both". Default: "png".

- plot_width:

  Numeric. Plot width in inches. Default: 10.

- plot_height_per_module:

  Numeric. Plot height per module in inches. Default: 3.

- save_rds:

  Logical. Save full R object as RDS. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible list of file paths created.

## Details

Creates the following files:

- \_combined.csv: Significance-filtered results in one table

- \_summary.csv: Summary counts per module/source

- \_by_module/.csv: ALL evaluated terms per module/cluster, unfiltered
  for significance (has a `significant` column) (default)

- \_by_source/.csv: Significance-filtered results split by source (if
  by_source = TRUE)

- *plots/dotplot*.png: Dotplot per source (if save_plots = TRUE)

- *plots/gostplot_full*.html: Interactive, unfiltered g:GOSt Manhattan
  plot per module (if save_full_gostplot = TRUE)

- .rds: Full gost_enrichment object (if save_rds = TRUE)

The "module" column contains whatever labels were used in the original
query (e.g., "C1", "C2" for PART clusters, "blue", "turquoise" for WGCNA
modules).

## Examples

``` r
if (FALSE) { # \dontrun{
enrich <- APOLLO_enrich_gost(module_genes)
ELEUTHIA_export_enrichment(enrich, "results/enrichment")

# Split by source instead of module
ELEUTHIA_export_enrichment(enrich, "results/enrichment",
                           by_module = FALSE, by_source = TRUE)

} # }
```
