# Export HOMER Motif Enrichment Results

Saves HOMER motif enrichment results to CSV files, a metadata summary,
and optional plots. Works with both single-set (`homer_motif`) and batch
(`homer_motif_batch`) objects.

## Usage

``` r
ELEUTHIA_export_homer_results(
  homer_result,
  output_dir,
  prefix = "homer",
  top_n = 15,
  q_thresh = 0.05,
  save_plots = TRUE,
  plot_format = "png",
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- homer_result:

  A `homer_motif` or `homer_motif_batch` object from
  [`APOLLO_homer_motif_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_homer_motif_enrichment.md),
  [`APOLLO_homer_motif_enrichment_batch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_homer_motif_enrichment_batch.md),
  or
  [`APOLLO_load_homer_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_load_homer_results.md)
  /
  [`APOLLO_load_homer_batch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_load_homer_batch.md).

- output_dir:

  Character. Directory to save results. Created if needed.

- prefix:

  Character. Prefix for output filenames. Default: "homer".

- top_n:

  Integer. Top N motifs for plots (passed to plot functions). Default:
  15.

- q_thresh:

  Numeric. q-value threshold for plot filtering. Default: 0.05.

- save_plots:

  Logical. Save diagnostic plots. Default: TRUE.

- plot_format:

  Character. Plot format: "png", "pdf", or "both". Default: "png".

- save_rds:

  Logical. Save full R object as RDS. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Exports the following files:

- **Known motifs CSV** — all known motif results
  (`{prefix}_known_motifs.csv`)

- **De novo motifs CSV** — de novo motifs if available
  (`{prefix}_denovo_motifs.csv`)

- **Metadata TXT** — run parameters and summary
  (`{prefix}_metadata.txt`)

- **Dotplot** — top enriched motifs (via `AETHER_plot_motif_enrichment`)

- **Heatmap** — motif enrichment across peak sets, batch only (via
  `AETHER_plot_motif_comparison`)

- **RDS** — full object for reloading

For batch objects, the known motifs CSV contains the combined table with
a `peak_set` column identifying each peak set.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single peak set
ELEUTHIA_export_homer_results(homer_res, "results/homer/")

# Batch
ELEUTHIA_export_homer_results(homer_batch, "results/homer/",
                              prefix = "atac_homer")
} # }
```
