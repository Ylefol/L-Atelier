# Export ATAC-seq Analysis Results

Orchestrating export function for the standard ATAC-seq analysis
pipeline. Delegates to dedicated exporters for each result type; any
`NULL` argument is silently skipped.

## Usage

``` r
ELEUTHIA_export_atac_results(
  output_dir,
  annotated_peaks = NULL,
  dea_result = NULL,
  enrichment = NULL,
  homer_result = NULL,
  sequence_composition = NULL,
  tss_enrichment = NULL,
  tss_color_by = "group",
  prefix = "atac",
  l2fc_thresh = 1,
  p_thresh = 0.05,
  sample_info = NULL,
  group_col = "group",
  save_plots = TRUE,
  plot_format = "png",
  save_sequences = FALSE,
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- output_dir:

  Character. Root directory to save all results. Created if needed.

- annotated_peaks:

  Named list of annotated peak data.frames (e.g., one entry per group:
  `list(WT = wt_peaks, KO = ko_peaks)`). Each data.frame is saved as
  `{name}_annotated_peaks.csv` in the `peaks/` subdirectory. Default:
  NULL.

- dea_result:

  DEA result from
  [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md).
  Passed to
  [`ELEUTHIA_export_dea_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_dea_results.md).
  Default: NULL.

- enrichment:

  Enrichment result from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md).
  Passed to
  [`ELEUTHIA_export_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_enrichment.md).
  Default: NULL.

- homer_result:

  A `homer_motif` or `homer_motif_batch` object. Passed to
  [`ELEUTHIA_export_homer_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_homer_results.md).
  Default: NULL.

- sequence_composition:

  A data.frame from
  [`APOLLO_sequence_composition()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_sequence_composition.md).
  Passed to
  [`ELEUTHIA_export_sequence_composition()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_sequence_composition.md).
  Default: NULL.

- tss_enrichment:

  A `hades_tss_enrichment` object from
  [`HADES_tss_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_tss_enrichment.md).
  Scores are saved as CSV; profile and score plots are saved via
  [`AETHER_plot_tss_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_tss_enrichment.md).
  Default: NULL.

- tss_color_by:

  Character. `color_by` argument passed to
  [`AETHER_plot_tss_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_tss_enrichment.md).
  Any sample sheet metadata column carried through by
  [`HADES_tss_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_tss_enrichment.md)
  can be used. Default: `"group"`.

- prefix:

  Character. Filename prefix for DEA, enrichment, HOMER, and composition
  exports. Default: "atac".

- l2fc_thresh:

  Numeric. log2FC threshold for DEA significance. Default: 1.0.

- p_thresh:

  Numeric. Adjusted p-value threshold for DEA significance. Default:
  0.05.

- sample_info:

  Data.frame with sample metadata for DEA heatmap. Default: NULL.

- group_col:

  Character. Group column in `sample_info`. Default: `"group"`.

- save_plots:

  Logical. Generate and save plots for all result types. Default: TRUE.

- plot_format:

  Character. "png", "pdf", or "both". Default: "png".

- save_sequences:

  Logical. Save raw sequences from sequence composition to a separate
  file. Default: FALSE.

- save_rds:

  Logical. Save full R objects as RDS for all result types. Default:
  TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible named list of character vectors, one element per result type
exported (e.g., `$peaks`, `$dea`, `$homer`). Each element contains the
file paths created by the corresponding exporter.

## Details

Output subdirectory layout:

    output_dir/
      peaks/                        # annotated peak CSVs + annotation barplot
      dea/                          # DEA results, plots, RDS
        plots/                      # volcano, MA, PART heatmap
      enrichment/                   # enrichment results (always top-level)
      homer/                        # HOMER motifs, plots, RDS
      composition/                  # sequence composition CSV + summary
      tss_enrichment/               # TSS enrichment scores CSV, plots, RDS

## Examples

``` r
if (FALSE) { # \dontrun{
ELEUTHIA_export_atac_results(
  output_dir        = "results/atac/",
  annotated_peaks   = list(WT = wt_peaks, KO = ko_peaks),
  dea_result        = dea_res,
  enrichment        = enrich_res,
  homer_result      = homer_batch,
  sequence_composition = comp,
  prefix            = "smug1_atac"
)
} # }
```
