# Export Olink QC Results

Writes all outputs from the Olink QC/exploration pipeline to a
structured directory: filtered data as RDS and CSVs, QC plots, a
PC-metadata association heatmap, per-variable PCA plots (coloured by the
two PCs most significantly associated with each variable), and a
plain-text QC summary.

## Usage

``` r
ELEUTHIA_export_olink_qc(
  olink_data,
  assoc,
  qc_plots,
  output_dir,
  top_n_vars = 5L,
  p_threshold = 0.05,
  heatmap_top_n = NULL,
  ntop_pca = NULL,
  batch_plots = NULL,
  qc_steps = NULL,
  verbose = TRUE
)
```

## Arguments

- olink_data:

  An `olink_data` object — the filtered, QC-cleaned object to export.

- assoc:

  An `artemis_pc_assoc` object from
  [`ARTEMIS_pc_metadata_association()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_pc_metadata_association.md).

- qc_plots:

  A named list from
  [`AETHER_plot_olink_qc()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_olink_qc.md).

- output_dir:

  Character. Root directory for all outputs. Created if it does not
  exist.

- top_n_vars:

  Integer. Number of top variables per PC (by p-value) to include in the
  PCA plot loop. Default: 5.

- p_threshold:

  Numeric. Significance threshold used in the association heatmap border
  annotation and in the summary file. Default: 0.05.

- heatmap_top_n:

  Integer or NULL. Passed to `AETHER_plot_pc_association(top_n)` to
  restrict the heatmap to the top N variables per PC. NULL shows all
  tested variables. Default: NULL.

- ntop_pca:

  Integer or NULL. Number of most variable proteins used for PCA. NULL
  uses all proteins. Default: NULL.

- batch_plots:

  Named list of ggplot objects for batch-check PCA plots. Each element
  is saved as `qc/pca_batch_{name}.png`. Supply one plot coloured by
  PlateID for a standard batch check, or a before/after pair when batch
  correction has been applied. Default: NULL (no plots saved).

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisibly returns NULL. All outputs are written to disk.

## Details

**Directory layout:**

    output_dir/
      data/
        olink_filtered.rds   -- filtered olink_data object
        sample_meta.csv      -- per-sample metadata (with outlier annotations)
        assay_meta.csv       -- per-protein metadata (warn_fraction etc.)
      qc/
        npx_distributions.png
        sample_qc.csv        -- SampleQC pass/fail table
        warn_proteins.png    -- only if present in qc_plots
        summary.txt          -- key QC numbers and parameter record
      association/
        pc_metadata_heatmap.png
      pca/
        <variable>.png       -- one per selected variable

Each PCA plot uses the two PCs most significantly associated with that
variable (lowest p-value in `assoc$pvalues`). The PC pair is shown in
the plot title.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- HADES_filter_olink(ol_raw)
ol <- HADES_filter_olink_proteins(ol)
ol <- HADES_detect_outliers_olink(ol, filter = TRUE)

qc_plots <- AETHER_plot_olink_qc(ol, plot_warn_proteins = TRUE)
assoc    <- ARTEMIS_pc_metadata_association(ol$wide, ol$sample_meta,
                                             sample_col = "SampleID")

ELEUTHIA_export_olink_qc(ol, assoc, qc_plots, output_dir = "results/Olink_QC")
} # }
```
