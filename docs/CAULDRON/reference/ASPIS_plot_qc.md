# QC metric violin plots

Visualises the distribution of per-cell quality control metrics computed
by
[`AEGIS_compute_qc_metrics`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_compute_qc_metrics.md).
Each metric is shown as a violin with an optional jitter overlay,
arranged in a faceted grid. Threshold lines can be drawn to inspect
prospective filtering cut-offs before calling
[`AEGIS_filter_cells`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_filter_cells.md).

## Usage

``` r
ASPIS_plot_qc(
  sce,
  metrics = c("sum", "detected", "subsets_mt_percent"),
  group_by = NULL,
  thresholds = NULL,
  log_scale = c("sum", "detected"),
  show_points = TRUE,
  max_points = 5000L,
  point_size = 0.3,
  point_alpha = 0.3,
  ncol = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with QC columns in `colData` (run
  [`AEGIS_compute_qc_metrics`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_compute_qc_metrics.md)
  first).

- metrics:

  Character vector of `colData` column names to plot. Default
  `c("sum", "detected", "subsets_mt_percent")`. Any numeric `colData`
  column is accepted.

- group_by:

  Character. A `colData` column used to split violins by group (e.g.
  `"sample_id"`, `"condition"`). `NULL` (default) shows a single violin
  per metric.

- thresholds:

  Named list of threshold values to draw as red dashed horizontal lines.
  Names must match entries in `metrics` (on the original, untransformed
  scale). Example: `list(subsets_mt_percent = 20, sum = 500)`. `NULL`
  (default) draws no lines.

- log_scale:

  Character vector of metric names to display on a log10 scale. Values
  are transformed as `log10(x + 1)` and axis labels updated accordingly.
  Default `c("sum", "detected")`.

- show_points:

  Logical. Overlay individual cell points as jitter. Automatically
  suppressed when the number of cells exceeds `max_points` to avoid
  overplotting. Default `TRUE`.

- max_points:

  Integer. Maximum number of cells for which jitter points are drawn. If
  `ncol(sce) > max_points`, a random subsample of `max_points` cells is
  shown. Default `5000`.

- point_size:

  Numeric. Jitter point size. Default `0.3`.

- point_alpha:

  Numeric. Jitter point transparency. Default `0.3`.

- ncol:

  Integer or `NULL`. Number of columns in the facet grid. `NULL`
  (default) uses `min(length(metrics), 3)`.

## Value

A `ggplot` object.

## Details

The following colData columns are produced by
[`AEGIS_compute_qc_metrics`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_compute_qc_metrics.md)
and can be passed to `metrics`:

- `sum`:

  Total UMI count per cell (library size). Low values indicate empty
  droplets or dead cells; extremely high values may indicate doublets.
  Typically inspected on a log scale.

- `detected`:

  Number of genes with at least one count. Follows a similar
  distribution to `sum` but is less sensitive to a few highly expressed
  genes. Typically inspected on a log scale.

- `subsets_mt_percent`:

  Percentage of counts from mitochondrial genes. High values suggest
  compromised cell membranes (cytoplasmic RNA lost, mitochondrial RNA
  retained). Inspected on a linear scale.

- `subsets_ribo_percent`:

  Percentage of counts from ribosomal protein genes. Unusually high
  values may indicate stressed or proliferating cells. Only present if
  ribosomal genes were detected.
