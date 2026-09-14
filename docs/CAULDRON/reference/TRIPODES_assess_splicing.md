# Assess spliced/unspliced read distribution across cells and genes

Computes per-cell and per-gene splicing metrics and stores them in the
SCE. This is the detailed QC step to run after
[`TRIPODES_check_velocity_ready`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_check_velocity_ready.md)
and before
[`TRIPODES_run_velocity`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_velocity.md).
Results are consumed by `ASPIS_plot_splicing` for visualisation.

## Usage

``` r
TRIPODES_assess_splicing(
  sce,
  spliced_assay = "spliced",
  unspliced_assay = "unspliced",
  min_unspliced_mean = 0.01,
  ratio_range = c(0.05, 0.99),
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`; must pass `TRIPODES_check_velocity_ready`.

- spliced_assay:

  Name of the spliced counts assay (default `"spliced"`).

- unspliced_assay:

  Name of the unspliced counts assay (default `"unspliced"`).

- min_unspliced_mean:

  Minimum mean unspliced counts for a gene to be flagged
  `velocity_reliable = TRUE` (default `0.01`).

- ratio_range:

  Numeric vector of length 2; genes whose `gene_splicing_ratio` falls
  outside this range are flagged as unreliable (default
  `c(0.05, 0.99)`). Genes with near-zero or near-total unspliced
  fractions have poorly constrained dynamics.

- verbose:

  Logical; print summary report (default `TRUE`).

## Value

The input SCE with per-cell metrics added to `colData` and per-gene
metrics added to `rowData`. `metadata(sce)$splicing_assessed` is set to
`TRUE`.

## Details

Per-cell metrics added to `colData`:

- splicing_ratio:

  Spliced UMI / (spliced + unspliced) UMI per cell.

- total_spliced:

  Total spliced UMI count per cell.

- total_unspliced:

  Total unspliced UMI count per cell.

Per-gene metrics added to `rowData`:

- mean_spliced:

  Mean spliced counts across cells.

- mean_unspliced:

  Mean unspliced counts across cells.

- gene_splicing_ratio:

  mean_spliced / (mean_spliced + mean_unspliced).

- velocity_reliable:

  Logical — gene has sufficient unspliced signal for reliable velocity
  estimation (see `min_unspliced_mean` and `ratio_range`).
