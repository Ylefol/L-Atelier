# Filter low-quality cells

Removes cells that fail QC thresholds, returning a filtered
`SingleCellExperiment`.

## Usage

``` r
AEGIS_filter_cells(
  sce,
  mode = c("adaptive", "fixed"),
  n_mads = 3,
  min_counts = NULL,
  max_counts = NULL,
  min_features = NULL,
  max_features = NULL,
  max_pct_mt = NULL,
  filter_mt = TRUE,
  remove_doublets = TRUE,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with QC columns in `colData`.

- mode:

  Character. `"adaptive"` (default) or `"fixed"`.

- n_mads:

  Numeric. Number of MADs for adaptive outlier detection. Default `3`.

- min_counts, max_counts:

  Numeric. Hard bounds on total counts (`sum`).

- min_features, max_features:

  Numeric. Hard bounds on detected genes (`detected`).

- max_pct_mt:

  Numeric. Hard ceiling on \\ Applied in both adaptive and fixed modes
  when provided, unless `filter_mt = FALSE`.

- filter_mt:

  Logical. Whether to apply mitochondrial filtering. When `FALSE`, MT\\
  and any `max_pct_mt` hard ceiling is ignored. Default `TRUE`. Set to
  `FALSE` for protocols where MT\\ indicator (e.g. PARSE Biosciences
  fixed-cell experiments).

- remove_doublets:

  Logical. Also remove cells classified as doublets by
  `scDblFinder.class`. Default `TRUE`.

- verbose:

  Logical. Print a removal summary. Default `TRUE`.

## Value

Filtered `SingleCellExperiment`.

## Details

Two modes are supported:

- `"adaptive"`:

  Outliers are identified per metric using median-absolute-deviation
  (MAD) via
  [`perCellQCFilters`](https://rdrr.io/pkg/scuttle/man/perCellQCFilters.html).
  Cells more than `n_mads` MADs below the median for `sum` or
  `detected`, or above the median for MT\\ adapts to dataset depth and
  is the recommended default.

- `"fixed"`:

  Hard thresholds supplied by the user. At least one of `min_counts`,
  `max_counts`, `min_features`, `max_features`, or `max_pct_mt` must be
  provided.

In either mode, `max_pct_mt` acts as an additional hard ceiling if
supplied alongside `mode = "adaptive"`, but only when
`filter_mt = TRUE`.

Requires
[`AEGIS_compute_qc_metrics`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_compute_qc_metrics.md)
to have been run first. If `remove_doublets = TRUE`, also requires
[`AEGIS_detect_doublets`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_detect_doublets.md).
