# Plot TSS Enrichment Profile and Per-sample Scores

Produces two complementary plots from a `hades_tss_enrichment` object
returned by
[`HADES_tss_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_tss_enrichment.md):

1.  **Profile plot**: normalised mean signal profile across all TSSs per
    sample (lines), with a vertical dashed line at the TSS and a
    horizontal dotted line at y = 1 (background level). A well-enriched
    library shows a sharp peak at the centre.

2.  **Score plot**: per-sample TSS enrichment score bar chart with an
    optional horizontal threshold line.

## Usage

``` r
AETHER_plot_tss_enrichment(
  tss_result,
  color_by = "condition",
  profile_alpha = 0.85,
  show_threshold = TRUE,
  title = NULL
)
```

## Arguments

- tss_result:

  A `hades_tss_enrichment` object from
  [`HADES_tss_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_tss_enrichment.md).

- color_by:

  Character. Column used to colour both plots. Can be `"sample"` (each
  sample a unique colour), `"pass"` (pass vs fail relative to the stored
  threshold), or the name of any metadata column carried through from
  the sample sheet (e.g. `"condition"`, `"batch"`). Falls back to
  `"sample"` with a message if the requested column is not present in
  `tss_result$scores`.

- profile_alpha:

  Numeric in (0, 1\]. Line transparency in the profile plot. Default
  `0.85`.

- show_threshold:

  Logical. Draw a horizontal reference line at the stored threshold on
  the score plot. Default `TRUE`.

- title:

  Character or `NULL`. Title for the profile plot. If `NULL` a default
  title is used.

## Value

A named list with two `ggplot` objects:

- `$profile`:

  Normalised signal profiles across the TSS window.

- `$score`:

  Per-sample TSS enrichment score bar chart.

Combine with `patchwork::wrap_plots(result, ncol = 1)` if desired.

## See also

[`HADES_tss_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_tss_enrichment.md)
