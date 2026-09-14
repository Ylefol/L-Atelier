# Plot a propeller cell-type proportion test result

Visualises the output of
[`KERAUNOS_propeller_proportions`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md):
one facet per cluster, showing each biological sample's proportion as a
point (coloured by group) with a crossbar at the group mean, plus the
test's FDR (or raw p-value) printed in the top-left corner of each
facet.

## Usage

``` r
ASPIS_plot_propeller(
  propeller_result,
  palette = NULL,
  ncol = NULL,
  point_size = 2.5,
  jitter_width = 0.08,
  show_stats = TRUE,
  stat_label = c("fdr", "pvalue"),
  order_by = c("significance", "cluster"),
  title = NULL
)
```

## Arguments

- propeller_result:

  A `keraunos_propeller` object from
  [`KERAUNOS_propeller_proportions`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md).

- palette:

  Character vector of colours, one per group level. Overrides any
  `group_colors` carried on `propeller_result`. Default `NULL`: uses
  `propeller_result$params$group_colors` if it was supplied to
  [`KERAUNOS_propeller_proportions`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md),
  otherwise falls back to the built-in 20-colour Tableau-inspired set
  (same as
  [`ASPIS_plot_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_umap.md)).

- ncol:

  Integer. Facets per row. Default `NULL` lets `facet_wrap` choose.

- point_size:

  Numeric. Sample point size. Default `2.5`.

- jitter_width:

  Numeric. Horizontal jitter applied to points so overlapping samples
  stay visible. Default `0.08`.

- show_stats:

  Logical. Annotate each facet with its test statistic. Default `TRUE`.

- stat_label:

  Character. `"fdr"` (default) or `"pvalue"` – which column from
  `propeller_result$results` to display.

- order_by:

  Character. `"significance"` (default) orders facets by ascending
  P.Value (most significant first, matching `propeller_result$results`);
  `"cluster"` orders alphabetically/ numerically by cluster label.

- title:

  Character. Plot title. Default `NULL`.

## Value

A ggplot object, faceted by cluster.

## Details

A boxplot is deliberately not used here – with as few as 2 samples per
group (as in designs like Nico_midbrain), a boxplot's quartiles are
degenerate and visually imply a distribution the data doesn't have.
Individual points plus a mean crossbar show exactly what was tested
without overstating it.
