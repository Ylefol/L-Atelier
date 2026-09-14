# Plot method agreement / variability

Shows mean activity with min-max range across methods for each source.
Visualizes which TFs/pathways have consistent results across methods vs
those where methods disagree.

## Usage

``` r
AETHER_plot_method_agreement(
  comparison,
  top_n = 30,
  order_by = "activity",
  scale = c("raw", "zscore"),
  title = "Activity Across Methods",
  point_size = 2,
  colors = c("#B2182B", "#2166AC")
)
```

## Arguments

- comparison:

  A decoupler_comparison object.

- top_n:

  Integer. Number of sources to show. Default: 30.

- order_by:

  Character. How to order sources: "activity" (absolute mean),
  "variance" (SD across methods), or "name". Default: "activity".

- scale:

  Character. Which per-source summary to plot: `"raw"` (each method's
  native activity scale) or `"zscore"` (each method's activity matrix
  z-scored before the mean/range are computed, i.e.
  `comparison$summary_zscored`). Default: `"raw"`, for backward
  compatibility. Use `"zscore"` whenever the compared methods differ in
  native output scale (e.g. wsum's unbounded weighted-sum scores vs.
  ulm/mlm's t-like statistics vs. viper's NES-like scores) – on `"raw"`,
  the method with the largest scale will dominate both the top_n ranking
  and the error-bar range regardless of actual cross-method agreement,
  which `"zscore"` corrects for.

- title:

  Character. Plot title. Default: "Activity Across Methods".

- point_size:

  Numeric. Size of mean activity points. Default: 2.

- colors:

  Character vector of length 2. Colors for positive and negative mean
  activities. Default: c("#B2182B", "#2166AC").

## Value

A ggplot object.

## Details

The plot shows:

- Point: mean activity across methods

- Error bars: range (min to max) across methods

- Color: whether mean activity is positive (red) or negative (blue)

Sources where methods agree will have small error bars. Sources where
methods disagree (or even disagree on direction) will have large bars
crossing zero.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_method_agreement(comparison)
AETHER_plot_method_agreement(comparison, order_by = "variance")

# Methods on very different native scales (e.g. wsum alongside ulm/mlm/viper)
AETHER_plot_method_agreement(comparison, scale = "zscore")

} # }
```
