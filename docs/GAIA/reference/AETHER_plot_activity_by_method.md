# Plot activity comparison across methods

Creates a faceted plot showing activity distributions for each method.

## Usage

``` r
AETHER_plot_activity_by_method(
  comparison,
  sources = NULL,
  top_n = 10,
  plot_type = "boxplot",
  title = "Activity by Method"
)
```

## Arguments

- comparison:

  A decoupler_comparison object.

- sources:

  Character vector. Specific sources to plot. If NULL, uses top sources
  by mean activity. Default: NULL.

- top_n:

  Integer. If sources is NULL, how many top sources to show. Default:
  10.

- plot_type:

  Character. "boxplot" or "violin". Default: "boxplot".

- title:

  Character. Plot title. Default: "Activity by Method".

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Top 10 sources across methods
AETHER_plot_activity_by_method(comparison)

# Specific sources
AETHER_plot_activity_by_method(comparison, sources = c("TP53", "MYC", "STAT3"))

} # }
```
