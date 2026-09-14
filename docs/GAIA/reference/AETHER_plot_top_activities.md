# Plot top activities as bar chart

Creates a bar chart showing top positive and negative activities. Useful
for visualizing which TFs/pathways are most active.

## Usage

``` r
AETHER_plot_top_activities(
  result,
  sample = NULL,
  top_n = 15,
  title = NULL,
  x_label = "Activity Score",
  fill_colors = c("#B2182B", "#2166AC")
)
```

## Arguments

- result:

  A decoupler_result object or activity matrix.

- sample:

  Character or integer. Which sample to plot. If NULL and multiple
  samples, uses mean across samples. Default: NULL.

- top_n:

  Integer. Number of top positive and negative to show. Default: 15.

- title:

  Character. Plot title. Default: auto-generated.

- x_label:

  Character. X-axis label. Default: "Activity Score".

- fill_colors:

  Character vector of length 2. Colors for positive and negative
  activities. Default: c("#B2182B", "#2166AC").

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mean across all samples
AETHER_plot_top_activities(tf_result)

# Specific sample
AETHER_plot_top_activities(tf_result, sample = "sample1")

} # }
```
