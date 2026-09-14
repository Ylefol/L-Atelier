# Volcano-Style Plot for Anchor Analysis

Creates a plot showing log2 fold change vs mean signal strength. Useful
for identifying whether high-signal or low-signal regions show different
behavior.

## Usage

``` r
AETHER_plot_fc_vs_signal(
  shift_result,
  fc_threshold = 0.5,
  n_label = 5,
  title = "Fold Change vs Signal Strength"
)
```

## Arguments

- shift_result:

  Result from ARTEMIS_global_shift_test().

- fc_threshold:

  Numeric. log2FC threshold for highlighting (default = 0.5).

- n_label:

  Integer. Number of extreme points to label (default = 5).

- title:

  Character. Plot title.

## Value

A ggplot object.
