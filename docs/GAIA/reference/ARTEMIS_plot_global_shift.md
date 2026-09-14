# Plot Global Shift Results

Visualizes the results of ARTEMIS_global_shift_test() with a density
plot of log2 fold changes and a paired dot plot.

## Usage

``` r
ARTEMIS_plot_global_shift(
  shift_result,
  plot_type = "both",
  highlight_threshold = 1,
  title = "Global Shift Analysis"
)
```

## Arguments

- shift_result:

  Result from ARTEMIS_global_shift_test().

- plot_type:

  Character. Type of plot: "density" (log2FC distribution), "paired"
  (paired dot plot), or "both" (default).

- highlight_threshold:

  Numeric. log2FC threshold for highlighting strongly changed regions
  (default = 1, i.e., 2-fold).

- title:

  Character. Plot title (default = "Global Shift Analysis").

## Value

A ggplot object (or list of ggplot objects if plot_type = "both").
