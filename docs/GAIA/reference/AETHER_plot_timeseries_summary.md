# Plot Time Series DEA Summary

Creates a bar chart showing the number of up- and down-regulated genes
per comparison from a time series DEA result. Up-regulated genes are
shown above the axis, down-regulated below.

## Usage

``` r
AETHER_plot_timeseries_summary(ts_de_result, colors = NULL, title = NULL)
```

## Arguments

- ts_de_result:

  An `artemis_ts_de` object from
  [`ARTEMIS_timeseries_conditional()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_conditional.md)
  or
  [`ARTEMIS_timeseries_temporal()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_timeseries_temporal.md).

- colors:

  Named character vector of length 2: c(up = ..., down = ...). Default:
  red for up, blue for down.

- title:

  Character or NULL. Plot title. Default: NULL (auto-generated).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
cond_de <- ARTEMIS_timeseries_conditional(counts, targets, "Ctrl", "Treat")
AETHER_plot_timeseries_summary(cond_de)

} # }
```
