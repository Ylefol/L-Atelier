# Plot soft threshold power selection diagnostics

Creates diagnostic plots for selecting the soft threshold power in
WGCNA. Shows scale-free topology fit and mean connectivity across
powers.

## Usage

``` r
AETHER_plot_wgcna_power(
  power_result,
  r2_cutoff = 0.8,
  highlight_power = NULL,
  return_plots = FALSE
)
```

## Arguments

- power_result:

  A wgcna_power object from ARTEMIS_wgcna_pick_power(), or a data.frame
  of fit indices with columns: Power, SFT.R.sq, slope, mean.k.

- r2_cutoff:

  R-squared threshold line to display. Default: 0.80.

- highlight_power:

  Power value to highlight. If NULL and power_result is a wgcna_power
  object, uses the selected power.

- return_plots:

  Logical. If TRUE, returns list of ggplot objects instead of plotting.
  Default: FALSE.

## Value

If return_plots = TRUE, returns list with scale_free and connectivity
plots. Otherwise, displays plots and returns invisible NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
AETHER_plot_wgcna_power(power_result)

} # }
```
