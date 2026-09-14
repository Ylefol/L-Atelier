# Select soft threshold power for WGCNA

Determines the optimal soft threshold power for network construction.
Can automatically select power or return diagnostics for manual
selection.

## Usage

``` r
ARTEMIS_wgcna_pick_power(
  wgcna_data,
  powers = NULL,
  power_auto = TRUE,
  r2_cutoff = 0.8,
  mean_k_cutoff = NULL,
  network_type = "unsigned",
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- wgcna_data:

  A wgcna_data object, wgcna_cluster object, or expression matrix
  (samples as rows).

- powers:

  Numeric vector of powers to test. Default: c(1:10, seq(12, 30, by =
  2)).

- power_auto:

  Logical. If TRUE, automatically select power. If FALSE, return
  diagnostics for manual selection. Default: TRUE.

- r2_cutoff:

  R-squared threshold for scale-free topology. Power is selected where
  R^2 first exceeds this value. Default: 0.80.

- mean_k_cutoff:

  Optional. If R^2 criterion yields high mean connectivity, may select
  power with lower connectivity. Default: NULL (disabled).

- network_type:

  Network type: "unsigned", "signed", or "signed hybrid". Default:
  "unsigned".

- plot:

  Logical. Generate diagnostic plots. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_power" containing:

- power:

  Selected power (or NA if auto-selection disabled/failed)

- power_estimate:

  WGCNA's built-in power estimate

- fit_indices:

  Data.frame of fit indices for all tested powers

- selection_method:

  How power was selected

- r2_at_power:

  R^2 value at selected power

- mean_k_at_power:

  Mean connectivity at selected power

- recommendation:

  Text recommendation

## Details

Power selection uses multiple criteria in order:

1.  WGCNA's built-in powerEstimate (if available and meets R^2
    threshold)

2.  First power where scale-free topology R^2 exceeds r2_cutoff

3.  Power with best R^2 if none exceed threshold (with warning)

The diagnostic plots show:

- Scale-free topology fit (R^2) vs power

- Mean connectivity vs power

## Examples

``` r
if (FALSE) { # \dontrun{
wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
# Use: power_result$power

# Manual selection (just show plots, don't auto-select)
power_result <- ARTEMIS_wgcna_pick_power(wgcna_data, power_auto = FALSE)

} # }
```
