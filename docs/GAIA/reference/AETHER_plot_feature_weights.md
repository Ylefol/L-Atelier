# Plot Feature Weights from Sparse CCA

Creates a stem plot visualization of canonical weights from sCCA,
highlighting non-zero features and optionally showing true signal
regions. Prints summary statistics to console.

## Usage

``` r
AETHER_plot_feature_weights(
  weights,
  method_name = "sCCA",
  true_features = NULL,
  threshold = 1e-06
)
```

## Arguments

- weights:

  Numeric vector of canonical weights (from sCCA result)

- method_name:

  Character string for plot title (default = "sCCA")

- true_features:

  Optional numeric vector c(start, end) indicating true signal region.
  If provided, region will be highlighted with blue shading. Useful for
  simulated data.

- threshold:

  Threshold for considering weights as non-zero (default = 1e-6)

## Value

ggplot2 object. Also prints summary statistics as side effect.

## Details

Creates a stem plot where:

- Non-zero weights (\|weight\| \> threshold) are shown in orange

- Zero/negligible weights are shown in gray

- True signal regions (if specified) are highlighted with blue
  background

- Prints: number of non-zero features, max weight, mean absolute weight

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate sCCA result
weights <- c(rnorm(20, mean = 0.5), rep(0, 80))

# Basic plot
p <- AETHER_plot_feature_weights(weights, method_name = "ConvCCA")

# Plot with true signal region highlighted
p <- AETHER_plot_feature_weights(weights, method_name = "ConvCCA",
                          true_features = c(1, 20))

} # }
```
