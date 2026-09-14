# Plot Feature Weights from Multi-Dataset Sparse CCA

Creates faceted stem plots showing canonical weights from multi-dataset
sCCA. Each dataset gets its own panel. Highlights non-zero features and
optionally shows true signal regions. Prints summary statistics for each
dataset.

## Usage

``` r
AETHER_plot_multi_feature_weights(
  weights_list,
  method_name = "sCCA",
  dataset_names = NULL,
  true_features_list = NULL,
  threshold = 1e-06
)
```

## Arguments

- weights_list:

  List of numeric vectors, each containing canonical weights for one
  dataset. Length of list = number of datasets integrated.

- method_name:

  Character string for plot title (default = "sCCA")

- dataset_names:

  Optional character vector of dataset names for panel labels. If NULL,
  uses "Dataset 1", "Dataset 2", etc.

- true_features_list:

  Optional list of numeric vectors, each c(start, end) indicating true
  signal region for that dataset. NULL elements mean no highlighting for
  that dataset. Useful for simulated data validation.

- threshold:

  Threshold for considering weights as non-zero (default = 1e-6)

## Value

ggplot2 object with faceted panels. Also prints summary statistics as
side effect.

## Details

Creates faceted stem plots (one panel per dataset) where:

- Non-zero weights (\|weight\| \> threshold) are shown in orange

- Zero/negligible weights are shown in gray

- True signal regions (if specified) are highlighted with blue
  background

- Each panel shows one dataset's canonical weights

- Prints for each dataset: number of non-zero features, max weight, mean
  absolute weight

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate multi-dataset sCCA results
weights1 <- c(rnorm(20, mean = 0.5), rep(0, 80))
weights2 <- c(rnorm(15, mean = 0.6), rep(0, 85))
weights3 <- c(rnorm(10, mean = 0.4), rep(0, 40))

# Basic plot
p <- AETHER_plot_multi_feature_weights(
  weights_list = list(weights1, weights2, weights3),
  method_name = "multi.convCCA",
  dataset_names = c("RNA-seq", "ATAC-seq", "CUT&TAG")
)

# Plot with true signal regions highlighted
p <- AETHER_plot_multi_feature_weights(
  weights_list = list(weights1, weights2, weights3),
  method_name = "multi.convCCA",
  dataset_names = c("RNA-seq", "ATAC-seq", "CUT&TAG"),
  true_features_list = list(c(1, 20), c(1, 15), c(1, 10))
)

} # }
```
