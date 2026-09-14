# Plot k Selection Metrics

Simple line/bar plot showing silhouette scores across k values. Useful
for determining optimal number of clusters.

## Usage

``` r
AETHER_plot_k_selection(
  k_eval,
  highlight_optimal = TRUE,
  show_sizes = FALSE,
  title = "Silhouette Score by k"
)
```

## Arguments

- k_eval:

  An artemis_k_evaluation object from ARTEMIS_evaluate_k_range().

- highlight_optimal:

  Logical. Highlight the optimal k. Default = TRUE.

- show_sizes:

  Logical. Show cluster sizes as secondary plot. Default = FALSE.

- title:

  Character. Plot title. Default = "Silhouette Score by k".

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
k_eval <- ARTEMIS_evaluate_k_range(my_data)
AETHER_plot_k_selection(k_eval)

} # }
```
