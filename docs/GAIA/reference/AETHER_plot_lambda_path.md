# Plot Lambda Selection Path

Visualizes the cross-validation error across lambda values, highlighting
lambda.min and lambda.1se choices.

## Usage

``` r
AETHER_plot_lambda_path(
  cv_result,
  show_1se = TRUE,
  title = "Cross-Validation for Lambda Selection"
)
```

## Arguments

- cv_result:

  Output from
  [`ARTEMIS_select_lambda()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_select_lambda.md).

- show_1se:

  Logical. Show lambda.1se in addition to lambda.min (default = TRUE).

- title:

  Character. Plot title.

## Value

A ggplot object.
