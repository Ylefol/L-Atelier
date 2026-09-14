# Plot activity vs trait correlation

Scatter plot showing correlation between activity scores and a sample
trait.

## Usage

``` r
AETHER_plot_activity_trait(
  result,
  trait,
  source = NULL,
  method = "pearson",
  show_labels = FALSE,
  label_size = 3,
  title = NULL
)
```

## Arguments

- result:

  A decoupler_result object.

- trait:

  Numeric vector of trait values (one per sample), OR a data.frame with
  a single numeric column.

- source:

  Character. Which source (TF/pathway) to plot. If NULL, shows the
  source with highest correlation. Default: NULL.

- method:

  Character. Correlation method: "pearson" or "spearman". Default:
  "pearson".

- show_labels:

  Logical. If TRUE, labels points with sample names using ggrepel for
  automatic non-overlapping positioning. Default: FALSE.

- label_size:

  Numeric. Size of point labels. Default: 3.

- title:

  Character. Plot title. Default: auto-generated.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot most correlated TF
AETHER_plot_activity_trait(tf_result, trait = sample_ages)

# Specific TF with labels
AETHER_plot_activity_trait(tf_result, trait = sample_ages, source = "TP53",
                           show_labels = TRUE)

# With non-overlapping labels (requires ggrepel)
AETHER_plot_activity_trait(tf_result, trait = sample_ages, show_labels = "repel")

} # }
```
