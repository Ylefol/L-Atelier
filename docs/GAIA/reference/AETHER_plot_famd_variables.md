# Plot FAMD Variables

Visualize variable contributions to FAMD dimensions. Shows quantitative
variables as arrows and qualitative variable categories as points.

## Usage

``` r
AETHER_plot_famd_variables(
  famd_result,
  dims = c(1, 2),
  var_type = "all",
  aggregate_quali = FALSE,
  color_by = "contrib",
  show_labels = TRUE,
  label_size = 3,
  arrow_size = 0.2,
  contrib_threshold = 0,
  title = "FAMD - Variables"
)
```

## Arguments

- famd_result:

  An artemis_famd object from ARTEMIS_famd()

- dims:

  Integer vector of length 2. Which dimensions to plot. Default = c(1,
  2).

- var_type:

  Character. Which variables to show: "all" (default), "quanti", or
  "quali".

- aggregate_quali:

  Logical. If TRUE, aggregate qualitative variables to variable level
  (one point per variable) instead of showing each category separately.
  Uses contribution-weighted centroid for positioning. Reduces clutter
  when many categories exist. Default = FALSE.

- color_by:

  Character. How to color variables: "type" (quanti vs quali), "contrib"
  (contribution strength). Default = "contrib".

- show_labels:

  Logical. Label variables. Default = TRUE.

- label_size:

  Numeric. Size of labels. Default = 3.

- arrow_size:

  Numeric. Size of arrow heads for quanti vars. Default = 0.2.

- contrib_threshold:

  Numeric. Only show variables with contribution above this threshold
  (as percentage). Default = 0 (show all).

- title:

  Character. Plot title. Default = "FAMD - Variables".

## Value

A ggplot object
