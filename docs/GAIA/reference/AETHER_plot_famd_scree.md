# Plot FAMD Scree Plot

Bar plot showing variance explained by each FAMD dimension, with
cumulative variance line.

## Usage

``` r
AETHER_plot_famd_scree(
  famd_result,
  n_dims = 10,
  show_cumulative = TRUE,
  bar_fill = "#4575b4",
  line_color = "#d73027",
  title = "FAMD - Scree Plot"
)
```

## Arguments

- famd_result:

  An artemis_famd object from ARTEMIS_famd()

- n_dims:

  Integer. Number of dimensions to show. Default = 10.

- show_cumulative:

  Logical. Show cumulative variance line. Default = TRUE.

- bar_fill:

  Character. Fill color for bars. Default = "#4575b4".

- line_color:

  Character. Color for cumulative line. Default = "#d73027".

- title:

  Character. Plot title. Default = "FAMD - Scree Plot".

## Value

A ggplot object
