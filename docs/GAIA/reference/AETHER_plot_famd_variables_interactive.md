# Interactive FAMD Variables Plot

Interactive plotly version of the FAMD variables plot. Useful for
exploring high-cardinality categorical variables where static plots
become unreadable.

## Usage

``` r
AETHER_plot_famd_variables_interactive(
  famd_result,
  dims = c(1, 2),
  var_type = "all",
  point_size = 8,
  show_arrows = TRUE,
  title = "FAMD - Variables (Interactive)",
  save_html = NULL
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

- point_size:

  Numeric. Size of points. Default = 8.

- show_arrows:

  Logical. Show arrows for quantitative variables. Default = TRUE.

- title:

  Character. Plot title. Default = "FAMD - Variables (Interactive)".

- save_html:

  Character. If provided, saves plot as standalone HTML file at this
  path. Default = NULL (no save).

## Value

A plotly object

## Details

This function creates an interactive scatter plot where you can:

- Hover over points to see variable/category names and contributions

- Zoom into regions of interest

- Pan around the plot

- Toggle variable types on/off via legend

The HTML output is self-contained and can be shared with others who have
a modern web browser (no R required to view).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ARTEMIS_famd(my_data)
AETHER_plot_famd_variables_interactive(result)

# Save as HTML for sharing
AETHER_plot_famd_variables_interactive(result, save_html = "famd_explore.html")

} # }
```
