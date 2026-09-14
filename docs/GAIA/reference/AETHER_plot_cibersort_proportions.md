# Plot cell type proportions as stacked bar chart

Creates a stacked bar plot showing the estimated cell type proportions
for each sample. Useful for comparing immune composition across samples.

## Usage

``` r
AETHER_plot_cibersort_proportions(
  result,
  group_by = NULL,
  colors = NULL,
  title = "Cell Type Proportions",
  text_size = 11,
  legend_position = "right",
  bar_width = 0.8,
  margins = c(5, 60, 5, 5)
)
```

## Arguments

- result:

  An artemis_cibersort object or a proportions matrix (samples x cell
  types).

- group_by:

  Optional factor or character vector of group labels per sample. Used
  to reorder samples by group and add faceting. Length must match number
  of samples. Default: NULL.

- colors:

  Named character vector of colors for cell types. Names should match
  cell type column names. If NULL, uses a default qualitative palette.

- title:

  Character. Plot title. Default: "Cell Type Proportions".

- text_size:

  Numeric. Base font size. Default: 11.

- legend_position:

  Character. Legend position. Default: "right".

- bar_width:

  Numeric. Width of bars. Default: 0.8.

- margins:

  Numeric vector of length 4. Plot margins in points (bottom, left, top,
  right). Default: c(80, 5, 5, 5). Increase bottom margin for long
  sample names.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic stacked bar
AETHER_plot_cibersort_proportions(cibersort_result)

# Grouped by condition
AETHER_plot_cibersort_proportions(cibersort_result, group_by = sample_groups)

} # }
```
