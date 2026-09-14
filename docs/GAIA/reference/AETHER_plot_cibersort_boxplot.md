# Plot cell type proportions as boxplots by condition

Creates faceted boxplots showing cell type proportions grouped by
condition. Each facet represents a cell type with boxplots per group.

## Usage

``` r
AETHER_plot_cibersort_boxplot(
  result,
  group_by,
  cell_types = NULL,
  top_n = NULL,
  colors = NULL,
  title = "Cell Type Proportions by Group",
  text_size = 11,
  point_size = 2,
  show_points = TRUE,
  ncol = NULL
)
```

## Arguments

- result:

  An artemis_cibersort object or a proportions matrix (samples x cell
  types).

- group_by:

  Character vector of condition labels per sample. Required. Length must
  match number of samples. If the vector is **named** (names = group
  labels, values = color hex codes), the names are used as group labels
  and the values as group colors, overriding the `colors` parameter.

- cell_types:

  Character vector of cell types to include. Default: NULL (all non-zero
  cell types).

- top_n:

  Integer. Show only the top N cell types by mean proportion. Default:
  NULL (show all).

- colors:

  Named character vector of colors per group. If NULL, uses a default
  palette. Ignored if `group_by` is a named vector.

- title:

  Character. Plot title. Default: "Cell Type Proportions by Group".

- text_size:

  Numeric. Base font size. Default: 11.

- point_size:

  Numeric. Size of jittered points. Default: 2.

- show_points:

  Logical. Overlay jittered data points. Default: TRUE.

- ncol:

  Integer. Number of facet columns. Default: NULL (auto).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_cibersort_boxplot(cibersort_result, group_by = sample_conditions)

# Top 10 cell types only
AETHER_plot_cibersort_boxplot(cibersort_result, group_by = groups, top_n = 10)

# Named vector: names = groups, values = colors
group_colors <- c(CTRL = "#1f77b4", UVC = "#d62728")
groups_named <- group_colors[sample_sheet$condition]
AETHER_plot_cibersort_boxplot(cibersort_result, group_by = groups_named)

} # }
```
