# Plot FAMD Individuals (Samples)

Scatter plot of individuals/samples in FAMD space, colored by a grouping
variable. Primary visualization for exploring sample clustering.

## Usage

``` r
AETHER_plot_famd_individuals(
  famd_result,
  dims = c(1, 2),
  color_by = NULL,
  point_size = 3,
  alpha = 0.8,
  show_ellipse = NULL,
  ellipse_level = 0.95,
  ellipse_alpha = 0.1,
  show_labels = FALSE,
  label_size = 3,
  colors = NULL,
  title = "FAMD - Individuals"
)
```

## Arguments

- famd_result:

  An artemis_famd object from ARTEMIS_famd()

- dims:

  Integer vector of length 2. Which dimensions to plot. Default = c(1,
  2).

- color_by:

  Character, vector, or factor. Can be:

  - Column name from the original data

  - Named vector (names must match FAMD row names - safest for external
    data)

  - Unnamed vector (must match FAMD sample order exactly)

  Should be categorical. If NULL, all points are same color.

- point_size:

  Numeric. Size of points. Default = 3.

- alpha:

  Numeric. Point transparency (0-1). Default = 0.8.

- show_ellipse:

  Logical. Draw confidence ellipses around groups. Default = TRUE if
  color_by is specified.

- ellipse_level:

  Numeric. Confidence level for ellipses. Default = 0.95.

- ellipse_alpha:

  Numeric. Ellipse fill transparency. Default = 0.1.

- show_labels:

  Logical. Label individual points. Default = FALSE.

- label_size:

  Numeric. Size of point labels. Default = 3.

- colors:

  Named character vector. Custom colors for groups. If NULL, uses
  default palette.

- title:

  Character. Plot title. Default = "FAMD - Individuals".

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ARTEMIS_famd(my_data)
AETHER_plot_famd_individuals(result, color_by = "treatment")
AETHER_plot_famd_individuals(result, color_by = "batch", dims = c(1, 3))

} # }
```
