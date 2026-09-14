# PC-Metadata Association Heatmap

Heatmap of -log10(p-values) from a PC-metadata association test. Rows
are principal components, columns are metadata variables. The colour
intensity indicates strength of association; a dashed line marks the
significance threshold.

## Usage

``` r
AETHER_plot_pc_association(
  assoc_result,
  p_threshold = 0.05,
  top_n = NULL,
  order_vars = TRUE,
  max_log10p = NULL,
  color_low = "white",
  color_high = "#b2182b",
  show_values = FALSE,
  font_size = 9,
  title = "PC-Metadata Association"
)
```

## Arguments

- assoc_result:

  An `artemis_pc_assoc` object from
  [`ARTEMIS_pc_metadata_association()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_pc_metadata_association.md).

- p_threshold:

  Numeric. Significance threshold for the cell border annotation.
  Default: 0.05.

- top_n:

  Integer or NULL. If set, retains only the `top_n` most significant
  variables (lowest p-value) per PC, then shows the union across all
  PCs. Useful when many variables are significant. NULL = show all
  tested variables. Default: NULL.

- order_vars:

  Logical. Order metadata variables (columns) by mean -log10(p)
  descending so the most associated variables appear first. Default:
  TRUE.

- max_log10p:

  Numeric or NULL. Cap the colour scale at this value. NULL = use the
  data maximum. Default: NULL.

- color_low:

  Character. Colour for low -log10(p) (no association). Default:
  "white".

- color_high:

  Character. Colour for high -log10(p) (strong association). Default:
  "#b2182b" (red).

- show_values:

  Logical. Overlay -log10(p) values as text in each cell. Default:
  FALSE.

- font_size:

  Numeric. Base font size. Default: 9.

- title:

  Character or NULL. Plot title. Default: "PC-Metadata Association".

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
assoc <- ARTEMIS_pc_metadata_association(ol$wide, ol$sample_meta,
                                          sample_col = "SampleID")
p <- AETHER_plot_pc_association(assoc)
print(p)
} # }
```
