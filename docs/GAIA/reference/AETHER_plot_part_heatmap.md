# Plot PART Clustering Heatmap

Creates a heatmap of PART clustered data using ComplexHeatmap. Matches
the TiSA package visualization style with:

- Row annotation showing cluster color blocks

- Top annotation with group blocks and (optionally) timepoint colors

- Column splitting by group

- Row splitting by cluster

- Custom legends for Z-score, clusters, groups, and timepoints

## Usage

``` r
AETHER_plot_part_heatmap(
  part_result,
  sample_info,
  sample_col = NULL,
  group_col = "group",
  time_col = NULL,
  group_colors = NULL,
  time_colors = NULL,
  show_row_names = FALSE,
  show_column_labels = TRUE,
  row_names_side = "right",
  save_path = NULL,
  width = 10,
  height = 12
)
```

## Arguments

- part_result:

  An `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).

- sample_info:

  Data.frame with sample metadata. Required column: group. Optional:
  timepoint (for time series data).

- sample_col:

  Character or NULL. Column in sample_info whose values are sample
  identifiers matching `colnames(part_result$data)`. When provided, sets
  `rownames(sample_info)` to that column before matching — use this when
  rownames are integer indices rather than SampleIDs (e.g. after
  subsetting with a logical mask). Default: NULL (rownames used as-is).

- group_col:

  Character. Column name in sample_info for group. Default: "group".

- time_col:

  Character or NULL. Column name in sample_info for timepoint. Default:
  NULL. If NULL or column not found, timepoint annotation is skipped
  (useful for non-time-series data).

- group_colors:

  Named vector of colors for groups. Names should match group levels.
  Default: NULL (auto-generated).

- time_colors:

  Named vector of colors for timepoints. Default: NULL (uses
  yellow-orange-red gradient).

- show_row_names:

  Logical. Show gene names. Default: FALSE.

- show_column_labels:

  Logical. Whether to show column (sample) labels on the heatmap.
  Default: TRUE.

- row_names_side:

  Character. Side for row names: "left" or "right". Default: "right".

- save_path:

  Character. If provided, saves the heatmap to this path (supports .png,
  .pdf, .svg). Default: NULL (returns plot object).

- width:

  Numeric. Width in inches for saved plot. Default: 10.

- height:

  Numeric. Height in inches for saved plot. Default: 12.

## Value

If save_path is NULL, returns the ComplexHeatmap draw object. Otherwise
saves to file and returns invisibly.

## Details

This function is ported from TiSA's PART_heat_map() function and uses
ComplexHeatmap for high-quality visualization. The heatmap displays:

- Genes as rows, ordered and split by PART cluster

- Samples as columns, split by group and ordered by timepoint (if
  available)

- Z-score color scale (blue-white-red)

- Cluster color blocks on the left

- Group and timepoint annotations on top (timepoint optional)

For non-time-series data (standard DEA), set `time_col = NULL` or leave
it as default. The function will create a heatmap with only group
annotations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with time series data
AETHER_plot_part_heatmap(part_result, sample_info, time_col = "timepoint")

# Non-time-series data (no timepoint annotation)
AETHER_plot_part_heatmap(part_result, sample_info)

# With custom colors
AETHER_plot_part_heatmap(
  part_result, sample_info,
  group_colors = c("Control" = "blue", "Treatment" = "red")
)

# Save to file
AETHER_plot_part_heatmap(part_result, sample_info, save_path = "heatmap.png")

} # }
```
