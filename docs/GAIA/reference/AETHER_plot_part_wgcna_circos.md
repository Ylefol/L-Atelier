# Plot PART-WGCNA comparison circos

Creates a circos plot showing the relationship between PART clustering
(DEA-based) and WGCNA modules. PART clusters occupy 270 degrees showing
per-gene L2FC values, WGCNA modules occupy 90 degrees showing trait
correlations. Chords connect clusters to modules based on shared genes.

## Usage

``` r
AETHER_plot_part_wgcna_circos(
  circos_data,
  output_file = "PART_WGCNA_circos.pdf",
  plot_width = 15,
  plot_height = 15,
  data_track_height = 0.05,
  label_track_height = 0.03,
  expression_limit = NULL,
  expression_colors = c("blue", "white", "red"),
  trait_colors = c("#276419", "white", "#8e0152"),
  na_color = "white",
  show_trait_values = FALSE,
  show_trait_significance = TRUE,
  show_sector_borders = TRUE,
  show_legend = TRUE,
  show_track_labels = FALSE,
  track_label_cex = 0.4,
  position_grouping_threshold = 0.25,
  min_chord_width_pct = 0.1,
  min_module_visual_size = 1000,
  part_degrees = 270,
  wgcna_degrees = 90,
  gap_degrees = 15,
  gap_after_track = NULL,
  gap_track_height = 0.02,
  verbose = TRUE
)
```

## Arguments

- circos_data:

  A `part_wgcna_circos` object from
  [`DEMETER_prepare_part_wgcna_circos()`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_prepare_part_wgcna_circos.md).

- output_file:

  Character. Output file path. Extension determines format (.pdf or
  .png). Default: "PART_WGCNA_circos.pdf".

- plot_width:

  Numeric. Plot width in inches. Default: 15.

- plot_height:

  Numeric. Plot height in inches. Default: 15.

- data_track_height:

  Numeric. Height of each data track (expression/trait). Default: 0.05.

- label_track_height:

  Numeric. Height of the innermost label track. Default: 0.03.

- expression_limit:

  Numeric or NULL. Symmetric cap for the L2FC color scale. If numeric
  (e.g., 3), the scale runs from -3 to 3 and values beyond are clamped
  to the extreme colors. If NULL (default), uses the full data range.

- expression_colors:

  Character vector of 2-3 colors for L2FC values (low, center, high).
  Default: c("blue", "white", "red").

- trait_colors:

  Character vector of 2-3 colors for trait correlations. Default:
  c("#276419", "white", "#8e0152").

- na_color:

  Character. Color for NA expression values. Default: "white".

- show_trait_values:

  Logical. Show correlation values as text in module sectors. Default:
  FALSE.

- show_trait_significance:

  Logical. Show significance markers in module sectors. Default: TRUE.

- show_sector_borders:

  Logical. Show borders on label track sectors. Default: TRUE.

- show_legend:

  Logical. Add a legend page to the output. Default: TRUE.

- show_track_labels:

  Logical. Add text labels in the gap identifying each track (sample
  names for PART, trait names for WGCNA). Experimental — may need
  adjustment depending on track count and label length. Default: FALSE.

- track_label_cex:

  Numeric. Text size for track labels. Default: 0.4.

- position_grouping_threshold:

  Numeric (0-1). Fraction of module sector within which same-cluster
  genes are grouped into a single chord. Default: 0.25.

- min_chord_width_pct:

  Numeric (0-1). Minimum chord width as fraction of module sector size.
  Default: 0.1.

- min_module_visual_size:

  Integer. Minimum visual size for small modules so they remain visible.
  Default: 1000.

- part_degrees:

  Numeric. Degrees allocated to PART section. Default: 270.

- wgcna_degrees:

  Numeric. Degrees allocated to WGCNA section. Default: 90.

- gap_degrees:

  Numeric. Gap between PART and WGCNA sections. Default: 15.

- gap_after_track:

  Integer or NULL. Insert an empty spacer track after this data track
  number (1-indexed). Affects all sectors equally. Useful for visually
  grouping tracks. Default: NULL (no gap).

- gap_track_height:

  Numeric. Height of the spacer track. Default: 0.02.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible NULL. Plot is saved to `output_file`.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_part_wgcna_circos(circos_data, output_file = "circos.pdf")

} # }
```
