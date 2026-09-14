# Dotplot of HOMER Motif Enrichment Results

Creates a dotplot showing enriched transcription factor motifs. For
single peak sets, shows one column of dots. For batch results (multiple
peak sets), shows motifs on the y-axis and peak sets on the x-axis.

## Usage

``` r
AETHER_plot_motif_enrichment(
  homer_result,
  top_n = 15,
  color_by = "q_value",
  size_by = "pct_target",
  sets = NULL,
  q_thresh = 0.05,
  title = NULL,
  font_size = 9,
  max_name_length = 40,
  color_low = "#2166ac",
  color_high = "#b2182b",
  dot_range = c(2, 8)
)
```

## Arguments

- homer_result:

  A homer_motif or homer_motif_batch object.

- top_n:

  Integer. Number of top motifs to show per peak set (default = 15).

- color_by:

  Character. Variable for color scale: "q_value" (default, shown as
  -log10) or "p_value".

- size_by:

  Character. Variable for dot size: "pct_target" (default) or
  "n_target".

- sets:

  Character vector. For batch results, specific peak sets to include. If
  NULL, all sets are shown.

- q_thresh:

  Numeric. Only show motifs significant in at least one set (default =
  0.05). Set to 1 to show all.

- title:

  Character. Plot title. If NULL, auto-generated.

- font_size:

  Numeric. Base font size (default = 9).

- max_name_length:

  Integer. Maximum characters for motif names (default = 40).

- color_low:

  Character. Color for low significance (default = "#2166ac" blue).

- color_high:

  Character. Color for high significance (default = "#b2182b" red).

- dot_range:

  Numeric vector of length 2. Min and max dot sizes (default = c(2, 8)).

## Value

A ggplot object.
