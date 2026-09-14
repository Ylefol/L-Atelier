# Atlas visualisation: metadata tracks + colour ring + dendrogram + scatter

Draws a circos-style figure with a central UMAP/tSNE scatter surrounded
by up to three concentric track groups (outside → inside):

1.  **Metadata tracks** (optional) — one ring per `meta_col` entry.
    Entries can be discrete `colData` columns, numeric `colData`
    columns, or gene names (expression from `meta_assay`). Track style
    is controlled by `meta_style`.

2.  **Colour ring** — one sector per category in `colour_by`, width
    proportional to cell/bin count.

3.  **Dendrogram track** (optional) — per-sector hierarchical clustering
    of cells or bin centroids.

## Usage

``` r
ASPIS_plot_atlas(
  sce,
  colour_by,
  dimred = c("UMAP", "tSNE"),
  style = c("points", "contour", "both"),
  point_size = 0.5,
  point_alpha = 0.5,
  ring_palette = NULL,
  label_clusters = FALSE,
  label_size = 3,
  label_split = FALSE,
  label_split_min_sil = 0.35,
  label_split_k_max = 5L,
  label_split_size_scale = TRUE,
  label_split_size_range = NULL,
  label_split_min_cells = 50L,
  label_split_nudge = 0,
  label_split_angle = TRUE,
  n_levels = 8L,
  contour_alpha = 0.7,
  contour_adjust = 1,
  boundary_pad = 0.1,
  inner_pad = 0.25,
  gap_degrees = 2,
  inner_scale = 0.95,
  show_label_track = FALSE,
  label_track_height = 0.12,
  label_track_cex = 0.55,
  label_track_col = "black",
  label_track_font = 1L,
  show_ring = TRUE,
  ring_height = 0.15,
  ring_labels = TRUE,
  ring_label_cex = 0.65,
  ring_label_col = "white",
  ring_gap = 0.01,
  ring_border_col = NA,
  min_seg_frac = NULL,
  bottom_gap_degrees = 15,
  show_meta_labels = FALSE,
  meta_label_names = NULL,
  meta_label_cex = 0.65,
  meta_label_col = "black",
  meta_label_lwd = 0.5,
  show_dendrogram = TRUE,
  dend_height = 0.15,
  dend_clip = 1,
  dend_power = 1,
  dend_n_pcs = 20L,
  dend_max_cells = 500L,
  dend_col = "grey30",
  dend_lwd = 0.5,
  n_bin = NULL,
  meta_col = NULL,
  meta_assay = "logcounts",
  meta_palette = NULL,
  meta_height = 0.04,
  meta_style = c("dot", "dot_border", "heatmap"),
  meta_cex = 0.3,
  meta_border_col = "white",
  meta_border_lwd = 0.3,
  meta_scale_q = c(0, 1),
  show_satellites = TRUE,
  sat_col = NULL,
  sat_custom = list(),
  sat_size = 0.22,
  sat_point_size = 0.3,
  sat_point_alpha = 0.5,
  show_sat_title = TRUE,
  sat_title_cex = 0.7,
  sat_title_col = "black",
  sat_legend = FALSE,
  sat_legend_max_cats = 10L,
  sat_legend_title_cex = 0.6,
  sat_legend_text_cex = 0.55,
  sat_legend_pad = 0.01,
  sat_legend_size = 1,
  output_file = NULL,
  width = 8,
  height = 8,
  res = 150
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a UMAP or tSNE embedding and, ideally, a
  PCA.

- colour_by:

  Character. Discrete `colData` column for the scatter colouring and
  ring sectors.

- dimred:

  Character. `"UMAP"` (default) or `"tSNE"`.

- style:

  Character. `"points"` (default), `"contour"`, or `"both"`.

- point_size:

  Numeric. Point size in the central scatter. Default `0.5`.

- point_alpha:

  Numeric. Point transparency (0–1). Default `0.5`.

- ring_palette:

  Character vector or `NULL`. Colours for `colour_by` categories —
  applied consistently to both the central scatter points and the colour
  ring segments so they always match. `NULL` uses the built-in 20-colour
  palette.

- label_clusters:

  Logical. Overlay centroid labels on the scatter. One label per
  `colour_by` category, placed at the overall centroid. Default `FALSE`.

- label_size:

  Numeric. Base centroid label size (ggplot size units). Default `3`.

- label_split:

  Logical. Enable density-aware sub-cluster labelling. When `TRUE`, each
  `colour_by` category is analysed independently in 2D embedding space:
  k-means splits (k = 2 to `label_split_k_max`) are tried and the best
  split is accepted if its mean silhouette score exceeds
  `label_split_min_sil`. Accepted splits place one label per
  sub-cluster, all showing the same category name. Overrides
  `label_clusters` when active. Default `FALSE`.

- label_split_min_sil:

  Numeric (0–1). Minimum mean silhouette score required to accept a k \>
  1 split. Higher values are stricter (require more clearly separated
  sub-clouds); lower values allow splits of weakly separated
  distributions. A value around `0.25` is permissive; `0.5` is strict.
  Default `0.35`.

- label_split_k_max:

  Integer. Maximum number of sub-clusters to try per category. Default
  `5L`.

- label_split_size_scale:

  Logical. Scale label text size proportionally to the square-root of
  the sub-cluster's fraction of the category's cells. The largest
  sub-cluster gets the upper bound of `label_split_size_range`; smaller
  sub-clusters are mapped linearly within the range. Default `TRUE`.

- label_split_size_range:

  Numeric vector of length 2: `c(min, max)` in ggplot text-size units.
  Controls the absolute size bounds when
  `label_split_size_scale = TRUE`. When `FALSE` all labels get the
  maximum value. `NULL` (default) derives the range automatically from
  `label_size`: `c(label_size * 0.5, label_size)`.

- label_split_min_cells:

  Integer. Minimum number of cells a category (or an individual
  sub-cluster after splitting) must contain for a label to be placed.
  Categories below this threshold are skipped entirely; sub-clusters
  below it are dropped and the remaining ones retain their labels. If
  all sub-clusters of a category are dropped the category falls back to
  a single centroid label (provided the category itself meets the
  threshold). Default `50L`.

- label_split_nudge:

  Numeric. Shift each label away from its centroid along the axis
  perpendicular to the text (PC2), in units of that sub-cluster's PC2
  standard deviation. Positive values push the label "upward" on screen
  (PC2 is always oriented toward positive y); negative values push it
  "downward". Useful for moving labels off the dense core of a cluster
  while keeping the angle unchanged. Default `0` (no shift).

- label_split_angle:

  Logical. Tilt each label along the principal axis (PC1) of the
  sub-cluster's point cloud. Text is constrained to the −90°–90° range
  so it remains readable. Default `TRUE`.

- n_levels:

  Integer or `"dynamic"`. Contour levels. Default `8`.

- contour_alpha:

  Numeric. Contour opacity. Default `0.7`.

- contour_adjust:

  Numeric. KDE bandwidth multiplier for contour silhouettes. Values
  above `1` produce smoother, more spread-out cluster outlines; values
  below `1` tighten the silhouettes to the data. Only applies when
  `style` is `"contour"` or `"both"`. Default `1` (automatic bandwidth).

- boundary_pad:

  Numeric. Fractional axis expansion for KDE. Default `0.10`.

- inner_pad:

  Numeric. Additional scatter viewport expansion beyond the data range.
  Default `0.25`.

- gap_degrees:

  Numeric. White-space degrees between sectors. Default `2`.

- inner_scale:

  Numeric (0–1). Fraction of the innermost ring boundary used for the
  scatter viewport. Default `0.95`.

- show_label_track:

  Logical. Draw an outermost text label track with one wrapped label per
  sector. Labels are the `colour_by` category names. Text is wrapped to
  avoid spilling beyond sector boundaries. Default `FALSE`.

- label_track_height:

  Numeric (0–1). Radial thickness of the label track. Increase to
  accommodate more lines of wrapped text. Default `0.12`.

- label_track_cex:

  Numeric. Label text size. Default `0.55`.

- label_track_col:

  Character. Label text colour. Default `"black"`.

- label_track_font:

  Integer. Font face: `1` plain (default), `2` bold, `3` italic.

- show_ring:

  Logical. Draw the colour ring track. Set to `FALSE` to hide it (useful
  when only the metadata tracks and/or dendrogram are needed). Default
  `TRUE`.

- ring_height:

  Numeric (0–1). Radial thickness of the colour ring. Default `0.15`.

- ring_labels:

  Logical. Draw category labels inside the colour ring. Default `TRUE`.

- ring_label_cex:

  Numeric. Ring label size. Default `0.65`.

- ring_label_col:

  Character. Ring label text colour. Default `"white"`.

- ring_gap:

  Numeric (0–1). Radial thickness of the empty spacer between the
  metadata tracks and the colour ring. Default `0.01`.

- ring_border_col:

  Character or `NA`. Border colour for the colour ring rectangles. When
  `NA` (default) and `meta_style = "heatmap"`, automatically inherits
  `meta_border_col`.

- min_seg_frac:

  Numeric (0–1) or `NULL`. Minimum arc fraction for any sector. When
  set, sectors smaller than this floor are expanded and larger sectors
  are proportionally compressed so the total arc remains 360°. The
  internal x-axis of each sector is unaffected, so dendrograms and
  metadata dots fill their sector correctly regardless. `NULL` (default)
  keeps sectors strictly proportional to cell count.

- bottom_gap_degrees:

  Numeric. Angular width (degrees) of the fixed gap placed at the visual
  bottom (6 o'clock) of the circos plot. This gap is always present and
  slightly wider than `gap_degrees` so that metadata track labels can be
  drawn without overlap. Default `15`. Note: fixing this gap at the
  bottom overrides the embedding-alignment rotation (sectors retain
  their relative order but not absolute orientation).

- show_meta_labels:

  Logical. Draw a short label for each metadata track in the bottom gap,
  at the radial level of that track. A short tick connects each label to
  its track. Default `FALSE`.

- meta_label_names:

  Named character vector or `NULL`. Override display names for metadata
  tracks. Names must match entries in `meta_col`; values are the text
  shown. Useful for shortening long column or gene names. `NULL` uses
  the original `meta_col` entries verbatim.

- meta_label_cex:

  Numeric. Text size for metadata track labels. Default `0.65`.

- meta_label_col:

  Character. Colour for metadata track label text and connector ticks.
  Default `"black"`.

- meta_label_lwd:

  Numeric. Line width of the connector ticks. Default `0.5`.

- show_dendrogram:

  Logical. Draw the per-sector dendrogram track. Default `TRUE`.

- dend_height:

  Numeric (0–1). Radial thickness of the dendrogram track. Default
  `0.15`.

- dend_clip:

  Numeric (0–1\]. Fraction of the global maximum dendrogram height to
  display. Values below 1 zoom into the lower portion of the height
  range: the track shows only heights 0–`dend_clip × max_h`, stretching
  that region to fill the full track while clipping the root at the
  boundary. Useful for inspecting fine leaf structure. Default `1` (no
  clipping).

- dend_power:

  Numeric (\> 0). Power exponent for a non-linear height transformation
  applied before drawing. Compresses the visual prominence of the root
  arch without clipping it. The transformation is
  `h_new = (h / max_h)^power × max_h`, which is concave for `power < 1`:
  small heights (leaf forks) expand radially while the gap between the
  last cluster-level fork and the root shrinks. A value of `0.5`
  (square-root) gives moderate root compression; `0.3` gives stronger
  compression. Default `1` (no transformation).

- dend_n_pcs:

  Integer. PCA dimensions used for cell/bin distances. Default `20`.

- dend_max_cells:

  Integer or `NULL`. Maximum cells per sector used for dendrogram
  construction when **not** binning (`n_bin = NULL`). Sectors exceeding
  this count are randomly subsampled before
  [`hclust()`](https://rdrr.io/r/stats/hclust.html) to keep memory and
  runtime manageable (full \\n \times n\\ distance matrix). Has no
  effect when `n_bin` is set, because in that case clustering operates
  on bin centroids rather than individual cells. `NULL` uses all cells.
  Default `500`.

- dend_col:

  Character. Dendrogram line colour. Default `"grey30"`.

- dend_lwd:

  Numeric. Dendrogram line width. Default `0.5`.

- n_bin:

  Integer or `NULL`. Cells per bin. `NULL` or `1` disables binning (one
  leaf per cell). When set, cells are sorted by PC1 and grouped into
  consecutive bins; dendrograms and metadata tracks operate on bin
  centroids / aggregated values. Default `NULL`.

- meta_col:

  Character vector or `NULL`. One or more track specifications drawn
  from outermost to innermost. Each entry can be:

  - A **discrete** `colData` column (character/factor) — per-bin
    majority vote, categorical palette.

  - A **numeric** `colData` column (integer/double) — per-bin mean,
    colour gradient.

  - A **gene name** present in `rownames(sce)` — per-bin mean expression
    from `meta_assay`, colour gradient.

  `NULL` (default) skips all metadata tracks.

- meta_assay:

  Character. Assay used for gene expression tracks. Default
  `"logcounts"`.

- meta_palette:

  Named list or `NULL`. Colour specification per track, keyed by
  column/gene name. For discrete tracks: a character vector. For
  continuous/gene tracks: a two-element vector `c(low, high)` for
  `colorRampPalette`. Missing entries use defaults (built-in categorical
  palette or `c("grey90", "#2171B5")`).

- meta_height:

  Numeric (0–1). Radial thickness of each metadata track. Default
  `0.04`.

- meta_style:

  Character. Visual style for all metadata tracks. `"dot"` (default) —
  filled circles; `"dot_border"` — filled circles with an outline;
  `"heatmap"` — filled colour tiles spanning the full bin width and
  track height.

- meta_cex:

  Numeric. Point size for `"dot"` and `"dot_border"` styles. Default
  `0.3`.

- meta_border_col:

  Character, `NA`, or `NULL`. Border colour applied to metadata track
  elements: outline of dots (`"dot_border"` style) or tile borders
  (`"heatmap"` style). `NULL` or `NA` draws no border. Also inherited by
  `ring_border_col` when that is `NA` and `meta_style = "heatmap"`.
  Default `"white"`.

- meta_border_lwd:

  Numeric. Border line width for `"dot_border"` style. Default `0.3`.

- meta_scale_q:

  Numeric vector of length 2. Quantiles used to define the low and high
  anchors of the colour gradient for continuous tracks (numeric
  `colData` columns or genes). Values outside this range are clamped to
  the nearest anchor colour. Default `c(0, 1)` (full range). For sparse
  genes expressed in only a subset of cells, try `c(0, 0.95)` or
  `c(0, 0.99)` to prevent a few zero-expressing cells from washing out
  the colour scale.

- show_satellites:

  Logical. Draw satellite mini-embedding panels in the device corners.
  Set to `FALSE` to suppress all satellites. Default `TRUE`.

- sat_col:

  Character vector or `NULL`. Which `meta_col` entries to display as
  satellite mini-embeddings in the four device corners (top-left,
  top-right, bottom-left, bottom-right). Must be a subset of `meta_col`;
  at most four entries are used. `NULL` (default) uses the first four
  entries of `meta_col`. If `meta_col` is `NULL` no satellites are drawn
  regardless of this parameter.

- sat_custom:

  Named list or `NULL`. Escape hatch for satellite panels that aren't a
  simple metadata/gene scatter – e.g. the output of
  [`ASPIS_plot_velocity_stream`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_velocity_stream.md).
  Names must be one or more of `"TL"`, `"TR"`, `"BL"`, `"BR"`; each
  value must be a `ggplot` object. A named corner is drawn exactly as
  supplied (structural chrome – panel border, plot margin, background –
  is stripped to match the other satellites, but the plot's own title
  and legend are left untouched, since there is no metadata track to
  derive a title from or an automatic legend for) instead of a
  `sat_col`-driven scatter for that corner. Any corner not named here
  still falls back to the usual `sat_col`/`meta_col` behaviour. Default
  [`list()`](https://rdrr.io/r/base/list.html) (no custom corners).

- sat_size:

  Numeric (0–1). Fraction of the device canvas occupied by each
  satellite panel. Default `0.22`.

- sat_point_size:

  Numeric. Point size in satellite scatters. Default `0.3`.

- sat_point_alpha:

  Numeric. Point transparency in satellite scatters. Default `0.5`.

- sat_title_cex:

  Numeric. Title text size for satellite panels, expressed as a
  multiplier on the base ggplot size (11 pt). Default `0.7` (≈ 7.7 pt).

- sat_title_col:

  Character. Satellite title colour. Default `"black"`.

- sat_legend:

  Logical or character vector. Controls which satellite panels display a
  legend. `FALSE` (default) shows no legends. `TRUE` shows a legend for
  every satellite. A character vector of `meta_col` names shows legends
  only for the named satellites. Discrete legends are automatically
  suppressed when the number of categories exceeds
  `sat_legend_max_cats`.

- sat_legend_max_cats:

  Integer. Maximum number of categories a discrete satellite track may
  have before its legend is automatically hidden. Continuous/gradient
  tracks are always shown regardless of this limit. Default `10L`.

- sat_legend_title_cex:

  Numeric. Legend title size multiplier (×11 pt). Default `0.6`.

- sat_legend_text_cex:

  Numeric. Legend item text size multiplier (×11 pt). Default `0.55`.

- sat_legend_pad:

  Numeric (0–1). Padding between each legend and its nearest device
  edge, in npc units. TL and TR legends are shifted inward from the left
  and right edges respectively; BL and BR legends are shifted upward
  from the bottom edge. Default `0.01`.

- sat_legend_size:

  Numeric. Size multiplier applied uniformly to all legend geometry:
  colourbar bar dimensions and discrete key size. Values above 1
  enlarge; below 1 shrink. Text sizes are unaffected (use
  `sat_legend_title_cex` / `sat_legend_text_cex` for those). Default
  `1`.

- output_file:

  Character or `NULL`. Path for automatic PNG device management. `NULL`
  (default) uses the current device.

- width:

  Numeric. Output width in inches. Default `8`.

- height:

  Numeric. Output height in inches. Default `8`.

- res:

  Integer. Output resolution in dpi. Default `150`.

## Value

`invisible(NULL)`. Output is drawn to the current device.
