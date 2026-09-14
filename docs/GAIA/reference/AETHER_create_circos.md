# Create Circos Plot

Creates a circos plot from prepared data, with signal tracks showing
genomic coverage across chromosomes. Supports multiple samples/groups
with customizable colors and track organization.

## Usage

``` r
AETHER_create_circos(
  circos_data,
  genome = "hg38",
  output = NULL,
  split_by = "omics",
  color_by = "track",
  invert_gradient = FALSE,
  colors = NULL,
  track_height = 0.12,
  track_margin = 0.01,
  line_width = 1,
  ideogram_height = 0.04,
  title = NULL,
  width = 10,
  height = 10,
  res = 300,
  highlights = NULL,
  highlight_padding = 0.005,
  annotation_regions = NULL,
  annotation_colors = NULL,
  show_chords = TRUE,
  chord_alpha = 0.5,
  annotation_sector_size = 2e+08,
  start_degree = 90,
  verbose = TRUE
)
```

## Arguments

- circos_data:

  List from AETHER_prepare_circos_data() containing bins, signal matrix,
  sample_info, and chr_sizes.

- genome:

  Character. Genome assembly name for cytoband data (default = "hg38").

- output:

  Character string. Output file path. Format is inferred from extension:
  ".png" for PNG, ".pdf" for PDF. If NULL, plot is displayed but not
  saved.

- split_by:

  Character vector specifying how to split tracks. Options:

  - "omics" - one track per omics type (always applied)

  - "group" - additionally split by group

  - "batch" - additionally split by batch

  Can combine: c("omics", "group") creates one track per omics+group.
  Default is "omics" only.

- color_by:

  Character string specifying color granularity. Options:

  - "track" (default) - one base color per track, gradients for samples

  - "group" - colors keyed by group name (e.g., "WT", "KO")

  - "omics_group" - colors keyed by omics_group (e.g., "ATACseq_WT")

  Allows finer color control than track splitting. For example,
  split_by="omics" with color_by="omics_group" puts WT and KO in same
  track but different colors. When multiple samples share the same color
  key, gradients are generated.

- invert_gradient:

  Logical. Controls gradient direction relative to signal:

  - FALSE (default) - highest signal = darkest color (plotted on top)

  - TRUE - lowest signal = darkest color (plotted on top)

- colors:

  Named character vector of colors. Names should match the color_by
  scheme:

  - color_by="track": names like "ATACseq", "CHIPseq"

  - color_by="group": names like "WT", "KO"

  - color_by="omics_group": names like "ATACseq_WT", "ATACseq_KO"

  If NULL, assigns colors automatically.

- track_height:

  Numeric. Height of each signal track as fraction of total radius
  (default = 0.12).

- track_margin:

  Numeric. Margin between tracks (default = 0.01).

- line_width:

  Numeric. Width of signal lines (default = 1.0).

- ideogram_height:

  Numeric. Height of chromosome ideogram (default = 0.04).

- title:

  Character string. Plot title displayed in center (default = NULL).

- width:

  Numeric. Output width in inches (default = 10).

- height:

  Numeric. Output height in inches (default = 10).

- res:

  Numeric. Resolution for PNG output in DPI (default = 300).

- highlights:

  Named list of highlight definitions. Each element should be a list
  with:

  - regions: data.frame with chr, start, end columns

  - tracks: character - track name(s) to highlight, or "all" for all
    signal tracks

  - color: hex color string (should include alpha for transparency,
    e.g., "#FF000033")

  Example: list("peaks" = list(regions = peaks_df, tracks = "ATACseq",
  color = "#FF000033"))

- highlight_padding:

  Numeric. Padding for highlights as fraction of chromosome length
  (default = 0.005, i.e., 0.5% on each side). Increase to make small
  regions more visible.

- annotation_regions:

  Data frame with annotated regions for chord connections. Must contain
  columns: chr, start, end, annotation. The annotation column should
  contain category labels (e.g., "Promoter", "Intron", "Exon").
  Compatible with APOLLO_annotate_peaks() output. User can simplify
  annotation labels before passing.

- annotation_colors:

  Named character vector of colors for annotation categories. Names
  should match values in annotation_regions\$annotation. If NULL,
  default colors are used.

- show_chords:

  Logical. If TRUE (default), draw chord connections from regions to
  annotation sector. If FALSE, only show annotation sector without
  chords.

- chord_alpha:

  Numeric. Transparency for chord fill (0-1, default = 0.5).

- annotation_sector_size:

  Numeric. Size of the annotation pseudo-chromosome in base pairs
  (default = 2e8). Adjust for visual balance.

- start_degree:

  Numeric. Starting angle in degrees for the first sector (default = 90,
  which is top/12 o'clock). Use smaller values to rotate clockwise. For
  annotation sector at center-right (3 o'clock), try values around -50
  to -70 depending on the number of chromosomes.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

Invisibly returns NULL. Side effect is creating the plot.

## Details

The function creates a circos plot with:

- Outer ring: Chromosome ideogram with labels and alternating colors

- Inner rings: Signal tracks, each independently scaled

Track organization examples:

- split_by = "omics": ATACseq track, CHIPseq track

- split_by = c("omics", "group"): ATACseq_WT, ATACseq_KO, CHIPseq_WT,
  CHIPseq_KO

- split_by = c("omics", "batch"): ATACseq_b1, ATACseq_b2, etc.

Each track is independently scaled (y-axis) to prevent squishing when
combining different omics types with different signal ranges.

Color assignment depends on the color_by parameter:

- color_by="track": One base color per track, with gradient shades for
  multiple samples (darkest = highest signal, plotted on top)

- color_by="group": Colors assigned by group name, allowing different
  colors within the same track

- color_by="omics_group": Colors assigned by omics+group combination,
  e.g., "ATACseq_WT" and "ATACseq_KO" can have different colors even in
  the same track when split_by="omics"

Samples are always plotted in order of mean signal (lowest first,
highest last on top) to preserve data visibility.

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare data
circos_data <- AETHER_prepare_circos_data(
  sample_sheet, chr_sizes, omics = c("ATACseq", "CHIPseq")
)

# One track per omics (default)
AETHER_create_circos(circos_data, output = "circos.png")

# Split by omics and group
AETHER_create_circos(circos_data, output = "circos.png",
                      split_by = c("omics", "group"))

# With centered title
AETHER_create_circos(circos_data, output = "circos.png",
                      title = "Multi-omics Overview")

# With region highlights
my_highlights <- list(
  "ATAC_peaks" = list(
    regions = atac_peaks,  # data.frame with chr, start, end
    tracks = "ATACseq",
    color = "#FF000033"    # red with transparency
  ),
  "shared_peaks" = list(
    regions = shared_peaks,
    tracks = "all",        # highlight on all signal tracks
    color = "#0000FF33"    # blue with transparency
  )
)
AETHER_create_circos(circos_data, output = "circos.png",
                      highlights = my_highlights)

# With annotation chords (Phase 2B)
# annotation_regions should have chr, start, end, annotation columns
# (compatible with APOLLO_annotate_peaks() output)
annotated_peaks <- APOLLO_annotate_peaks(my_peaks, txdb)

# Optional: simplify annotation labels
annotated_peaks$annotation <- gsub("Distal Intergenic", "Intergenic",
                                    annotated_peaks$annotation)

AETHER_create_circos(circos_data, output = "circos_with_chords.png",
                      annotation_regions = annotated_peaks,
                      show_chords = TRUE,
                      chord_alpha = 0.4)

# Custom annotation colors
my_anno_colors <- c(
  "Promoter" = "#E78AC3",
  "Intron" = "#66C2A5",
  "Exon" = "#8DA0CB",
  "Intergenic" = "#FC8D62"
)
AETHER_create_circos(circos_data, output = "circos_custom.png",
                      annotation_regions = annotated_peaks,
                      annotation_colors = my_anno_colors)
} # }
```
