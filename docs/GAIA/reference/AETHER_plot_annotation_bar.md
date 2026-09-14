# Genomic Annotation Distribution Bar Plot

Creates a stacked horizontal bar plot showing the distribution of
genomic features (Promoter, Exon, Intron, Intergenic, etc.) for one or
more sets of annotated peaks. Useful for comparing annotation profiles
between different datasets (e.g., ATAC vs ChIP peaks).

## Usage

``` r
AETHER_plot_annotation_bar(
  annotated_list,
  use_simple = TRUE,
  show_percentage = TRUE,
  colors = NULL,
  title = "Genomic Feature Distribution",
  bar_height = 0.7,
  show_labels = TRUE,
  min_label_pct = 5
)
```

## Arguments

- annotated_list:

  A named list of data.frames from APOLLO_annotate_peaks(), or a single
  data.frame. Each data.frame must have an 'annotation_simple' or
  'annotation' column.

- use_simple:

  Logical. If TRUE (default), uses 'annotation_simple' column. If FALSE,
  uses full 'annotation' column.

- show_percentage:

  Logical. If TRUE (default), shows percentages. If FALSE, shows raw
  counts.

- colors:

  Named character vector of colors for each annotation category. If
  NULL, uses a default color scheme.

- title:

  Character. Plot title (default = "Genomic Feature Distribution").

- bar_height:

  Numeric. Height of each bar (default = 0.7).

- show_labels:

  Logical. If TRUE (default), shows percentage labels on bars.

- min_label_pct:

  Numeric. Minimum percentage to show label (default = 5). Prevents
  cluttering small segments with labels.

## Value

A ggplot object.

## Details

The plot shows one horizontal stacked bar per dataset in annotated_list.
Each segment's width is proportional to the count/percentage of peaks in
that genomic category.

Default annotation categories (from annotation_simple):

- Promoter: regions within TSS window

- 5' UTR, 3' UTR: untranslated regions

- Exon: exonic regions

- Intron: intronic regions

- Downstream: downstream of genes

- Intergenic: between genes

## Examples

``` r
if (FALSE) { # \dontrun{
# Single dataset
p <- AETHER_plot_annotation_bar(atac_annotated)

# Compare ATAC vs ChIP annotation profiles
annotated_list <- list(
  "ATAC peaks" = atac_annotated,
  "ChIP peaks" = chip_annotated
)
p <- AETHER_plot_annotation_bar(annotated_list)

# Custom colors
my_colors <- c(Promoter = "#e41a1c", Exon = "#377eb8", Intron = "#4daf4a")
p <- AETHER_plot_annotation_bar(annotated_list, colors = my_colors)

} # }
```
