# Heatmap of Motif Enrichment Across Peak Sets

Creates a heatmap showing -log10(q-value) for motifs across multiple
peak sets. Useful for comparing which TF motifs are enriched in
different conditions or clusters.

## Usage

``` r
AETHER_plot_motif_comparison(
  homer_batch,
  top_n = 20,
  motifs = NULL,
  cluster_sets = FALSE,
  cluster_motifs = TRUE,
  q_thresh = 0.05,
  colors = NULL,
  title = "Motif Enrichment Comparison",
  show_values = FALSE,
  ...
)
```

## Arguments

- homer_batch:

  A homer_motif_batch object.

- top_n:

  Integer. Number of top motifs to include (default = 20). Selected by
  best q-value across all sets.

- motifs:

  Character vector. Specific motif family names to include. Overrides
  top_n if provided.

- cluster_sets:

  Logical. Hierarchically cluster columns (peak sets). Default FALSE
  (preserves input order).

- cluster_motifs:

  Logical. Hierarchically cluster rows (motifs). Default TRUE.

- q_thresh:

  Numeric. Significance threshold for motif selection (default = 0.05).
  Motifs must be significant in at least one set.

- colors:

  Character vector of length 3. Low, mid, high colors for the heatmap.
  Default: white to red.

- title:

  Character. Plot title (default = "Motif Enrichment Comparison").

- show_values:

  Logical. Show -log10(q) values in cells. Default FALSE.

- ...:

  Additional arguments passed to pheatmap.

## Value

A pheatmap object (invisibly).
