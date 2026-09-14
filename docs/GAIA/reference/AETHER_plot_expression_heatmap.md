# Heatmap of Normalized Expression for a Feature Subset

Draws a pheatmap of normalized counts for an arbitrary set of features
(genes, TEs, or any other row ID present in the count matrix), features
as rows and samples as columns. Unlike
[`AETHER_plot_pca()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_pca.md)
(which summarizes the whole matrix), this is meant for zooming into a
small, caller-chosen subset – e.g. all members of one TE family, or one
gene module – to see whether they move together across samples. Works
the same way regardless of subset size: rows just compress for larger
sets, unlike a per-feature facet grid which stops being readable past a
few dozen features (see
[`AETHER_plot_expression_boxplot()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_expression_boxplot.md)
for that case).

## Usage

``` r
AETHER_plot_expression_heatmap(
  counts,
  features,
  sample_info = NULL,
  sample_col = NULL,
  group_col = "group",
  log_transform = TRUE,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color_palette = NULL,
  title = "Expression Heatmap",
  show_rownames = TRUE,
  show_colnames = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- counts:

  An `artemis_norm` object from
  [`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md),
  or a normalized count matrix (features x samples). If an
  `artemis_norm` object, `sample_info` and `group_col` default to its
  `$targets` and `$parameters$group_col` (but can be overridden).

- features:

  Character vector of feature IDs to include (must be a subset of
  `rownames(counts)`). IDs not found are dropped (reported when
  `verbose = TRUE`).

- sample_info:

  Data.frame of sample metadata, rownames matching `colnames(counts)`.
  Required when `counts` is a plain matrix and `group_col` is not NULL.

- sample_col:

  Character or NULL. Column in `sample_info` holding sample IDs, used as
  rownames for matching if provided. Default: NULL.

- group_col:

  Character or NULL. Column in `sample_info` to show as a column
  annotation strip. Set to NULL to omit annotation. Default: "group".

- log_transform:

  Logical. Log2(x + 1)-transform counts before scaling. Default: TRUE.

- scale:

  Character. Passed to
  [`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html):
  "row", "column", or "none". Default: "row" (z-score each feature
  across samples, the usual choice for comparing features on different
  scales).

- cluster_rows, cluster_cols:

  Logical. Cluster features / samples. Default: TRUE / FALSE (samples
  usually kept in a meaningful order, e.g. grouped by condition, rather
  than re-clustered).

- color_palette:

  Character vector of colors for the heatmap gradient. Default: NULL
  (diverging blue-white-red, appropriate for z-scored data).

- title:

  Character. Plot title. Default: "Expression Heatmap".

- show_rownames, show_colnames:

  Logical. Default: TRUE.

- verbose:

  Logical. Print notes about dropped/missing features and samples.
  Default: TRUE.

- ...:

  Additional arguments passed to
  [`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html).

## Value

A pheatmap object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
# All members of one TE family, from an artemis_norm object
AETHER_plot_expression_heatmap(norm_data, features = l2_family_ids)

# From a plain matrix + sample sheet
AETHER_plot_expression_heatmap(norm_counts, features = my_genes,
                                sample_info = targets, group_col = "group")
} # }
```
