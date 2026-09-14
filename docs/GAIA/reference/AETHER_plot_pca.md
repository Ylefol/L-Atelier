# PCA plot of normalized counts

Performs PCA on normalized count data and creates a scatter plot colored
by group with optional shape mapping for a second variable (e.g.,
timepoint).

## Usage

``` r
AETHER_plot_pca(
  counts,
  sample_info = NULL,
  sample_col = NULL,
  group_col = "group",
  var_type = "categorical",
  shape_col = NULL,
  dims = c(1, 2),
  ntop = 500,
  log_transform = TRUE,
  colors = NULL,
  point_size = 3,
  label_samples = FALSE,
  title = "PCA",
  verbose = TRUE
)
```

## Arguments

- counts:

  Normalized count matrix (genes x samples), an artemis_norm object, or
  an artemis_ts_norm object. Genes as rows, samples as columns.

- sample_info:

  Data.frame with sample metadata. Rownames or a column must match the
  column names of the count matrix. If counts is an
  artemis_norm/artemis_ts_norm object, this is extracted automatically
  (but can be overridden).

- sample_col:

  Character or NULL. Column in `sample_info` holding sample IDs matching
  `colnames(counts)`. When provided, this column is used as rownames for
  matching instead of the existing rownames. Required when `sample_info`
  has sequential integer rownames (e.g. from `$sample_meta` of an
  `olink_data` object). Default: NULL.

- group_col:

  Character. Column in sample_info to use for point color. Default:
  "group".

- var_type:

  Character. How to treat `group_col`: `"categorical"` uses a discrete
  colour scale (values coerced to factor); `"continuous"` uses a viridis
  gradient (values kept numeric). Default: `"categorical"`.

- shape_col:

  Character or NULL. Column in sample_info to use for point shape.
  Default: NULL (all points same shape).

- dims:

  Integer vector of length 2. Which PCs to plot. Default: c(1, 2).

- ntop:

  Integer. Number of most variable genes to use for PCA. Set to NULL to
  use all genes. Default: 500.

- log_transform:

  Logical. Log2-transform counts before PCA (recommended for raw
  normalized counts). Adds a pseudocount of 1. Default: TRUE.

- colors:

  Named character vector of colors for groups, or NULL for default
  ggplot2 colors. Names should match levels in group_col.

- point_size:

  Numeric. Size of points. Default: 3.

- label_samples:

  Logical or "repel". TRUE for direct labels, "repel" for
  non-overlapping labels (requires ggrepel). Default: FALSE.

- title:

  Character or NULL. Plot title. Default: "PCA".

- verbose:

  Logical. Print PCA summary. Default: TRUE.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
# From artemis_norm object
p <- AETHER_plot_pca(norm_data, group_col = "group")

# With timepoint shapes
p <- AETHER_plot_pca(norm_data, group_col = "group", shape_col = "timepoint")

# From matrix + sample sheet
p <- AETHER_plot_pca(norm_counts, sample_info = targets,
                      group_col = "condition", shape_col = "batch")

} # }
```
