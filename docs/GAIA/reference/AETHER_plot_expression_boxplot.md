# Faceted Boxplots of Normalized Expression by Condition

Draws faceted boxplots of normalized expression for an arbitrary set of
features, one facet per feature, grouped by condition on the x-axis.
When the feature set is large, this is capped to the top N features by a
chosen DEA statistic (`sort_by`) so the same call works whether
`features` has 10 elements or 300: for small sets the cap has no effect,
for large sets it narrows to the most notable members rather than
producing an unreadable facet grid.

## Usage

``` r
AETHER_plot_expression_boxplot(
  counts,
  features,
  sample_info = NULL,
  sample_col = NULL,
  group_col = "group",
  dea_result = NULL,
  sort_by = "padj",
  top_n = 12,
  log_transform = TRUE,
  colors = NULL,
  title = "Expression by Group",
  text_size = 11,
  point_size = 2,
  show_points = TRUE,
  ncol = NULL,
  verbose = TRUE
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
  `rownames(counts)`).

- sample_info:

  Data.frame of sample metadata, rownames matching `colnames(counts)`.
  Required when `counts` is a plain matrix.

- sample_col:

  Character or NULL. Column in `sample_info` holding sample IDs, used as
  rownames for matching if provided. Default: NULL.

- group_col:

  Character. Column in `sample_info` to group by (x-axis / fill).
  Default: "group".

- dea_result:

  DEA result from
  [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
  (or a data.frame with `feature_id` and `sort_by` columns). Required
  when `length(features) > top_n`, used to rank features for capping.
  Default: NULL.

- sort_by:

  Character. Column in `dea_result` used to rank features before capping
  to `top_n`. Default: "padj" (most statistically confident first; NAs
  sort last). Use "log2FoldChange" to rank by effect size instead
  (ranked by absolute value).

- top_n:

  Integer or NULL. Maximum number of features to show. Default: 12. Set
  to NULL to disable capping and show every feature in `features` (only
  advisable for small sets – see
  [`AETHER_plot_expression_heatmap()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_expression_heatmap.md)
  for a view that scales to large sets without capping).

- log_transform:

  Logical. Log2(x + 1)-transform counts. Default: TRUE.

- colors:

  Named character vector of colors per group. Default: NULL (default
  palette).

- title:

  Character. Plot title. Default: "Expression by Group".

- text_size:

  Numeric. Base font size. Default: 11.

- point_size:

  Numeric. Size of jittered points. Default: 2.

- show_points:

  Logical. Overlay jittered replicate points. Default: TRUE.

- ncol:

  Integer or NULL. Number of facet columns. Default: NULL (auto).

- verbose:

  Logical. Print notes about dropped features, samples, and capping.
  Default: TRUE.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Small family (n <= top_n): every member shown
AETHER_plot_expression_boxplot(norm_data, features = l2_family_ids)

# Large family: capped to the 12 most significant members
AETHER_plot_expression_boxplot(norm_data, features = l1_family_ids,
                                dea_result = dea_result, top_n = 12)
} # }
```
