# Artemis - Anchor-Based Analysis Functions

Functions for statistical analysis of signal at anchor regions,
comparing conditions (e.g., WT vs KO) at predefined genomic loci. Global
Shift Analysis at Anchor Regions

Tests whether signal at anchor regions is systematically different
between two conditions (e.g., WT vs KO). Useful for answering "Is there
an overall effect at these sites?"

## Usage

``` r
ARTEMIS_global_shift_test(
  quant_result,
  group_col = "group",
  group_a = "WT",
  group_b = "SMUG1_KO",
  summary_method = "mean",
  test = "wilcoxon",
  pseudocount = 1,
  verbose = TRUE
)
```

## Arguments

- quant_result:

  A quantification result from ELEUTHIA_quantify_bed(), containing
  counts matrix, annotation, and targets with group column.

- group_col:

  Character. Column name in targets containing group labels (default =
  "group").

- group_a:

  Character. First group label (default = "WT"). This is the
  reference/baseline group.

- group_b:

  Character. Second group label (default = "SMUG1_KO"). This is the
  comparison group.

- summary_method:

  Character. How to summarize signal per region per group: "mean"
  (default) or "median".

- test:

  Character. Statistical test to use: "wilcoxon" (default, paired
  Wilcoxon signed-rank test) or "t.test" (paired t-test).

- pseudocount:

  Numeric. Added to counts before log transformation to avoid log(0)
  (default = 1).

- verbose:

  Logical. Print results summary (default = TRUE).

## Value

A list containing:

- per_region:

  Data.frame with per-region statistics: region_id, mean_a, mean_b,
  log2FC, raw_diff

- global_test:

  List with test statistic, p.value, and test name

- summary:

  List with summary statistics: n_regions, n_up, n_down, median_log2FC,
  mean_log2FC

- groups:

  Character vector of group names used (c(group_a, group_b))

## Details

The function:

1.  Calculates mean/median signal per region for each group

2.  Computes log2 fold change (group_b / group_a) per region

3.  Tests for systematic shift using paired test across all regions

4.  Returns per-region statistics and global test result

A significant p-value indicates that signal at anchor regions is
systematically higher or lower in group_b compared to group_a.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare ATAC signal at ChIP anchors between WT and KO
result <- ARTEMIS_global_shift_test(atac_at_anchors,
                                     group_a = "WT",
                                     group_b = "SMUG1_KO")

# View summary
result$summary

# Get per-region fold changes
head(result$per_region)

} # }
```
