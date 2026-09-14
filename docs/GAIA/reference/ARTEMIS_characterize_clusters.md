# Characterize Clusters by Variable Discriminative Power

Determines which variables best distinguish between clusters by testing
each variable's association with cluster membership and calculating
effect sizes.

## Usage

``` r
ARTEMIS_characterize_clusters(
  cluster_result,
  data = NULL,
  p_adjust = "BH",
  verbose = TRUE
)
```

## Arguments

- cluster_result:

  An artemis_cluster object from ARTEMIS_cluster_mixed().

- data:

  Optional data frame. If NULL (default), uses data stored in
  cluster_result. Provide if you want to test additional variables not
  used in clustering.

- p_adjust:

  Character. Method for p-value adjustment. One of "BH"
  (Benjamini-Hochberg), "bonferroni", "holm", "none". Default = "BH".

- verbose:

  Logical. Print progress and summary. Default = TRUE.

## Value

A list with class "artemis_characterization" containing:

- results:

  Data frame with per-variable statistics, sorted by effect size

- top_variables:

  Names of variables with large effect sizes

- cluster_profiles:

  Summary statistics per cluster for top variables

- n_clusters:

  Number of clusters

- n_variables:

  Number of variables tested

## Details

For quantitative variables:

- Test: Kruskal-Wallis rank sum test (non-parametric ANOVA)

- Effect size: Eta-squared (η²) = H / (n-1), where H is the test
  statistic

- Interpretation: 0.01 = small, 0.06 = medium, 0.14 = large

For qualitative variables:

- Test: Chi-squared test of independence

- Effect size: Cramér's V = sqrt(χ² / (n \* min(r-1, c-1)))

- Interpretation: 0.1 = small, 0.3 = medium, 0.5 = large

Variables are ranked by effect size. Those with high effect sizes are
the most discriminative - they differ substantially between clusters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Characterize clusters
char <- ARTEMIS_characterize_clusters(cluster_result)

# View top discriminating variables
head(char$results, 10)

# Get variable names with large effects
char$top_variables

} # }
```
