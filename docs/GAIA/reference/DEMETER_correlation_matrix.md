# Calculate correlation matrix with p-values

Computes pairwise correlations between all columns of a
data.frame/matrix, along with p-values for each correlation.

## Usage

``` r
DEMETER_correlation_matrix(
  data,
  method = "spearman",
  alternative = "two.sided",
  adjust = "none",
  use = "pairwise.complete.obs"
)
```

## Arguments

- data:

  Data.frame or matrix with variables as columns.

- method:

  Correlation method: "spearman", "pearson", or "kendall". Default:
  "spearman".

- alternative:

  Alternative hypothesis: "two.sided", "less", "greater". Default:
  "two.sided".

- adjust:

  P-value adjustment method: "BH", "bonferroni", "none", etc. Default:
  "none".

- use:

  How to handle missing values: "pairwise.complete.obs" (default),
  "complete.obs", "everything".

## Value

A list with:

- correlation:

  Correlation matrix

- pvalue:

  P-value matrix

- padj:

  Adjusted p-value matrix (if adjust != "none")

- n:

  Sample size matrix (pairwise)

- method:

  Correlation method used

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
result <- DEMETER_correlation_matrix(data)
result$correlation
result$pvalue

} # }
```
