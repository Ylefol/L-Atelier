# Create lower/upper triangular p-value matrix for plotting

Formats a p-value matrix for use with correlation plots, where only the
upper or lower triangle is shown.

## Usage

``` r
DEMETER_format_pvalue_triangle(
  pvalue_matrix,
  triangle = "upper",
  format = "scientific",
  sig_threshold = 0.05,
  hide_ns = TRUE
)
```

## Arguments

- pvalue_matrix:

  Square matrix of p-values.

- triangle:

  Which triangle to keep: "upper" or "lower". Default: "upper".

- format:

  Format for p-values (passed to DEMETER_format_pvalues).

- sig_threshold:

  Significance threshold.

- hide_ns:

  Hide non-significant values.

## Value

Character matrix with formatted p-values in specified triangle, empty
strings elsewhere.
