# Assign four-category significance labels to an isoform switch data.frame

Adds a .cat column: "up", "down", "low_dIF", "non_sig". "up"/"down"
refer to increased/decreased isoform usage in the experiment group. NA
q-values are treated as non-significant.

## Usage

``` r
.iso_assign_categories(df, q_col, alpha, dIF_cutoff)
```
