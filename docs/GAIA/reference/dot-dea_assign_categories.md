# Assign four-category significance labels to a DEA data.frame

Adds a .cat column: "up", "down", "low_reg", "non_sig". NA filter values
are treated as non-significant.

## Usage

``` r
.dea_assign_categories(df, filter_choice, p_thresh, l2fc_thresh)
```
