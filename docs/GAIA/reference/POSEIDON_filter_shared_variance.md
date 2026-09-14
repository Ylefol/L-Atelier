# Filter Features by Shared Variance Across Datasets

For multi-dataset integration, removes genes that have zero variance in
ANY dataset. Ensures all datasets have valid features.

## Usage

``` r
POSEIDON_filter_shared_variance(
  X_list,
  by_row = TRUE,
  variance_threshold = 0,
  verbose = TRUE
)
```

## Arguments

- X_list:

  List of data matrices (all must have same features in rows or columns)

- by_row:

  Logical, are features in rows? (default = TRUE for gene expression)

- variance_threshold:

  Minimum variance (default = 0)

- verbose:

  Print filtering summary (default = TRUE)

## Value

List of filtered matrices with only features having variance \> 0 in ALL
datasets
