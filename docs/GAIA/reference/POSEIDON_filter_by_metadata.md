# Filter Samples by Metadata Value

General function to filter samples by any metadata column. More flexible
than filter_by_group.

## Usage

``` r
POSEIDON_filter_by_metadata(quant_result, column, values, verbose = TRUE)
```

## Arguments

- quant_result:

  A list containing counts and targets.

- column:

  Character string. Column name in targets to filter by.

- values:

  Vector. Value(s) to keep.

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

Filtered quant_result.

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter by batch
batch1_data <- POSEIDON_filter_by_metadata(data, column = "batch", values = 1)

# Filter by biological replicate
bio1_data <- POSEIDON_filter_by_metadata(data, column = "bio_rep", values = c(1, 2))

} # }
```
