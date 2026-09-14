# Filter Samples by Group

Filters a count dataset to include only samples belonging to specified
group(s). Works with the standard quant_result structure.

## Usage

``` r
POSEIDON_filter_by_group(
  quant_result,
  group_col = "group",
  groups,
  verbose = TRUE
)
```

## Arguments

- quant_result:

  A list containing:

  - counts: Matrix of counts (features x samples)

  - targets: Data.frame with sample metadata

- group_col:

  Character string. Column name in targets containing group information
  (default = "group").

- groups:

  Character vector. Group value(s) to keep.

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

Filtered quant_result with only samples from specified group(s).

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter for WT samples only
wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")

# Filter for multiple groups
subset_data <- POSEIDON_filter_by_group(data, groups = c("WT", "Control"))

} # }
```
